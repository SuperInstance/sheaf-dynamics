//! Sheaf dynamics: cellular sheaves on graphs with Laplacian-driven diffusion.
//! The sheaf Laplacian generalizes the graph Laplacian to structured data on stalks.

/// A stalk: the data type attached to a node (vector space dimension)
pub type StalkDim = usize;

/// Assignment: maps edge data to node data (restriction map as a matrix)
#[derive(Clone)]
pub struct RestrictionMap {
    pub from_node: usize,
    pub to_node: usize,
    pub matrix: Vec<Vec<f64>>,
}

impl RestrictionMap {
    pub fn identity(node: usize, dim: usize) -> Self {
        Self {
            from_node: node,
            to_node: node,
            matrix: (0..dim).map(|i| (0..dim).map(|j| if i == j { 1.0 } else { 0.0 }).collect()).collect(),
        }
    }

    pub fn apply(&self, v: &[f64]) -> Vec<f64> {
        let rows = self.matrix.len();
        let cols = if rows > 0 { self.matrix[0].len() } else { 0 };
        let n = v.len().min(cols);
        (0..rows).map(|i| {
            (0..n).map(|j| self.matrix[i][j] * v[j]).sum::<f64>()
        }).collect()
    }
}

/// A cellular sheaf on a graph
pub struct Sheaf {
    pub n_nodes: usize,
    pub stalk_dims: Vec<StalkDim>,
    pub restrictions: Vec<RestrictionMap>,
    pub adj: Vec<Vec<f64>>, // underlying graph adjacency
}

impl Sheaf {
    pub fn new(adj: Vec<Vec<f64>>, stalk_dims: Vec<StalkDim>) -> Self {
        let n = adj.len();
        let mut restrictions = Vec::new();
        for i in 0..n {
            for j in (i+1)..n {
                if adj[i][j] > 0.0 {
                    // Default: identity restriction maps
                    restrictions.push(RestrictionMap {
                        from_node: i, to_node: j,
                        matrix: (0..stalk_dims[j]).map(|r|
                            (0..stalk_dims[i]).map(|c| if r == c { 1.0 } else { 0.0 }).collect()
                        ).collect(),
                    });
                    restrictions.push(RestrictionMap {
                        from_node: j, to_node: i,
                        matrix: (0..stalk_dims[i]).map(|r|
                            (0..stalk_dims[j]).map(|c| if r == c { 1.0 } else { 0.0 }).collect()
                        ).collect(),
                    });
                }
            }
        }
        Self { n_nodes: n, stalk_dims, restrictions, adj }
    }

    /// Total dimension of the sheaf (sum of all stalk dimensions)
    pub fn total_dim(&self) -> usize {
        self.stalk_dims.iter().sum()
    }

    /// Build the sheaf Laplacian as a block matrix
    /// L = L_low + L_high where L_low is the graph Laplacian weighted by stalk dims
    /// and L_high encodes the restriction map compatibility
    pub fn sheaf_laplacian(&self) -> Vec<Vec<f64>> {
        let n = self.total_dim();
        let mut lap = vec![vec![0.0_f64; n]; n];

        let offsets = self.compute_offsets();

        for rm in &self.restrictions {
            let di = self.stalk_dims[rm.from_node];
            let dj = self.stalk_dims[rm.to_node];
            let oi = offsets[rm.from_node];
            let oj = offsets[rm.to_node];

            // F^T F on diagonal block of from_node
            for r in 0..di {
                for s in 0..di {
                    let val: f64 = (0..dj).map(|k| {
                        rm.matrix[k].get(r).copied().unwrap_or(0.0) * rm.matrix[k].get(s).copied().unwrap_or(0.0)
                    }).sum();
                    lap[oi + r][oi + s] += val;
                }
            }

            // -F^T on off-diagonal block (from_node, to_node)
            for r in 0..di {
                for s in 0..dj {
                    lap[oi + r][oj + s] -= rm.matrix[s].get(r).copied().unwrap_or(0.0);
                }
            }

            // -F on off-diagonal block (to_node, from_node)
            for r in 0..dj {
                for s in 0..di {
                    lap[oj + r][oi + s] -= rm.matrix[r].get(s).copied().unwrap_or(0.0);
                }
            }

            // FF^T on diagonal block of to_node
            for r in 0..dj {
                for s in 0..dj {
                    let val: f64 = (0..di).map(|k| {
                        rm.matrix[r].get(k).copied().unwrap_or(0.0) * rm.matrix[s].get(k).copied().unwrap_or(0.0)
                    }).sum();
                    lap[oj + r][oj + s] += val;
                }
            }
        }

        lap
    }

    fn compute_offsets(&self) -> Vec<usize> {
        let mut offsets = vec![0usize; self.n_nodes];
        for i in 1..self.n_nodes {
            offsets[i] = offsets[i-1] + self.stalk_dims[i-1];
        }
        offsets
    }

    /// Eigenvalues of the sheaf Laplacian
    pub fn eigenvalues(&self) -> Vec<f64> {
        let mut lap = self.sheaf_laplacian();
        jacobi(&mut lap)
    }

    /// Sheaf CR: λ₂/λₙ of the sheaf Laplacian
    pub fn cr(&self) -> f64 {
        let eigs = self.eigenvalues();
        if eigs.len() < 2 { return 0.0; }
        let l2 = eigs[1];
        let ln = *eigs.last().unwrap_or(&1.0);
        if ln <= 0.0 { 0.0 } else { l2 / ln }
    }

    /// Diffuse a signal on the sheaf using the sheaf Laplacian
    /// dx/dt = -L_s x, Euler method
    pub fn diffuse(&self, signal: &[f64], dt: f64, steps: usize) -> Vec<f64> {
        let lap = self.sheaf_laplacian();
        let n = signal.len();
        let mut x = signal.to_vec();
        for _ in 0..steps {
            let mut dx = vec![0.0; n];
            for i in 0..n {
                for j in 0..n {
                    dx[i] -= lap[i][j] * x[j];
                }
            }
            for i in 0..n {
                x[i] += dt * dx[i];
            }
        }
        x
    }

    /// Global sections: the null space of the sheaf Laplacian (dimension)
    pub fn global_section_dim(&self) -> usize {
        let eigs = self.eigenvalues();
        eigs.iter().filter(|&&e| e < 1e-10).count()
    }
}

fn jacobi(a: &mut Vec<Vec<f64>>) -> Vec<f64> {
    let n = a.len();
    if n == 0 { return vec![]; }
    for _ in 0..100 * n * n {
        let (mut p, mut q) = (0, 1);
        let mut max_val = 0.0_f64;
        for i in 0..n { for j in (i+1)..n { if a[i][j].abs() > max_val { max_val = a[i][j].abs(); p = i; q = j; } } }
        if max_val < 1e-14 { break; }
        let app = a[p][p]; let aqq = a[q][q]; let apq = a[p][q];
        let theta = if (app - aqq).abs() < 1e-30 { std::f64::consts::FRAC_PI_4 }
                     else { 0.5 * (2.0 * apq / (app - aqq)).atan() };
        let (c, s) = (theta.cos(), theta.sin());
        for i in 0..n { if i != p && i != q { let aip = a[i][p]; let aiq = a[i][q]; a[i][p] = c*aip+s*aiq; a[p][i]=a[i][p]; a[i][q]=-s*aip+c*aiq; a[q][i]=a[i][q]; } }
        a[p][p] = c*c*app+2.0*s*c*apq+s*s*aqq;
        a[q][q] = s*s*app-2.0*s*c*apq+c*c*aqq;
        a[p][q] = 0.0; a[q][p] = 0.0;
    }
    let mut eigs: Vec<f64> = (0..n).map(|i| a[i][i]).collect();
    eigs.sort_by(|a, b| a.partial_cmp(b).unwrap());
    eigs
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn total_dim_correct() {
        let s = Sheaf::new(
            vec![vec![0.0, 1.0], vec![1.0, 0.0]],
            vec![3, 2],
        );
        assert_eq!(s.total_dim(), 5);
    }

    #[test]
    fn sheaf_laplacian_positive_semidefinite() {
        let s = Sheaf::new(
            vec![vec![0.0, 1.0, 0.0], vec![1.0, 0.0, 1.0], vec![0.0, 1.0, 0.0]],
            vec![2, 2, 2],
        );
        let eigs = s.eigenvalues();
        for &e in &eigs {
            assert!(e >= -1e-10, "Eigenvalue {} is negative", e);
        }
    }

    #[test]
    fn sheaf_has_zero_eigenvalue() {
        let s = Sheaf::new(
            vec![vec![0.0, 1.0], vec![1.0, 0.0]],
            vec![2, 2],
        );
        let eigs = s.eigenvalues();
        assert!(eigs[0].abs() < 1e-10, "First eigenvalue should be ~0: {}", eigs[0]);
    }

    #[test]
    fn cr_in_range() {
        let s = Sheaf::new(
            vec![vec![0.0, 1.0, 1.0], vec![1.0, 0.0, 1.0], vec![1.0, 1.0, 0.0]],
            vec![2, 2, 2],
        );
        let cr = s.cr();
        assert!(cr >= 0.0 && cr <= 1.0, "CR = {} out of range", cr);
    }

    #[test]
    fn diffusion_smooths_signal() {
        let s = Sheaf::new(
            vec![vec![0.0, 1.0], vec![1.0, 0.0]],
            vec![2, 2],
        );
        let signal = vec![10.0, 0.0, -10.0, 0.0]; // sharp difference
        let diffused = s.diffuse(&signal, 0.01, 100);
        // After diffusion, signal should be smoother (less variance)
        let var_before: f64 = signal.iter().map(|x| x * x).sum::<f64>() / signal.len() as f64;
        let var_after: f64 = diffused.iter().map(|x| x * x).sum::<f64>() / diffused.len() as f64;
        assert!(var_after < var_before, "Variance should decrease: {} vs {}", var_after, var_before);
    }

    #[test]
    fn global_sections_constant_sheaf() {
        // Constant sheaf on connected graph: 1 global section
        let s = Sheaf::new(
            vec![vec![0.0, 1.0, 1.0], vec![1.0, 0.0, 1.0], vec![1.0, 1.0, 0.0]],
            vec![1, 1, 1],
        );
        let dim = s.global_section_dim();
        assert!(dim >= 1, "Connected graph should have at least 1 global section, got {}", dim);
    }

    #[test]
    fn restriction_map_identity() {
        let rm = RestrictionMap::identity(0, 3);
        let v = vec![1.0, 2.0, 3.0];
        let result = rm.apply(&v);
        assert_eq!(result, v);
    }

    #[test]
    fn restriction_map_applies() {
        let rm = RestrictionMap {
            from_node: 0, to_node: 1,
            matrix: vec![vec![1.0, 0.0], vec![0.0, 1.0]],
        };
        let v = vec![3.0, 4.0];
        let result = rm.apply(&v);
        assert!((result[0] - 3.0).abs() < 1e-10);
        assert!((result[1] - 4.0).abs() < 1e-10);
    }
}
