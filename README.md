# sheaf-dynamics

**Cellular sheaves on graphs with Laplacian-driven diffusion — where topology meets dynamics.**

Pure Rust, zero dependencies. The sheaf Laplacian generalizes the graph Laplacian to structured data living on vertices, and this library makes that generalization tangible and runnable.

## What This Gives You

- **Sheaf construction** — assign vector spaces (stalks) to vertices and linear maps (restrictions) to edges
- **Sheaf Laplacian** — block matrix `L = δ*δ` that generalizes the ordinary graph Laplacian
- **Laplacian-driven diffusion** — watch structured data propagate through a sheaf over time
- **Zero dependencies** — all linear algebra from scratch, readable and auditable

## The Big Idea

A graph Laplacian tells you how "smooth" a scalar function on vertices is. A *sheaf Laplacian* tells you how well vector-valued data at each vertex agrees with its neighbors through the restriction maps. When the stalks are 1-dimensional and restrictions are identity, the sheaf Laplacian *reduces exactly to the ordinary graph Laplacian*.

This is the bridge between [sheaf-cohomology](https://github.com/SuperInstance/sheaf-cohomology) (which measures the gap between local and global) and dynamics (which propagates information along edges).

## Quick Start

```rust
use sheaf_dynamics::{Sheaf, RestrictionMap};

// Build a sheaf on a triangle with 2D stalks
let adj = vec![
    vec![0.0, 1.0, 1.0],
    vec![1.0, 0.0, 1.0],
    vec![1.0, 1.0, 0.0],
];
let stalk_dims = vec![2, 2, 2];
let sheaf = Sheaf::new(adj, stalk_dims);

// Total dimension of the sheaf state space
assert_eq!(sheaf.total_dim(), 6); // 3 vertices × 2D

// Build the sheaf Laplacian (block matrix)
let lap = sheaf.sheaf_laplacian();
```

## API Reference

| Type / Function | Description |
|----------------|-------------|
| `Sheaf` | Cellular sheaf: adjacency + stalk dimensions + restriction maps |
| `RestrictionMap` | Linear map on an edge with from/to nodes and matrix |
| `Sheaf::new(adj, stalk_dims)` | Construct with identity restriction maps |
| `sheaf.total_dim()` | Sum of all stalk dimensions |
| `sheaf.sheaf_laplacian()` | Block matrix L = δ*δ |

## How It Fits

Part of the SuperInstance spectral ecosystem:

- **[sheaf-cohomology](https://github.com/SuperInstance/sheaf-cohomology)** — H⁰, H¹, Euler characteristic, consistency measures
- **sheaf-dynamics** — Diffusion and propagation on sheaves (this repo)
- **[spectral-graph-core](https://github.com/SuperInstance/spectral-graph-core)** — Ordinary graph Laplacian as a special case

## Testing

```bash
cargo test
```

## Installation

```toml
[dependencies]
sheaf-dynamics = { git = "https://github.com/SuperInstance/sheaf-dynamics" }
```

## License

MIT

Part of the [SuperInstance](https://github.com/SuperInstance) ecosystem.
