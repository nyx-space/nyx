# Tech Stack

## Languages
- **Rust (edition 2021):** The core engine of Nyx, providing memory safety, high performance, and robust type-checking.
- **Python:** Used for examples, plotting scripts, and data analysis.

## Core Libraries
- **nalgebra:** Used for precise linear algebra computations (vectors, matrices).
- **hifitime:** Provides high-fidelity time management and conversions.
- **anise:** The underlying engine for ephemerides and frames.
- **serde:** For flexible serialization and deserialization of data structures.
- **rayon:** Enables data-parallelism for performance-critical tasks (e.g., Monte Carlo simulations).
- **parquet / arrow:** High-performance columnar storage formats for large datasets (e.g., simulation results).

## Infrastructure & Tooling
- **GitHub Actions:** Automates the project's CI/CD pipeline, including linting, testing, and documentation deployment.
- **Firebase:** Used for specific configuration and hosting needs.
- **Cargo:** The standard Rust package manager and build tool.