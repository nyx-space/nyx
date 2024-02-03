# Change Log
All notable changes to Nyx starting with the version 1.0.0 will be documented here.

This project adheres to [Semantic Versioning](https://semver.org/), unless a given technical breaking change is unlikely to be one.

## 1.0.1
### Unlikely breaking changes
- NyxError enum no longer has `OutOfInterpolationWindow` or `TrajectoryCreationError`. These are now part of the more detailed `TrajError` error enum.