
# Space-Filling Designs: MATLAB Implementation

## Overview
This repository provides MATLAB implementations for:
- The incremental construction of space-filling designs, based on the algorithms presented in the book:
  "T. Karvonen, L. Pronzato, A. Zhigljavsky: Space-Filling Design and Kernels: Theory and Algorithms, Springer, 2026."
- The evaluation of various space-filling performance criteria for a given design.
  

---

## Examples
The repository includes two main scripts to demonstrate the usage of the provided functions:
- `Examples_generation_cube.m`: Examples of design constructions in the hypercube [0,1]^d.
- `Examples_performance_cube.m`: Examples of performance evaluation for designs in [0,1]^d.

---

## Extending to Other Domains
The provided MATLAB functions can be adapted to generate and evaluate designs in other sets (e.g., balls or simplices) by modifying the files in the following directories:
- `design_generation/`: MATLAB functions for generating designs.
- `design_performance/`: MATLAB functions for evaluating the performance of a given design.
- `auxiliary_functions/`: Auxiliary MATLAB functions used by the main scripts.

---

## Disclaimer
- These MATLAB functions are provided as-is, in the hope that they will be useful, but **without any warranty or guarantee of fitness for a particular purpose**.
- The implementations of the algorithms from "Space-Filling Design and Kernels: Theory and Algorithms" (2026) have not been optimized for efficiency or numerical accuracy. For example, using a Cholesky decomposition for the kernel matrices involved could improve computational efficiency and numerical stability.
- All functions have been written without the assistance of AI tools, and there is likely room for improvement.

---
## Acknowledgments
If you use these functions in your research, please cite the following book:
"T. Karvonen, L. Pronzato, A. Zhigljavsky: Space-Filling Design and Kernels: Theory and Algorithms, Springer, 2026."
