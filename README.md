# Space-Filling Designs: MATLAB Implementation

## Overview
This repository provides MATLAB implementations for:
- The incremental construction of space-filling designs, based on the algorithms presented in Chaps. 6 and 7 of the book:
  "T. Karvonen, L. Pronzato, A. Zhigljavsky: Space-Filling Design and Kernels: Theory and Algorithms, Springer, 2026."
- The evaluation of various space-filling performance criteria for a given design.

---

## Examples in a cube [0,1]^d (Chap. 6)
The repository includes two main scripts to demonstrate the usage of the provided functions:
- `Examples_generation_cube.m`: Examples of design constructions in the hypercube [0,1]^d;
- `Examples_performance_cube.m`: Examples of performance evaluation for designs in [0,1]^d.

The directory `design_generation/` contains MATLAB functions for generating designs;
the directory `design_performance/` contains MATLAB functions for evaluating the performance of a given design.
Auxiliary MATLAB functions are provided in the directory `auxiliary_functions/` 

---

## Extending to Other Domains
The provided MATLAB functions can be adapted to generate and evaluate designs in other sets (e.g., balls or simplices) by modifying the files in `design_generation/` and `design_performance/`. In particular, Appendix A.6.2 of [Karvonen, Pronzato and Zhigljavsky, Springer, 2026] indicates how to generate random and quasi random (Sobol') points in a ball, on a sphere, or in a simplex.  

---

## High dimension (Chap. 7)
The directory `high_d/` contains examples and MATLAB functions for constructing specific types of designs in high-dimensional spaces:

- **Fractional Factorial Design:**
  - `Example_fractional_factorial.m`: Example of construction of a fractional factorial design with **large** (if not maximum) **packing radius**.

- **Random Designs with small L_s-Mean Quantization Error:**
  - `Example_random_designs_ball_and_cube`: Example of construction of n-point random designs **with small L_s-mean quantization error** of the uniform measure in a ball or a cube.
    - In the **unit ball** B_d(0,1): the design points are uniformly distributed in a ball B_d(0,R), or on a sphere S_{d-1}(0,r) (Section 7.1.4 of [Karvonen, Pronzato and Zhigljavsky, Springer, 2026]).
    - In the **cube** [-1,1]^d: the design points are uniformly distributed in a cube [-delta_1,delta_1]^d, or are sampled from the vertices of a (smaller) cube [-delta_0,delta_0]^d (Section 7.2 of [Karvonen, Pronzato and Zhigljavsky, Springer, 2026]).

All `.m` files required to determine the **optimal values** of R, r, delta_1, and delta_0 as functions of d, n, and s are contained in the directory `high_d/`.
  
---

## Disclaimer
- These MATLAB functions are provided as-is, in the hope that they will be useful, but **without any warranty or guarantee of fitness for a particular purpose**.
- The implementations of the algorithms from "Space-Filling Design and Kernels: Theory and Algorithms" (2026) have not been optimized for efficiency or numerical accuracy. For example, using a Cholesky decomposition for the kernel matrices involved could improve computational efficiency and numerical stability.
- All functions have been written without the assistance of AI tools, and there is likely room for improvement.

---
## Acknowledgments
If you use these functions in your research, please cite the following book:
"T. Karvonen, L. Pronzato, A. Zhigljavsky: Space-Filling Design and Kernels: Theory and Algorithms, Springer, 2026."
