# Space-Filling Designs: MATLAB Implementation

## Overview
This repository provides MATLAB implementations for:
- The incremental construction of space-filling designs, based on the algorithms presented in Chap. 6 of the book:
  "T. Karvonen, L. Pronzato, A. Zhigljavsky: Space-Filling Design and Kernels: Theory and Algorithms, Springer, 2026."
- The evaluation of various space-filling performance criteria for a given design.

---

## Examples in a cube [0,1]^d
The repository includes two main scripts to demonstrate the usage of the provided functions:
- `Examples_generation_cube.m`: Examples of design constructions in the hypercube [0,1]^d.
- `Examples_performance_cube.m`: Examples of performance evaluation for designs in [0,1]^d.

---

## Extending to Other Domains
The provided MATLAB functions can be adapted to generate and evaluate designs in other sets (e.g., balls or simplices) by modifying the files in the following directories:
- `design_generation/`: MATLAB functions for generating designs.
- `design_performance/`: MATLAB functions for evaluating the performance of a given design.
- `auxiliary_functions/`: Auxiliary MATLAB functions used by the main scripts.
In particular, Appendix A.6.2 of [Karvonen, Pronzato and Zhigljavsky, Springer, 2026] indicates how to generate random and quasi random (Sobol') points in a ball, on a sphere, or in a simplex.  

---

## Directory: `high_d`
This directory contains examples and MATLAB files for constructing specific types of designs in high-dimensional spaces:

- **Fractional Factorial Design:**
  - Example of construction of a fractional factorial design with a **large (if not maximum) packing radius**.

- **Optimal Random Designs for Ls-Mean Quantization Error:**
  - Example of construction of optimal random designs for the **Ls-mean quantization error** in:
    - The **unit ball** \( B_d(0,1) \), using uniform designs in a ball \( B_d(0,R) \), or a sphere \( S_{d-1}(0,r) \) (Section 7.1.4 of [Karvonen, Pronzato and Zhigljavsky, Springer, 2026]).
    - The **cube** \([-1,1]^d\), using uniform designs in a cube \([-delta_1,delta_1]^d\), or on the vertices of a (smaller) cube \([-delta_0,delta_0]^d\) (Section 7.2 of [Karvonen, Pronzato and Zhigljavsky, Springer, 2026]).

- **Included Files:**
  - All required `.m` files to determine the **optimal values** of \( R \), \( r \), \( delta_1 \), and \( delta_0 \) as functions of \( d \), \( n \), and \( s \).
---

## Disclaimer
- These MATLAB functions are provided as-is, in the hope that they will be useful, but **without any warranty or guarantee of fitness for a particular purpose**.
- The implementations of the algorithms from "Space-Filling Design and Kernels: Theory and Algorithms" (2026) have not been optimized for efficiency or numerical accuracy. For example, using a Cholesky decomposition for the kernel matrices involved could improve computational efficiency and numerical stability.
- All functions have been written without the assistance of AI tools, and there is likely room for improvement.

---
## Acknowledgments
If you use these functions in your research, please cite the following book:
"T. Karvonen, L. Pronzato, A. Zhigljavsky: Space-Filling Design and Kernels: Theory and Algorithms, Springer, 2026."
