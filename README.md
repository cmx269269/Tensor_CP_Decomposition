# Tensor_CP_Decomposition

[Introduction]
 version 0: 
A MATLAB software for tensor CP decomposition
Copyright (c) 2026 by Minru Bai and Mingxin Chen.

This is software package for solving a group sparse regularized least squares problem of the form: 
(P) min A,B,C {f(A,B,C) + lambda* ||C||_{2,0}: A,B have unit column},
where f(A,B,C) = norm(X - [A,B,C],‘fro’)^2.

The main algorithm is a block coordinate proximal gradient algorithm with extrapolation to (P); details can be found in the following 
reference: 







This repository provides the MATLAB implementation of the method proposed in the paper "GROUP SPARSE-BASED CP DECOMPOSITION FOR INCOMPLETE TENSORS VIA ADAPTIVE PARAMETER SELECTION"

The proposed adaptive parameter selection method automatically estimates the regularization parameter from the input data. The estimated parameter generally yields satisfactory results in most cases. However, slight adjustments to the estimated parameter may be necessary to achieve better numerical performance in certain cases.

The main script is V_test_OP_sc.m. Run this script in MATLAB to reproduce the numerical experiments.
