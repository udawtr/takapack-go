# Takapack - Sparse Matrix Solver

## Overview

Implementation of sparse LU decomposition and CG method (conjugate gradient method).
Original source code was written by (Takashi Ijiri)[http://takashiijiri.com/study/miscs/takapack.html]

## Features

* Compatible with [Sparse matrix formats module](https://pkg.go.dev/github.com/james-bowman/sparse)'s CSR struct (Compressed Sparse Row).
* Solve Ax=b with LU decomposition
* Solve Ax=b with CG method

## Usage

Tha CSR struct of james-bowman's sparse module is compatible with this package. First create DOK or other fomart matrix and convert to CSR format.
Second just Call takapack's function. Just it.


```
// Construct a new 5x5 DOK (Dictionary Of Keys) matrix
mat_dok := sparse.NewDOK(5, 5)

// Populate it with some non-zero values
//       2  3  0  0  0        8
//       3  0  4  0  6       45
// A =   0 -1 -3  2  0    b =-3
//       0  0  1  0  0        3
//       0  4  2  0  1       19
mat_dok.Set(0, 0, 2)
mat_dok.Set(0, 1, 3)
mat_dok.Set(1, 0, 3)
mat_dok.Set(1, 2, 4)
mat_dok.Set(1, 4, 6)
mat_dok.Set(2, 1, -1)
mat_dok.Set(2, 2, -3)
mat_dok.Set(2, 3, 2)
mat_dok.Set(3, 2, 1)
mat_dok.Set(4, 1, 4)
mat_dok.Set(4, 2, 2)
mat_dok.Set(4, 4, 1)

// Convert to CSR format
mat_csr := mat_dok.ToCSR()

// Show A
takapack.TraceMat(mat_csr)

// Solve Ax=b
b := []float64{8, 45, -3, 3, 19}
lu := takapack.LUFactorization(mat_csr)
x := LUSolve(lu, b)
fmt.Println("Solution:", x)

// If you want to use CG
// threshold := 0.00000001
// x := takapack.CGSparseSolve(mat_csr, b, threshold)
```

## Installation

With Go installed, package installation is performed using go get.

```
go get -u github.com/udawtr/takapack-go/takapack
```

## Acknowledgements

- [takapack(original)](http://takashiijiri.com/study/miscs/takapack.html)
- Barrett, Richard et al. (1994). Section 2.3.1 Conjugate Gradient Method (CG). In Templates for the Solution of Linear Systems: Building Blocks for Iterative Methods (2nd ed.) (pp. 12-15). Philadelphia, PA: SIAM. Retrieved from http://www.netlib.org/templates/templates.pdf
- Hestenes, M., and Stiefel, E. (1952). Methods of conjugate gradients for solving linear systems. Journal of Research of the National Bureau of Standards, 49(6), 409. doi:10.6028/jres.049.044
- Málek, J. and Strakoš, Z. (2015). Preconditioning and the Conjugate Gradient Method in the Context of Solving PDEs. Philadelphia, PA: SIAM.

## See Also

* [james-bowman/sparse](https://pkg.go.dev/github.com/james-bowman/sparse)


## License

MIT