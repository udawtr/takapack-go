package takapack

import (
	"fmt"
	"testing"

	"github.com/james-bowman/sparse"
)

// takapackTest demonstrates the use of LU factorization and CG methods
func Test_takapack(t *testing.T) {
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
	// TraceMat(mat_csr)

	b := []float64{8, 45, -3, 3, 19}

	// LU factorization and solving
	Lu := LUFactorization(mat_csr)
	TraceMat(mat_csr)
	x1 := LUSolve(Lu, b)
	fmt.Printf("%f %f %f %f %f\n", x1[0], x1[1], x1[2], x1[3], x1[4])

	// CG sparse solve
	x2 := CGSparseSolve(mat_csr, b, 0.00000001)
	fmt.Printf("%f %f %f %f %f\n", x2[0], x2[1], x2[2], x2[3], x2[4])
}
