package takapack

import (
	"fmt"
	"testing"
)

// takapackTest demonstrates the use of LU factorization and CG methods
func Test_takapack(t *testing.T) {

	/*----------------------------
		2  3  0  0  0        8
		3  0  4  0  6       45
	A =	0 -1 -3  2  0    b =-3
		0  0  1  0  0        3
		0  4  2  0  1       19
	----------------------------*/

	// Matrix A and vector b are defined here
	//(umfpackは complessed columnだけど、実装の都合上 Compressed Rowを採用 ())
	// umfpackと併用するときは umfpack_solve の第一引数に UMFPACK_At を食わせて転置すればOK
	//compressed row form of A 	N := 5
	N := 5
	Ap := []int{0, 2, 5, 8, 9, 12}
	Ai := []int{0, 1, 0, 2, 4, 1, 2, 3, 2, 1, 2, 4}
	Ax := []float64{2.0, 3.0, 3, 4, 6, -1, -3, 2, 1, 4, 2, 1}

	TraceMat(N, Ap, Ai, Ax)

	b := []float64{8, 45, -3, 3, 19}
	x1 := make([]float64, 5)
	x2 := make([]float64, 5)

	// LU factorization and solving
	LUp, LUi, LUx, LU_rowflip := LUFactorization(N, Ap, Ai, Ax)
	TraceMat(N, LUp, LUi, LUx)
	x1 = LUSolve(N, LUp, LUi, LUx, LU_rowflip, b)
	fmt.Printf("%f %f %f %f %f\n", x1[0], x1[1], x1[2], x1[3], x1[4])

	// CG sparse solve
	x2 = CGSparseSolve(N, Ap, Ai, Ax, b, 0.00000001)
	fmt.Printf("%f %f %f %f %f\n", x2[0], x2[1], x2[2], x2[3], x2[4])
}
