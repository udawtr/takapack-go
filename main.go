package main

import (
	"fmt"

	"github.com/james-bowman/sparse"
	"github.com/udawtr/takapack_go/takapack"
)

func main() {

	// Example usage
	Ap := []int{0, 2, 5, 8, 9, 12}
	Ai := []int{0, 1, 0, 2, 4, 1, 2, 3, 2, 1, 2, 4}
	Ax := []float64{2.0, 3.0, 3, 4, 6, -1, -3, 2, 1, 4, 2, 1}

	mat := sparse.NewCSR(5, 5, Ap, Ai, Ax)

	takapack.TraceMat(mat)

	b := []float64{8, 45, -3, 3, 19}
	threshold := 0.00000001
	x := takapack.CGSparseSolve(mat, b, threshold)

	fmt.Println("Solution:", x)
}
