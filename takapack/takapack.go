package takapack

import (
	"fmt"
	"math"

	"github.com/james-bowman/sparse"
)

// takapackCGSparseSolve solves a sparse linear system using the Conjugate Gradient method
func CGSparseSolve(mat *sparse.CSR, b []float64, threshold float64) []float64 {
	rawMat := mat.RawMatrix()
	r, c := rawMat.I, rawMat.J
	Ap := rawMat.Indptr
	Ai := rawMat.Ind
	Ax := rawMat.Data

	if r != c {
		panic("r != c")
	}
	N := r

	result := make([]float64, N)
	cgR := make([]float64, N)
	cgD := make([]float64, N)
	cgQ := make([]float64, N)

	// Initialize r = b - Ax (A is symmetric)
	copy(cgR, b)
	for i := 0; i < N; i++ {
		for j := Ap[i]; j < Ap[i+1]; j++ {
			cgR[i] -= Ax[j] * result[Ai[j]]
		}
	}

	// d = r
	copy(cgD, cgR)

	// deltaNew = r^T * r
	deltaNew := 0.0
	for i := range cgR {
		deltaNew += cgR[i] * cgR[i]
	}

	//deltaZero := deltaNew
	iteration := 0

	for iteration < N && deltaNew > threshold {
		// q = Ad
		for i := range cgQ {
			cgQ[i] = 0
		}
		for i := 0; i < N; i++ {
			for j := Ap[i]; j < Ap[i+1]; j++ {
				cgQ[i] += Ax[j] * cgD[Ai[j]]
			}
		}

		// alpha = deltaNew / (d^T * q)
		alpha := 0.0
		for i := range cgD {
			alpha += cgD[i] * cgQ[i]
		}
		alpha = deltaNew / alpha

		// x = x + alpha * d
		for i := range result {
			result[i] += alpha * cgD[i]
		}

		// r = r - alpha * q
		for i := range cgR {
			cgR[i] -= alpha * cgQ[i]
		}

		// deltaOld = deltaNew
		deltaOld := deltaNew

		// deltaNew = r^T * r
		deltaNew = 0.0
		for i := range cgR {
			deltaNew += cgR[i] * cgR[i]
		}

		// beta = deltaNew / deltaOld
		beta := deltaNew / deltaOld

		// d = r + beta * d
		for i := range cgD {
			cgD[i] = cgR[i] + beta*cgD[i]
		}

		iteration++
	}

	return result
}

// pair is a utility structure to hold a pair of values
type pair struct {
	first  int
	second float64
}

func LUFactorization(mat *sparse.CSR) *SparseLU {
	rawMat := mat.RawMatrix()
	r, c := rawMat.I, rawMat.J
	Ap := rawMat.Indptr
	Ai := rawMat.Ind
	Ax := rawMat.Data

	if r != c {
		panic("r != c")
	}
	N := r

	//compressed row form用の LU分解
	LUrowFlip := make([]int, N)
	for i := range LUrowFlip {
		LUrowFlip[i] = i
	}

	nonZeroEntryNum := 0
	Row := make([][]pair, N)
	MyiTmp := make([]float64, N)
	Myi := make([]float64, N)

	// for each column(各列(縦)について)
	for I := 0; I < N; I++ {
		// 1)-----------I列 (縦)を全て集める : Myi -------------------------------------------------------------
		for y := 0; y < N; y++ {
			MyiTmp[y] = 0
			for kk := Ap[y]; kk < Ap[y+1]; kk++ {
				if Ai[kk] == I {
					MyiTmp[y] = Ax[kk]
					break
				}
			}
		}
		// row flipを適用
		for y := 0; y < N; y++ {
			Myi[y] = MyiTmp[LUrowFlip[y]]
		}

		//2)-----------Myiの値を計算 (この時Aiにアクセスしてはダメ (計算結果は入っていないから))---------------
		for y := 0; y <= I; y++ {
			v := Myi[y]
			for _, pair := range Row[y] {
				if pair.first < y {
					v -= pair.second * Myi[pair.first]
				}
			}
			Myi[y] = v
		}

		big := math.Abs(Myi[I])
		piv := I
		for y := I + 1; y < N; y++ {
			v := Myi[y]
			for _, pair := range Row[y] {
				if pair.first < I {
					v -= pair.second * Myi[pair.first]
				}
			}
			Myi[y] = v
			if tmp := math.Abs(Myi[y]); tmp > big {
				big = tmp
				piv = y
			}
		}

		//3)------------係数掛ける-----------------------
		pivV := Myi[piv]
		coef := 1.0 / pivV
		for y := I; y < N; y++ {
			Myi[y] *= coef
		}
		Myi[piv] = pivV

		//4)------------軸選択(I<-->piv)----------------
		if piv != I {
			Myi[I], Myi[piv] = Myi[piv], Myi[I]
			Row[I], Row[piv] = Row[piv], Row[I]
			LUrowFlip[I], LUrowFlip[piv] = LUrowFlip[piv], LUrowFlip[I]
		}

		//5) -----------値を代入 LU <-- Myi-------------
		for y := 0; y < N; y++ {
			if Myi[y] != 0.0 {
				nonZeroEntryNum++
				Row[y] = append(Row[y], pair{I, Myi[y]})
			}
		}
	}

	//最後に Row を Compressed row formにまとめる
	LUp := make([]int, N+1)
	LUx := make([]float64, nonZeroEntryNum)
	LUi := make([]int, nonZeroEntryNum)

	index := 0
	for i := 0; i < N; i++ {
		LUp[i] = index
		for _, pair := range Row[i] {
			LUi[index] = pair.first
			LUx[index] = pair.second
			index++
		}
	}
	LUp[N] = index

	return &SparseLU{
		N:         N,
		LUp:       LUp,
		LUi:       LUi,
		LUx:       LUx,
		LUrowFlip: LUrowFlip,
	}
}

// takapackTraceMat prints the contents of a sparse matrix
func TraceMat(mat *sparse.CSR) {
	rawMat := mat.RawMatrix()
	r, c := rawMat.I, rawMat.J
	Ap := rawMat.Indptr
	Ai := rawMat.Ind
	Ax := rawMat.Data

	for y := 0; y < r; y++ {
		idx := Ap[y]
		for x := 0; x < c; x++ {
			if idx < Ap[y+1] && x == Ai[idx] {
				fmt.Printf("%.4f  ", Ax[idx])
				idx++
			} else {
				fmt.Printf("*.****  ")
			}
		}
		fmt.Println()
	}
}

type SparseLU struct {
	N         int
	LUp, LUi  []int
	LUx       []float64
	LUrowFlip []int
}

func LUSolve(mat *SparseLU, b []float64) []float64 {
	N := mat.N
	LUrowFlip := mat.LUrowFlip
	LUp := mat.LUp
	LUi := mat.LUi
	LUx := mat.LUx

	// Initialize flipped B and Y
	flippedB := make([]float64, N)
	flippedY := make([]float64, N)
	for i, v := range LUrowFlip {
		flippedB[i] = b[v]
	}

	// Forward substitution to solve L*flippedY = flippedB
	for y := 0; y < N; y++ {
		val := flippedB[y]
		for kk := LUp[y]; kk < LUp[y+1] && LUi[kk] < y; kk++ {
			val -= LUx[kk] * flippedY[LUi[kk]]
		}
		flippedY[y] = val
	}

	// Initialize result X
	result := make([]float64, N)

	// Backward substitution to solve U*X = flippedY
	for y := N - 1; y >= 0; y-- {
		centerIdx := 0
		for kk := LUp[y]; kk < LUp[y+1]; kk++ {
			if LUi[kk] == y {
				centerIdx = kk
				break
			}
		}

		val := flippedY[y]
		for kk := centerIdx + 1; kk < LUp[y+1]; kk++ {
			val -= LUx[kk] * result[LUi[kk]]
		}
		result[y] = val / LUx[centerIdx]
	}

	return result
}
