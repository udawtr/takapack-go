package takapack

import (
	"fmt"
	"math"
)

// takapackCGSparseSolve solves a sparse linear system using the Conjugate Gradient method
func CGSparseSolve(N int, Ap, Ai []int, Ax, b []float64, threshold float64) []float64 {
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

// Pair is a utility structure to hold a pair of values
type Pair struct {
	first  int
	second float64
}

func LUFactorization(N int, Ap, Ai []int, Ax []float64) (LUp, LUi []int, LUx []float64, LUrowFlip []int) {
	//compressed row form用の LU分解
	LUrowFlip = make([]int, N)
	for i := range LUrowFlip {
		LUrowFlip[i] = i
	}

	nonZeroEntryNum := 0
	Row := make([][]Pair, N)
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
				Row[y] = append(Row[y], Pair{I, Myi[y]})
			}
		}
	}

	//最後に Row を Compressed row formにまとめる
	LUp = make([]int, N+1)
	LUx = make([]float64, nonZeroEntryNum)
	LUi = make([]int, nonZeroEntryNum)

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

	return LUp, LUi, LUx, LUrowFlip
}

// takapackTraceMat prints the contents of a sparse matrix
func TraceMat(N int, Ap, Ai []int, Ax []float64) {
	fmt.Println("\ntakapack_tracemat")

	for y := 0; y < N; y++ {
		idx := Ap[y]
		for x := 0; x < N; x++ {
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

func LUSolve(N int, LUp, LUi []int, LUx []float64, LUrowFlip []int, b []float64) []float64 {
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
