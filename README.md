# takapack(go)

http://takashiijiri.com/study/miscs/takapack.html より、

　
- 疎な連立方程式を解く
  -  与える行列は, umfpackと同様の comressed row form (umfpackではcompressed column formだった)
  - CG法による　　連立方程式 Ax = b の 求解をサポート
  - LU分解による 連立方程式 Ax = b の 求解をサポート
　　
- NYSL（Version 0.9982）ライセンスなので、好きに使ってください（大学のレポートとかの参考になれば良いなと思っています.）．

## 使い方
　
### 1) 解きたい連立方程式 Ax = b の行列Aをcompressed row formで作る。
　　　　配列bもつくる。

```
/*----------------------------
    2  3  0  0  0        8
    3  0  4  0  6       45
A =   0 -1 -3  2  0    b =-3
    0  0  1  0  0        3
    0  4  2  0  1       19
----------------------------*/

    //(umfpackは complessed columnだけど、実装の都合上 Compressed Rowを採用 ())
    // umfpackと併用するときは umfpack_solve の第一引数に UMFPACK_At を食わせて転置すればOK
    //compressed row form of A 
	N := 5
	Ap := []int{0, 2, 5, 8, 9, 12}
	Ai := []int{0, 1, 0, 2, 4, 1, 2, 3, 2, 1, 2, 4}
	Ax := []float64{2.0, 3.0, 3, 4, 6, -1, -3, 2, 1, 4, 2, 1}
    b := []float64{8, 45, -3, 3, 19}
```

### 2) CG法で解くなら CG法の関数を呼ぶ

```
    x := CGSparseSolve(N, Ap, Ai, Ax, b, 0.00000001);
    fmt.Println("Solution:", x)
```

### 3)LU分解で解くなら i)LU分解関数 ii)解く関数 iii)解放する関数 を順に呼ぶ

```
	LUp, LUi, LUx, LU_rowflip := LUFactorization(N, Ap, Ai, Ax)
    x := LUSolve(N, LUp, LUi, LUx, LU_rowflip, b)
    fmt.Println("Solution:", x)
``` 

## 雑な動作確認方法

```
$ go run main.go

takapack_tracemat
2.0000  3.0000  *.****  *.****  *.****  
3.0000  *.****  4.0000  *.****  6.0000  
*.****  -1.0000  -3.0000  2.0000  *.****  
*.****  *.****  1.0000  *.****  *.****  
*.****  4.0000  2.0000  *.****  1.0000  
Solution: [-5.293829947609467 1.2166710160167535 -3.3297217218117185 5.597773881583578 16.722980224896652]
```