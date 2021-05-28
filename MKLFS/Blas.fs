module Blas

#nowarn "9" // Disable warning when using pointers.

open EnumTypes
open System.Numerics
open System.Security
open System.Runtime.InteropServices
open System.Runtime.CompilerServices

//TODO Use correct filename depending on OS.
[<Literal>]
let DLL = "mkl_rt.1.dll"

//LEVEL 1 DOUBLE

[<DllImport(DLL, CallingConvention = CallingConvention.Cdecl, ExactSpelling = true)>]
extern float cblas_dasum(int N, float* X, int incX)
let dasum N (X:float array) incX iniX = 
    use xp = fixed & (X.[iniX]) 
    cblas_dasum(N, xp, incX)

[<DllImport(DLL, CallingConvention = CallingConvention.Cdecl, ExactSpelling = true)>]
extern void cblas_daxpy(int N, float a, float* X, int incX, float* Y, int incY)
let daxpy N a (X:float array) incX iniX (Y:float array) incY iniY = 
    use xp = fixed &X.[iniX]
    use yp = fixed &Y.[iniY]
    cblas_daxpy(N, a, xp, incX, yp, incY)

[<DllImport(DLL, CallingConvention = CallingConvention.Cdecl, ExactSpelling = true)>]
extern float cblass_dcabs1(Complex* Z)
let dcabs1 (Z:Complex array) iniZ = 
    use zp = fixed &Z.[iniZ]
    cblass_dcabs1(zp)

[<DllImport(DLL, CallingConvention = CallingConvention.Cdecl, ExactSpelling = true)>]
extern void cblas_dcopy(int N, float* X, int incX, float* Y, int incY)
let dcopy N (X:float array) incX iniX (Y:float array) incY iniY = 
    use xp = fixed &X.[iniX]
    use yp = fixed &Y.[iniY]
    cblas_dcopy(N, xp, incX, yp, incY)

[<DllImport(DLL, CallingConvention = CallingConvention.Cdecl, ExactSpelling = true)>]
extern float cblas_ddot(int N, float* X, int incX, float* Y, int incY)
let ddot N (X:float array) incX iniX (Y:float array) incY iniY = 
    use xp = fixed &X.[iniX]
    use yp = fixed &Y.[iniY]
    cblas_ddot(N, xp, incX, yp, incY)

[<DllImport(DLL, CallingConvention = CallingConvention.Cdecl, ExactSpelling = true)>]
extern float cblas_dnrm2(int N, float* X, int incX)
let dnrm2 N (X:float array) incX iniX = 
    use xp = fixed &X.[iniX]
    cblas_dnrm2(N, xp, incX)

[<DllImport(DLL, CallingConvention = CallingConvention.Cdecl, ExactSpelling = true)>]
extern void cblas_drot(int N, float* X, int incX, float* Y, int incY, float C, float S)
let drot N (X:float array) incX iniX (Y:float array) incY iniY C S = 
    use xp = fixed &X.[iniX]
    use yp = fixed &Y.[iniY]
    cblas_drot(N, xp, incX, yp, incY, C, S)

//TODO drot
//[<DllImport(DLL, CallingConvention = CallingConvention.Cdecl)>]
//extern void cblas_drotg(float& a, float& b, float& c, float& s)
//let drotg (a:byref<float>) (b:byref<float>) (c:byref<float>) (s:byref<float>) = cblas_drotg((a:byref<float>), (b:byref<float>), (c:byref<float>), (s:byref<float>))

[<DllImport(DLL, CallingConvention = CallingConvention.Cdecl, ExactSpelling = true)>]
extern void cblas_drotm(int N, float* X, int incX, float* Y, int incY, float* param)
let drotm N (X:float array) incX iniX (Y:float array) incY iniY param = 
    use xp = fixed &X.[iniX]
    use yp = fixed &Y.[iniY]
    cblas_drotm(N, xp, incX, yp, incY, param)

//TODO drotmg
//[<DllImport(DLL, CallingConvention = CallingConvention.Cdecl, ExactSpelling = true)>]
//extern void cblas_drotmg(float& d1, float& d2, float& x1, float y1, float* param)

[<DllImport(DLL, CallingConvention = CallingConvention.Cdecl, ExactSpelling = true)>]
extern void cblas_dscal(int N, float a, float* X, int incX)
let dscal N a (X:float array) incX iniX = 
    use xp = fixed &X.[iniX]
    cblas_dscal(N, a, xp, incX)

[<DllImport(DLL, CallingConvention = CallingConvention.Cdecl, ExactSpelling = true)>]
extern float cblas_dsdot(int N, float* X, int incX, float* Y, int incY)
let dsdot N (X:float array) incX iniX (Y:float array) incY iniY= 
    use xp = fixed &X.[iniX]
    use yp = fixed &Y.[iniY]
    cblas_dsdot(N, xp, incX, yp, incY)

[<DllImport(DLL, CallingConvention = CallingConvention.Cdecl, ExactSpelling = true)>]
extern void cblas_dswap(int N, float* X, int incX, float* Y, int incY);
let dswap N (X:float array) incX iniX (Y:float array) incY iniY = 
    use xp = fixed &X.[iniX]
    use yp = fixed &Y.[iniY]
    cblas_dswap(N, xp, incX, yp, incY)

[<DllImport(DLL, CallingConvention = CallingConvention.Cdecl, ExactSpelling = true)>]
extern float cblass_dzasum(int N, Complex* X, int incX)
let dzasum N X incX = cblass_dzasum(N, X, incX)

[<DllImport(DLL, CallingConvention = CallingConvention.Cdecl, ExactSpelling = true)>]
extern float cblass_dznrm2(int N, Complex* X, int incX)
let dznrm2 N X incX = cblass_dznrm2(N, X, incX)

//LEVEL 2 DOUBLE

[<DllImport(DLL, CallingConvention = CallingConvention.Cdecl, ExactSpelling = true)>]
extern void cblas_dgbmv(Layout Layout,
    Trans TransA, int M, int N,
    int KL, int KU, float alpha,
    float* A, int lda, float* X,
    int incX, float beta, float* Y, int incY)
let dgmbv layout transA M N KL KU alpha (A:float array) lda (X:float array) incX iniX beta (Y:float array) incY iniY =
    use ap = fixed &A.[0]
    use xp = fixed &X.[iniX]
    use yp = fixed &Y.[iniY]
    cblas_dgbmv(layout, transA, M, N, KL, KU, alpha, ap, lda, xp, incX, beta, yp, incY)

[<DllImport(DLL, CallingConvention = CallingConvention.Cdecl, ExactSpelling = true)>]
extern void cblas_dgemv(Layout Layout,
    Trans TransA, int M, int N,
    float alpha, float* A, int lda,
    float* X, int incX, float beta,
    float* Y, int incY)
let dgemv layout transA M N alpha (A:float array) lda (X:float array) incX iniX beta (Y:float array) incY iniY =
    use ap = fixed &A.[0]
    use xp = fixed &X.[iniX]
    use yp = fixed &Y.[iniY]
    cblas_dgemv(layout, transA, M, N, alpha, ap, lda, xp, incX, beta, yp, incY)

[<DllImport(DLL, CallingConvention = CallingConvention.Cdecl, ExactSpelling = true)>]
extern void cblas_dger(Layout Layout, int M, int N,
    float alpha, float* X, int incX,
    float* Y, int incY, float* A, int lda)
let dger layout M N alpha (X:float array) incX iniX (Y:float array) incY iniY (A:float array) lda =
    use ap = fixed &A.[0]
    use xp = fixed &X.[iniX]
    use yp = fixed &Y.[iniY]
    cblas_dger(layout, M, N, alpha, xp, incX, yp, incY, ap, lda)

[<DllImport(DLL, CallingConvention = CallingConvention.Cdecl, ExactSpelling = true)>]
extern void cblas_dsbmv(Layout Layout, UpLo UpLo,
    int N, int K, float alpha, float* A,
    int lda, float* X, int incX,
    float beta, float* Y, int incY)
let dsbmv layout uplo N K alpha (A:float array) lda (X:float array) incX iniX beta (Y:float array) incY iniY =
    use ap = fixed &A.[0]
    use xp = fixed &X.[iniX]
    use yp = fixed &Y.[iniY]
    cblas_dsbmv(layout, uplo, N, K, alpha, ap, lda, xp, incX, beta, yp, incY)

[<DllImport(DLL, CallingConvention = CallingConvention.Cdecl, ExactSpelling = true)>]
extern void cblas_dspmv(Layout Layout, UpLo UpLo,
    int N, float alpha, float* Ap,
    float* X, int incX,
    float beta, float* Y, int incY)
let dspmv layout uplo N alpha (A:float array) (X:float array) incX iniX beta (Y:float array) incY iniY =
    use ap = fixed &A.[0]
    use xp = fixed &X.[iniX]
    use yp = fixed &Y.[iniY]
    cblas_dspmv(layout, uplo, N, alpha, ap, xp, incX, beta, yp, incY)

[<DllImport(DLL, CallingConvention = CallingConvention.Cdecl, ExactSpelling = true)>]
extern void cblas_dspr(Layout Layout, UpLo UpLo,
    int N, float alpha, float* X,
    int incX, float* Ap)
let dspr layout uplo N alpha (X:float array) incX iniX (Ap:float array) =
    use ap = fixed &Ap.[0]
    use xp = fixed &X.[iniX]
    cblas_dspr(layout, uplo, N, alpha, xp, incX, ap)

[<DllImport(DLL, CallingConvention = CallingConvention.Cdecl, ExactSpelling = true)>]
extern void cblas_dspr2(Layout Layout, UpLo UPLO,
    int N, float alpha, float* X,
    int incX, float* Y, int incY, float* A)
let dspr2 layout uplo N alpha (X:float array) incX iniX (Y:float array) incY iniY (A:float array) =
    use ap = fixed &A.[0]
    use xp = fixed &X.[iniX]
    use yp = fixed &Y.[iniY]
    cblas_dspr2(layout, uplo, N, alpha, xp, incX, yp, incY, ap)

[<DllImport(DLL, CallingConvention = CallingConvention.Cdecl, ExactSpelling = true)>]
extern void cblas_dsymv(Layout Layout, UpLo UpLo,
    int N, float alpha, float* A,
    int lda, float* X, int incX,
    float beta, float* Y, int incY)
let dsymv layout uplo N alpha (A:float array) lda (X:float array) incX iniX beta (Y:float array) incY iniY = 
    use ap = fixed &A.[0]
    use xp = fixed &X.[iniX]
    use yp = fixed &Y.[iniY]
    cblas_dsymv(layout, uplo, N, alpha, ap, lda, xp, incX, beta, yp, incY)

[<DllImport(DLL, CallingConvention = CallingConvention.Cdecl, ExactSpelling = true)>]
extern void cblas_dsyr(Layout Layout, UpLo UpLo,
        int N, float alpha, float* X,
        int incX, float* A, int lda)
let dsyr layout uplo N alpha (X:float array) incX iniX (A:float array) lda =
    use ap = fixed &A.[0]
    use xp = fixed &X.[iniX]
    cblas_dsyr(layout, uplo, N, alpha, xp, incX, ap, lda)

[<DllImport(DLL, CallingConvention = CallingConvention.Cdecl, ExactSpelling = true)>]
extern void cblas_dsyr2(Layout Layout, UpLo UpLo,
    int N, float alpha, float* X,
    int incX, float* Y, int incY, float* A,
    int lda)
let dsyr2 layout uplo N alpha (X:float array) incX iniX (Y:float array) incY iniY (A:float array) lda =
    use ap = fixed &A.[0]
    use xp = fixed &X.[iniX]
    use yp = fixed &Y.[iniY]
    cblas_dsyr2(layout, uplo, N, alpha, xp, incX, yp, incY, ap, lda)

[<DllImport(DLL, CallingConvention = CallingConvention.Cdecl, ExactSpelling = true)>]
extern void cblas_dtbmv(Layout Layout, UpLo UpLo,
    Trans TransA, Diag Diag,
    int N, int K, float* A, int lda,
    float* X, int incX)
let dtbmv layout uplo transA diag N K (A:float array) lda (X:float array) incX iniX = 
    use ap = fixed &A.[0]
    use xp = fixed &X.[iniX]
    cblas_dtbmv(layout, uplo, transA, diag, N, K, ap, lda, xp, incX)

[<DllImport(DLL, CallingConvention = CallingConvention.Cdecl, ExactSpelling = true)>]
extern void cblas_dtbsv(Layout Layout, UpLo UpLo,
    Trans TransA, Diag Diag,
    int N, int K, float* A, int lda,
    float* X, int incX)
let dtbsv layout uplo transA diag N K (A:float array) lda (X:float array) incX iniX = 
    use ap = fixed &A.[0]
    use xp = fixed &X.[iniX]
    cblas_dtbsv(layout, uplo, transA, diag, N, K, ap, lda, xp, incX)

[<DllImport(DLL, CallingConvention = CallingConvention.Cdecl, ExactSpelling = true)>]
extern void cblas_dtpmv(Layout Layout, UpLo UpLo,
    Trans TransA, Diag Diag,
    int N, float* Ap, float* X, int incX)
let dtpmv layout uplo transA diag N (Ap:float array) (X:float array) incX iniX =
    use ap = fixed &Ap.[0]
    use xp = fixed &X.[iniX]
    cblas_dtpmv(layout, uplo, transA, diag, N, ap, xp, incX)

[<DllImport(DLL, CallingConvention = CallingConvention.Cdecl, ExactSpelling = true)>]
extern void cblas_dtpsv(Layout Layout, UpLo UpLo,
    Trans TransA, Diag Diag,
    int N, float* Ap, float* X, int incX)
let dtpsv layout uplo transA diag N (Ap:float array) (X:float array) incX iniX =
    use ap = fixed &Ap.[0]
    use xp = fixed &X.[iniX]
    cblas_dtpsv(layout, uplo, transA, diag, N, ap, xp, incX)

[<DllImport(DLL, CallingConvention = CallingConvention.Cdecl, ExactSpelling = true)>]
extern void cblas_dtrmv(Layout Layout, UpLo UpLo,
    Trans TransA, Diag Diag,
    int N, float* A, int lda,
    float* X, int incX)
let dtrmv layout uplo transA diag N (A:float array) lda (X:float array) incX iniX =
    use ap = fixed &A.[0]
    use xp = fixed &X.[iniX]
    cblas_dtrmv(layout, uplo, transA, diag, N, ap, lda, xp, incX)

[<DllImport(DLL, CallingConvention = CallingConvention.Cdecl, ExactSpelling = true)>]
extern void cblas_dtrsv(Layout Layout, UpLo UpLo, Trans Trans, Diag Diag, int N, float* A, int lda, float* X, int incX)
let dtrsv layout uplo trans diag N (A:float array) lda (X:float array) incX iniX = 
    use ap = fixed &A.[0]
    use xp = fixed &A.[iniX]
    cblas_dtrsv(layout, uplo, trans, diag, N, ap, lda, xp, incX)

// LEVEL 3 DOUBLE

[<DllImport(DLL, CallingConvention = CallingConvention.Cdecl, ExactSpelling = true)>]
extern void cblas_dgemm(Layout Layout, Trans TransA, Trans TransB, int M, int N, int K, float alpha, float* A, int lda, float* B, int ldb, float beta, float* C, int ldc)
let dgemm layout transA transB M N K alpha (A:float array) lda (B:float array) ldb beta (C:float array) ldc =
    use ap = fixed &A.[0]
    use bp = fixed &B.[0]
    use cp = fixed &C.[0]
    cblas_dgemm(layout, transA, transB, M, N, K, alpha, ap, lda, bp, ldb, beta, cp, ldc)

[<DllImport(DLL, CallingConvention = CallingConvention.Cdecl, ExactSpelling = true)>]
extern void cblas_dsymm(Layout Layout, Side Side,
    UpLo UpLo, int M, int N,
    float alpha, float* A, int lda,
    float* B, int ldb, float beta,
    float* C, int ldc)
let dsymm layout side uplo M N alpha (A:float array) lda (B:float array) ldb beta (C:float array) ldc =
    use ap = fixed &A.[0]
    use bp = fixed &B.[0]
    use cp = fixed &C.[0]
    cblas_dsymm(layout, side, uplo, M, N, alpha, ap, lda, bp, ldb, beta, cp, ldc)

[<DllImport(DLL, CallingConvention = CallingConvention.Cdecl, ExactSpelling = true)>]
extern void cblas_dsyrk(Layout Layout, UpLo UpLo,
    Trans Trans, int N, int K,
    float alpha, float* A, int lda,
    float beta, float* C, int ldc)
let dsyrk layout uplo trans N K alpha (A:float array) lda beta (C:float array) ldc =
    use ap = fixed &A.[0]
    use cp = fixed &C.[0]
    cblas_dsyrk(layout, uplo, trans, N, K, alpha, ap, lda, beta, cp, ldc)

[<DllImport(DLL, CallingConvention = CallingConvention.Cdecl, ExactSpelling = true)>]
extern void cblas_dsyr2k(Layout Layout, UpLo UpLo,
    Trans Trans, int N, int K,
    float alpha, float* A, int lda,
    float* B, int ldb, float beta,
    float* C, int ldc)
let dsyr2k layout uplo trans N K alpha (A:float array) lda (B:float array) ldb beta (C:float array) ldc =
    use ap = fixed &A.[0]
    use bp = fixed &B.[0]
    use cp = fixed &C.[0]
    cblas_dsyr2k(layout, uplo, trans, N, K, alpha, ap, lda, bp, ldb, beta, cp, ldc)

[<DllImport(DLL, CallingConvention = CallingConvention.Cdecl, ExactSpelling = true)>]
extern void cblas_dtrmm(Layout Layout, Side Side,
    UpLo UpLo, Trans TransA,
    Diag Diag, int M, int N,
    float alpha, float* A, int lda,
    float* B, int ldb)
let dtrmm layout side uplo transA diag M N alpha (A:float array) lda (B:float array) ldb =
    use ap = fixed &A.[0]
    use bp = fixed &B.[0]
    cblas_dtrmm(layout, side, uplo, transA, diag, M, N, alpha, ap, lda, bp, ldb)

[<DllImport(DLL, CallingConvention = CallingConvention.Cdecl, ExactSpelling = true)>]
extern void cblas_dtrsm(Layout Layout, Side Side,
    UpLo UpLo, Trans TransA,
    Diag Diag, int M, int N,
    float alpha, float* A, int lda,
    float* B, int ldb)
let dtrsm layout side uplo transA diag M N alpha (A:float array) lda (B:float array) ldb =
    use ap = fixed &A.[0]
    use bp = fixed &B.[0]
    cblas_dtrsm(layout, side, uplo, transA, diag, M, N, alpha, ap, lda, bp, ldb)