module rec LinearAlgebra.Double

open EnumTypes
open Blas

exception DimensionMismatch of string

type Vector(data:float array) =
    member this.data = data
    member this.length = this.data.Length
    
    static member sum (V:Vector) = dasum V.length V.data 1

    static member copy (V:Vector) = 
        let ans = Array.zeroCreate V.length
        dcopy V.length V.data 1 0 ans 1 0 
        Vector(ans)

    static member (+) (Va:Vector, Vb:Vector) =
        if Va.length = Vb.length then
            let ans = Vector.copy Vb
            axpy' 1.0 Va ans
            ans
        else raise (DimensionMismatch $"The length of Va ({Va.length}) is not the same as the length of Vb ({Vb.length}).")

    static member (-) (Va:Vector, Vb:Vector) =
        if Va.length = Vb.length then
            let ans = Vector.copy Va
            axpy' -1.0 Vb ans
            ans
        else raise (DimensionMismatch $"The length of Va ({Va.length}) is not the same as the length of Vb ({Vb.length}).")

    static member (*) (scalar:float, V:Vector) = 
        let ans = Vector.copy V
        scal' scalar ans
        ans
        
type Matrix(data:float array, rows:int, cols:int, ?layout:Layout) =
    member this.layout = defaultArg layout Layout.ColMajor
    member this.data = data
    member this.rows = rows
    member this.cols = cols

    member this.MVOp ((alpha:float, Va:Vector, beta:float, Vb:Vector)) =
        if this.rows = Va.length && Va.length = Vb.length then
            let ld = // Get leading dimension of M depending on layout. 
                if this.layout = Layout.ColMajor 
                then this.rows else this.cols
            dgemv this.layout Trans.NoTrans this.rows this.cols alpha this.data ld Va.data 1 0 beta Vb.data 1 0
        else raise (DimensionMismatch $"The numbers of cols of M ({this.cols}), the length of Va ({Va.length}) and the length of Vb ({Vb.length}) must be the same.")



type Transposed(data:float array, length1:int, length2:int, ?layout:Layout) =
    inherit Matrix(data, length1, length2)
    member this.layout = defaultArg layout Layout.ColMajor

    member this.MVOp ((alpha:float, Va:Vector, beta:float, Vb:Vector)) =
        if this.rows = Va.length && Va.length = Vb.length then 
            let ld = // Get leading dimension of M depending on layout. 
                if this.layout = Layout.ColMajor 
                then this.rows else this.cols
            dgemv this.layout Trans.Trans this.rows this.cols alpha this.data ld Va.data 1 0 beta Vb.data 1 0
        else raise (DimensionMismatch $"The numbers of cols of M ({this.cols}), the length of Va ({Va.length}) and the length of Vb ({Vb.length}) must be the same.")

type Symmetric(data:float array, length1:int, length2:int, uplo:UpLo, ?layout:Layout) =
    inherit Matrix(data, length1, length2)
    member this.layout = defaultArg layout Layout.ColMajor
    member this.uplo = uplo
    
    member this.MVOp ((alpha:float, Va:Vector, beta:float, Vb:Vector, ld:int)) =
        if this.rows = Va.length && Va.length = Vb.length then 
            let ld = // Get leading dimension of M depending on layout. 
                if this.layout = Layout.ColMajor 
                then this.rows else this.cols
            dsymv this.layout this.uplo this.rows alpha this.data ld Va.data 1 0 beta Vb.data 1 0
        else raise (DimensionMismatch $"The numbers of cols of M ({this.cols}), the length of Va ({Va.length}) and the length of Vb ({Vb.length}) must be the same.")

// Vb <- (alpha * Va) + Vb
let axpy' (alpha:float) (Va:Vector) (Vb:Vector) = daxpy Va.length alpha Va.data 1 0 Vb.data 1 0

let dot (Va:Vector) (Vb:Vector) = 
    if Va.length = Vb.length then ddot Va.length Va.data 0 1 Vb.data 0 1
    else raise (DimensionMismatch $"The length of Va ({Va.length}) is not the same as the length of Vb ({Vb.length}).")

// V <- alpha * V
let scal' (alpha:float) (V:Vector) = dscal V.length alpha V.data 1 0

// Vb <- alpha * (M * Va) + (beta * Vb)
let mv' (alpha:float) (M:Matrix) (Va:Vector) (beta:float) (Vb:Vector) = M.MVOp(alpha, Va, beta, Vb)
    //if M.rows = Va.length && Va.length = Vb.length then 
    //    let ld = // Get leading dimension of M depending on layout. 
    //        if M.layout = Layout.ColMajor 
    //        then M.rows else M.cols
    //    M.MVOp(alpha, M, Va, beta, Vb)
    //else raise (DimensionMismatch $"The numbers of cols of M ({M.cols}), the length of Va ({Va.length}) and the length of Vb ({Vb.length}) must be the same.")

// Vb <- alpha * (M * Va) + (beta * Vb)
type private MV = MV with
    static member op (alpha:float, M:Matrix, Va:Vector, beta:float, Vb:Vector, ld:int) = // Case Matrix
        dgemv M.layout Trans.NoTrans M.rows M.cols alpha M.data ld Va.data 1 0 beta Vb.data 1 0
    static member op (alpha:float, M:Transposed, Va:Vector, beta:float, Vb:Vector, ld:int) = // Case Transposed
        dgemv M.layout Trans.Trans M.rows M.cols alpha M.data ld Va.data 1 0 beta Vb.data 1 0
    static member op (alpha:float, M:Symmetric, Va:Vector, beta:float, Vb:Vector, ld:int) = // Case Symmetric
        dsymv M.layout M.uplo M.rows alpha M.data ld Va.data 1 0 beta Vb.data 1 0

let mv2' (alpha:float) (M:Matrix) (Va:Vector) (beta:float) (Vb:Vector) = 
    if M.rows = Va.length && Va.length = Vb.length then 
        let ld = // Get leading dimension of M depending on layout. 
            if M.layout = Layout.ColMajor 
            then M.rows else M.cols
        MV.op(alpha, M, Va, beta, Vb, ld)
    else raise (DimensionMismatch $"The numbers of cols of M ({M.cols}), the length of Va ({Va.length}) and the length of Vb ({Vb.length}) must be the same.")

