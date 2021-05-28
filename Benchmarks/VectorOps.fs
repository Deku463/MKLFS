module VectorOps

open System
open BenchmarkDotNet.Attributes
open LinearAlgebra.Double

type VectorOpsBench () = 
    let rng = Random()
    let v1 = Vector(Array.init 1000 (fun x -> rng.NextDouble()))
    let v2 = Vector(Array.init 1000 (fun x -> rng.NextDouble()))

    [<Benchmark>]
    member this.Addition () = v1 + v2

    [<Benchmark>]
    member this.Dot () = dot v1 v2
    
    [<Benchmark>]
    member this.Scale () = 2.0 * v1

    [<Benchmark>]
    member this.Sum () = Vector.sum v1
    

