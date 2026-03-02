#load "../components.fsx"

// -----------------------------------------------------------------------------
// Supporting script: approximate shrub allometry with a shifted power-law
//
// Goal:
//   Replace the expensive mechanistic biomass–radius allometry with a fast,
//   closed-form surrogate:
//
//       R(B) = (B / a)^(1/p) + c
//
//   where R is stem radius (mm) and B is biomass (g).
//
// Usage:
//   1. Run this script.
//   2. Copy the fitted (a, p, c) and the Bristlecone expression back into
//      components.fsx.
// -----------------------------------------------------------------------------

open Bristlecone
open Bristlecone.Language

// -----------------------------------------------------------------------------
// 1. Build a radius–biomass grid from the mechanistic allometry
// -----------------------------------------------------------------------------

// Dummy radius parameter (only used to drive the allometry).
let radiusParam =
    parameter
        "radius"
        NoConstraints
        0.01<Dendro.Units.millimetre>
        0.24<Dendro.Units.millimetre>

let fakePool =
    Parameter.Pool.fromList
        [ radiusParam.ParamId.Inner,
          Parameter.Pool.boxParam radiusParam.ParamId.Inner.Value radiusParam.Parameter ]

// Compile the mechanistic biomass allometry: radius -> biomass (g)
let bioFn =
    ExpressionCompiler.compileMeasure
        fakePool
        (ShrubModel.Allometry.Proxies.toBiomassMM (P radiusParam))

// Dummy state tensor (not used by the allometry itself).
let dummyTensor = Tensors.Typed.ofScalar 2.1<ModelSystem.state>

// Build a fine grid of radii (1–100 mm) and evaluate biomass at each.
let grid =
    [ 1.53<Dendro.Units.millimetre> .. 0.01<Dendro.Units.millimetre> .. 100.<Dendro.Units.millimetre> ]
    |> List.map (fun r ->
        let pool =
            Tensors.Typed.ofVector
                [| r * 1.<parameter / Dendro.Units.millimetre> |]

        let biomass =
            bioFn pool Map.empty dummyTensor 1
            |> Tensors.Typed.toFloatScalar

        r, biomass)

// -----------------------------------------------------------------------------
// 2. Fit a shifted power-law surrogate: R(B) = (B / a)^(1/p) + c
// -----------------------------------------------------------------------------

#r "nuget: MathNet.Numerics.FSharp"

open MathNet.Numerics
open MathNet.Numerics.Optimization
open MathNet.Numerics.LinearAlgebra

/// Shifted power-law surrogate: radius (mm) as a function of biomass (g).
let model (a: float) (p: float) (c: float) (biomass: float) =
    (biomass / a) ** (1.0 / p) + c

/// Sum of squared errors in radius space.
let objective (pars: Vector<float>) =
    let a = pars.[0]
    let p = pars.[1]
    let c = pars.[2]

    grid
    |> List.sumBy (fun (rObs, bObs) ->
        let rPred = model a p c (float bObs)
        let err = rPred - float rObs
        err * err)

/// Initial guess for (a, p, c)
let initial = vector [| 0.5; 3.0; 0.1 |]

/// Fit using Nelder–Mead simplex.
let solver = NelderMeadSimplex(1e-8, 2000)
let result = solver.FindMinimum(ObjectiveFunction.Value objective, initial)

let a_fit = result.MinimizingPoint.[0]
let p_fit = result.MinimizingPoint.[1]
let c_fit = result.MinimizingPoint.[2]

printfn "Fitted shifted power-law parameters:"
printfn "  a = %.9f" a_fit
printfn "  p = %.9f" p_fit
printfn "  c = %.9f" c_fit

// -----------------------------------------------------------------------------
// 3. Diagnostics: residuals, RMSE, max error
// -----------------------------------------------------------------------------

let predictRadius (biomass: float) =
    model a_fit p_fit c_fit biomass

/// (radius_obs_mm, residual_mm)
let residuals =
    grid
    |> List.map (fun (rObs, bObs) ->
        let rPred = predictRadius (float bObs)
        let err = rPred - float rObs
        rObs, err)

let rmse =
    residuals
    |> List.averageBy (fun (_, err) -> err * err)
    |> sqrt

let maxAbsRadius, maxAbsErr =
    residuals
    |> List.maxBy (fun (_, err) -> abs err)
    |> fun (r, e) -> r, e

printfn "Residual diagnostics:"
printfn "  RMSE       = %.4f mm" rmse
printfn "  Max |err|  = %.4f mm at radius %.1f mm" (abs maxAbsErr) (float maxAbsRadius)
