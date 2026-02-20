module ShrubModel

#load "units.fsx"
#load "constants.fsx"

open Bristlecone.Language
open Bristlecone.Dendro.Units
open Units

[<Measure>] type stems

let private c = Constant

module Allometry =

    let pi = c System.Math.PI

    let maxExpr (a: ModelExpression<'u>) (b: ModelExpression<'u>) =
        Conditional (a .> b) a b

    module NiklasAndSpatz =

        let basalRadius (k5: ModelExpression<1>) (k6: ModelExpression<m>)
                        (stemLength: ModelExpression<cm>) : ModelExpression<cm> =
            let Lm = toMetres stemLength
            let term = maxExpr ((Lm + k6) / k5) (c 0.0<m>)
            let r = term ** c 1.5
            let rCm = toCentimetres r
            let rRadius = rCm / c 2.0
            maxExpr rRadius (c 1e-6<Units.cm>)

        let stemLength (k5: ModelExpression<1>) (k6: ModelExpression<m>)
                    (radius: ModelExpression<cm>) : ModelExpression<cm> =
            let diameter = (radius * c 2.0) |> toMetres
            // Dimensionless factor (dM / 1 m)^(2/3)
            let factor = (diameter / c 1.0<m>) ** c (2.0/3.0)
            let h = (k5 * factor * c 1.<m> - k6) |> toCentimetres
            maxExpr h (c 1e-6<cm>)

    module Gotmark2016 =

        /// Volume of a truncated cone (cm³)
        let coneVolume (h: ModelExpression<cm>) (r1: ModelExpression<cm>) (r2: ModelExpression<cm>) =
            pi * h * (r1 * r1 + r1 * r2 + r2 * r2) / c 3.0

        /// Volume of a cylinder (cm³)
        let cylinderVolume (h: ModelExpression<cm>) (r: ModelExpression<cm>) =
            pi * h * (r * r)

        /// Main stem volume
        let mainStemVolume n rtip (radius: ModelExpression<cm>) h =
            Conditional (radius .> rtip)
                (n * coneVolume h radius rtip)
                (n * cylinderVolume h rtip)

        let private safeLen x = maxExpr x (c 1e-6<cm>)

        /// b = dimensionless branching coefficient, not biomass.
        let branchVolume b a rtip p lmin n (radius: ModelExpression<cm> -> ModelExpression<cm>)
                         h (kf: ModelExpression<1>) =
            let pk = (p ** kf) * h
            let factor = n * a * ((a + c 1.0) ** kf)

            // Long-branch regime: pk >= lmin
            let longTerm =
                let pkLong = (p ** (kf + c 1.0)) * h
                Conditional (radius pkLong .> rtip)
                    (factor * coneVolume pkLong (radius pkLong) rtip)
                    (factor * cylinderVolume pkLong rtip)

            // Short-branch regime: pk < lmin
            let shortTerm =
                let pkShort = safeLen <| c 3.0 * p * (pk - c (2.0/3.0) * lmin)
                Conditional (b * pkShort .> rtip)
                    (factor * coneVolume pkShort (radius pkShort) rtip)
                    (factor * cylinderVolume pkShort rtip)

            // Regime selection without recursion
            let regime =
                Conditional (pk .< lmin) shortTerm longTerm

            // While-condition gate: only active generations contribute
            Conditional (pk .> c (2.0/3.0) * lmin) regime (c 0.0<cm^3>)

        /// Because Bristlecone's typed expression trees don't allow recursion,
        /// we can't just stop when we reach the end of branching. We have to define
        /// a maximum amount of realistic branching and work to that.
        let maxBranchingGenerations = 7 // Safe below around 7m stem length

        /// Shrub volume with expression-based parameters
        let shrubVolume b a rtip p lmin k5 k6 n h =
            let radius = NiklasAndSpatz.basalRadius k5 k6
            let main = mainStemVolume n rtip (radius h) h
            let kIndices : ModelExpression<1> list =
                [0 .. maxBranchingGenerations] |> List.map (fun i -> c (float i))
            let branches = kIndices |> List.map (branchVolume b a rtip p lmin n radius h)
            let branchTotal = branches |> List.fold (+) (c 0.0<cm^3>)
            Conditional (List.last branches .> c 1e-6<cm^3>) Invalid (main + branchTotal)


    module Allometrics =

        let mass woodDensity volume = volume * woodDensity

        let massToVolume woodDensity mass = mass / woodDensity

        let shrubBiomass b a rtip p lmin k5 k6 n woodDensity (radius: ModelExpression<cm>) =
            radius
            |> NiklasAndSpatz.stemLength k5 k6
            |> Gotmark2016.shrubVolume b a rtip p lmin k5 k6 n
            |> mass woodDensity

        let shrubRadius b a rtip p lmin k5 k6 n woodDensity mass : ModelExpression<cm> =
            let findRadius volume =
                let v x =
                    x
                    |> NiklasAndSpatz.stemLength k5 k6
                    |> Gotmark2016.shrubVolume b a rtip p lmin k5 k6 n
                Inverse (fun x -> v x) volume (c 0.01<cm>) (c 20.0<cm>)

            mass |> massToVolume woodDensity |> findRadius
 
        let shrubHeight k5 k6 radius =
            radius |> NiklasAndSpatz.stemLength k5 k6


    module Proxies =

        /// Radius in millimetres
        let toBiomassMM (radius: ModelExpression<millimetre>) =
            let radiusCm = radius |> mmToCm
            Allometrics.shrubBiomass
                Constants.Allometrics.b
                Constants.Allometrics.a
                Constants.Allometrics.rtip
                Constants.Allometrics.p
                Constants.Allometrics.lmin
                Constants.Allometrics.k5
                Constants.Allometrics.k6
                Constants.Allometrics.numberOfStems
                Constants.Allometrics.salixWoodDensity
                radiusCm

        /// Biomass in grams.
        let toRadiusMM (biomassGrams: ModelExpression<g>) : ModelExpression<mm> =
            biomassGrams
            |> Allometrics.shrubRadius
                Constants.Allometrics.b
                Constants.Allometrics.a
                Constants.Allometrics.rtip
                Constants.Allometrics.p
                Constants.Allometrics.lmin
                Constants.Allometrics.k5
                Constants.Allometrics.k6
                Constants.Allometrics.numberOfStems
                Constants.Allometrics.salixWoodDensity
            |> Units.cmToMm


module GrowthLimitation =

    // N availability range of our data is approx. -0.10 to 1.72
    let nMin, nMax = Constant -0.10, Constant 1.72
    
    /// Returns a penalty if the function is flat over the realistic range.
    let tooFlatPenalty uptakeMax uptakeMin =
        Conditional (uptakeMax - uptakeMin .< Constant 1e-6</g/year>)
            Invalid
            (Constant 0.</g/year>)

    /// Returns a penalty if the function is saturated over the realistic range.
    let tooSaturatedPenalty uptakeMax uptakeMin =
        Conditional ((uptakeMax / uptakeMin) .> Constant 50.0)
            Invalid
            (Constant 0.</g/year>)    

    /// <summary>A function that represents the efficiency of two simultaneous
    /// processes. The combined processes have a single 'handling time',
    /// which represents the rate at which the processes occur.</summary>
    /// <param name="a">efficiency of process A</param>
    /// <param name="b">efficiency of process B</param>
    /// <param name="h">integrated handling time / rate of process A+B (per item)</param>
    let hollingDiscModelDual (a:ModelExpression<1/g/year>) (b:ModelExpression<g/1>) (h: ModelExpression<year>) : ModelExpression<1> -> ModelExpression</g/year> =
        let model r : ModelExpression</g/year> = (a * r) / (Constant 1. + (a * b * h * r))
        let uptakeMax, uptakeMin = model nMax, model nMin
        fun resource ->
            model resource + tooSaturatedPenalty uptakeMax uptakeMin + tooFlatPenalty uptakeMax uptakeMin

    /// Monod model once saturation has been reached
    let linear (a:ModelExpression</g/year>) =
        fun (resource: ModelExpression<1>) ->
            Conditional (a * nMax .< Constant 1e-12</g/year>) Invalid (a * resource)

    /// <summary>A monotonically increasing function of a resource `r`.</summary>
    ///  <param name="h">soil resource concentration required for growth at half the maximum rate</param>
    ///  <param name="r">the current resource concentration</param>**Parameters**
    let michaelisMenten h r : ModelExpression<1> =
        r / (h + r)

    /// From Jabot and Pottier 2012
    let monod k r : ModelExpression<1> =
        r / (k + r)


module AbioticResource =

    /// **Description**
    /// A standard chemostat-type model for replenshment of an abiotic resource.
    /// **Parameters**
    ///   * `d` - a rate constant
    ///   * `s` - resource concentration of the inflow
    ///   * `n` - the current resource concentration
    let chemostat d s (n:float) =
        d * (s - n)
