module ShrubModel

#r "/Users/andrewmartin/Documents/GitHub Projects/bristlecone/src/Bristlecone/bin/Debug/net5.0/Bristlecone.dll"
#r "/Users/andrewmartin/Documents/GitHub Projects/bristlecone/src/Bristlecone.Dendro/bin/Debug/net5.0/Bristlecone.Dendro.dll"
#load "units.fsx"
#load "constants.fsx"

open Bristlecone
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
            let term = (Lm + k6) / k5
            let r = c 50.0 * (term ** c 1.5) |> toCentimetres
            maxExpr r (c 1e-6<cm>)

        let stemLength (k5: ModelExpression<1>) (k6: ModelExpression<m>)
                    (radius: ModelExpression<cm>) : ModelExpression<cm> =
            let rm = toMetres radius
            let term = k5 * ((rm / c 50.0) ** c (2.0/3.0))
            let Lm = term - k6
            maxExpr (toCentimetres Lm) (c 1e-6<cm>)


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
                let pkShort = c 3.0 * p * (pk - c (2.0/3.0) * lmin)
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
        let maxBranchingGenerations = 20

        /// Shrub volume with expression-based parameters
        let shrubVolume b a rtip p lmin k5 k6 n h =
            let radius = NiklasAndSpatz.basalRadius k5 k6
            let main = mainStemVolume n rtip (radius h) h
            let kIndices : ModelExpression<1> list =
                [0 .. maxBranchingGenerations] |> List.map (fun i -> c (float i))
            let branches =
                kIndices
                |> List.map (branchVolume b a rtip p lmin n radius h)
                |> List.fold (+) (c 0.0<cm^3>)
            main + branches

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
                Inverse (fun x -> v x) volume (c 0.01<cm>) (c 100.0<cm>)

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
                (Constants.Allometrics.k6 |> toMetres)
                Constants.Allometrics.numberOfStems
                Constants.Allometrics.salixWoodDensity
                radiusCm

        /// Biomass in grams.
        let toRadiusMM (biomassGrams: ModelExpression<g>) : ModelExpression<mm> =
            Conditional (Bool.isFinite biomassGrams)
                (
                    biomassGrams
                    |> Allometrics.shrubRadius
                        Constants.Allometrics.b
                        Constants.Allometrics.a
                        Constants.Allometrics.rtip
                        Constants.Allometrics.p
                        Constants.Allometrics.lmin
                        Constants.Allometrics.k5
                        (toMetres Constants.Allometrics.k6)
                        Constants.Allometrics.numberOfStems
                        Constants.Allometrics.salixWoodDensity
                    |> Units.cmToMm
                ) Invalid


module GrowthLimitation =

    /// <summary>A function that represents the efficiency of two simultaneous
    /// processes. The combined processes have a single 'handling time',
    /// which represents the rate at which the processes occur.</summary>
    /// <param name="a">efficiency of process A</param>
    /// <param name="b">efficiency of process B</param>
    /// <param name="h">integrated handling time / rate of process A+B (per item)</param>
    /// <param name="min">a level of resource at which the resultant process is `> 1e-12`</param>
    let hollingDiscModelDual (a:ModelExpression<1/g/year>) (b:ModelExpression<g/1>) (h: ModelExpression<year>) min : ModelExpression<1> -> ModelExpression</g/year> =
        let model r : ModelExpression</g/year> = (a * r) / (Constant 1. + (a * b * h * r))
        fun resource ->
            Conditional (model min .< Constant 1e-12</g/year>) Invalid (model resource)

    /// Monod model once saturation has been reached
    let linear (a:ModelExpression</g/year>) min =
        fun (resource: ModelExpression<1>) ->
            Conditional (a * min .< Constant 1e-12</g/year>) Invalid (a * resource)

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
