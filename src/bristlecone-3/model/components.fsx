module ShrubModel

#r "/Users/andrewmartin/Documents/GitHub Projects/bristlecone/src/Bristlecone/bin/Debug/net5.0/Bristlecone.dll"
#r "/Users/andrewmartin/Documents/GitHub Projects/bristlecone/src/Bristlecone.Dendro/bin/Debug/net5.0/Bristlecone.Dendro.dll"
#load "constants.fsx"
#load "units.fsx"

open Bristlecone
open Bristlecone.Language
open Bristlecone.Dendro.Units
open Units

[<Measure>] type stems

module Allometry =

    let c = Constant
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

        /// Shrub volume with expression-based parameters
        let shrubVolume b a rtip p lmin k5 k6 n h =
            let radius = NiklasAndSpatz.basalRadius k5 k6
            let main = mainStemVolume n rtip (radius h) h
            let kIndices : ModelExpression<1> list =
                [0..15] |> List.map (fun i -> c (float i))
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

        let shrubRadius b a rtip p lmin k5 k6 n woodDensity mass =
            let findRadius volume =
                let v x =
                    x
                    |> NiklasAndSpatz.stemLength k5 k6
                    |> Gotmark2016.shrubVolume b a rtip p lmin k5 k6 n

                let f = (fun x -> (v x) - volume)
                Statistics.RootFinding.bisect 0 200 f 0.01 100.00 1e-8 // Assumption that shrub radius is between 0.01 and 100.0cm.

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
            if
                System.Double.IsNaN biomassGrams
                || System.Double.IsInfinity biomassGrams
                || System.Double.IsNegativeInfinity biomassGrams
            then
                Invalid
            else
                let radiusCm =
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

                radiusCm * c 10.




module Allometric =

    let pi = System.Math.PI

    module NiklasAndSpatz_Allometry =

        let nthroot n A =
            let rec f x =
                let m = n - 1.
                let x' = (m * x + A / x ** m) / n

                match abs (x' - x) with
                | t when t < abs (x * 1e-9) -> x'
                | _ -> f x'

            f (A / double n)

        /// Gives the basal radius in centimetres of a stem/branch given its length in centimetres. Function from Niklas and Spatz (2004).
        /// The basal radius is always positive.
        let basalRadius k5 k6 stemLength =
            max (100. * ((0.01 * stemLength + k6) / k5) ** (3. / 2.) / 2.) 1e-06

        /// Inverse equation of basalRadius.
        let stemLength k5 k6 radius =
            max (2. * ((nthroot 3. 2.) * 5. ** (2. / 3.) * k5 * radius ** (2. / 3.) - 50. * k6)) 1e-06

    module Götmark2016_ShrubModel =

        /// Total shrub volume given height and number of stems
        let shrubVolume b a rtip p lmin k5 k6 n h =

            let radius = NiklasAndSpatz_Allometry.basalRadius k5 k6

            let mainStemVolume =
                match radius h with
                | r when r > rtip -> n * pi * h * ((radius h) ** 2. + (radius h) * rtip + rtip ** 2.) / 3.
                | _ -> n * pi * h * rtip ** 2.

            let mutable volume = mainStemVolume
            let mutable k = 0.

            while (p ** k * h > lmin * 2. / 3.) do
                let volToAdd =
                    match (p ** k * h < lmin) with
                    | true ->
                        match (b * 3. * p * (p ** k * h - 2. * lmin / 3.) > rtip) with
                        | true ->
                            n
                            * a
                            * (a + 1.) ** (float k)
                            * pi
                            * 3.
                            * p
                            * (p ** k * h - 2. * lmin / 3.)
                            * ((radius (3. * p * (p ** k * h - 2. * lmin / 3.)))
                               * (radius (3. * p * (p ** k * h - 2. * lmin / 3.)) * rtip + rtip ** 2.))
                            / 3.
                        | false ->
                            n
                            * a
                            * (a + 1.) ** (float k)
                            * 3.
                            * p
                            * (p ** k * h - 2. * lmin / 3.)
                            * pi
                            * rtip ** 2.
                    | false ->
                        match (radius (p ** (k + 1.) * h) > rtip) with
                        | true ->
                            n
                            * a
                            * (a + 1.) ** (float k)
                            * pi
                            * p ** (k + 1.)
                            * h
                            * ((radius (p ** (k + 1.) * h)) ** 2.
                               + (radius (p ** (k + 1.) * h)) * rtip
                               + rtip ** 2.)
                            / 3.
                        | false -> n * a * (a + 1.) ** (float k) * p ** (k + 1.) * h * pi * rtip ** 2.

                volume <- volume + volToAdd
                k <- k + 1.

            k, volume


    module Allometrics =

        open Götmark2016_ShrubModel
        open Bristlecone.Dendro

        let private removeUnit (x: float<_>) = float x

        let mass woodDensity volume = volume * woodDensity

        let massToVolume woodDensity mass = mass / woodDensity

        let shrubBiomass b a rtip p lmin k5 k6 n woodDensity (radius: float<millimetre>) =
            radius
            |> removeUnit
            |> NiklasAndSpatz_Allometry.stemLength k5 k6
            |> shrubVolume b a rtip p lmin k5 k6 n
            |> snd
            |> mass woodDensity

        let shrubRadius b a rtip p lmin k5 k6 n woodDensity mass =
            let findRadius volume =
                let v x =
                    x
                    |> NiklasAndSpatz_Allometry.stemLength k5 k6
                    |> shrubVolume b a rtip p lmin k5 k6 n
                    |> snd

                let f = (fun x -> (v x) - volume)
                Statistics.RootFinding.bisect 0 200 f 0.01 100.00 1e-8 // Assumption that shrub radius is between 0.01 and 100.0cm.

            mass |> massToVolume woodDensity |> findRadius

        let shrubHeight k5 k6 radius =
            radius |> NiklasAndSpatz_Allometry.stemLength k5 k6


module GrowthLimitation =

    /// <summary>A function that represents the efficiency of two simultaneous
    /// processes. The combined processes have a single 'handling time',
    /// which represents the rate at which the processes occur.</summary>
    /// <param name="a">efficiency of process A</param>
    /// <param name="b">efficiency of process B</param>
    /// <param name="h">integrated handling time / rate of process A+B (per item)</param>
    /// <param name="min">a level of resource at which the resultant process is `> 1e-12`</param>
    let hollingDiscModelDual (a:ModelExpression<1/g/year>) (b:ModelExpression<g>) (h: ModelExpression</year>) min : ModelExpression<1> -> ModelExpression</g/year> =
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
    let michaelisMenten h r : ModelExpression<_> =
        r / (h + r)

    /// From Jabot and Pottier 2012
    let monod k r : ModelExpression<_> =
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

module Proxies =

    open Bristlecone.Dendro

    /// Radius in millimetres
    let toBiomassMM (radius: float<millimetre>) =
        radius / 10.
        |> Allometric.Allometrics.shrubBiomass
            Constants.Allometrics.b
            Constants.Allometrics.a
            Constants.Allometrics.rtip
            Constants.Allometrics.p
            Constants.Allometrics.lmin
            Constants.Allometrics.k5
            Constants.Allometrics.k6
            Constants.Allometrics.numberOfStems
            Constants.Allometrics.salixWoodDensity

    /// Biomass in grams.
    let toRadiusMM (biomassGrams: ModelExpression<g>) : ModelExpression<mm> =
        if
            System.Double.IsNaN biomassGrams
            || System.Double.IsInfinity biomassGrams
            || System.Double.IsNegativeInfinity biomassGrams
        then
            nan
        else
            let radiusCm =
                biomassGrams
                |> Allometric.Allometrics.shrubRadius
                    Constants.Allometrics.b
                    Constants.Allometrics.a
                    Constants.Allometrics.rtip
                    Constants.Allometrics.p
                    Constants.Allometrics.lmin
                    Constants.Allometrics.k5
                    Constants.Allometrics.k6
                    Constants.Allometrics.numberOfStems
                    Constants.Allometrics.salixWoodDensity

            radiusCm * 10.

    /// d15N to "Standardised N Availability". From Craine 2009, as shown in Craine 2015 (Plant and Soil).
    /// Assuming d15N is a linear index of N availability, the minimum supported value of d15N is -3.09, as 0 N availability.
    let d15NtoAvailability d15N =
        (100. * d15N + 309.) / 359.

    let shrubHeightCm radiusMM =
        radiusMM / 10. |> Allometric.Allometrics.shrubHeight Constants.Allometrics.k5 Constants.Allometrics.k6


module GeometricConstraint = 

    /// Linear growth rate in dM/dt form.
    let none _ = Constant 1.

    /// A dM/dt form of the von Bertalanffy monomollecular growth function, where M = mass.
    let vonBertalanffy k m : ModelExpression<_> = 
        (k / m)

    /// The dM/dt form of the Chapman-Richards growth function, where M = mass.
    let chapmanRichards k m =
        (Constant 1. - (m / k))


module FeedbackToSoil =

    let none b = Constant 0. * b
    let withBiomassLoss alpha gammab (b: ModelExpression<g>) : ModelExpression<_> = alpha * b * gammab
