module ModelComponents

#r "nuget: Bristlecone.Dendro,2.0.0"
#load "constants.fsx"

open Bristlecone
open Bristlecone.Language

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

        let shrubBiomass b a rtip p lmin k5 k6 n woodDensity (radius: float<mm>) =
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
    /// <param name="h">integrated handling time / rate of process A+B</param>
    /// <param name="min">a level of resource at which the resultant process is `> 1e-12`</param>
    let hollingDiscModelDual a b h min =
        let model r = (a * r) / (Constant 1. + (a * b * h * r))
        fun resource ->
            Conditional
            <| fun compute ->
                if compute (model min) < 1e-12 then
                    Invalid
                else
                    model resource

    /// Monod model once saturation has been reached
    let linear a min =
        fun resource ->
            Conditional
            <| fun compute ->
                if compute (a * min) < 1e-12 then
                    Invalid
                else
                    a * resource

    /// <summary>A monotonically increasing function of a resource `r`.</summary>
    ///  <param name="h">soil resource concentration required for growth at half the maximum rate</param>
    ///  <param name="r">the current resource concentration</param>**Parameters**
    let michaelisMenten h (r:ModelExpression) =
        r / (h + r)

    /// From Jabot and Pottier 2012
    let monod k (r:ModelExpression) =
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
    let toBiomassMM (radius: float<mm>) =
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
    let toRadiusMM biomassGrams =
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

    /// d15N to N availability. From Craine 2009, as shown in Craine 2015 (Plant and Soil).
    /// Assuming d15N is a linear index of N availability, the minimum supported value of d15N is -3.09, as 0 N availability.
    let d15NtoAvailability d15N =
        (100. * d15N + 309.) / 359.

    let shrubHeightCm radiusMM =
        radiusMM / 10. |> Allometric.Allometrics.shrubHeight Constants.Allometrics.k5 Constants.Allometrics.k6


module GeometricConstraint = 

    /// Linear growth rate in dM/dt form.
    let none _ = Constant 1.

    /// A dM/dt form of the von Bertalanffy monomollecular growth function, where M = mass.
    let vonBertalanffy k m : ModelExpression = 
        (k / m)

    /// The dM/dt form of the Chapman-Richards growth function, where M = mass.
    let chapmanRichards k m =
        (Constant 1. - (m / k))


module FeedbackToSoil =

    let none b = Constant 0. * b
    let withBiomassLoss alpha gammab b : float = alpha * b * gammab
