module ModelComponents

#r "../../packages/Bristlecone/lib/netstandard2.0/bristlecone.dll"
#load "constants.fsx"

open Bristlecone

let pi = System.Math.PI

module NiklasAndSpatz_Allometry =

    let nthroot n A =
        let rec f x =
            let m = n - 1.
            let x' = (m * x + A/x**m) / n
            match abs(x' - x) with
            | t when t < abs(x * 1e-9) -> x'
            | _ -> f x'
        f (A / double n)

    /// Gives the basal radius in centimeters of a stem/branch given its length in centimeters. Function from Niklas and Spatz (2004). 
    let basalRadius k5 k6 stemLength =
        100. * (( 0.01 * stemLength + k6) / k5) ** (3. / 2.) / 2.

    /// Inverse equation of basalRadius, rearranged using Wolfram Alpha
    /// http://www.wolframalpha.com/input/?i=solve+r+%3D+100*((0.01*h%2Bk_6)%2Fk_5)%5E(3%2F2)%2F2+for+h
    let stemLength k5 k6 radius =
        2. * ((nthroot 3. 2.) * 5. ** (2./3.) * k5 * radius ** (2./3.) - 50. * k6)

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

        while (p ** k * h > lmin * 2./3.) do
            let volToAdd =
                match (p ** k * h < lmin) with
                | true ->
                    match (b * 3. * p * (p ** k * h - 2. * lmin / 3.) > rtip) with
                    | true ->
                        n * a * (a + 1.) ** (float k) * pi * 3. * p * (p ** k * h - 2. * lmin / 3.) * ((radius (3. * p * (p ** k * h - 2. * lmin / 3.))) * (radius (3. * p * (p ** k * h - 2. * lmin / 3.)) * rtip + rtip ** 2.)) / 3.
                    | false ->
                        n * a * (a + 1.) ** (float k) * 3. * p * (p ** k * h - 2. * lmin / 3.) * pi * rtip ** 2.
                | false ->
                    match (radius (p ** (k + 1.) * h) > rtip) with
                    | true ->
                        n * a * (a + 1.) ** (float k) * pi * p ** (k + 1.) * h * ((radius (p ** (k+1.) * h)) ** 2. + (radius (p ** (k + 1.) * h)) * rtip + rtip ** 2.) / 3.
                    | false ->
                        n * a * (a + 1.) ** (float k) * p ** (k + 1.) * h * pi * rtip ** 2.
            
            volume <- volume + volToAdd
            k <- k + 1.

        k, volume


module Allometrics =

    open Götmark2016_ShrubModel

    let mass woodDensity volume =
        volume * woodDensity

    let massToVolume woodDensity mass =
        mass / woodDensity

    let shrubBiomass b a rtip p lmin k5 k6 n woodDensity radius =
        radius
        |> NiklasAndSpatz_Allometry.stemLength k5 k6
        |> shrubVolume b a rtip p lmin k5 k6 n |> snd
        |> mass woodDensity

    let shrubRadius b a rtip p lmin k5 k6 n woodDensity mass =
        let findRadius volume =
            let v x = x |> NiklasAndSpatz_Allometry.stemLength k5 k6 |> shrubVolume b a rtip p lmin k5 k6 n |> snd
            let f = (fun x -> (v x) - volume )
            Statistics.RootFinding.bisect 0 200 f 0.01 1000.00 1e-8 // Assumption that shrub radius is between 0.01 and 100.0cm.
        mass
        |> massToVolume woodDensity
        |> findRadius

    let shrubHeight k5 k6 radius =
        radius |> NiklasAndSpatz_Allometry.stemLength k5 k6


module GrowthLimitation =

    /// A rearranged version of a Monod model
    let hollingDiscModel a b h =
        Some <| fun r -> (a * r) / (1. + (a * b * h * r))

    /// TEST: An integrated supply and use model
    let saturatingSupplySaturatingGrowth r a b h rootMass =
        Some <| fun resource -> (a * r * rootMass * resource) / (1. + a * b * r * rootMass * resource + a * h * resource)

    /// Monod model once saturation has been reached
    let linear (a:float) =
        Some <| fun r -> a * r

    /// The resource enforces no limitation on growth, and is negated
    let none = None

    /// **Description**
    /// A monotonically increasing function of a resource `r`.
    /// **Parameters**
    ///   * `h` - soil resource concentration required for growth at half the maximum rate
    ///   * `r` - the current resource concentration
    let michaelisMenten h (r:float) =
        r / (h + r)

    /// From Jabot and Pottier 2012
    let monod k (r:float) =
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

    /// Radius in millimetres
    let toBiomassMM radiusMM = 
        radiusMM / 10. |> Allometrics.shrubBiomass Constants.Allometrics.b Constants.Allometrics.a Constants.Allometrics.rtip Constants.Allometrics.p Constants.Allometrics.lmin Constants.Allometrics.k5 Constants.Allometrics.k6 Constants.Allometrics.numberOfStems Constants.Allometrics.salixWoodDensity

    /// Biomass in grams
    let toRadiusMM biomassGrams = 
        let radiusCm = biomassGrams |> Allometrics.shrubRadius Constants.Allometrics.b Constants.Allometrics.a Constants.Allometrics.rtip Constants.Allometrics.p Constants.Allometrics.lmin Constants.Allometrics.k5 Constants.Allometrics.k6 Constants.Allometrics.numberOfStems Constants.Allometrics.salixWoodDensity
        radiusCm * 10.

    /// d15N to N availability. From Craine 2009, as shown in Craine 2015 (Plant and Soil).
    /// Assuming d15N is a linear index of N availability, the minimum supported value of d15N is -3.09, as 0 N availability.
    let d15NtoAvailability d15N =
        (100. * d15N + 309.) / 359.

    let shrubHeightCm radiusMM =
        radiusMM / 10. |> Allometrics.shrubHeight Constants.Allometrics.k5 Constants.Allometrics.k6


module GeometricConstraint = 

    /// Linear growth rate in dM/dt form.
    let none _ = 1.

    /// A dM/dt form of the von Bertalanffy monomollecular growth function, where M = mass.
    let vonBertalanffy k m = 
        (k / m)

    /// The dM/dt form of the Chapman-Richards growth function, where M = mass.
    let chapmanRichards k m =
        (1. - (m / k))


module FeedbackToSoil =

    let none b : float = 0. * b
    let withBiomassLoss alpha gammab b : float = alpha * b * gammab
