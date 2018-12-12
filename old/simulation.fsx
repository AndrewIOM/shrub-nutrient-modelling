#load "../src/Bristlecone/bristlecone.fsx"

////////////////////////////////////////////////////
/// Yamal Salix lanata Shrub - Nitrogen Interactions
////////////////////////////////////////////////////

(* Shrub ring width modelled with a single
   resource limitation. The models presented are based on
   Tilman (1990 / 1988). *)

open DendroFit
open Types.ParameterEstimation
open Types
open Time



module LikelihoodTest =

    let private getData s (predictions:CodedMap<PredictedSeries>) = predictions.Item (ShortCode.create s)

    let bivariateGaussian' (p:ParameterPool) obsx obsy expx expy = 
        let diffx = obsx - expx
        let diffy = obsy - expy
        let sigmax = p |> Pool.getEstimate "sigmax"
        let sigmay = p |> Pool.getEstimate "sigmay"
        let rho = p |> Pool.getEstimate "rho"
        let zta1 = diffx ** 2. / sigmax ** 2.
        let zta2 = 2. * rho * ((diffx * diffy) / sigmax * sigmay)
        let zta3 = diffy ** 2. / sigmay ** 2.
        let z = zta1 - zta2 + zta3
        printfn "X: %f %f - Y: %f %f" obsx expx obsy expy
        printfn "sigmax = %f; sigma y = %f; rho = %f; L = %f (%f)" sigmax sigmay rho (-log((1./(2.*System.Math.PI*sigmax*sigmay*sqrt(1.-(rho*rho)))) * exp (-z / (2. * (1. - (rho * rho)))))) (((1./(2.*System.Math.PI*sigmax*sigmay*sqrt(1.-(rho*rho)))) * exp (-z / (2. * (1. - (rho * rho))))))
        -log((1./(2.*System.Math.PI*sigmax*sigmay*sqrt(1.-(rho*rho)))) * exp (-z / (2. * (1. - (rho * rho)))))


    let bivarGaussTest' sigmax sigmay rho obsx obsy expx expy =



        2. 

    /// <summary> 
    /// Log likelihood function for dual simultaneous system, assuming Gaussian error for both x and y.
    /// </summary> 
    let bivariateGaussian key1 key2 p data = 
        let x = data |> getData key1
        let y = data |> getData key2
        [1 .. (Array.length x.Observed) - 1] 
        |> List.sumBy (fun i -> (bivariateGaussian' p x.Observed.[i] y.Observed.[i] x.Expected.[i] y.Expected.[i])) 


// 0. Configure Options
// ----------------------------

module Options =

    let resolution = Annual
    let iterations = 10000
    let testSeriesLength = 50


module Constants =

    // Empirically-derived parameters:
    let k5 = 80.200334 // Allometric fit to Yamal shrub BD-length data #1
    let k6 = 9.918091 // Allometric fit to Yamal shrub BD-length data #2

    // Constants from the literature:
    let a = 2. // the number of child branches added to previous branches (including the tops of the stems) for a shrub
    let p = 0.5 // the length of a child branch as a proportion of its parent branch/stem
    let lmin = 20. //cm. the length at which a stem or branch gets child branches
    let rtip = 0.1 //cm. the radius of the outermost tip of a stem or branch
    let b = 0.0075 // the ratio of the basal radius of a stem or branch and its length
    let salixWoodDensity = 0.5 // g / cm3 (from internet)
    let numberOfStems = 2.2


// 1. Setup model components
// ----------------------------

module ModelComponents =

    let pi = System.Math.PI

    module Götmark2016_ShrubModel =

        /// Gives the basal radius in centimeters of a stem/branch given its length in centimeters. Function from Niklas and Spatz (2004). 
        let basalRadius k5 k6 stemLength =
            100. * (( 0.01 * stemLength + k6) / k5) ** (3. / 2.) / 2.

        /// Total shrub volume given height and number of stems
        let shrubVolume b a rtip p lmin k5 k6 n h =

            let radius = basalRadius k5 k6
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

        let stemLength k5 k6 radius =
            k5 * (2. * radius) ** (2. / 3.) - k6

        let mass woodDensity volume =
            volume * woodDensity

        let massToVolume woodDensity mass =
            mass / woodDensity

        let shrubBiomass b a rtip p lmin k5 k6 n woodDensity radius =
            radius
            |> stemLength k5 k6
            |> shrubVolume b a rtip p lmin k5 k6 n |> snd
            |> mass woodDensity

        let shrubRadius b a rtip p lmin k5 k6 n woodDensity mass =
            let findRadius volume =
                let v x = x |> stemLength k5 k6 |> shrubVolume b a rtip p lmin k5 k6 n |> snd
                let f = (fun x -> 
                    //printfn "~ %f / %f" (v x) x
                    (v x) - volume )
                // printfn "Finding solution to inverse allometric - volume is %f" volume
                // MathNet.Numerics.RootFinding.Bisection.FindRoot(System.Func<float,float> f, 0.01, 1000., 1e-8, 1000)
                ODE.RootFinding.bisect 0 200 f 0.01 1000.00 1e-8 // Assumption that shrub radius is between 0.01 and 20.0cm.
            mass
            |> massToVolume woodDensity
            |> findRadius

    module GrowthRate =

        /// **Description**
        /// A monotonically increasing function of a resource `r`.
        /// **Parameters**
        ///   * `fMax` - maximum specific growth rate
        ///   * `h` - soil resource concentration required for growth at half the maximum rate
        ///   * `r` - the current resource concentration
        let michaelisMenten fMax h (r:float) =
            fMax * (r / (h + r))

    module AbioticResource =

        /// **Description**
        /// A standard chemostat-type model for replenshment of an abiotic resource.
        /// **Parameters**
        ///   * `d` - a rate constant
        ///   * `s` - resource concentration of the inflow
        ///   * `n` - the current resource concentration
        let chemostat d s (n:float) =
            d * (s - n)


// 2. Create Hypotheses
// ----------------------------

module Model =

    /// Radius in millimetres
    let toBiomass radius = 
        radius / 10. |> ModelComponents.Allometrics.shrubBiomass Constants.b Constants.a Constants.rtip Constants.p Constants.lmin Constants.k5 Constants.k6 Constants.numberOfStems Constants.salixWoodDensity

    /// Biomass in grams
    let toRadius biomass = 
        let radiusCm = biomass |> ModelComponents.Allometrics.shrubRadius Constants.b Constants.a Constants.rtip Constants.p Constants.lmin Constants.k5 Constants.k6 Constants.numberOfStems Constants.salixWoodDensity
        radiusCm * 10.

    /// Conversion of d15N proxy into bioavailable N - explicit fractionation
    let isotopeToBioavailableN fracScale d15N : float =
        d15N * fracScale

    /// Uptake rate (modelled as 'effective nutrient availability')
    let effectiveNutrient n rootMass leafMass : float = 
        n * rootMass / leafMass

    /// Resource-dependent growth function
    /// AKA Specific growth rate
    /// Maximum rate of photosynthesis per unit of photosynthetic tissue under optimal conditions
    let f = ModelComponents.GrowthRate.michaelisMenten

    /// Nitrogen replenishment rate.
    let y = ModelComponents.AbioticResource.chemostat
    
    /// Cumulative stem biomass [bs].
    let dbsdt' bs n fMax h al ar r m = 
        let bl = bs * al
        let br = bs * ar
        let wholePlantGrowth = bl * (f fMax h (effectiveNutrient n br bl)) - r * (bs + bl + br) - m * (bs + bl + br)
        // printfn "Start biomass is (Stem %f) (Leaf %f) (Root %f). Total growth during this time is %fg" bs bl br wholePlantGrowth
        printfn "B: Biomass growth of %f" wholePlantGrowth
        wholePlantGrowth * (1. / (1. + al + ar ))

    /// Bioavailable soil nitrogen [N]
    let dndt' (bs:float) (n:float) d s q fMax h al ar r : float =
        let bl = bs * al
        let br = bs * ar
        let wholePlantGrowth = bl * (f fMax h (effectiveNutrient n br bl)) - r * (bs + bl + br)
        // printfn "N: Biomass growth of %f" wholePlantGrowth
        (y d s n) - q * wholePlantGrowth

    /// Measurement variable: stem radius [rw].
    /// If there has been biomass accumulation during the time interval, wood accumulation occurs, according to shrub allometric relationships. 
    let drwdt' bs n r k al ar respR eb = 
        let biomassStemChange = dbsdt' bs n r k al ar respR eb
        printfn "W: Stem Biomass growth of %f" biomassStemChange
        if biomassStemChange > 0.
            then 
                let oldRadius = bs |> toRadius
                let newRadius = (bs + biomassStemChange) |> toRadius
                printfn "Positive stem biomass: %f, %f, %f, %f" biomassStemChange newRadius oldRadius (newRadius - oldRadius)
                newRadius - oldRadius
            else 
                printfn "Negative stem biomass (%f)" biomassStemChange
                0.

    /// DendroFit function for dBs/dt
    let dbsdt p _ bs (e:Environment) =
        dbsdt' bs (e.[ShortCode.create "N"]) (p |> Pool.getEstimate "fMax") (p |> Pool.getEstimate "h")
            (p |> Pool.getEstimate "al") (p |> Pool.getEstimate "ar") (p |> Pool.getEstimate "r") (p |> Pool.getEstimate "m")

    /// DendroFit function for dN/dt
    let dndt p _ n (e:Environment) =
        dndt' (e.[ShortCode.create "x"]) n (p |> Pool.getEstimate "d") (p |> Pool.getEstimate "s") (p |> Pool.getEstimate "q") (p |> Pool.getEstimate "fMax")
            (p |> Pool.getEstimate "h") (p |> Pool.getEstimate "al") (p |> Pool.getEstimate "ar") (p |> Pool.getEstimate "r")

    /// DendroFit function for dr/dt
    let drwdt p _ _ (e:Environment) =
        drwdt' (e.[ShortCode.create "bs"]) (e.[ShortCode.create "N"]) (p |> Pool.getEstimate "fMax") (p |> Pool.getEstimate "h")
            (p |> Pool.getEstimate "al") (p |> Pool.getEstimate "ar") (p |> Pool.getEstimate "r") (p |> Pool.getEstimate "m")

    let system =
        { Equations  = [ ShortCode.create "x",      drwdt
                         ShortCode.create "bs",     dbsdt
                         ShortCode.create "N",      dndt ] |> Map.ofList
          Parameters = [ // for growth function
                         ShortCode.create "fMax",   Parameter.create PositiveOnly   0.01 1.00   // Maximum rate of resource-saturated growth per unit biomass
                         ShortCode.create "h",      Parameter.create PositiveOnly   0.01 1.00   // Soil resource concentration for growth at half of maximum rate
                         // for nitrogen replenishment
                         ShortCode.create "d",      Parameter.create PositiveOnly   0.01 1.00   // Rate of nitrogen replenishment
                         ShortCode.create "s",      Parameter.create PositiveOnly   0.01 1.00   // Nitrogen concentration of the environmental inflow
                         // for shrub allocation and physiology
                         ShortCode.create "al",     Parameter.create PositiveOnly   0.01 0.50   // Allocation fraction to leaves, relative to stem
                         ShortCode.create "ar",     Parameter.create PositiveOnly   0.01 5.00   // Allocation fraction to roots, relative to stem
                         ShortCode.create "r",      Parameter.create PositiveOnly   0.0001 1.00   // Respiration rate per unit biomass
                         ShortCode.create "m",      Parameter.create PositiveOnly   0.01 1.00   // Biomass loss rate (e.g. from herbivory)
                         ShortCode.create "q",      Parameter.create PositiveOnly   0.00000001 0.75   // Nutrient concentration per unit biomass
                         // for likelihood function
                         ShortCode.create "sigmax", Parameter.create Unconstrained -0.25 0.25   // Standard deviation of x (biomass)
                         ShortCode.create "sigmay", Parameter.create Unconstrained -0.25 0.25   // Standard deviation of y (nitrogen)
                         ShortCode.create "rho",    Parameter.create Unconstrained -0.15 0.15   // Covariance between growth and nitrogen
                        ] |> Map.ofList
          Likelihood = LikelihoodTest.bivariateGaussian "x" "N" }


// Parameters - playground
open MathNet.Numerics.LinearAlgebra

let simulation =
    // This set works:
    // let fMax = 0.50
    // let h = 1.0
    // let al = 0.2
    // let ar = 1.00
    // let r = 0.0001
    // let m = 0.02
    // let d = 0.10
    // let s = 0.30
    // let q = 0.0004 // OK
    let fMax = 0.50
    let h = 0.8
    let al = 0.3
    let ar = 0.5
    let r = 0.0001
    let m = 0.02
    let d = 0.003
    let s = 0.25
    let q = 0.00000001 // OK

    let initial = [
        6.21 |> Model.toBiomass
        4.3
        6.21 ]

    let xEq t x y z = Model.dbsdt' x y fMax h al ar r m
    let yEq t x y z = Model.dndt' x y d s q fMax h al ar r
    let zEq t x y z = Model.drwdt' x y fMax h al ar r m

    let tInitial = 1.
    let tEnd = 35.
    let tStep = 1.

    let rp t (x:Vector<float>) = 
        [ xEq; yEq; zEq ]
        |> List.map (fun m -> m t x.[0] x.[1] x.[2])
        |> vector

    let n = (tEnd - tInitial + 1.) / tStep |> int

    let f = System.Func<float, Vector<float>, Vector<float>> rp
    let result = MathNet.Numerics.OdeSolvers.RungeKutta.FourthOrder(initial |> vector, tInitial, tEnd, n, f)
    
    [0 .. 2]
    |> Seq.mapi(fun i _ -> result |> Array.map(fun x -> x.[i]))
    |> Seq.toList


// Plot simulation

#load "plotting.fsx"
open RProvider
open RProvider.graphics
open RProvider.ggplot2

R.par  (namedParams [ "mfrow", [1;3]] ) |> ignore
R.plot (namedParams [ "x",box (simulation.[0]); "type", box "l"; "xlab", box "Year"; "ylab", box "Stem Biomass (g)" ]) |> ignore
R.plot (namedParams [ "x",box (simulation.[1]); "type", box "l"; "xlab", box "Year"; "ylab", box "Soil Nitrogen Concentration (d15N)" ]) |> ignore
R.plot (namedParams [ "x",box (simulation.[2]); "type", box "l"; "xlab", box "Year"; "ylab", box "Stem Radius (mm)" ]) |> ignore

// Notes
// ---------
// > Total / stem biomass should not be able to go negative - how to fix?
//      - Was caused by choosing an initial size below 1mm. Allometrics don't work at low size
// > Nitrogen concentration should not be able to go negative - how to fix?
//      - Fixed. Was an issue with biomass calculation in N equation (bs and n were mixed up, and was not multipled by leaf mass)

let fakePlantData =
    [
        ShortCode.create "x", simulation.[0] |> TimeSeries.create (System.DateTime(1990,01,01)) (System.TimeSpan.FromDays(365.))
        ShortCode.create "N", simulation.[1] |> TimeSeries.create (System.DateTime(1990,01,01)) (System.TimeSpan.FromDays(365.))
        ShortCode.create "bs",simulation.[2] |> TimeSeries.create (System.DateTime(1990,01,01)) (System.TimeSpan.FromDays(365.))
    ] |> Map.ofList
let estimate = DendroFit.estimate' Annual 1 Model.system StartingValues.FirstDataItem fakePlantData