#load "../src/Bristlecone/bristlecone.fsx"

////////////////////////////////////////////////////
/// Yamal Salix lanata Shrub - Nitrogen Interactions
////////////////////////////////////////////////////

(* Shrub ring width modelled with a single
   resource limitation. The models presented are based on
   Tilman (1990 / 1988). *)

open Bristlecone
open Types.ParameterEstimation
open Types
open Time

// 0. Configure Options
// ----------------------------

module Options =

    let resolution = Annual
    let iterations = 5000


module Constants =

    // Empirically-derived parameters (centimetres - Niklas and Spatz):
    let k5 =  19.98237 // Allometric fit to Yamal shrub radius-length data #1
    let k6 = 0.42091 // Allometric fit to Yamal shrub radius-length data #2

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

            let radius h = basalRadius k5 k6 h

            let mainStemVolume =
                match radius h with
                | r when r > rtip -> n * pi * h * ((radius h) ** 2. + (radius h) * rtip + rtip ** 2.) / 3.
                | _ -> n * pi * h * rtip ** 2.

            let mutable volume = mainStemVolume
            let mutable k = 0.

            // while ((p ** k) * h > lmin * 2./3.) do
            //     let volToAdd =
            //         match ((p ** k) * h < lmin) with
            //         | true ->
            //             // If the branch to be addded is less than lmin, then it is scaled to make the function continuous
            //             match (b * 3. * p * ((p ** k) * h - 2. * lmin / 3.) > rtip) with
            //             | true ->
            //                 // Child branch = truncated cone
            //                 n * a * (a + 1.) ** k * pi * 3. * p * (p ** k * h - 2. * lmin / 3.) * ((radius (3. * p * (p ** k * h - 2. * lmin / 3.))) * (radius (3. * p * (p ** k * h - 2. * lmin / 3.)) * rtip + rtip ** 2.)) / 3.
            //             | false ->
            //                 // Child branch = cylinder
            //                 n * a * (a + 1.) ** k * 3. * p * (p ** k * h - 2. * lmin / 3.) * pi * rtip ** 2.
            //         | false ->
            //             match (radius (p ** (k + 1.) * h) > rtip) with
            //             | true ->
            //                 // Child branch = truncated cone
            //                 n * a * (a + 1.) ** k * pi * p ** (k + 1.) * h * ((radius (p ** (k+1.) * h)) ** 2. + (radius (p ** (k + 1.) * h)) * rtip + rtip ** 2.) / 3.
            //             | false ->
            //                 // Child branch = cylinder
            //                 n * a * (a + 1.) ** k * p ** (k + 1.) * h * pi * rtip ** 2.
            //     volume <- volume + volToAdd
            //     k <- k + 1.
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


// 2. Create Hypotheses
// ----------------------------

let ``nutrient-dependent growth`` =

    /// d15N to N availability. From Craine 2009, as shown in Craine 2015 (Plant and Soil).
    /// Assuming d15N is a linear index of N availability, the minimum supported value of d15N is -3.09, as 0 N availability.
    let d15NtoAvailability d15N =
        (100. * d15N + 309.) / 359.

    /// Radius in millimetres
    let toBiomass radius = 
        radius / 10. |> ModelComponents.Allometrics.shrubBiomass Constants.b Constants.a Constants.rtip Constants.p Constants.lmin Constants.k5 Constants.k6 Constants.numberOfStems Constants.salixWoodDensity

    /// Biomass in grams
    let toRadius biomass = 
        let radiusCm = biomass |> ModelComponents.Allometrics.shrubRadius Constants.b Constants.a Constants.rtip Constants.p Constants.lmin Constants.k5 Constants.k6 Constants.numberOfStems Constants.salixWoodDensity
        radiusCm * 10.

    // /// Conversion of d15N proxy into bioavailable N - explicit fractionation
    // let isotopeToBioavailableN fractionationEffect proxyError d15N : float =
    //     d15N + fractionationEffect + proxyError

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
    let dbsdt' bs n fMax h al ar m r = 
        let bl = bs * al
        let br = bs * ar
        let photosyntheticCapacity = bl * (f fMax h (effectiveNutrient n br bl)) - r * (bs + bl + br)
        // printfn "B = %f (Stem %f) (Leaf %f) (Root %f). Photosynthetic capacity = %fg" (bs + bl + br) bs bl br photosyntheticCapacity
        // printfn "B Change = %f" ((1. / (1. + al + ar )) * (photosyntheticCapacity - m * (bs + bl + br)))
        (1. / (1. + al + ar )) * (photosyntheticCapacity - m * (bs + bl + br))

    /// Bioavailable soil nitrogen [N]
    let dndt' (bs:float) (n:float) lambda gamman q fMax h al ar r : float =
        let bl = bs * al
        let br = bs * ar
        let photosyntheticCapacity = bl * (f fMax h (effectiveNutrient n br bl)) - r * (bs + bl + br)
        // printfn "N = %f; Env gain = %f; Env loss = %f; Plant use = %f" n lambda (gamman * n) (q * photosyntheticCapacity)
        // printfn "N Change = %f" (lambda - q * photosyntheticCapacity - gamman * n)
        lambda - q * photosyntheticCapacity - gamman * n

    /// Measurement variable: stem radius [rw].
    /// If there has been biomass accumulation during the time interval, wood accumulation occurs, according to shrub allometric relationships. 
    let drwdt' bs n fMax h al ar m r = 
        let biomassStemChange = dbsdt' bs n fMax h al ar m r
        // printfn "W: Stem Biomass growth of %f" biomassStemChange
        if biomassStemChange > 0.
            then 
                let oldRadius = bs |> toRadius
                let newRadius = (bs + biomassStemChange) |> toRadius
                // printfn "Biomass stem change is positive: %f, %f, %f" biomassStemChange newRadius oldRadius
                newRadius - oldRadius
            else 0.

    /// Bristlecone function for dBs/dt
    let dbsdt p _ bs (e:Environment) =
        dbsdt' bs ((e.[ShortCode.create "N"]) |> d15NtoAvailability) (p |> Pool.getEstimate "fMax") (p |> Pool.getEstimate "h")
            (p |> Pool.getEstimate "al") (p |> Pool.getEstimate "ar") (p |> Pool.getEstimate "gammab") (p |> Pool.getEstimate "r")

    /// Bristlecone function for dN/dt
    let dndt p _ n (e:Environment) =
        dndt' (e.[ShortCode.create "x"]) (n |> d15NtoAvailability) (p |> Pool.getEstimate "lambda") (p |> Pool.getEstimate "gamman") 
            (p |> Pool.getEstimate "q") (p |> Pool.getEstimate "fMax") (p |> Pool.getEstimate "h") (p |> Pool.getEstimate "al") (p |> Pool.getEstimate "ar") (p |> Pool.getEstimate "r")

    /// Bristlecone function for dr/dt
    let drwdt p _ _ (e:Environment) =
        drwdt' (e.[ShortCode.create "bs"]) ((e.[ShortCode.create "N"]) |> d15NtoAvailability) (p |> Pool.getEstimate "fMax") 
            (p |> Pool.getEstimate "h") (p |> Pool.getEstimate "al") (p |> Pool.getEstimate "ar") (p |> Pool.getEstimate "gammab") (p |> Pool.getEstimate "r")

    { Equations  = [ ShortCode.create "x",      drwdt
                     ShortCode.create "bs",     dbsdt
                     ShortCode.create "N",      dndt ] |> Map.ofList
      Parameters = [ // for growth function
                     ShortCode.create "fMax",   Parameter.create Unconstrained   1.000 1.001   // Maximum rate of resource-saturated growth per unit biomass
                     ShortCode.create "h",      Parameter.create Unconstrained   0.001 0.002   // Soil resource availability for growth at half of maximum rate
                     // for nitrogen replenishment
                     ShortCode.create "lambda", Parameter.create Unconstrained   0.010 0.020   // Rate of nitrogen replenishment
                     ShortCode.create "gamman", Parameter.create Unconstrained   0.001 0.005   // Loss rate of nitrogen
                     // for shrub allocation and physiology
                     ShortCode.create "al",     Parameter.create Unconstrained   1.000 1.000   // Allocation fraction to leaves, relative to stem
                     ShortCode.create "ar",     Parameter.create Unconstrained   1.000 1.000   // Allocation fraction to roots, relative to stem
                     ShortCode.create "r",      Parameter.create Unconstrained   0.001 0.005   // Respiration cost per unit biomass
                     ShortCode.create "q",      Parameter.create Unconstrained   0.001 0.005   // Nutrient requirement per unit biomass
                     ShortCode.create "gammab", Parameter.create Unconstrained   0.001 0.005   // Loss rate of biomass
                     // for likelihood function
                    //  ShortCode.create "sigmax", Parameter.create Unconstrained -0.50 0.50   // Standard deviation of x (biomass)
                    //  ShortCode.create "sigmay", Parameter.create Unconstrained -0.50 0.50   // Standard deviation of y (nitrogen)
                    //  ShortCode.create "rho",    Parameter.create Unconstrained -0.25 0.25   // Covariance between growth and nitrogen
                    ] |> Map.ofList
      Likelihood = ModelLibrary.Likelihood.sumOfSquares ["x"; "N"] }


// // 2. Test Hypotheses Work
// // ----------------------------

// // let startValues = [ ShortCode.create "x", 0.23; ShortCode.create "N", 4.64; ShortCode.create "bs", 471.5475542] |> Map.ofList

// // ``nutrient-dependent growth``
// // |> Bristlecone.test' Options.resolution 1000 Options.testSeriesLength startValues


// // 3. Load Real Data and Estimate
// // ----------------------------

// let shrubs = 
//     let yuribei = DataAccess.Shrub.loadRingWidths (__SOURCE_DIRECTORY__ + "/data/yuribei-rw.csv")
//     let d15N = DataAccess.Shrub.loadLocalEnvironmentVariable (__SOURCE_DIRECTORY__ + "/data/yuribei-d15N-imputed.csv")
//     yuribei
//     |> Seq.map (fun s -> s.Identifier.Value, s)
//     |> Seq.keyMatch d15N
//     |> Seq.map (fun (_,plant,d15N) -> PlantIndividual.zipEnv (ShortCode.create "N") plant d15N)
//     |> Seq.toList

// // Start value configuration
// let getStartValues startDate (plant:Types.PlantIndividual.PlantIndividual) =
//     let initialRadius =
//         match plant.Growth with
//         | Types.PlantIndividual.PlantGrowth.RingWidth s -> 
//             match s with
//             | GrowthSeries.Absolute c -> c.Head |> fst |> removeUnit
//             | GrowthSeries.Cumulative c -> 6.3 //c.Head |> fst |> removeUnit
//             | _ -> invalidOp "Not implemented"
//         | _ -> invalidOp "Not implemented 2"
//     let initialMass = initialRadius |> removeUnit //|> toBiomass
//     let initialNitrogen = plant.Environment.[ShortCode.create "N"].Head |> fst
//     [ ShortCode.create "x", initialRadius
//       ShortCode.create "N", initialNitrogen 
//       ShortCode.create "bs", initialMass ] |> Map.ofList

// // Match environment and growth lengths
// // If n-1 time point has a value in the original set, append to new set

// // Run 4 chains at once to assess convergence

// let estimated =
//     [| 1 .. 8 |]
//     |> Array.Parallel.map(fun _ ->
//         [shrubs.[0]] |> List.map (fun s ->
//             let shrub = s |> PlantIndividual.toCumulativeGrowth
//             let common = shrub |> PlantIndividual.keepCommonYears
//             let startDate = common.Environment.[ShortCode.create "N"] |> TimeSeries.start
//             let startConditions = getStartValues startDate shrub
//             common
//             |> Bristlecone.estimatePlant Options.resolution Options.iterations ``nutrient-dependent growth`` (Custom startConditions))
//     )

// // 4. Plot Results (using R)
// // ----------------------------

#load "plotting.fsx"
open RProvider
open RProvider.graphics
open RProvider.ggplot2

// // i. Plot likelihood by iteration (for MCMC)
// R.plot( namedParams [ "x", box (estimated.[0].[0].Trace |> fst |> List.rev |> List.skip 10) ; "type", box "l"; "xlab", box "Iteration"; "ylab", box "-log likelihood" ]) |> ignore

// // i. Estimated versus Observed Series
// R.par(namedParams [ "mfrow", [2;2]] ) |> ignore
// R.plot (namedParams [ "x",box estimated.[0].[0].Series.[ShortCode.create "x"].Observed; "type", box "l"; "xlab", box "Year"; "ylab", box "Stem Radius (Observed)" ]) |> ignore
// R.plot (namedParams [ "x",box estimated.[0].[0].Series.[ShortCode.create "x"].Expected; "type", box "l"; "xlab", box "Year"; "ylab", box "Stem Radius (Model)" ]) |> ignore
// R.plot (namedParams [ "x",box estimated.[0].[0].Series.[ShortCode.create "N"].Observed; "type", box "l"; "xlab", box "Year"; "ylab", box "Nitrogen (Observed)" ]) |> ignore
// R.plot (namedParams [ "x",box estimated.[0].[0].Series.[ShortCode.create "N"].Expected; "type", box "l"; "xlab", box "Year"; "ylab", box "Nitrogen (Model)" ]) |> ignore

// // iii. For each parameter, find mean and standard deviation, and plot graph with histogram behind it

// let parameterValues = 
//     [0 .. ``nutrient-dependent growth``.Parameters.Count - 1]
//     |> List.map(fun i ->
//             estimated.[0].Trace
//             |> snd
//             |> List.map (fun x -> x.[i])
//             |> List.toArray )

// let parameters = ``nutrient-dependent growth``.Parameters |> Map.toList

// let everyNth n seq = 
//     seq |> Seq.mapi (fun i el -> el, i)              // Add index to element
//         |> Seq.filter (fun (el, i) -> i % n = n - 1) // Take every nth element
//         |> Seq.map fst                               // Drop index from the result

// R.par(namedParams [ "mfrow", [4;4]] ) |> ignore
// [0 .. ``nutrient-dependent growth``.Parameters.Count - 1]
// |> List.map (fun i -> R.plot( namedParams [ "x", box (parameterValues.[i] |> Array.rev |> Array.skip 500 |> everyNth 100) ; "type", box "l"; "xlab", box "Iteration"; "ylab", box ((parameters.[i] |> fst).Value) ]) |> ignore)

// // Alternative: kernel density plot for each parameter
// open RProvider.stats

// R.par(namedParams [ "mfrow", [4;4]] ) |> ignore
// [0 .. ``nutrient-dependent growth``.Parameters.Count - 1]
// |> List.map (fun i -> 
//     let d = parameterValues.[i] |> Array.rev |> Array.skip 100 |> R.density
//     R.plot(d) )

// // Correlation between parameters
// let correlation = (MathNet.Numerics.Statistics.Correlation.PearsonMatrix parameterValues).ToArray()
// R.image(correlation)


// ///// Use GGMCMC package to generate outputs
// /// 
// open RProvider.ggmcmc

// // 1. Convert data into a data frame with four variables:
// // > Iteration Number of iteration.
// // > Chain Number of the chain.
// // > Parameter Name of the parameter.
// // > value value sampled.

// let chainsDataFrame =
//     estimated
//     |> Array.mapi (fun chainNumber chain -> 
//         chain.[0].Trace
//         |> snd
//         |> List.rev
//         |> List.skip 0 // Skip burn in period
//         |> everyNth 5 // Thinning by this amount
//         |> Seq.mapi (fun iterationNumber values ->
//             chain.[0].Parameters
//             |> Map.toList
//             |> List.mapi(fun i (name,p) -> 
//                 iterationNumber + 1,
//                 chainNumber + 1,
//                 name.Value,
//                 values.[i] )
//         )
//         |> List.concat
//         |> Seq.toList
//     )
//     |> Array.toList
//     |> List.concat

// let df =
//     namedParams [
//         "Iteration", chainsDataFrame |> List.map(fun (a,_,_,_) -> a) |> box 
//         "Chain", chainsDataFrame |> List.map(fun (_,a,_,_) -> a) |> box
//         "Parameter", chainsDataFrame |> List.map(fun (_,_,a,_) -> a) |> box
//         "value", chainsDataFrame |> List.map(fun (_,_,_,a) -> a) |> box
//     ]
//     |> R.data_frame

// // Save out R data frame of MCMC results
// R.assign("chains", df)
// R.save(list = ["chains"], file = "dphil-shrub-chain.rdata")


// Plots of Allometric Fit
// __________________________

let (++) (plot1:RDotNet.SymbolicExpression) (plot2:RDotNet.SymbolicExpression) = 
    R.``+``(plot1, plot2)

(* A. Shrub Height - Radius Relationship *)
let shrubHeightVersusRadius =
    let df = 
        namedParams [
            "Radius", [1. .. 400. ] |> List.map (ModelComponents.Götmark2016_ShrubModel.basalRadius Constants.k5 Constants.k6) |> box
            "Height", [ 1. .. 400. ] |> box
        ] |> R.data_frame
    R.ggplot(df, R.aes__string(namedParams [ "x", box "Radius"; "y", box "Height" ])) 
    ++ R.labs(namedParams ["x", box "Basal Radius (cm)"; "y", box "Stem Length (cm)" ])
    ++ R.geom__line()


(* B. Shrub Height - Volume Relationship *)
let shrubHeightVolume =
    let volume height =
        height |> ModelComponents.Götmark2016_ShrubModel.shrubVolume Constants.b Constants.a Constants.rtip Constants.p Constants.lmin Constants.k5 Constants.k6 Constants.numberOfStems
    let df3 =
        namedParams [
            "StemVolume", [1. .. 180. ] |> List.map (volume >> snd) |> box
            "BranchGenerations", [1. .. 180. ] |> List.map (volume >> fst) |> box
            "Height", [ 1. .. 180. ] |> box
        ] |> R.data_frame
    R.ggplot(df3, R.aes__string(namedParams [ "x", box "StemVolume"; "y", box "Height" ])) 
    ++ R.geom__line()


(* C. Shrub Radius - Stem Volume Relationship *)
let shrubRadiusVolume =

    let toBiomass stems radiusMM = 
        radiusMM
        |> ModelComponents.Allometrics.shrubBiomass Constants.b Constants.a Constants.rtip Constants.p Constants.lmin Constants.k5 Constants.k6 stems Constants.salixWoodDensity

    let toRadius biomassGrams = 
        biomassGrams
        |> ModelComponents.Allometrics.shrubRadius Constants.b Constants.a Constants.rtip Constants.p Constants.lmin Constants.k5 Constants.k6 Constants.numberOfStems 1.
        |> (/) 1. // Convert cm back to mm

    let r,b,s =
        [ 2.; 3.; 5. ]
        |> List.collect(fun stems ->
            [ 1. .. 50. ] // Radius in mm
            |> List.map(fun r -> r, toBiomass (float stems) r, float stems ) )
        |> List.unzip3

    let df = 
        namedParams [
            "StemRadius", r |> box
            "Volume", b |> box
            // "Branches", br |> box
            "StemCount", s |> box ]
        |> R.data_frame

    R.ggplot(df, R.aes__string(namedParams [ "x", box "StemRadius"; "y", box "Volume"; "group", box "StemCount"; "color", box "StemCount" ])) 
    ++ R.geom__line()


