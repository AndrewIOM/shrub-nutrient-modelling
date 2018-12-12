
#load "../src/bristlecone.fsx"

////////////////////////////////////////////////////
/// Yamal Salix lanata Shrub - Nitrogen Interactions
////////////////////////////////////////////////////

// Shrub ring width modelled with a single
// resource limitation. The models presented are based on
// Tilman (1990 / 1988).

open Bristlecone
open Bristlecone.ModelSystem
open Bristlecone.PlantIndividual

// 0. Configure Options
// ----------------------------

module Options =

    let resolution = Annual
    let iterations = 1
    let burn = 12000
    let chains = 3
    let testSeriesLength = 40


module Constants =

    // Empirically-derived parameters:
    let k5 = 19.98239 // Allometric fit to Yamal shrub BD-length data #1 (in centimetres)
    let k6 = 0.42092 // Allometric fit to Yamal shrub BD-length data #2 (in centimetres)

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
                Optimisation.RootFinding.bisect 0 200 f 0.01 1000000.00 1e-8 // Assumption that shrub radius is between 0.01 and 100.0cm.
            mass
            |> massToVolume woodDensity
            |> findRadius

    module GrowthLimitation =

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

module Proxies =

    /// Radius in millimetres
    let toBiomassMM radiusMM = 
        radiusMM / 10. |> ModelComponents.Allometrics.shrubBiomass Constants.b Constants.a Constants.rtip Constants.p Constants.lmin Constants.k5 Constants.k6 Constants.numberOfStems Constants.salixWoodDensity

    /// Biomass in grams
    let toRadiusMM biomassGrams = 
        let radiusCm = biomassGrams |> ModelComponents.Allometrics.shrubRadius Constants.b Constants.a Constants.rtip Constants.p Constants.lmin Constants.k5 Constants.k6 Constants.numberOfStems Constants.salixWoodDensity
        radiusCm * 10.

    /// d15N to N availability. From Craine 2009, as shown in Craine 2015 (Plant and Soil).
    /// Assuming d15N is a linear index of N availability, the minimum supported value of d15N is -3.09, as 0 N availability.
    let d15NtoAvailability d15N =
        (100. * d15N + 309.) / 359.


module Allocation =

    type PlantBiomass = {
        Leaf: float
        Root: float
        Stem: float }

    /// All biomass is assumed to be photosynthetic, with no allocation based on tissue function. 
    let none b = { Stem = b; Root = 0.; Leaf = 0. }

    /// Roots, leaves, and stem are in equal 1:1:1 ratio. 
    let equal bs = { Stem = bs; Root = bs; Leaf = bs }

    /// Floating allocation. 
    let floating al ar bs = { Stem = bs; Leaf = (bs * al); Root = (bs * ar) }


module GrowthRate = 

    let linear x = x
    let logistic x k = 1. - (x / k)



module FeedbackToSoil =

    let none b : float = 0. * b
    let withBiomassLoss alpha gammab b : float = alpha * b * gammab


module BaseEquations =
    
    open Allocation

    /// **Description**
    /// Cumulative stem biomass [db/dt].
    /// If bs>0 and bl and br are set at zero, allocation is disabled and bs represents total plant biomass.
    /// If bl and br are >0, allocation enabled and equation returns stem biomass production.
    /// 
    /// **Parameters**
    ///   * `bs` - Current stem biomass
    ///   * `bl` - Current leaf biomass
    ///   * `br` - Current root biomass
    ///   * `fMax` - Maximum rate of resource-saturated growth
    ///   * `m` - Loss rate of biomass due to disturbance
    ///   * `r` - Respiration rate per unit biomass
    ///
    /// **Output Type**
    ///   * `float` - Change in [stem] biomass during this time point
    let plantBiomass (b:PlantBiomass) fMax m r = 
        (b.Stem * fMax - r * (b.Stem + b.Leaf + b.Root)) - m * (b.Stem + b.Leaf + b.Root)

    let soilResource (b:PlantBiomass) (n:float) lambda gamman q fMax r feedback : float =
        lambda - q * (b.Stem * fMax - r * (b.Stem * b.Leaf + b.Root)) - gamman * n + feedback (b.Stem + b.Leaf + b.Root)


let ``base model`` fMax allocationMode nitrogenFeedback additionalParameters =

    /// Cumulative stem biomass [bs].
    let dbsdt' (bs:float) n gammab r fMax allocationMode = 
        BaseEquations.plantBiomass (bs |> allocationMode) (fMax n) gammab r

    /// Bioavailable soil nitrogen [N]
    let dndt' bs n lambda gamman q r fMax allocationMode feedback = 
        BaseEquations.soilResource (bs |> allocationMode) n lambda gamman q (fMax n) r feedback

    /// Measurement variable: stem radius [rw].
    let drwdt' bs n m r fMax allocationMode = 
        let biomassStemChange = dbsdt' bs n m r fMax allocationMode
        if biomassStemChange > 0.
            then
                let oldRadius = bs |> Proxies.toRadiusMM
                let newRadius = (bs + biomassStemChange) |> Proxies.toRadiusMM
                newRadius - oldRadius
            else 0.

    /// Bristlecone function for dBs/dt
    let dbsdt p _ bs (e:Environment) =
        dbsdt' bs ((e.[ShortCode.create "N"]) |> Proxies.d15NtoAvailability)
            (p |> Pool.getEstimate "gammab") (p |> Pool.getEstimate "r") (fMax p) (allocationMode p)

    /// Bristlecone function for dN/dt
    let dndt p _ n (e:Environment) =
        dndt' (e.[ShortCode.create "x"]) (n |> Proxies.d15NtoAvailability) (p |> Pool.getEstimate "lambda") (p |> Pool.getEstimate "gamman") 
            (p |> Pool.getEstimate "q") (p |> Pool.getEstimate "r") (fMax p) (allocationMode p) (nitrogenFeedback p)

    /// Bristlecone function for dr/dt
    let drwdt p _ _ (e:Environment) =
        drwdt' (e.[ShortCode.create "bs"]) ((e.[ShortCode.create "N"]) |> Proxies.d15NtoAvailability)
            (p |> Pool.getEstimate "gammab") (p |> Pool.getEstimate "r") (fMax p) (allocationMode p)

    { Equations  = [ ShortCode.create "x",      drwdt
                     ShortCode.create "bs",     dbsdt
                     ShortCode.create "N",      dndt ] |> Map.ofList
      Parameters = [ // for nitrogen replenishment
                     ShortCode.create "lambda", Parameter.create PositiveOnly   0.010 0.350   // Rate of nitrogen replenishment
                     ShortCode.create "gamman", Parameter.create PositiveOnly   0.001 0.400   // Loss rate of nitrogen
                     // for shrub allocation and physiology
                     ShortCode.create "r",      Parameter.create PositiveOnly   0.001 0.050   // Respiration cost per unit biomass
                     ShortCode.create "q",      Parameter.create PositiveOnly   0.100 0.100   // Nutrient requirement per unit biomass
                     ShortCode.create "gammab", Parameter.create PositiveOnly   0.001 0.002   // Loss rate of biomass
                     // for likelihood function
                     ShortCode.create "rho",    Parameter.create Unconstrained  -0.50 0.25    // Covariance between growth and nitrogen
                     ShortCode.create "sigmax", Parameter.create PositiveOnly    0.50 0.55    // Standard deviation of x (biomass)
                     ShortCode.create "sigmay", Parameter.create PositiveOnly    0.40 0.50    // Standard deviation of y (nitrogen)
                    ] |> List.append additionalParameters |> Map.ofList
      Likelihood = ModelLibrary.Likelihood.bivariateGaussian "x" "N" }

let combine4 xs ys zs bs = [
    for x in xs do
    for y in ys do
    for z in zs do
    for b in bs do
    yield x, y, z, b ]

let hypotheses =

    let allocationModes =
        [ (fun _ -> Allocation.none), []
          (fun _ -> Allocation.equal), []
          (fun p -> Allocation.floating (p |> Pool.getEstimate "al") (p |> Pool.getEstimate "ar")),
          [ ShortCode.create "al",     Parameter.create PositiveOnly   0.100 0.200 
            ShortCode.create "ar",     Parameter.create PositiveOnly   0.900 1.000 ] ]

    let limitationModes =
        [ (fun p -> ModelComponents.GrowthLimitation.michaelisMenten (p |> Pool.getEstimate "fMax") (p |> Pool.getEstimate "h")),
          [ ShortCode.create "fMax",   Parameter.create PositiveOnly   0.200 0.600   // Maximum rate of resource-saturated growth per unit biomass
            ShortCode.create "h",      Parameter.create PositiveOnly   0.001 0.400   // Soil resource availability for growth at half of maximum rate 
          ] ]

    let feedbackModes =
        [ (fun p -> FeedbackToSoil.none), []
          (fun p -> FeedbackToSoil.withBiomassLoss (p |> Pool.getEstimate "alpha") (p |> Pool.getEstimate "gammab") ),
          [ ShortCode.create "alpha",  Parameter.create PositiveOnly   0.001 0.005 ] ]

    let growthRateModes = 
        [ fun _ -> id ]

    combine4 allocationModes limitationModes feedbackModes growthRateModes
    |> List.map (fun ((a,ap),(l,lp),(f,fp),(g)) -> ``base model`` l a f (List.concat [lp; ap; fp]))



// 3. Test Engine and Model
// ----------------------------
// Running a full test is strongly recommended. The test will demonstrate if the current
// configuration can find known parameters for a model. If this step fails, there is an
// issue with either your model, or the Bristlecone configuration.

let startValues =
    [ ShortCode.create "x",  0.23
      ShortCode.create "N",  4.64
      ShortCode.create "bs", 4.64 |> Proxies.toBiomassMM ] 
    |> Map.ofList

hypotheses.[0] |> Bristlecone.testModel (Bristlecone.mkContinuous |> Bristlecone.withConditioning (Custom startValues)) Options.testSeriesLength startValues Options.iterations Options.burn



// 3. Load Real Data and Estimate
// ----------------------------

let shrubs = 
    let yuribei = DataAccess.Shrub.loadRingWidths (__SOURCE_DIRECTORY__ + "/data/yuribei-rw.csv")
    let d15N = DataAccess.Shrub.loadLocalEnvironmentVariable (__SOURCE_DIRECTORY__ + "/data/yuribei-d15N-imputed.csv")
    yuribei
    |> Seq.map (fun s -> s.Identifier.Value, s)
    |> Seq.keyMatch d15N
    |> Seq.map (fun (_,plant,d15N) -> PlantIndividual.zipEnv (ShortCode.create "N") plant d15N)
    |> Seq.toList

let getStartValues startDate (plant:PlantIndividual) =
    let initialRadius =
        match plant.Growth with
        | PlantIndividual.PlantGrowth.RingWidth s -> 
            match s with
            | GrowthSeries.Absolute c -> c.Head |> fst |> removeUnit
            | GrowthSeries.Cumulative c -> 6.3 //c.Head |> fst |> removeUnit
            | _ -> invalidOp "Not implemented"
        | _ -> invalidOp "Not implemented 2"
    let initialMass = initialRadius |> removeUnit |> Proxies.toBiomassMM
    let initialNitrogen = plant.Environment.[ShortCode.create "N"].Head |> fst
    [ ShortCode.create "x", initialRadius
      ShortCode.create "N", initialNitrogen 
      ShortCode.create "bs", initialMass ] |> Map.ofList

// Generate results for all models, chains, and plants
let estimates =
    [hypotheses.[0] ]
    |> List.toArray
    |> Array.Parallel.map(fun h ->
        [| 1 .. Options.chains |]
        |> Array.Parallel.map(fun _ ->
            [shrubs.[0]] |> List.map (fun s ->
                let shrub = s |> PlantIndividual.toCumulativeGrowth
                let common = shrub |> PlantIndividual.keepCommonYears
                let startDate = common.Environment.[ShortCode.create "N"] |> TimeSeries.start
                let startConditions = getStartValues startDate shrub
                s.Identifier, h, common |> Bristlecone.PlantIndividual.fit (Bristlecone.mkContinuous |> Bristlecone.withConditioning (Custom startConditions)) Options.iterations Options.burn h )))



// 4. Save Results to R data file
// ----------------------------
let burnin = 0

let everyNth n seq = 
    seq |> Seq.mapi (fun i el -> el, i)              // Add index to element
        |> Seq.filter (fun (el, i) -> i % n = n - 1) // Take every nth element
        |> Seq.map fst                               // Drop index from the result

// PlantCode, Hypothesis, Iteration, Chain, Parameter, Likelihood
let chainsDataFrame =
    estimates
    |> Array.mapi (fun hypothesisNumber hypothesis -> 
        hypothesis
        |> Array.mapi (fun chainNumber chain ->
            chain
            |> List.map (fun (plantCode,model,result) ->
                result.Trace
                |> List.unzip
                |> snd
                |> List.rev
                |> List.skip burnin // Skip burn in period
                |> everyNth 5 // Thinning by this amount
                |> Seq.mapi (fun iterationNumber values ->
                    result.Parameters
                    |> Map.toList
                    |> List.mapi(fun i (name,p) -> 
                        plantCode.Value,
                        hypothesisNumber,
                        iterationNumber + 1,
                        chainNumber + 1,
                        name.Value,
                        result.Likelihood,
                        values.[i] )
                )
                |> List.concat
                |> Seq.toList
            )
            |> List.concat )
        |> Array.toList
        |> List.concat )
    |> Array.toList
    |> List.concat

// Save as CSV file
open FSharp.Data
type BristleconeResult = CsvProvider<Sample = "PlantCode (string), Hypothesis (int), Iteration (int), Chain (int), Parameter (string), Likelihood (float), value (float)">
let buildRowFromObject = fun (a,b,c,d,e,f,g) -> BristleconeResult.Row(a,b,c,d,e,f,g)
let buildTableFromObjects = (Seq.map buildRowFromObject) >> Seq.toList >> BristleconeResult
let myCsv = chainsDataFrame |> buildTableFromObjects
myCsv.Save("dphil-shrub-output.csv")

// 5. Plot results using RProvider
// ----------------------------

#load "../packages/RProvider/RProvider.fsx"
//#load "plotting.fsx"
open RProvider

open RProvider.grDevices
open RProvider.graphics
R.x11()

// i. Plot of likelihoods for every chain, and every individual?
// One plot for each hypothesis?
module Plot =

    let (++) (plot1:RDotNet.SymbolicExpression) (plot2:RDotNet.SymbolicExpression) = 
        R.``+``(plot1, plot2)
    
    // let ggplotLikelihood df =
    //     R.ggplot(df, namedParams [ "x", box "value" ])
    //     ++ R.geom__line()
    //     ++ R.facet__wrap("PlantCode")

    // ggplotLikelihood df

    // Display plot of likelihood:
    // - For each hypothesis
    // - Every plant and chain


    // let likelihoodPlots =
    //     R.par(namedParams [ "mfrow", [(series |> List.length) / 3;(series |> List.length) / 3]] ) |> ignore
    //     series |> List.map(fun (s,h,c,k) -> R.plot( namedParams [ "x", box (k.Trace |> List.map fst |> List.rev |> List.skip 2500) ; "type", box "l"; "xlab", box "Iteration"; "ylab", box "-log likelihood" ]) |> ignore )

    // let predictionPlots =
    //     R.par(namedParams [ "mfrow", [(series |> List.length) / 2; 5]] ) |> ignore
    //     series |> List.map(fun (s,h,c,k) -> R.plot( namedParams [ "x", box (k.Series.[ShortCode.create "x"].Expected) ; "type", box "l"; "xlab", box "Iteration"; "ylab", box "Stem Radius (Model)" ]) |> ignore ) |> ignore
    //     series |> List.map(fun (s,h,c,k) -> R.plot( namedParams [ "x", box (k.Series.[ShortCode.create "N"].Expected) ; "type", box "l"; "xlab", box "Iteration"; "ylab", box "Nitrogen (Model)" ]) |> ignore )


// iii. For each parameter, find mean and standard deviation, and plot graph with histogram behind it



// // Test out new stuff
// open MathNet.Numerics

// let random = MathNet.Numerics.Random.MersenneTwister.Default

// let p =
//     [|    Parameter.create PositiveOnly   1000. 2000.   
//           Parameter.create PositiveOnly   1000. 2000.   
//           Parameter.create PositiveOnly   1000. 2000.   
//           Parameter.create PositiveOnly   1000. 2000.   
//           Parameter.create PositiveOnly   1000. 2000.   
//           Parameter.create Unconstrained  1000. 2000.   
//           Parameter.create PositiveOnly   1000. 2000.   
//           Parameter.create PositiveOnly   1000. 2000. |]

// let cov = Optimisation.MonteCarlo.TuningMode.defaultCovariance p.Length
// let sample = Optimisation.Distribution.MutlivariateNormal.sample cov random
// sample()
// let zeroMean : LinearAlgebra.Matrix<float> = (MathNet.Numerics.LinearAlgebra.DenseMatrix.identity p.Length) * 0.

// let x = Distributions.MatrixNormal(zeroMean, cov, cov, random)
// x.Sample()


// let cov2 =
//     [| 1 .. 10000 |]
//     |> Array.map(fun _ ->
//         p
//         |> Array.map Parameter.bounds
//         |> Array.map (fun (low,high) -> Optimisation.MonteCarlo.draw random 0. ((high - low) / 4.) () ) )
//         |> Optimisation.MonteCarlo.TuningMode.samplesToMatrix
//         |> Optimisation.MonteCarlo.TuningMode.computeCovariance

// let y = Distributions.MatrixNormal(zeroMean, cov2, cov2, random)
// y.Sample().Diagonal().ToArray()

// let sample2 = Optimisation.Distribution.MutlivariateNormal.sample cov2 random
// sample2().ToArray()