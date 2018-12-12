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
open FSharp.Data

type BristleconeResult = CsvProvider<Sample = "PlantCode (string), Hypothesis (int), Iteration (int), Chain (int), Parameter (string), Likelihood (float), value (float)", CacheRows = true>

// 0. Configure Options
// ----------------------------

module Options =
    let resultsDirectory = "/Volumes/Macbook Time Machine/Shrub Results 250k/paper1-no-n-limitation-models/"
    let iterations = 10
    let chains = 3
    let engine =
        Bristlecone.mkContinuous 
        |> Bristlecone.withContinuousTime Integration.MathNet.integrate
        //|> Bristlecone.withContinuousTime (Integration.MsftOslo.integrateWithErrorHandling (Integration.MsftOslo.Options.custom 1e-03 1e-03))
        |> Bristlecone.withTunedMCMC [ Optimisation.MonteCarlo.TuneMethod.Scale, 2000, 10000
                                       Optimisation.MonteCarlo.TuneMethod.CovarianceWithScale 0.250, 500, 100000 ]


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
                Optimisation.RootFinding.bisect 0 200 f 0.01 1000.00 1e-8 // Assumption that shrub radius is between 0.01 and 100.0cm.
            mass
            |> massToVolume woodDensity
            |> findRadius

    module GrowthLimitation =

        /// A rearranged version of a Monod model
        let hollingDiscModel a b h =
            Some <| fun r -> (a * r) / (1. + (a * b * h * r))

        /// TEST: An integrated supply and use model
        let saturatingSupplySaturatingGrowth r a b h rootMass =
            Some <| fun resource -> (a * r * rootMass * resource) / (1. + a * b * r * rootMass * resource + a * h * resource)

        /// Monod model once saturation has been reached
        let linear a =
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
            radiusMM / 10. |> Allometrics.shrubBiomass Constants.b Constants.a Constants.rtip Constants.p Constants.lmin Constants.k5 Constants.k6 Constants.numberOfStems Constants.salixWoodDensity

        /// Biomass in grams
        let toRadiusMM biomassGrams = 
            let radiusCm = biomassGrams |> Allometrics.shrubRadius Constants.b Constants.a Constants.rtip Constants.p Constants.lmin Constants.k5 Constants.k6 Constants.numberOfStems Constants.salixWoodDensity
            radiusCm * 10.

        /// d15N to N availability. From Craine 2009, as shown in Craine 2015 (Plant and Soil).
        /// Assuming d15N is a linear index of N availability, the minimum supported value of d15N is -3.09, as 0 N availability.
        let d15NtoAvailability d15N =
            (100. * d15N + 309.) / 359.


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


    module BaseEquations =

        /// Cumulative stem biomass [dBs/dt]
        let biomass b n r gammab geom f : float =
            b * r * (f n) * geom(b) - gammab * b

        let soilNitrogen n b gamman y geom f feedback : float =
            y - (geom(b) * b * (f n)) - gamman * n + feedback(b)

        let soilNitrogenNoUptake n b gamman y feedback : float =
            y - gamman * n + feedback(b)


// 2. Create Hypotheses
// ----------------------------

let ``base model`` maxGrowthRate nLimitation nitrogenFeedback additionalParameters =

    let randomLog environment parameters =
        if System.Random().Next(0,500000) = 1 then
            printfn "RANDOM LOG: Environment = %A; Parameters = %A" environment parameters

    /// Cumulative stem biomass [b].
    let dbsdt' (b:float) n gammab r maxGrowthRate limit =
        match limit with
        | Some l -> ModelComponents.BaseEquations.biomass b n r gammab maxGrowthRate l
        | None -> ModelComponents.BaseEquations.biomass b n r gammab maxGrowthRate (fun _ -> 1.)

    /// Bioavailable soil nitrogen [N]
    let dndt' bs n lambda gamman maxGrowthRate feedback limit = 
        match limit with
        | Some l -> ModelComponents.BaseEquations.soilNitrogen n bs gamman lambda maxGrowthRate l feedback
        | None -> ModelComponents.BaseEquations.soilNitrogenNoUptake n bs gamman lambda feedback

    /// Measurement variable: stem radius [rw].
    let drwdt' bs n r gammab maxGrowthRate limit =
        let biomassStemChange = dbsdt' bs n gammab r maxGrowthRate limit
        if biomassStemChange > 0.
            then
                let oldRadius = bs |> ModelComponents.Proxies.toRadiusMM
                let newRadius = (bs + biomassStemChange) |> ModelComponents.Proxies.toRadiusMM
                newRadius - oldRadius
            else 0.

    /// Bristlecone function for dBs/dt
    let dbsdt p _ bs (e:Environment) =
        randomLog e p
        dbsdt' bs ((e.[ShortCode.create "N"]) |> ModelComponents.Proxies.d15NtoAvailability)
            (p |> Pool.getEstimate "gamma[b]") (p |> Pool.getEstimate "r") (maxGrowthRate p) (nLimitation p)

    /// Bristlecone function for dN/dt
    let dndt p _ n (e:Environment) =
        dndt' (e.[ShortCode.create "bs"]) (n |> ModelComponents.Proxies.d15NtoAvailability) 
            (p |> Pool.getEstimate "lambda") (p |> Pool.getEstimate "gamma[n]") (maxGrowthRate p) (nitrogenFeedback p) (nLimitation p)

    /// Bristlecone function for dr/dt
    let drwdt p _ _ (e:Environment) =
        drwdt' (e.[ShortCode.create "bs"]) ((e.[ShortCode.create "N"]) |> ModelComponents.Proxies.d15NtoAvailability) 
            (p |> Pool.getEstimate "r") (p |> Pool.getEstimate "gamma[b]") (maxGrowthRate p) (nLimitation p)

    { Equations  = [ code "x",         drwdt
                     code "bs",        dbsdt
                     code "N",         dndt ] |> Map.ofList
      Parameters = [ // for nitrogen dynamics
                     code "lambda",    parameter PositiveOnly   0.001 0.500   // Rate of nitrogen replenishment
                     code "gamma[n]",  parameter PositiveOnly   0.001 0.200   // Loss rate of nitrogen
                     // for shrub physiology
                     code "gamma[b]",  parameter PositiveOnly   0.001 0.200   // Loss rate of biomass
                     //code "r",         parameter PositiveOnly   0.500 1.000   // Photosynthetic efficiency (either N-limited or N-unlimited)
                     // for likelihood function
                     code "rho",       parameter Unconstrained  -0.50 0.500   // Covariance between growth and nitrogen
                     code "sigma[x]",  parameter PositiveOnly   0.100 1.200   // Standard deviation of x (biomass)
                     code "sigma[y]",  parameter PositiveOnly   0.250 0.750   // Standard deviation of y (nitrogen)
                    ] |> List.append additionalParameters |> Map.ofList
      Likelihood = ModelLibrary.Likelihood.bivariateGaussian "x" "N" }

let hypotheses =

    // [A] N may limited growth via combined N-limitations on (a) photosynthetic and (b) uptake rates
    let limitationModes =
        [ (fun p -> ModelComponents.GrowthLimitation.hollingDiscModel (p |> Pool.getEstimate "a") (p |> Pool.getEstimate "r") (p |> Pool.getEstimate "h")),
           [ code "a",      parameter PositiveOnly   0.001 0.010          // N-uptake efficiency
             code "h",      parameter PositiveOnly   0.001 0.010
             code "r",      parameter PositiveOnly   100.0 500.0 ]      // N-handling time (including uptake and incorporation)
          (fun p -> ModelComponents.GrowthLimitation.linear (p |> Pool.getEstimate "a")), 
           [ code "a",      parameter PositiveOnly   0.001 0.010
             code "r",      parameter PositiveOnly   100.0 500.0 ]        // N-uptake efficiency
          (fun _ -> ModelComponents.GrowthLimitation.none), 
          [ code "r",         parameter PositiveOnly   0.001 1.000 ] ]

    // [B] Loss of plant material may feedback into the soil pool of available nitrogen (instant)
    let feedbackModes =
        [ (fun _ -> ModelComponents.FeedbackToSoil.none), []
          (fun p -> ModelComponents.FeedbackToSoil.withBiomassLoss (p |> Pool.getEstimate "alpha") (p |> Pool.getEstimate "gamma[b]") ),
          [ ShortCode.create "alpha",  Parameter.create PositiveOnly   0.0001 0.0010 ] ]    // N-recycling efficiency

    // [C] A plant may be subject to mechanical constraints on its maximum size
    let geometricModes = 
        [  (fun p -> ModelComponents.GeometricConstraint.none), []
           (fun p -> ModelComponents.GeometricConstraint.chapmanRichards (p |> Pool.getEstimate "k")),
           [ ShortCode.create "k",  Parameter.create PositiveOnly   3000.00 5000.00 ] ]     // Asymptotic biomass (grams)

    List.combine3 geometricModes limitationModes feedbackModes
    |> List.map (fun ((growth,gp),(limit,lp),(feedback,fp)) -> 
        ``base model`` growth limit feedback (List.concat [lp; fp; gp]))


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

let getStartValues (startDate:System.DateTime) (plant:PlantIndividual) =
    let initialRadius =
        match plant.Growth with
        | PlantIndividual.PlantGrowth.RingWidth s -> 
            match s with
            | GrowthSeries.Absolute c -> c.Head |> fst |> removeUnit
            | GrowthSeries.Cumulative c -> 
                let start = (c |> TimeSeries.trimStart (startDate - System.TimeSpan.FromDays(366.))).Values |> Array.head |> removeUnit
                printfn "Start cumulative growth = %f" start
                start
            | GrowthSeries.Relative _ -> invalidOp "Not implemented"
        | _ -> invalidOp "Not implemented 2"
    let initialMass = initialRadius |> removeUnit |> ModelComponents.Proxies.toBiomassMM
    let initialNitrogen = plant.Environment.[ShortCode.create "N"].Head |> fst
    [ ShortCode.create "x", initialRadius
      ShortCode.create "N", initialNitrogen 
      ShortCode.create "bs", initialMass ] |> Map.ofList

let everyNth n seq = 
    seq |> Seq.mapi (fun i el -> el, i)              // Add index to element
        |> Seq.filter (fun (el, i) -> i % n = n - 1) // Take every nth element
        |> Seq.map fst                               // Drop index from the result

let workPackages shrubs (hypotheses:ModelSystem list) engine saveDirectory =
    seq {
        for s in shrubs do
            for hi in [ 1.. hypotheses.Length ] do
                    let estimate() =
                        printfn "Starting estimate for shrub %s (H%i)" s.Identifier.Value hi
                        let estimates =
                            [hypotheses.[hi-1]]
                            |> List.toArray
                            |> Array.map(fun h ->
                                [| 1 .. Options.chains |]
                                |> Array.Parallel.map(fun _ ->
                                    [s] |> List.map (fun s ->
                                        let shrub = s |> PlantIndividual.toCumulativeGrowth
                                        let common = shrub |> PlantIndividual.keepCommonYears
                                        let startDate = common.Environment.[ShortCode.create "N"] |> TimeSeries.start
                                        let startConditions = getStartValues startDate shrub
                                        try 
                                            let e = engine |> Bristlecone.withConditioning (Custom startConditions)
                                            let result = common |> Bristlecone.PlantIndividual.fit e Options.iterations h
                                            printfn "Result %A" result
                                            Some (s.Identifier, h, result)
                                        with
                                        | e ->
                                            printfn "A chain failed with expection: %A" e
                                            None
                                    )))
                        printfn "Done with estimates for %s. Saving..." s.Identifier.Value

                        // PlantCode, Hypothesis, Iteration, Chain, Parameter, Likelihood
                        let chainsDataFrame =
                            estimates
                            |> Array.map (fun hypothesis -> 
                                hypothesis
                                |> Array.mapi (fun chainNumber chain ->
                                    chain
                                    |> List.choose id
                                    |> List.map (fun (plantCode,model,result) ->
                                        result.Trace
                                        |> List.rev
                                        |> Seq.mapi (fun iterationNumber (likelihood,values) ->
                                            result.Parameters
                                            |> Map.toList
                                            |> List.mapi(fun i (name,p) -> 
                                                plantCode.Value,
                                                hi,
                                                iterationNumber + 1,
                                                chainNumber + 1,
                                                name.Value,
                                                likelihood,
                                                values.[i] )
                                        )
                                        |> everyNth 5 // Thinning by this amount
                                        |> List.concat
                                        |> Seq.toList
                                    )
                                    |> List.concat )
                                |> Array.toList
                                |> List.concat )
                            |> Array.toList
                            |> List.concat

                        // Save as CSV file
                        let buildRowFromObject = fun (a,b,c,d,e,f,g) -> BristleconeResult.Row(a,b,c,d,e,f,g)
                        let buildTableFromObjects = (Seq.map buildRowFromObject) >> Seq.toList >> BristleconeResult
                        let myCsv = chainsDataFrame |> buildTableFromObjects
                        myCsv.Save(sprintf "%sdphil-shrub-output-%s-H%i.csv" saveDirectory s.Identifier.Value hi)
                    yield estimate
            }

(shrubs |> List.map (fun s -> printfn "%s" s.Identifier.Value))

workPackages (shrubs |> List.skip 1) hypotheses Options.engine Options.resultsDirectory
|> Seq.chunkBySize 3
|> Seq.toArray
|> Array.map (fun f -> f |> Array.Parallel.map(fun g -> g()))




/////////////////////////////
/// PLOTTING BELOW
/// //////////////////////////

// A. Load in pre-computed results
// ________________________________

let data =
    System.IO.Directory.GetFiles(Options.resultsDirectory, "*.csv")

let bestParameterPool (r: BristleconeResult) =
    let lowestLikelihood = (r.Rows |> Seq.minBy (fun x -> x.Likelihood)).Likelihood
    printfn "Lowest likelihood: %f" lowestLikelihood
    r.Rows
    |> Seq.sortBy(fun x -> x.Likelihood)
    |> Seq.takeWhile(fun x -> x.Likelihood = lowestLikelihood)
    |> Seq.map(fun row -> ShortCode.create row.Parameter, Parameter.create Unconstrained row.Value row.Value)
    |> Map.ofSeq

let hypothesisNumberFromFilename (name:string) =
    name.Split('H').[1].Split('.').[0] |> int

let fit s h =
    let shrub = s |> PlantIndividual.toCumulativeGrowth
    let common = shrub |> PlantIndividual.keepCommonYears
    let startDate = common.Environment.[ShortCode.create "N"] |> TimeSeries.start
    let startConditions = getStartValues startDate shrub
    let e = Bristlecone.mkContinuous |> Bristlecone.withContinuousTime Integration.MathNet.integrate |> Bristlecone.withConditioning (Custom startConditions)
    common |> Bristlecone.PlantIndividual.fit e 1 h

let fitData =
    data
    |> Array.rev
    |> Array.choose(fun fileName ->
        let data = BristleconeResult.Load fileName
        let headRow = data.Rows |> Seq.head
        let hypothesisNumber = hypothesisNumberFromFilename fileName
        let realShrub = shrubs |> List.find (fun s -> s.Identifier.Value = headRow.PlantCode)
        let hypothesis =  { hypotheses.[hypothesisNumber-1] with Parameters = bestParameterPool data}
        printfn "Hypothesis = %A" hypothesis
        try
            let a = fit realShrub hypothesis
            (headRow.PlantCode, hypothesisNumber, a) |> Some
        with
        | _ -> None )
    |> Seq.groupBy(fun (plant,_,_) -> plant)
    |> Seq.collect(fun (plant,data) -> 
        let weights = data |> Seq.map (fun (x,y,z) -> z) |> ModelSelection.Akaike.akaikeWeights
        Seq.zip weights data )
    |> Seq.collect(fun ((e,w),(shrubId, h, estimate)) ->
        estimate.Series 
        |> Map.toArray
        |> Array.collect (fun (s,v) ->
            v.Expected |> Array.mapi (fun i x -> shrubId, h, s.Value, i, x, estimate.Likelihood, w ) ) )
    |> Seq.toArray

type Prediction = CsvProvider<Sample = "PlantCode (string), Hypothesis (int), Variable (string), Time (int), Expected (float), Likelihood (float), AkaikeWeight (float)">

// Save as CSV file
let buildRowFromObject = fun (a,b,c,d,e,f,g) -> Prediction.Row(a,b,c,d,e,f,g)
let buildTableFromObjects = (Seq.map buildRowFromObject) >> Seq.toList >> Prediction
let myCsv = fitData |> buildTableFromObjects
myCsv.Save(sprintf "%sdphil-shrub-predictions.csv" "/Users/andrewmartin/Desktop/")
