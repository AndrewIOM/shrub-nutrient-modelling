#r "../packages/NETStandard.Library.NETFramework/build/net461/lib/netstandard.dll"
#load "../packages/Bristlecone/bristlecone.fsx"
// #load "../packages/Bristlecone/charts.fsx"
#load "components/components.fsx"
#r "../packages/Bristlecone.Dendro/lib/netstandard2.0/Bristlecone.Dendro.dll"
#load "components/temperature.fsx"

////////////////////////////////////////////////////
/// Long-Term N-Shrub relations in Yamal, Russia
////////////////////////////////////////////////////

// This script completes model-fitting and model-selection
// for 24 shrub individuals. Models for nitrogen limitation,
// temperature limitation, and size-related asymptotes are
// included.

open Bristlecone
open Bristlecone.ModelSystem
open Bristlecone.Dendro
open Bristlecone.Dendro.PlantIndividual
open Bristlecone.Workflow.Orchestration

// Temporary logger

module CustomLog =

    open System.Threading
    open Bristlecone.Logging

    let print threadId (x:LogEvent) = 
        match x with
        | Bristlecone.Logging.OptimisationEvent e ->
            if e.Iteration % 5000 = 0
            then printfn "##%i## At iteration %i (-logL = %f) %A" threadId e.Iteration e.Likelihood e.Theta
        | _ -> printfn "##%i## %A" threadId x

    let logger () =

        let agent = MailboxProcessor.Start(fun inbox -> 
            let rec messageLoop () = async {
                let! threadId,msg = inbox.Receive()
                print threadId msg
                return! messageLoop ()
                }
            messageLoop () )

        agent.Error.Add(fun e -> printfn "Error = %A" e)

        fun msg -> 
            let threadId = Thread.CurrentThread.ManagedThreadId
            agent.Post (threadId, msg)


// 1. Configure Options
// ----------------------------

module Options =
    let resultsDirectory = "/Users/andrewmartin/Desktop/Bristlecone Results/Paper3-Repeated-AM-CleanMCMC/"
    let chains = 1
    let endWhen = Optimisation.EndConditions.afterIteration 100000
    let logger = CustomLog.logger() //Logging.RealTimeTrace.graphWithConsole 60. 10000
    let filzbachOptions : Optimisation.MonteCarlo.Filzbach.FilzbachSettings<float> = {
        TuneAfterChanges = 20
        MaxScaleChange = 100.00
        MinScaleChange = 0.0010
        BurnLength = Optimisation.EndConditions.afterIteration 100000 }
    let engine =
        Bristlecone.mkContinuous 
        |> Bristlecone.withContinuousTime Integration.MathNet.integrate
        |> Bristlecone.withOutput logger
        |> Bristlecone.withTunedMCMC [ Optimisation.MonteCarlo.TuneMethod.CovarianceWithScale 0.100, 1000, Optimisation.EndConditions.afterIteration 100000 ]
        // |> Bristlecone.withCustomOptimisation (Optimisation.MonteCarlo.Filzbach.filzbach filzbachOptions)
        // |> Bristlecone.withCustomOptimisation (Optimisation.MonteCarlo.adaptiveMetropolis 0.250 500)
        // |> Bristlecone.withCustomOptimisation (Optimisation.MonteCarlo.SimulatedAnnealing.fastSimulatedAnnealing 0.01 false
        //     { Optimisation.MonteCarlo.SimulatedAnnealing.AnnealSettings<float>.Default with 
        //         BoilingAcceptanceRate = 0.85
        //         HeatRamp = (fun t -> t + sqrt t); TemperatureCeiling = Some 100.
        //         HeatStepLength = Optimisation.EndConditions.afterIteration 1000
        //         AnnealStepLength = (fun x -> Optimisation.MonteCarlo.SimulatedAnnealing.EndConditions.improvementCount 5000 250 x || Optimisation.EndConditions.afterIteration 10000 x) }) //(Optimisation.EndConditions.afterIteration 10000) })

    let orchestrator = OrchestrationAgent(logger, System.Environment.ProcessorCount, false)//System.Environment.ProcessorCount, false)


// 2. Create Hypotheses
// ----------------------------

module BaseEquations =

    /// Cumulative stem biomass [dBs/dt]
    let biomass b n r gammab geom f tempEffect sigmaz : float =
        if r > 100000. then nan
        else b * r * (f n) * geom(b) * tempEffect - gammab * b //+ sigmaz

    let soilNitrogen n b gamman y geom f feedback tempEffect : float =
        y - (geom(b) * b * (f n) * tempEffect) - gamman * n + feedback(b)

    let soilNitrogenNoUptake n b gamman y feedback : float =
        y - gamman * n + feedback(b)


let ``base model`` maxGrowthRate nLimitation nitrogenFeedback tempEffect additionalParameters =

    /// Cumulative stem biomass [b].
    let dbsdt' (b:float) n gammab r maxGrowthRate limit tLimit = //sigmaz =
        match limit with
        | Some l -> BaseEquations.biomass b n (r * 1000.) gammab maxGrowthRate l tLimit 0.//sigmaz
        | None -> BaseEquations.biomass b n r gammab maxGrowthRate (fun _ -> 1.) tLimit 0.//sigmaz

    /// Bioavailable soil nitrogen [N]
    let dndt' bs n lambda gamman maxGrowthRate feedback limit tLimit = 
        match limit with
        | Some l -> BaseEquations.soilNitrogen n bs gamman lambda maxGrowthRate l feedback tLimit
        | None -> BaseEquations.soilNitrogenNoUptake n bs gamman lambda feedback

    /// Bristlecone function for dBs/dt
    let dbsdt p t bs (e:Environment) =
        //printfn "T at %f is %f (tEffect = %f)" t (lookup e "T[max]" - 273.15) (tempEffect p e)
        dbsdt' bs ((e.[ShortCode.create "N"]) |> ModelComponents.Proxies.d15NtoAvailability)
            (p |> Pool.getEstimate "gamma[b]") ((p |> Pool.getEstimate "r")) (maxGrowthRate p) (nLimitation p) (tempEffect p e) //(p |> Pool.getEstimate "sigma[z]")

    /// Bristlecone function for dN/dt
    let dndt p _ n (e:Environment) =
        dndt' (e.[ShortCode.create "bs"]) (n |> ModelComponents.Proxies.d15NtoAvailability) 
            (p |> Pool.getEstimate "lambda") (p |> Pool.getEstimate "gamma[n]") (maxGrowthRate p) (nitrogenFeedback p) (nLimitation p) (tempEffect p e)

    /// Measurement (Size) variable: stem radius
    let stemRadius lastRadius lastEnv env =
        let oldCumulativeMass = lookup lastEnv "bs"
        let newCumulativeMass = lookup env "bs"
        if (newCumulativeMass - oldCumulativeMass) > 0.
        then 
            // This clause mandates that stem radius can never decrease.
            let newRadius = newCumulativeMass |> ModelComponents.Proxies.toRadiusMM
            if newRadius > lastRadius then newRadius else lastRadius
        else lastRadius

    { Equations  = [ code "bs",        dbsdt
                     code "N",         dndt ] |> Map.ofList
      Measures   = [ code "x",         stemRadius ] |> Map.ofList
      Parameters = [ // for nitrogen dynamics
                     code "lambda",    parameter PositiveOnly   0.100 2.000   // Rate of nitrogen replenishment
                     code "gamma[n]",  parameter PositiveOnly   0.001 5.000   // Loss rate of nitrogen
                     // for shrub physiology
                     code "gamma[b]",  parameter PositiveOnly   0.001 0.200   // Loss rate of biomass
                     // for likelihood function
                     code "rho",       parameter Unconstrained  -0.50 0.500   // Covariance between growth and nitrogen
                     code "sigma[x]",  parameter PositiveOnly   0.100 1.200   // Standard deviation of x (biomass)
                     code "sigma[y]",  parameter PositiveOnly   0.250 0.750   // Standard deviation of y (nitrogen)
                     //code "sigma[z]",  parameter Unconstrained  0.250 0.750   // Inter-shrub comparison term
                    ] |> List.append additionalParameters |> Map.ofList
      Likelihood = ModelLibrary.Likelihood.bivariateGaussian "x" "N" }


let hypotheses =

    // [A] N may limited growth via combined N-limitations on (a) photosynthetic and (b) uptake rates
    let limitationModes =
        [ (fun p -> ModelComponents.GrowthLimitation.hollingDiscModelDual ((p |> Pool.getEstimate "a") / 1000.) ((p |> Pool.getEstimate "r") * 1000.) ((p |> Pool.getEstimate "h") / 1.) 5.00),
           [ code "a",      parameter PositiveOnly   0.100 5.000          // N-uptake efficiency
             code "h",      parameter PositiveOnly   0.001 0.250
             code "r",      parameter PositiveOnly   1.000 10.00 ]      // N-handling time (including uptake and incorporation)
          (fun p -> ModelComponents.GrowthLimitation.linear ((p |> Pool.getEstimate "a") / 1000.)), 
           [ code "a",      parameter PositiveOnly   0.001 0.010
             code "r",      parameter PositiveOnly   1.000 10.00 ]        // N-uptake efficiency
          (fun _ -> ModelComponents.GrowthLimitation.none), 
           [ code "r",       parameter PositiveOnly  1.000 10.00 ] ]

    // [B] Loss of plant material may feedback into the soil pool of available nitrogen (instant)
    let feedbackModes =
        [ (fun _ -> ModelComponents.FeedbackToSoil.none), []
          (fun p -> ModelComponents.FeedbackToSoil.withBiomassLoss ((p |> Pool.getEstimate "alpha") / 100.) (p |> Pool.getEstimate "gamma[b]") ),
          [ code "alpha",  parameter PositiveOnly   0.001 0.01 ] ]    // N-recycling efficiency

    // [C] A plant may be subject to mechanical constraints on its maximum size
    let geometricModes = 
        [  (fun p -> ModelComponents.GeometricConstraint.none), []
           (fun p -> ModelComponents.GeometricConstraint.chapmanRichards ((p |> Pool.getEstimate "k") * 1000.)),
           [ code "k",  parameter PositiveOnly   3.00 5.00 ] ]     // Asymptotic biomass (grams)

    /// The universal gas constant in J mol−1 K−1
    let gasConstant = 8.314

    /// An Arrhenius function to represent temperature limitation on growth.
    /// Form of equation from paper: https://pubag.nal.usda.gov/download/13565/PDF
    let temperatureLimitation preExp activationEnergy temperature =
        //System.Math.E ** (- ((activationEnergy * 1000.) / (8.314 * temperature)))
        preExp * System.Math.E ** ((1000. * activationEnergy * (temperature - 298.)) / (298. * gasConstant * temperature))

    // [D] Net photosynthetic rate is temperature-dependent
    let temperature =
        [ (fun _ -> ModelComponents.Temperature.none), []
          (fun p e -> temperatureLimitation 1. (p |> Pool.getEstimate "Ea") (lookup e "T[max]")),
           [ code "Ea",             parameter PositiveOnly 0.10 0.50  ] ]

    List.combine4 temperature geometricModes limitationModes feedbackModes
    |> List.map (fun ((temperature,tp),(growth,gp),(limit,lp),(feedback,fp)) -> 
        ``base model`` growth limit feedback temperature (List.concat [lp; fp; gp; tp]))


// 3. Load Real Data and Estimate
// ----------------------------

open FSharp.Data

// 1. Load in ring widths
let rw = Bristlecone.Data.PlantIndividual.loadRingWidths (__SOURCE_DIRECTORY__ + "/../data/yamal-rw.csv")

// 2. Load in temperature data
type TemperatureData = CsvProvider<"/Users/andrewmartin/Projects/GitHub-Projects/shrub-nutrient-modelling/data/yamal-mean-temperatures.csv">

let summerTemperature =
    let maxTemperatures = TemperatureData.Load "/Users/andrewmartin/Projects/GitHub-Projects/shrub-nutrient-modelling/data/yamal-mean-temperatures.csv"
    maxTemperatures.Rows
    |> Seq.groupBy(fun r -> r.Station)
    |> Seq.map(fun (station,r) -> station, r |> Seq.map(fun r -> ( float r.JJA + 273.15, System.DateTime.Create 31 12 r.Year)) |> TimeSeries.fromObservations)
    |> Seq.toList

(summerTemperature.Tail.Head |> snd).Resolution

(summerTemperature.Head |> snd) |> TimeSeries.toObservations |> Seq.toList |> List.map(fun (v,d) -> printfn "%i" d.Year)

// 2. Define common timeline for three year bins
// Read in yuribei annual data and regional 3-yearly data
type Regional15N = CsvProvider<"/Users/andrewmartin/Projects/GitHub-Projects/shrub-nutrient-modelling/data/yamal-d15N-lowres.csv">
type YuribeiD15N = CsvProvider<"/Users/andrewmartin/Projects/GitHub-Projects/shrub-nutrient-modelling/data/yuribei-d15N-imputed.csv">

let bins = seq { 1900 .. 3 .. 2018 }
let regionalD15N = Regional15N.Load (__SOURCE_DIRECTORY__ + "/../data/yamal-d15N-lowres.csv")
let yuribeiD15N = YuribeiD15N.Load (__SOURCE_DIRECTORY__ + "/../data/yuribei-d15N-imputed.csv")

// 3. Lower resolution of data where it is higher resolution than bins.
// - Determine the period that each reading represents
// - Where time-points are repeated, find the average
// - Calculate the weighted average of each period (mixture)

/// Converts input data into 3-year binned d15N and increment data
let lowerResolution (plant:PlantIndividual) =

    printfn "Plant is %A" plant

    let growth =
        match plant.Growth with
        | PlantIndividual.PlantGrowth.RingWidth x ->
            match x with
            | Absolute rw -> rw

    let isotopeByBin = 
        bins
        |> Seq.choose(fun binStart ->
            let lowResBins =
                regionalD15N.Rows
                |> Seq.where(fun n -> n.Shrub = plant.Identifier.Value)
                |> Seq.where(fun n -> n.BinLatest <= (binStart + 2) && n.BinOldest >= binStart)
                |> Seq.map(fun r -> (r.BinOldest, r.BinLatest, r.D15N))
                |> Seq.toList

            // These are complete (pre-interpolated) TODO do this in code instead
            let highResBins =
                yuribeiD15N.Rows
                |> Seq.where(fun n -> n.``Plant Code`` = plant.Identifier.Value)
                |> Seq.where(fun n -> n.Date.Year <= (binStart + 2) && n.Date.Year >= binStart)
                |> Seq.map(fun r -> (r.Date.Year, r.Date.Year, (float r.``Predictor 2``) ))
                |> Seq.toList

            let allData = 
                highResBins 
                |> Seq.append lowResBins 
                |> Seq.sortBy (fun (s,_,_) -> s)
                |> Seq.filter(fun (_,_,x) -> not <| System.Double.IsNaN x)

            printfn "All data is %A" allData
            printfn "Low data is %A" lowResBins
            printfn "High data is %A" highResBins

            if allData |> Seq.isEmpty then None
            else
                // Weighted average is:
                // Get the ring increment for the whole bin and times isotope value
                // Average by these values
                try
                    let weightedAverage =
                        let weights = 
                            let increments = allData |> Seq.map (fun (start,en,_) -> [ start .. en ] |> List.sumBy(fun y -> growth |> TimeSeries.findExact (System.DateTime(y, 12, 31)) |> fst |> removeUnit ))
                            let totalIncrement = increments |> Seq.sum
                            increments |> Seq.map(fun i -> i / totalIncrement)
                        printfn "Weights are %A" weights
                        allData
                        |> Seq.zip weights
                        |> Seq.averageBy(fun (weight,(s,e,d15n)) -> d15n * weight )
                    Some (binStart + 2, weightedAverage)
                with | _ -> None )
        |> Seq.toList
    isotopeByBin



/// Converts input data into 3-year binned d15N and increment data
let threeYearToAnnual (plant:PlantIndividual) =

    printfn "Plant is %A" plant

    let growth =
        match plant.Growth with
        | PlantIndividual.PlantGrowth.RingWidth x ->
            match x with
            | Absolute rw -> rw

    let (isotopeByBin:(int*float) list list) = 
        bins
        |> Seq.choose(fun binStart ->
            let lowResBins =
                regionalD15N.Rows
                |> Seq.where(fun n -> n.Shrub = plant.Identifier.Value)
                |> Seq.where(fun n -> n.BinLatest <= (binStart + 2) && n.BinOldest >= binStart)
                |> Seq.map(fun r -> (r.BinOldest, r.BinLatest, r.D15N))

            let highResBins =
                yuribeiD15N.Rows
                |> Seq.where(fun n -> n.``Plant Code`` = plant.Identifier.Value)
                |> Seq.where(fun n -> n.Date.Year <= (binStart + 2) && n.Date.Year >= binStart)
                |> Seq.map(fun r -> (r.Date.Year, r.Date.Year, (float r.``Predictor 2``) ))

            let allData = 
                highResBins 
                |> Seq.append lowResBins 
                |> Seq.sortBy (fun (s,_,_) -> s)
                |> Seq.filter(fun (_,_,x) -> not <| System.Double.IsNaN x)

            if allData |> Seq.isEmpty then None
            else
                match highResBins |> Seq.length with
                | i when i >= 0 && i < 3 ->
                    try
                        let weightedAverage =
                            let weights = 
                                let increments = allData |> Seq.map (fun (start,en,_) -> [ start .. en ] |> List.sumBy(fun y -> growth |> TimeSeries.findExact (System.DateTime(y, 12, 31)) |> fst |> removeUnit ))
                                let totalIncrement = increments |> Seq.sum
                                increments |> Seq.map(fun i -> i / totalIncrement)
                            printfn "Weights are %A" weights
                            allData
                            |> Seq.zip weights
                            |> Seq.averageBy(fun (weight,(s,e,d15n)) -> d15n * weight )
                        Some <| ([ binStart .. binStart + 2 ] |> List.map (fun y -> y, weightedAverage))
                    with | _ -> None
                | _ ->
                    highResBins 
                    |> Seq.map(fun (x,y,z) -> (x, z) ) 
                    |> Seq.toList
                    |> Some )
        |> Seq.toList
    List.concat isotopeByBin


let shrubs =
    rw
    |> List.filter(fun s -> regionalD15N.Rows |> Seq.exists(fun x -> x.Shrub = s.Identifier.Value ) || yuribeiD15N.Rows |> Seq.exists(fun x -> x.``Plant Code`` = s.Identifier.Value ))
    |> List.map(fun plant ->
        let isotope3Year = 
            threeYearToAnnual plant
            //lowerResolution plant
            |> List.map (fun (x,y) -> (y, System.DateTime(x, 12, 31)))
            |> TimeSeries.fromObservations
        plant |> PlantIndividual.zipEnv (code "N") isotope3Year )

shrubs
|> List.map(fun x -> printfn "[%s] Res = %A" x.Identifier.Value x.Environment.[code "N"].Resolution)



let getStartValues (startDate:System.DateTime) (plant:PlantIndividual) =
    let initialRadius =
        match plant.Growth with
        | PlantIndividual.PlantGrowth.RingWidth s -> 
            match s with
            | GrowthSeries.Absolute c -> c.Head |> fst |> removeUnit
            | GrowthSeries.Cumulative c -> 
                let start = (c |> TimeSeries.trimStart (startDate - System.TimeSpan.FromDays(366.))).Values |> Seq.head |> removeUnit
                printfn "Start cumulative growth = %f" start
                start
            | GrowthSeries.Relative _ -> invalidOp "Not implemented"
        | _ -> invalidOp "Not implemented 2"
    let initialMass = initialRadius |> removeUnit |> ModelComponents.Proxies.toBiomassMM
    let initialNitrogen = plant.Environment.[ShortCode.create "N"].Head |> fst
    let initialT = plant.Environment.[code "T[max]"] |> TimeSeries.findExact startDate |> fst
    [ (ShortCode.create "x", initialRadius)
      (ShortCode.create "N", initialNitrogen)
      (ShortCode.create "T[max]", initialT)
      (ShortCode.create "bs", initialMass) ] |> Map.ofList

let enforceMinimumRadius minSize plant =
    let bounded = 
        match plant.Growth with
        | RingWidth rw ->
            match rw with
            | Cumulative g -> 
                g 
                |> TimeSeries.toObservations 
                |> Seq.skipWhile (fun (x,_) -> x < minSize) 
                |> TimeSeries.fromObservations 
                |> Cumulative
            | _ -> invalidOp "Not implemented"
        | _ -> invalidOp "Not implemented"
    { plant with Growth = bounded |> RingWidth }

let tailGrowth plant =
    let bounded = 
        match plant.Growth with
        | RingWidth rw ->
            match rw with
            | Cumulative g -> 
                g 
                |> TimeSeries.toObservations 
                |> Seq.tail
                |> TimeSeries.fromObservations 
                |> Cumulative
            | _ -> invalidOp "Not implemented"
        | _ -> invalidOp "Not implemented"
    { plant with Growth = bounded |> RingWidth }


// The `EstimationEngine` must be configured for each shrub:
// - Each shrub requires a custom starting condition based on its
//   time-series. 
// - Each shrub requires clipping to a custom timebound before fitting.
let fit engine endCondition s hypothesis =
    let tMax = 
        if s.Identifier.Value.Contains "S8"
        then summerTemperature |> Seq.find(fun (s,_) -> s = "Marre Sale") // A Varandei shrub
        else summerTemperature |> Seq.find(fun (s,_) -> s = "Hoseda Hard") // A Yamal shrub
    let shrub = s |> PlantIndividual.toCumulativeGrowth |> enforceMinimumRadius 1.53<mm> |> PlantIndividual.zipEnv (code "T[max]") (snd tMax)
    printfn "Shrub = %A" shrub
    let common = shrub |> tailGrowth |> PlantIndividual.keepCommonYears |> PlantIndividual.zipEnv (code "T[max]") (snd tMax)
    printfn "Common = %A" common
    let startDate = (common.Environment.[ShortCode.create "N"]).StartDate |> snd
    let startConditions = getStartValues startDate shrub
    let e = engine |> Bristlecone.withConditioning (Custom startConditions)
    Bristlecone.PlantIndividual.fit e endCondition hypothesis common

let workPackages shrubs (hypotheses:ModelSystem list) engine saveDirectory =
    seq {
        for s in shrubs do
            for h in [ 1 .. hypotheses |> List.length ] do //hypotheses |> List.length ] do // Skip first (non-temperature) hypotheses
                for _ in [ 1 .. Options.chains ] do
                    yield async {
                            let result = fit engine Options.endWhen s hypotheses.[h-1] |> fst
                            Bristlecone.Data.EstimationResult.saveAll saveDirectory s.Identifier.Value h 25 result
                            return result }
    }

// Orchestrate the analyses
let work = workPackages (shrubs |> Seq.where(fun x -> x.Environment.[code "N"].Resolution <> TemporalResolution.Variable)) hypotheses Options.engine Options.resultsDirectory
let run () =
    work |> Seq.take 100 |> Seq.iter (OrchestrationMessage.StartWorkPackage >> Options.orchestrator.Post)



// Post-Processing (Model Selection, Confidence Intervals)
// _______________________________________________________

open Bristlecone.Optimisation.ConfidenceInterval
open Bristlecone.Data

// A. Load best MLE estimates
let results =
    List.allPairs shrubs (hypotheses |> List.mapi(fun i v -> (i+1,v)))
    |> List.map (fun (s,(hi,h)) ->
        let r = 
            [ //EstimationResult.loadAll Options.resultsDirectory s.Identifier.Value h hi
              EstimationResult.loadAll "/Users/andrewmartin/Desktop/Bristlecone Results/LongTermWithRepeatedNValues-TempIndependent" s.Identifier.Value h hi ] |> Seq.concat
        if r |> Seq.isEmpty
        then (s, h, hi, None )
        else (s, h, hi, r |> Seq.minBy(fun x -> x.Likelihood) |> Some ))

// // B. Calculate Akaike weights
// let weights = 
//     results
//     |> Seq.groupBy(fun (s,_,_,_) -> s.Identifier)
//     |> Seq.collect(fun (_,r) ->
//         let weights =
//             r 
//             |> Seq.choose(fun (s,h,hi,r) -> r) 
//             |> ModelSelection.Akaike.akaikeWeights 
//         r
//         |> Seq.zip weights
//         |> Seq.map (fun (w,(s,h,hi,r)) -> s.Identifier.Value, hi, fst w, snd w)
//         |> Seq.toList )
//     |> Seq.toList

// Bristlecone.Data.ModelSelection.save Options.resultsDirectory weights

// // TODO:
// // Save Akaike weights: SubjectId * HypothesisId * -logL * ParamCount * AIC * AICc * AkaikeWeight

// // C. Calculate confidence intervals
// let confidenceIntervals() =
//     results
//     |> Seq.toArray
//     |> Array.choose(fun (s,h,hi,r) ->
//         r |> Option.map (fun res -> 
//             let a = ProfileLikelihood.profile fit Options.engine s h 10 res
//             Bristlecone.Data.Confidence.save Options.resultsDirectory s.Identifier.Value hi res.ResultId a
//             a))


// // Test
// // One-step-ahead analysis
let pretransform (data:CodedMap<TimeSeries<float>>) =
    data
    |> Map.toList
    |> List.collect(fun (k,v) ->
        if k = code "x"
        then [ (k, v); (code "bs", v |> TimeSeries.map(fun (x,_) -> x |> ModelComponents.Proxies.toBiomassMM)) ]
        else [ (k, v)] )
    |> Map.ofList

let oneStepAhead' engine hypothesis s pretransform parameters =
    let tMax = 
        if s.Identifier.Value.Contains "S8"
        then summerTemperature |> Seq.find(fun (s,_) -> s = "Marre Sale") // A Varandei shrub
        else summerTemperature |> Seq.find(fun (s,_) -> s = "Hoseda Hard") // A Yamal shrub
    let shrub = s |> PlantIndividual.toCumulativeGrowth |> enforceMinimumRadius 1.53<mm> |> PlantIndividual.zipEnv (code "T[max]") (snd tMax)
    printfn "Shrub = %A" shrub
    let common = shrub |> tailGrowth |> PlantIndividual.keepCommonYears
    printfn "Common = %A" common
    let startDate = (common.Environment.[ShortCode.create "N"]).StartDate |> snd
    let startConditions = getStartValues startDate shrub
    let e = engine |> Bristlecone.withConditioning (Custom startConditions)
    Bristlecone.PlantIndividual.predictAhead e hypothesis common pretransform parameters

let oneStepAhead = 
    results
    |> Seq.toArray
    |> Array.choose(fun (s,h,hi,r) ->
        r |> Option.map (fun res -> 
            oneStepAhead' Options.engine h s pretransform res.Parameters |> List.map (fun (r,x) -> r.Series) ))

type SOSRow = {
    Year: int
    Variable: string
    Observation: float
    OneStepPrediction: float
}

oneStepAhead
|> Array.collect(fun x -> x |> List.toArray |> Array.map(fun x -> x |> Map.map(fun k v -> v |> TimeSeries.toObservations |> Seq.head)))
|> Array.collect(fun x -> x |> Map.toArray )
|> Array.map(fun (k,o) -> {
    Year = (snd o).Year
    Variable = k.Value
    Observation = (fst o).Obs
    OneStepPrediction = (fst o).Fit
})
|> Array.map (fun x -> sprintf "%i,%s,%f,%f" x.Year x.Variable x.Observation x.OneStepPrediction)
|> String.concat "\n"
|> fun x -> System.IO.File.WriteAllText("/Users/andrewmartin/Desktop/onestepahead2.csv", x)