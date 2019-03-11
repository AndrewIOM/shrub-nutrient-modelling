#r "../packages/NETStandard.Library.NETFramework/build/net461/lib/netstandard.dll"
#load "../packages/Bristlecone/bristlecone.fsx"
#load "../packages/Bristlecone/charts.fsx"
#load "components/components.fsx"
#r "../packages/Bristlecone.Dendro/lib/netstandard2.0/Bristlecone.Dendro.dll"

////////////////////////////////////////////////////
/// Long-Term N-Shrub relations in Yamal, Russia
////////////////////////////////////////////////////

// Individual-based modelling of plant growth on a three-yearly resolution.

open Bristlecone
open Bristlecone.ModelSystem
open Bristlecone.PlantIndividual
open Bristlecone.Workflow.Orchestration

// 1. Configure Options
// ----------------------------

module Options =
    let resultsDirectory = "/Users/andrewmartin/Desktop/Bristlecone Results/LongTermWithRepeatedNValues/"
    let chains = 6
    let endWhen = Optimisation.EndConditions.afterIteration 1000000
    let logger = Logging.RealTimeTrace.graphWithConsole 30. 10000
    let engine =
        Bristlecone.mkContinuous 
        |> Bristlecone.withContinuousTime Integration.MathNet.integrate
        |> Bristlecone.withOutput logger
        |> Bristlecone.withCustomOptimisation (Optimisation.MonteCarlo.Filzbach.filzbach 0.01 (Optimisation.EndConditions.afterIteration 200000))
        // |> Bristlecone.withCustomOptimisation (Optimisation.MonteCarlo.SimulatedAnnealing.fastSimulatedAnnealing 0.0001 
        //     { Optimisation.MonteCarlo.SimulatedAnnealing.AnnealSettings<float>.Default with 
        //         BoilingAcceptanceRate = 0.60
        //         HeatRamp = (fun t -> t + sqrt t); TemperatureCeiling = Some 500.
        //         HeatStepLength = Optimisation.EndConditions.afterIteration 1000
        //         AnnealStepLength = Optimisation.EndConditions.afterIteration 5000
        //         InitialTemperature = 1. }) //(Optimisation.EndConditions.afterIteration 10000) })

    let orchestrator = OrchestrationAgent(logger, 6)//System.Environment.ProcessorCount)

// TEMP - Likelihood function with shrub comparison term
let private getData s (predictions:CodedMap<PredictedSeries>) = predictions.Item (ShortCode.create s)
let private pi = System.Math.PI

/// Negative log likelihood for a bivariate normal distribution.
/// For two random variables with bivariate normal N(u1,u2,sigma1,sigma2,rho).
let bivariateGaussian' (p:ParameterPool) obsx obsy expx expy = 
    let diffx = obsx - expx
    let diffy = obsy - expy
    let sigmax = p |> Pool.getEstimate "sigma[x]"
    let sigmay = p |> Pool.getEstimate "sigma[y]"
    // let sigmaz = (p |> Pool.getEstimate "sigma[z]") * 100.
    let rho = p |> Pool.getEstimate "rho"
    let zta1 = (diffx / sigmax) ** 2.
    let zta2 = 2. * rho * ((diffx / sigmax) ** 1.) * ((diffy / sigmay) ** 1.)
    let zta3 = (diffy / sigmay) ** 2.
    // let zta4 = (diffx / sigmaz) ** 2.
    let vNegLog = 2. * pi * sigmax * sigmay * sqrt (1. - rho ** 2.)
    let q = (1. / (1. - rho ** 2.)) * (zta1 + zta3 - zta2)
    vNegLog + (1./2.) * q

/// <summary>
/// Log likelihood function for dual simultaneous system, assuming Gaussian error for both x and y.
/// </summary> 
let bivariateGaussian key1 key2 p data = 
    let x = data |> getData key1
    let y = data |> getData key2
    [1 .. (Array.length x.Observed) - 1] 
    |> List.sumBy (fun i -> (bivariateGaussian' p x.Observed.[i] y.Observed.[i] x.Expected.[i] y.Expected.[i])) 




// 2. Create Hypotheses
// ----------------------------

module BaseEquations =

    /// Cumulative stem biomass [dBs/dt]
    let biomass b n r gammab geom f : float =
        b * r * (f n) * geom(b) - gammab * b

    let soilNitrogen n b gamman y geom f feedback : float =
        y - (geom(b) * b * (f n)) - gamman * n + feedback(b)

    let soilNitrogenNoUptake n b gamman y feedback : float =
        y - gamman * n + feedback(b)


let ``base model`` maxGrowthRate nLimitation nitrogenFeedback additionalParameters =

    /// Cumulative stem biomass [b].
    let dbsdt' (b:float) n gammab r maxGrowthRate limit =
        match limit with
        | Some l -> BaseEquations.biomass b n (r * 1000.) gammab maxGrowthRate l
        | None -> BaseEquations.biomass b n r gammab maxGrowthRate (fun _ -> 1.)

    /// Bioavailable soil nitrogen [N]
    let dndt' bs n lambda gamman maxGrowthRate feedback limit = 
        match limit with
        | Some l -> BaseEquations.soilNitrogen n bs gamman lambda maxGrowthRate l feedback
        | None -> BaseEquations.soilNitrogenNoUptake n bs gamman lambda feedback

    /// Bristlecone function for dBs/dt
    let dbsdt p _ bs (e:Environment) =
        dbsdt' bs ((e.[ShortCode.create "N"]) |> ModelComponents.Proxies.d15NtoAvailability)
            (p |> Pool.getEstimate "gamma[b]") ((p |> Pool.getEstimate "r")) (maxGrowthRate p) (nLimitation p)

    /// Bristlecone function for dN/dt
    let dndt p _ n (e:Environment) =
        dndt' (e.[ShortCode.create "bs"]) (n |> ModelComponents.Proxies.d15NtoAvailability) 
            (p |> Pool.getEstimate "lambda") (p |> Pool.getEstimate "gamma[n]") (maxGrowthRate p) (nitrogenFeedback p) (nLimitation p)

    /// Measurement (Size) variable: stem radius
    let stemRadius lastRadius lastEnv env =
        let oldCumulativeMass = lookup lastEnv "bs"
        let newCumulativeMass = lookup env "bs"
        // printfn "R = %f; New mass = %f; Old mass = %f" lastRadius newCumulativeMass oldCumulativeMass
        if (newCumulativeMass - oldCumulativeMass) > 0.
        then newCumulativeMass |> ModelComponents.Proxies.toRadiusMM
        else lastRadius

    { Equations  = [ code "bs",        dbsdt
                     code "N",         dndt ] |> Map.ofList
      Measures   = [ code "x",         stemRadius ] |> Map.ofList
      Parameters = [ // for nitrogen dynamics
                     code "lambda",    parameter PositiveOnly   0.500 1.500   // Rate of nitrogen replenishment
                     code "gamma[n]",  parameter PositiveOnly   0.001 0.200   // Loss rate of nitrogen
                     // for shrub physiology
                     code "gamma[b]",  parameter PositiveOnly   0.001 0.200   // Loss rate of biomass
                     // for likelihood function
                     code "rho",       parameter Unconstrained  -0.50 0.500   // Covariance between growth and nitrogen
                     code "sigma[x]",  parameter PositiveOnly   0.100 1.200   // Standard deviation of x (biomass)
                     code "sigma[y]",  parameter PositiveOnly   0.250 0.750   // Standard deviation of y (nitrogen)
                     //code "sigma[z]",  parameter Unconstrained  0.250 0.750   // Inter-shrub comparison term
                    ] |> List.append additionalParameters |> Map.ofList
      Likelihood = bivariateGaussian "x" "N" }

let hypotheses =

    // [A] N may limited growth via combined N-limitations on (a) photosynthetic and (b) uptake rates
    let limitationModes =
        [ (fun p -> ModelComponents.GrowthLimitation.hollingDiscModel ((p |> Pool.getEstimate "a") / 1000.) ((p |> Pool.getEstimate "r") * 1000.) (p |> Pool.getEstimate "h")),
           [ code "a",      parameter PositiveOnly   0.100 4.000          // N-uptake efficiency
             code "h",      parameter PositiveOnly   0.100 4.000
             code "r",      parameter PositiveOnly   0.500 1.000 ]      // N-handling time (including uptake and incorporation)
          (fun p -> ModelComponents.GrowthLimitation.linear ((p |> Pool.getEstimate "a") / 1000.)), 
           [ code "a",      parameter PositiveOnly   0.001 0.010
             code "r",      parameter PositiveOnly   0.500 1.000 ]        // N-uptake efficiency
          (fun _ -> ModelComponents.GrowthLimitation.none), 
          [ code "r",       parameter PositiveOnly   0.500 1.000 ] ]

    // [B] Loss of plant material may feedback into the soil pool of available nitrogen (instant)
    let feedbackModes =
        [ (fun _ -> ModelComponents.FeedbackToSoil.none), []
          (fun p -> ModelComponents.FeedbackToSoil.withBiomassLoss ((p |> Pool.getEstimate "alpha") / 100.) (p |> Pool.getEstimate "gamma[b]") ),
          [ code "alpha",  parameter PositiveOnly   0.01 1.00 ] ]    // N-recycling efficiency

    // [C] A plant may be subject to mechanical constraints on its maximum size
    let geometricModes = 
        [  (fun p -> ModelComponents.GeometricConstraint.none), []
           (fun p -> ModelComponents.GeometricConstraint.chapmanRichards ((p |> Pool.getEstimate "k") * 1000.)),
           [ code "k",  parameter PositiveOnly   3.00 5.00 ] ]     // Asymptotic biomass (grams)

    List.combine3 geometricModes limitationModes feedbackModes
    |> List.map (fun ((growth,gp),(limit,lp),(feedback,fp)) -> 
        ``base model`` growth limit feedback (List.concat [lp; fp; gp]))


// 3. Load Real Data and Estimate
// ----------------------------

open FSharp.Data

// 1. Load in ring widths
let rw = Bristlecone.Data.PlantIndividual.loadRingWidths (__SOURCE_DIRECTORY__ + "/../data/yamal-rw.csv")

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
    //|> List.skip 10 |> List.take 1

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
    [ (ShortCode.create "x", initialRadius)
      (ShortCode.create "N", initialNitrogen)
      (ShortCode.create "bs", initialMass) ] |> Map.ofList

let workPackages shrubs hypotheses engine saveDirectory =
    seq {
        for s in shrubs do

            // 1. Arrange the subject and settings
            let shrub = s |> PlantIndividual.toCumulativeGrowth
            printfn "Shrub = %A" shrub
            let common = shrub |> PlantIndividual.keepCommonYears
            printfn "Common = %A" common
            let startDate = (common.Environment.[ShortCode.create "N"]).StartDate |> snd
            let startConditions = getStartValues startDate shrub
            let e = engine |> Bristlecone.withConditioning (Custom startConditions)

            printfn "Start conditions: %A" startConditions

            // 2. Setup batches of dependent analyses
            for h in [ 1 .. hypotheses |> List.length ] do
                for i in [ 1 .. Options.chains ] do
                    let jobId = System.Guid.NewGuid()
                    yield async {
                            // A. Compute result
                            let result = Bristlecone.PlantIndividual.fit e Options.endWhen hypotheses.[h-1] common
                            // B. Save to file
                            Bristlecone.Data.Cache.saveTrace saveDirectory s.Identifier.Value h jobId [ result ]
                            return result }
    }

// Orchestrate the analyses
let work = workPackages (shrubs |> Seq.where(fun x -> x.Identifier.Value = "S8N0115A") |> Seq.where(fun x -> x.Environment.[code "N"].Resolution <> TemporalResolution.Variable)) hypotheses Options.engine Options.resultsDirectory
work  |> Seq.skip 0 |> Seq.take 6 |> Seq.iter (OrchestrationMessage.StartWorkPackage >> Options.orchestrator.Post)

// work |> Seq.head |> Async.RunSynchronously
