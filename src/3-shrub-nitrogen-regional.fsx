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
    let resultsDirectory = "/Users/andrewmartin/Desktop/Bristlecone Results/Paper3-Repeated-Final/"
    let chains = 1
    let endWhen = Optimisation.EndConditions.afterIteration 100000
    let logger = CustomLog.logger() //Logging.RealTimeTrace.graphWithConsole 60. 10000
    let filzbachOptions : Optimisation.MonteCarlo.Filzbach.FilzbachSettings<float> = {
        TuneAfterChanges = 20
        MaxScaleChange = 100.00
        MinScaleChange = 0.0010
        BurnLength = Optimisation.EndConditions.afterIteration 500000 }
    let engine =
        Bristlecone.mkContinuous 
        |> Bristlecone.withContinuousTime Integration.MathNet.integrate
        |> Bristlecone.withOutput logger
        // |> Bristlecone.withTunedMCMC [ Optimisation.MonteCarlo.TuneMethod.CovarianceWithScale 0.100, 1000, Optimisation.EndConditions.afterIteration 100000 ]
        // |> Bristlecone.withCustomOptimisation (Optimisation.MonteCarlo.Filzbach.filzbach filzbachOptions)
        // |> Bristlecone.withCustomOptimisation (Optimisation.MonteCarlo.adaptiveMetropolis 0.250 500)
        |> Bristlecone.withCustomOptimisation (Optimisation.MonteCarlo.SimulatedAnnealing.fastSimulatedAnnealing 0.01 true
            { Optimisation.MonteCarlo.SimulatedAnnealing.AnnealSettings<float>.Default with 
                InitialTemperature = 100.
                TemperatureCeiling = Some 100.
                HeatRamp = (fun t -> t + 5.00)
                BoilingAcceptanceRate = 0.85
                HeatStepLength = Optimisation.EndConditions.afterIteration 1000
                TuneLength = 1000
                AnnealStepLength = (fun x -> (*Optimisation.MonteCarlo.SimulatedAnnealing.EndConditions.improvementCount 5000 250 x ||*) Optimisation.EndConditions.afterIteration 10000 x) }) //(Optimisation.EndConditions.afterIteration 10000) })

    let orchestrator = OrchestrationAgent(logger, System.Environment.ProcessorCount, false)


type ComponentLogger(turnOff) =

    let mutable (data:Map<string,Map<float,float>>) = [] |> Map.ofList
    with
        member __.StoreValue(componentId, t, v) =
            if not turnOff then
                let existing = data |> Map.tryFind componentId
                match existing with
                | Some e -> data <- data |> Map.add componentId (e |> Map.add t v)
                | None -> data <- data |> Map.add componentId ([t,v] |> Map.ofList)
            v

        member __.GetAll() = data


// 2. Create Hypotheses
// ----------------------------

module BaseEquations =

    /// Cumulative stem biomass [dBs/dt]
    let biomass b n r gammab geom f tempEffect cLog : float =
        if r > 100000. then nan
        else 
            cLog "growthRate" (b * r * (f n) * geom(b) * tempEffect)
            cLog "biomassLossRate" (gammab * b)
            cLog "nUptake" (f n)
            cLog "geometryEffect" (geom(b))
            cLog "tLimitation" (tempEffect)
            cLog "nLimitation" (r * (f n))
            cLog "nAndtempLimitation" (r * (f n) * tempEffect)
            b * r * (f n) * geom(b) * tempEffect - gammab * b //+ sigmaz

    let soilNitrogen n b gamman y geom f feedback tempEffect cLog : float =
        cLog "plantSoilFeedback" (feedback(b))
        cLog "abioticReplenishment" y
        cLog "nEnvironmentalLoss" (-gamman * n)
        cLog "nUseByPlant" ((geom(b) * b * (f n) * tempEffect))
        y - (geom(b) * b * (f n) * tempEffect) - gamman * n + feedback(b)

    let soilNitrogenNoUptake n b gamman y feedback cLog : float =
        cLog "plantSoilFeedback" (feedback(b))
        cLog "abioticReplenishment" y
        cLog "nEnvironmentalLoss" (-gamman * n)
        cLog "nUseByPlant" 0.
        y - gamman * n + feedback(b)


let ``base model`` maxGrowthRate nLimitation nitrogenFeedback tempEffect additionalParameters (cLog:ComponentLogger) =

    /// Cumulative stem biomass [b].
    let dbsdt' t (b:float) n gammab r maxGrowthRate limit tLimit =
        match limit with
        | Some l -> BaseEquations.biomass b n (r * 1000.) gammab maxGrowthRate l tLimit (fun name l -> cLog.StoreValue(name,t,l) |> ignore)
        | None -> BaseEquations.biomass b n r gammab maxGrowthRate (fun _ -> 1.) tLimit (fun name l -> cLog.StoreValue(name,t,l) |> ignore)

    /// Bioavailable soil nitrogen [N]
    let dndt' t bs n lambda gamman maxGrowthRate feedback limit tLimit = 
        match limit with
        | Some l -> BaseEquations.soilNitrogen n bs gamman lambda maxGrowthRate l feedback tLimit (fun name l -> cLog.StoreValue(name,t,l) |> ignore)
        | None -> BaseEquations.soilNitrogenNoUptake n bs gamman lambda feedback (fun name l -> cLog.StoreValue(name,t,l) |> ignore)

    /// Bristlecone function for dBs/dt
    let dbsdt p t bs (e:Environment) =
        dbsdt' t bs ((e.[ShortCode.create "N"]) |> ModelComponents.Proxies.d15NtoAvailability)
            (p |> Pool.getEstimate "gamma[b]") ((p |> Pool.getEstimate "r")) (maxGrowthRate p) (nLimitation p) (tempEffect p e) //(p |> Pool.getEstimate "sigma[z]")

    /// Bristlecone function for dN/dt
    let dndt p t n (e:Environment) =
        dndt' t (e.[ShortCode.create "bs"]) (n |> ModelComponents.Proxies.d15NtoAvailability) 
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
          (fun p -> ModelComponents.GrowthLimitation.linear ((p |> Pool.getEstimate "a") / 1000.) 5.00), 
           [ code "a",      parameter PositiveOnly   0.001 0.010
             code "r",      parameter PositiveOnly   1.000 10.00 ]        // N-uptake efficiency
          (fun _ -> ModelComponents.GrowthLimitation.none), 
           [ code "r",       parameter PositiveOnly  1.000 10.00 ] ]

    // [B] Loss of plant material may feedback into the soil pool of available nitrogen (instant)
    let feedbackModes =
        [ (fun _ -> ModelComponents.FeedbackToSoil.none), []
          (fun p -> ModelComponents.FeedbackToSoil.withBiomassLoss ((p |> Pool.getEstimate "alpha") / 100.) (p |> Pool.getEstimate "gamma[b]") ),
          [ code "alpha",  parameter PositiveOnly   0.001 0.100 ] ]    // N-recycling efficiency

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
           [ code "Ea",             parameter PositiveOnly 10.00 30.00  ] ]

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
            | Cumulative(_) -> failwith "Not Implemented"
            | Relative(_) -> failwith "Not Implemented"
        | _ -> failwith "Not implemented"

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
            | Cumulative(_) -> failwith "Not Implemented"
            | Relative(_) -> failwith "Not Implemented"
        | _ -> failwith "not implemented"

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
                // printfn "Start cumulative growth = %f" start
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
    // printfn "Shrub = %A" shrub
    let common = shrub |> tailGrowth |> PlantIndividual.keepCommonYears |> PlantIndividual.zipEnv (code "T[max]") (snd tMax)
    // printfn "Common = %A" common
    let startDate = (common.Environment.[ShortCode.create "N"]).StartDate |> snd
    let startConditions = getStartValues startDate shrub
    let e = engine |> Bristlecone.withConditioning (Custom startConditions)
    Bristlecone.PlantIndividual.fit e endCondition hypothesis common

let saveLog (cLog:ComponentLogger) (h:int) plantCode (guid:System.Guid) =
    let filePath = (sprintf "%sbristlecone-%s-%i-components-%s.csv" Options.resultsDirectory plantCode h (guid.ToString()))
    if System.IO.File.Exists filePath then System.IO.File.Delete filePath
    use csv = System.IO.File.AppendText filePath
    let x = cLog.GetAll()
    for m in x do
        for ts in m.Value do
            csv.WriteLine (sprintf "%s,%f,%.16e" m.Key ts.Key ts.Value)

let workPackages shrubs hypotheses engine saveDirectory =
    seq {
        for s in shrubs do
            for h in [ 1 .. hypotheses |> List.length ] do //hypotheses |> List.length ] do // Skip first (non-temperature) hypotheses
                for _ in [ 1 .. Options.chains ] do
                    yield async {
                            let cLog = ComponentLogger(true)
                            let result = fit engine Options.endWhen s (hypotheses.[h-1] cLog) |> fst
                            saveLog cLog h s.Identifier.Value (result).ResultId
                            Bristlecone.Data.EstimationResult.saveAll saveDirectory s.Identifier.Value h 100 result
                            return result }
    }

// Orchestrate the analyses
let work = workPackages (shrubs |> Seq.where(fun x -> x.Environment.[code "N"].Resolution <> TemporalResolution.Variable)) hypotheses Options.engine Options.resultsDirectory
let run () =
    work |> Seq.iter (OrchestrationMessage.StartWorkPackage >> Options.orchestrator.Post)



// Post-Processing (Model Selection, Confidence Intervals)
// _______________________________________________________

open Bristlecone.Optimisation.ConfidenceInterval
open Bristlecone.Data

let loadAll directory subject (modelSystem:ModelSystem) modelId =
    let mles = MLE.load directory subject modelId |> Seq.map(fun (k,v) -> k.ToString(), v)
    let series = Series.load directory subject modelId |> Seq.map(fun (k,v) -> k.ToString(), v)
    mles
    |> Seq.keyMatch series
    |> Seq.map(fun (k,v1,v2) -> (k, (v1, v2)))
    |> Seq.choose(fun (k,(s,(l,p))) ->
        try
            let parameters = 
                // Convert 'b' into 'h'
                modelSystem.Parameters |> Map.map(fun k v -> 
                    if Option.isNone (p |> Map.tryFind k)
                    then
                        printfn "Found a 'b' model..."
                        // h = b / (r * a)
                        let h = (p |> Map.find (code "b")) / ((p |> Map.find (code "r")) * (p |> Map.find (code "a")))
                        Parameter.setEstimate v h
                    else Parameter.setEstimate v (p |> Map.find k))
            { ResultId = k |> System.Guid.Parse
              Likelihood = l
              Parameters = parameters
              Series = s
              Trace = [] } |> Some
        with _ -> None )

// A. Load best MLE estimates
let results =
    List.allPairs shrubs (hypotheses |> List.mapi(fun i v -> (i+1,v)))
    |> List.map (fun (s,(hi,h)) ->
        let r = 
            [ loadAll Options.resultsDirectory s.Identifier.Value (h (ComponentLogger(false))) hi
              loadAll "/Users/andrewmartin/Desktop/Bristlecone Results/LongTermWithRepeatedNValues-TempIndependent" s.Identifier.Value (h (ComponentLogger(false))) hi ] |> Seq.concat
        printfn "[%s %i] %A" s.Identifier.Value hi (r |> Seq.map(fun r -> r.Likelihood))
        if r |> Seq.isEmpty
        then (s, h, hi, None )
        else 
            let r' = r |> Seq.filter(fun x -> not (System.Double.IsNaN(x.Likelihood)))
            if Seq.isEmpty r'
            then (s, h, hi, None )
            else 
                printfn "[%s %i] %A" s.Identifier.Value hi (r' |> Seq.map(fun r -> r.Likelihood))
                (s, h, hi, r' |> Seq.minBy(fun x -> x.Likelihood) |> Some ))

// A1. Save out clog for each using its parameters.
let logOutComponents () =
    let engine = Options.engine |> Bristlecone.withCustomOptimisation (Optimisation.None.passThrough)
    results
    |> List.iter(fun (s,h,hi,e) -> 
        let cLog = ComponentLogger(false)
        let p = (e.Value.Parameters |> Map.map(fun k v -> parameter Unconstrained (v |> Parameter.getEstimate) (v |> Parameter.getEstimate)))
        // printfn "Parameters are %A" p
        let hypothesis = { (h cLog) with Parameters = p }
        let result = fit engine (Optimisation.EndConditions.afterIteration 0 ) s hypothesis |> fst
        saveLog cLog hi s.Identifier.Value (e).Value.ResultId )


// B. Calculate Akaike weights
let weights = 
    results
    |> Seq.groupBy(fun (s,_,_,_) -> s.Identifier)
    |> Seq.collect(fun (_,r) ->
        let weights =
            r 
            |> Seq.choose(fun (s,h,hi,r) -> r) 
            |> fun x -> printfn "Seq is %A" x; x
            |> ModelSelection.Akaike.akaikeWeights 
        r
        |> Seq.zip weights
        |> Seq.map (fun (w,(s,h,hi,r)) -> s, h, hi, fst w, snd w)
        |> Seq.toList )
    |> Seq.toList

results 
|> Seq.filter(fun (p,m,a,b) -> a = 13 && p.Identifier.Value = "YUSL39A")
|> Seq.map(fun (p,m,a,b) -> b)

//Bristlecone.Data.ModelSelection.save Options.resultsDirectory (weights |> List.map(fun (x,y,z,a,b) -> x.Identifier.Value,z,a,b))

// // TODO:
// // Save Akaike weights: SubjectId * HypothesisId * -logL * ParamCount * AIC * AICc * AkaikeWeight

// C. Calculate confidence intervals
let confidenceIntervals() =
    weights
    |> Seq.toArray
    |> Array.rev
    |> Array.splitInto 12
    |> Array.Parallel.map(fun x -> x |> Array.map (fun (s,h,hi,res,w) ->
        let a = ProfileLikelihood.profile (fun x y z a -> fit x y z a |> fst) Options.engine s (h (ComponentLogger(true))) 10 res
        Bristlecone.Data.Confidence.save Options.resultsDirectory s.Identifier.Value hi res.ResultId a
        a))


module ProfileWithLogging =

    open Bristlecone.EstimationEngine
    open Bristlecone.Logging
    open Bristlecone.Optimisation.ConfidenceInterval.ProfileLikelihood

    let logExists (h:int) plantCode (guid:System.Guid) =
        let filePath = (sprintf "%sbristlecone-%s-%i-componentsci-%s.csv" Options.resultsDirectory plantCode h (guid.ToString()))
        System.IO.File.Exists filePath

    let saveLog (results:seq<string>) (h:int) plantCode (guid:System.Guid) =
        let filePath = (sprintf "%sbristlecone-%s-%i-componentsci-%s.csv" Options.resultsDirectory plantCode h (guid.ToString()))
        if System.IO.File.Exists filePath then System.IO.File.Delete filePath
        use csv = System.IO.File.AppendText filePath
        csv.WriteLine "Component,Time,CI,MinValue,MaxValue"
        for m in results do csv.WriteLine m

    /// The profile likelihood method samples the likelihood space
    /// around the Maximum Likelihood Estimate 
    let profile iterations fit engine subject (timeSeries:System.DateTime seq) (hypothesis:ComponentLogger-> ModelSystem) n (result:EstimationResult) =

        // 1. Set estimation bounds to the MLE
        let mleToBounds mlePool = mlePool |> Map.map(fun k v -> Parameter.create (Parameter.detatchConstraint v |> snd) (v |> Parameter.getEstimate) (v |> Parameter.getEstimate))
        let hypothesisMle a = { hypothesis a with Parameters = mleToBounds result.Parameters }

        // 2. Generate a trace of at least n samples that deviate in L less than 2.0
        let results = OptimOutput()
        let mle = result.Likelihood

        let customFit = fit { engine with OptimiseWith = CustomOptimisationMethod.classic {CustomOptimisationMethod.TuneSettings.Default with MinimumSampleSize = iterations; KMax = iterations}; LogTo = fun e -> engine.LogTo e; results.SaveEvent e }
        let rec fit' currentTrace =
            let a = customFit (Optimisation.EndConditions.afterIteration n) subject (hypothesisMle (ComponentLogger(true)))
            let validTrace = a.Trace |> List.filter(fun (l,_) -> (l - mle) < 2.00 && (l - mle) > 0.00) |> List.distinct
            engine.LogTo <| GeneralEvent (sprintf "Profiling efficiency: %f/1.0." ((validTrace |> List.length |> float) / (float n)))
            // if currentTrace |> List.append validTrace |> List.length < n
            // then currentTrace |> List.append validTrace |> fit'
            // else 
            currentTrace |> List.append validTrace
        fit' [] |> ignore

        let trace = 
            results.GetAll()
            |> List.map(fun s -> (s.Likelihood, s.Theta |> Seq.toArray))
        engine.LogTo <| GeneralEvent (sprintf "Actual trace was %A" trace)

        printfn "Getting components for %i combinations" (trace |> List.distinct |> List.filter(fun x -> (x |> fst) - mle < Bounds.upperBound) |> List.length)

        let customFit2 = fit (engine |> Bristlecone.withCustomOptimisation  Optimisation.None.passThrough |> Bristlecone.withOutput ignore)
        let mleToBounds2 mlePool = mlePool |> Map.map(fun k v -> Parameter.create Unconstrained (v |> Parameter.getEstimate) (v |> Parameter.getEstimate))
        let components =
            let hypothesisMle a = { hypothesis a with Parameters = mleToBounds2 result.Parameters }
            trace
            |> List.distinct
            |> List.filter(fun x -> (x |> fst) - mle < Bounds.upperBound)
            |> List.map(fun t -> 
                let p = Optimisation.ParameterPool.fromPoint ((hypothesisMle (ComponentLogger(false))).Parameters) (t |> snd)
                let log = ComponentLogger(false)
                let h = { (hypothesisMle log) with Parameters = p |> mleToBounds2 }
                let result = customFit2 (Optimisation.EndConditions.afterIteration n) subject h
                let logs = log.GetAll()
                // printfn "Logs are %A" logs
                fst t, logs )

        printfn "Calculating at mle"

        let mleLogs =
            let hypothesisMle a = { hypothesis a with Parameters = mleToBounds2 result.Parameters }
            let log = ComponentLogger(false)
            let result = customFit2 (Optimisation.EndConditions.afterIteration n) subject (hypothesisMle log)
            (result.Likelihood, log.GetAll())

        let ts = timeSeries |> Seq.toArray
        let results bound ci = 
            components
            |> Seq.append [mleLogs]
            |> Seq.distinct
            |> Seq.filter(fun x -> (x |> fst) - mle <= bound)
            |> Seq.collect(fun (l,m) ->
                m
                |> Seq.collect(fun kv ->
                    kv.Value
                    |> Seq.map(fun ts -> kv.Key, ts.Key, ts.Value, l )))
            |> Seq.groupBy(fun (a,b,c,d) -> (a,b))
            |> Seq.map(fun (g,s) -> 
                (fst g), (snd g), s |> Seq.map(fun (a,b,c,d) -> c) |> Seq.min, s |> Seq.map(fun (a,b,c,d) -> c) |> Seq.max)
            |> Seq.filter(fun (_,b,_,_) -> b % 1.0 = 0. && b >= 0.)
            |> Seq.map(fun (a,b,c,d) -> sprintf "%s,%s,%s,%.16e,%.16e" a (ts.[int b].ToString()) ci c d)

        // 3. Calculate min and max at the specified limit for each parameter
        let lowerInterval = trace |> interval (hypothesisMle (ComponentLogger(true))).Parameters.Count mle Bounds.lowerBound
        let upperInterval = trace |> interval (hypothesisMle (ComponentLogger(true))).Parameters.Count mle Bounds.upperBound

        let paramResult = 
            result.Parameters
            |> Seq.zip3 lowerInterval upperInterval
            |> Seq.map(fun ((l1,l2),(u1,u2),p) -> 
                p.Key, {
                    Estimate = p.Value |> Parameter.getEstimate
                    ``68%`` = { Lower = l1; Upper = l2 }
                    ``95%`` = { Lower = u1; Upper = u2 }})
            |> Map.ofSeq

        let componentResult = 
            results Bounds.upperBound "95"
            |> Seq.append (results Bounds.lowerBound "68")
            |> Seq.append (results 0. "mle")

        paramResult, componentResult

    let confidenceIntervals () =
        weights
        |> Seq.sortByDescending(fun (a,b,c,d,e) -> e.Weight)
        |> Seq.where(fun (a,b,c,d,e) -> e.Weight > 0.05)
        |> Seq.toArray
        |> Array.splitInto System.Environment.ProcessorCount
        |> Array.Parallel.map(fun x -> x |> Array.map (fun (s,h,hi,res,w) ->
            if not <| logExists hi s.Identifier.Value res.ResultId then
                let ts = s.Environment.[ShortCode.create "N"] |> TimeSeries.dates
                let (pr,ar) = profile 50000 (fun x y z a -> fit x y z a |> fst) Options.engine s ts h 10 res
                saveLog ar hi s.Identifier.Value res.ResultId
                Bristlecone.Data.Confidence.save Options.resultsDirectory s.Identifier.Value hi res.ResultId pr
            ))

// ProfileWithLogging.confidenceIntervals()


// // // Test
// // // One-step-ahead analysis
// let pretransform (data:CodedMap<TimeSeries<float>>) =
//     data
//     |> Map.toList
//     |> List.collect(fun (k,v) ->
//         if k = code "x"
//         then [ (k, v); (code "bs", v |> TimeSeries.map(fun (x,_) -> x |> ModelComponents.Proxies.toBiomassMM)) ]
//         else [ (k, v)] )
//     |> Map.ofList

// let oneStepAhead' engine hypothesis s pretransform parameters =
//     let tMax = 
//         if s.Identifier.Value.Contains "S8"
//         then summerTemperature |> Seq.find(fun (s,_) -> s = "Marre Sale") // A Varandei shrub
//         else summerTemperature |> Seq.find(fun (s,_) -> s = "Hoseda Hard") // A Yamal shrub
//     let shrub = s |> PlantIndividual.toCumulativeGrowth |> enforceMinimumRadius 1.53<mm> |> PlantIndividual.zipEnv (code "T[max]") (snd tMax)
//     // printfn "Shrub = %A" shrub
//     let common = shrub |> tailGrowth |> PlantIndividual.keepCommonYears
//     // printfn "Common = %A" common
//     let startDate = (common.Environment.[ShortCode.create "N"]).StartDate |> snd
//     let startConditions = getStartValues startDate shrub
//     let e = engine |> Bristlecone.withConditioning (Custom startConditions)
//     Bristlecone.PlantIndividual.predictAhead e hypothesis common pretransform parameters

// let oneStepAhead = 
//     results
//     |> Seq.toArray
//     |> Array.choose(fun (s,h,hi,r) ->
//         r |> Option.map (fun res -> 
//             oneStepAhead' Options.engine h s pretransform res.Parameters |> List.map (fun (r,x) -> r.Series) ))

// type SOSRow = {
//     Year: int
//     Variable: string
//     Observation: float
//     OneStepPrediction: float
// }

// oneStepAhead
// |> Array.collect(fun x -> x |> List.toArray |> Array.map(fun x -> x |> Map.map(fun k v -> v |> TimeSeries.toObservations |> Seq.head)))
// |> Array.collect(fun x -> x |> Map.toArray )
// |> Array.map(fun (k,o) -> {
//     Year = (snd o).Year
//     Variable = k.Value
//     Observation = (fst o).Obs
//     OneStepPrediction = (fst o).Fit
// })
// |> Array.map (fun x -> sprintf "%i,%s,%f,%f" x.Year x.Variable x.Observation x.OneStepPrediction)
// |> String.concat "\n"
// |> fun x -> System.IO.File.WriteAllText("/Users/andrewmartin/Desktop/onestepahead2.csv", x)




// Start chains that have already finished and run again..

let workPackage (shrub:PlantIndividual) hi hypothesis engine saveDirectory =
    async {
        let cLog = ComponentLogger(true)
        let result = fit engine Options.endWhen shrub (hypothesis cLog) |> fst
        saveLog cLog hi shrub.Identifier.Value (result).ResultId
        Bristlecone.Data.EstimationResult.saveAll saveDirectory shrub.Identifier.Value hi 100 result
        return result
    }

let work2 =
    results
    |> List.collect(fun (plant,h,hi,r) ->
        match r with
        | Some r ->
            let mleToBounds mlePool = mlePool |> Map.map(fun k v -> Parameter.create (Parameter.detatchConstraint v |> snd) (v |> Parameter.getEstimate) (v |> Parameter.getEstimate))
            let hypothesisMle a = { h a with Parameters = mleToBounds r.Parameters }
            [ workPackage plant hi hypothesisMle Options.engine Options.resultsDirectory ]
        | None -> [] )

let runAgain () =
    work2 |> Seq.skip 204 |> Seq.iter (OrchestrationMessage.StartWorkPackage >> Options.orchestrator.Post)




/// Stepwise solving: use MLE from less complex models in more complex models.
/// A. Non-temperature 12 vs temperature 12
let runAugmented () =

    let work3 =
        results
        |> List.collect(fun (plant,h,hi,r) ->
            match r with
            | Some r ->
                if hi > 12 then []
                else
                    // Is NOT a temperature hypothesis.
                    // Get the temperature-eqivalent hypothesis and add Ea to this parameter pool.
                    let h = hypotheses.[hi + 12 - 1]
                    let mleToBounds mlePool = 
                        mlePool 
                        |> Map.map(fun k v -> Parameter.create (Parameter.detatchConstraint v |> snd) (v |> Parameter.getEstimate) ((v |> Parameter.getEstimate) + 0.01))
                        |> Map.add (code "Ea") (parameter PositiveOnly 0.001 0.002 )
                    let hypothesisMle a = { h a with Parameters = mleToBounds r.Parameters }
                    [ workPackage plant hi hypothesisMle Options.engine Options.resultsDirectory ]
            | None -> [] )

    work3 |> Seq.iter (OrchestrationMessage.StartWorkPackage >> Options.orchestrator.Post)

