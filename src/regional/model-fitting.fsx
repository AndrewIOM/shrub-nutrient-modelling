(**
# Long-Term N-Shrub relations in Yamal, Russia

The following code demonstrates model-fitting 
and model-selection for the mechanisms controlling nutrient 
limitation to shrub individuals in the Arctic tundra.

This script completes model-fitting and model-selection
for 24 shrub individuals. Models for nitrogen limitation,
temperature limitation, and size-related asymptotes are
included.

The hypotheses are tested on a per-shrub basis. For each shrub,
there are four components that may vary in the model, which
when multiplied together results in 12 possible combinations
of mechanisms for each shrub.

First, we must load Bristlecone and the seperate script that
defines the allometric equations that we will use later.
We can simply load the fsx file that we previously defined the
model system within; this contains the reference to Bristlecone.
*)

#load "model/model.fsx"

open Bristlecone // Opens Bristlecone core library and estimation engine
open Bristlecone.Language // Open the language for writing Bristlecone models
open Bristlecone.Time

(**
### Loading in the model

*)

(**
### Setting up a *Bristlecone engine*

A bristlecone engine provides a fixed setup for estimating parameters from data.
We use the same engine for all model fits within a single study.

Here, we scaffold an engine from `Bristlecone.mkContinuous`, as we are working
with continuous-time models.
*)

let output = Logging.Console.logger 25<iteration>

let engine: EstimationEngine.EstimationEngine<DatingMethods.Annual,int<year>,year,1> =
    Bristlecone.mkContinuous ()
    |> Bristlecone.withContinuousTime Integration.RungeKutta.rk4
    |> Bristlecone.withTimeConversion DateMode.Conversion.Annual.toYears
    |> Bristlecone.withOutput output
    |> Bristlecone.withConditioning Conditioning.RepeatFirstDataPoint
    |> Bristlecone.withCustomOptimisation (Optimisation.MonteCarlo.SimulatedAnnealing.fastSimulatedAnnealing 0.01<``optim-space``> false
                { HeatStepLength = Optimisation.EndConditions.Profiles.SimulatedAnnealing.heating
                  HeatRamp = fun t -> t * 1.10
                  BoilingAcceptanceRate = 0.85
                  TemperatureCeiling = Some 200.
                  TemperatureFloor = Some 0.75
                  InitialTemperature = 1.00
                  PreTuneEnd = Optimisation.EndConditions.Profiles.SimulatedAnnealing.preTuning
                  Tuning = { MinTuneLength = 300<iteration>; MaxTuneLength = 10000<iteration>; RequiredStableCount = Some 5; TuneN = 20<iteration> }
                  AnnealStepEnd = Optimisation.EndConditions.Profiles.SimulatedAnnealing.annealing })

(**
### Read in ring width, isotope, and weather station data

Here, we are using the Bristlecone.Dendro package to 
read in dendroecological data. For other problems, any
method to wrangle the data into a `TimeSeries` is acceptable.
*)

open Bristlecone.Dendro

let ringWidthDataset =
    Data.PlantIndividual.Csv.loadRingWidths (__SOURCE_DIRECTORY__ + "/../../data/yamal-rw.csv")

(**
#### Read in air temperature data
*)

open FSharp.Data
open FSharp.Data.UnitSystems.SI.UnitSymbols

type TemperatureData = CsvProvider<"../../data/yamal-mean-temperatures.csv">

let jjaTemperaturesByStation =
    let maxTemperatures = TemperatureData.Load "../../data/yamal-mean-temperatures.csv"
    maxTemperatures.Rows
    |> Seq.groupBy(fun r -> r.Station)
    |> Seq.map(fun (station,r) -> station, r |> Seq.map(fun r -> (float r.JJA + 273.15) * 1.<K>, DatingMethods.Annual <| r.Year * 1<year>) |> TimeSeries.fromObservations DateMode.annualDateMode)
    |> Seq.toList

(**
#### Read in two isotope datasets at three-yearly resolution

For this analysis, we conduct model-fitting on a common three-year timeline
across all shrubs ending at [2009 - 2012]. The below code ensures that data
resolution is lowered when it is higher than the three-year bins. For example,
some shrubs have higher resolution (annual) d15N data from 1980.
*)

type Regional15N = CsvProvider<"../../data/yamal-d15N-lowres.csv">
type YuribeiD15N = CsvProvider<"../../data/yuribei-d15N-imputed.csv">

let bins = seq { 1900<year> .. 3<year> .. 2018<year> }
let regionalD15N = Regional15N.Load (__SOURCE_DIRECTORY__ + "/../../data/yamal-d15N-lowres.csv")
let yuribeiD15N = YuribeiD15N.Load (__SOURCE_DIRECTORY__ + "/../../data/yuribei-d15N-imputed.csv")

/// Converts input data into 3-year binned d15N and increment data
let threeYearToAnnual (plant:PlantIndividual.PlantIndividual<Units.millimetre, DatingMethods.Annual,int<year>,int<year>>) =

    let growth = GrowthSeries.absolute plant.Growth

    let isotopeByBin : (int<year>*float) list list = 
        bins
        |> Seq.choose(fun binStart ->
            let lowResBins =
                regionalD15N.Rows
                |> Seq.where(fun n -> n.Shrub = plant.Identifier.Value)
                |> Seq.where(fun n -> n.BinLatest * 1<year> <= (binStart + 2<year>) && n.BinOldest * 1<year> >= binStart)
                |> Seq.map(fun r -> (r.BinOldest * 1<year>, r.BinLatest * 1<year>, r.D15N))

            let highResBins =
                yuribeiD15N.Rows
                |> Seq.where(fun n -> n.``Plant Code`` = plant.Identifier.Value)
                |> Seq.where(fun n -> n.Date.Year * 1<year> <= (binStart + 2<year>) && n.Date.Year * 1<year> >= binStart)
                |> Seq.map(fun r -> (r.Date.Year * 1<year>, r.Date.Year * 1<year>, (float r.``Predictor 2``) ))

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
                                let increments = allData |> Seq.map (fun (start,en,_) -> [ start .. 1<year> .. en ] |> List.sumBy(fun y -> growth |> TimeSeries.findExact (DatingMethods.Annual y) |> fst ))
                                let totalIncrement = increments |> Seq.sum
                                increments |> Seq.map(fun i -> i / totalIncrement)
                            allData
                            |> Seq.zip weights
                            |> Seq.averageBy(fun (weight,(s,e,d15n)) -> d15n * weight )
                        Some <| ([ binStart .. 1<year> .. binStart + 2<year> ] |> List.map (fun y -> y, weightedAverage))
                    with | _ -> None
                | _ ->
                    highResBins 
                    |> Seq.map(fun (x,y,z) -> (x, z) ) 
                    |> Seq.toList
                    |> Some )
        |> Seq.toList
    List.concat isotopeByBin

let dataset =
    ringWidthDataset
    |> Seq.choose(fun plant ->
        let isotope3Year = threeYearToAnnual plant
        if isotope3Year.Length = 0 then None
        else
            let isoTs =
                isotope3Year
                |> List.map (fun (x,y) -> y, DatingMethods.Annual x)
                |> TimeSeries.fromObservations DateMode.annualDateMode
            plant |> PlantIndividual.zipEnvironment Model.N isoTs |> Some )
    |> Seq.toList

(**
The final dataset in `dataset` includes each plant individual with
its associated isotope data ('N'), both at three-yearly resoluion.

### Model-fitting to real data

In this step, we will apply our `EstimationEngine` to the real
shrub data and model hypotheses.

Because there are so many hypotheses (x12) and individuals (x10),
and we want to run replicate analyses (x3), it makes sense to
run these 360 workloads in parallel. To do this, we can make use of
the workflow tools in the `Bristlecone.Workflow` namespace.

An orchestration agent is an agent that queues and runs work items
in parallel. Below, we demonstrate setting up an agent:
*)

open Bristlecone.Workflow

let orchestrator: Orchestration.OrchestrationAgent<DatingMethods.Annual,int<year>,int<year>> =
    Orchestration.OrchestrationAgent(
        writeOut = output,
        maxSimultaneous = System.Environment.ProcessorCount,
        retainResults = false
    )

(**

Before we can run the models, there are some specific considerations
required for this problem.

#### Setting the start values

Given time-series data, a common modelling problem is the question of
what to set t = 0 as. One strategy is to repeat the first data point (t1)
as t0. In this instance, out isotope time-series are shorter than the
ring-width time-series, so we generally know the real value of one
time-series but not the other. Because of this, we use a custom start
point for each shrub individual.
*)

let startValues (startDate:DatingMethods.Annual) (plant: PlantIndividual.PlantIndividual<Units.millimetre,DatingMethods.Annual,'timeunit,'timespan>) =
    let removeUnit (x: float<_>) = float x

    let initialRadius =
        match plant.Growth with
        | GrowthSeries.Cumulative c ->
            let trimmed = c |> TimeSeries.trimStart (startDate.Value - 1<year> |> DatingMethods.Annual)
            match trimmed with
            | Some t -> t.Values |> Seq.head
            | None -> failwith "Could not get t0 from ring-width series"
        | _ -> invalidOp "Not applicable"

    let initialMass = Constant initialRadius |> ShrubModel.Allometry.Proxies.toBiomassMM |> ExpressionCompiler.compileSimple
    let initialNitrogen = plant.Environment.[Model.N.Code].Head |> fst

    [ Model.SR.Code, initialRadius |> removeUnit
      Model.B.Code,  initialMass |> removeUnit
      Model.N.Code,  initialNitrogen ]
    |> Map.ofList

(**
Next, we set up some configuration settings, then wrap up our hypotheses and shrubs
as work packages for the orchestrator, where a work package is a `Async<ModelSystem.EstimationResult>`.

The function to make the work packages first sets up the common time period and custom
start values for each shrub, then creates per-hypothesis work packages for that shrub.
*)

module Config =

    let numberOfReplicates = 3
    let resultsDirectory = "~/Desktop/Bristlecone-3.0/Regional/"
    let thinTrace = Some 100
    let endWhen = Optimisation.EndConditions.atIteration 100000<iteration>

let s = dataset.Head
let h = Model.hypotheses.Head

open Constants.Allometrics

// Function to scaffold work packages
let workPackages (shrubs: PlantIndividual.PlantIndividual<Units.millimetre,DatingMethods.Annual, int<year>, int<year>> seq) (hypotheses: Hypotheses.Hypothesis<'timeindex> list) engine saveDirectory =
    seq {
        for s in shrubs do

            let tMax = 
                if s.Identifier.Value.Contains "S8"
                then jjaTemperaturesByStation |> Seq.find(fun (s,_) -> s = "Marre Sale") // A Varandei shrub
                else jjaTemperaturesByStation |> Seq.find(fun (s,_) -> s = "Hoseda Hard") // A Yamal shrub

            // 1. Arrange the subject and settings
            let shrub = 
                s
                |> PlantIndividual.enforceMinimumSize 1.53<Units.millimetre>
                |> PlantIndividual.zipEnvironment Model.TMax (snd tMax |> TimeSeries.map(fun (v,_) -> v))
            let common =
                shrub 
                |> PlantIndividual.tailGrowth
                |> PlantIndividual.commonTimeline
                |> PlantIndividual.zipEnvironment Model.TMax (snd tMax |> TimeSeries.map(fun (v,_) -> v))
            let startDate = common.Environment.[(code "N").Value].StartDate |> snd
            let startConditions = startValues startDate shrub
            let e = engine |> Bristlecone.withConditioning (Conditioning.Custom startConditions)

            // 2. Setup batches of dependent analyses
            for h in hypotheses do
                for _ in [ 1 .. Config.numberOfReplicates ] do
                    yield
                        async {
                            // A. Compute result
                            let result =
                                Bristlecone.tryFitDendro e Config.endWhen h.Model
                                    Bristlecone.FittingMethod.CumulativeGrowth Model.SR.Code common
                            // B. Save to file
                            match result with
                            | Ok r ->
                                Bristlecone.Data.EstimationResult.saveAll
                                    (fun (s:DatingMethods.Annual) -> sprintf "%i" s.Value)
                                    saveDirectory
                                    s.Identifier.Value
                                    h.ReferenceCode
                                    Config.thinTrace
                                    r

                                return r
                            | Error e -> return failwithf "Error in package: %s" e
                        }
    }

let work = workPackages dataset Model.hypotheses engine Config.resultsDirectory

work
|> Seq.iter (Orchestration.OrchestrationMessage.StartWorkPackage >> orchestrator.Post)

System.Console.ReadLine()