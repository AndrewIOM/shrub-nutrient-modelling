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
*)

#r "nuget: Bristlecone.Dendro,2.0.0"

open Bristlecone // Opens Bristlecone core library and estimation engine
open Bristlecone.Language // Open the language for writing Bristlecone models
open Bristlecone.Time

(**
### Loading in the model

We can simply load the fsx file that we previously defined the
model system within.
*)

#load "model/model.fsx"


(**
### Setting up a *Bristlecone engine*

A bristlecone engine provides a fixed setup for estimating parameters from data.
We use the same engine for all model fits within a single study.

Here, we scaffold an engine from `Bristlecone.mkContinuous`, as we are working
with continuous-time models.
*)

let output = Logging.Console.logger 1000

let engine =
    Bristlecone.mkContinuous
    |> Bristlecone.withContinuousTime Integration.MathNet.integrate
    |> Bristlecone.withOutput output
    |> Bristlecone.withConditioning Conditioning.RepeatFirstDataPoint
    |> Bristlecone.withCustomOptimisation (
        Optimisation.MonteCarlo.SimulatedAnnealing.fastSimulatedAnnealing
            0.01
            true
            { Optimisation.MonteCarlo.SimulatedAnnealing.AnnealSettings<float>.Default with
                InitialTemperature = 100.
                TemperatureCeiling = Some 100.
                HeatRamp = (fun t -> t + 5.00)
                BoilingAcceptanceRate = 0.85
                HeatStepLength = Optimisation.EndConditions.afterIteration 1000
                TuneLength = 1000
                AnnealStepLength =
                    (fun x n ->
                        Optimisation.MonteCarlo.SimulatedAnnealing.EndConditions.improvementCount 5000 250 x n
                        || Optimisation.EndConditions.afterIteration 10000 x n) }
    )

(**
### Load in the real dataset

Here, we are using the Bristlecone.Dendro package to 
read in dendroecological data. For other problems, any
method to wrangle the data into a `TimeSeries` is acceptable.

*)

open Bristlecone.Dendro

let dataset =
    let plants =
        Data.PlantIndividual.loadRingWidths (__SOURCE_DIRECTORY__ + "/../data/yamal-rw.csv")

    let isotopeData =
        Data.PlantIndividual.loadLocalEnvironmentVariable (__SOURCE_DIRECTORY__ + "/../data/yuribei-d15N-imputed.csv")

    plants |> PlantIndividual.zipEnvMany "N" isotopeData

(**
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

let orchestrator =
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

let startValues (startDate: System.DateTime) (plant: PlantIndividual.PlantIndividual) =
    let removeUnit (x: float<_>) = float x

    let initialRadius =
        match plant.Growth with
        | PlantIndividual.PlantGrowth.RingWidth s ->
            match s with
            | GrowthSeries.Cumulative c ->
                let trimmed = c |> TimeSeries.trimStart (startDate - System.TimeSpan.FromDays(366.))

                match trimmed with
                | Some t -> t.Values |> Seq.head
                | None -> failwith "Could not get t0 from ring-width series"
            | _ -> invalidOp "Not applicable"
        | _ -> invalidOp "Not applicable"

    let initialMass = initialRadius |> Allometric.Proxies.toBiomassMM
    let initialNitrogen = plant.Environment.[(code "N").Value].Head |> fst

    [ ((code "x").Value, initialRadius |> removeUnit)
      ((code "N").Value, initialNitrogen)
      ((code "bs").Value, initialMass) ]
    |> Map.ofList

(**
Next, we set up some configuration settings, then wrap up our hypotheses and shrubs
as work packages for the orchestrator, where a work package is a `Async<ModelSystem.EstimationResult>`.

The function to make the work packages first sets up the common time period and custom
start values for each shrub, then creates per-hypothesis work packages for that shrub.
*)

module Config =

    let numberOfReplicates = 3
    let resultsDirectory = "~/Desktop/Bristlecone Results Clean/Regional/"
    let thinTrace = Some 100
    let endWhen = Optimisation.EndConditions.afterIteration 100000


// Function to scaffold work packages
let workPackages shrubs (hypotheses: Hypotheses.Hypothesis list) engine saveDirectory =
    seq {
        for s in shrubs do

            // 1. Arrange the subject and settings
            let shrub = s |> PlantIndividual.toCumulativeGrowth
            let common = shrub |> PlantIndividual.keepCommonYears
            let startDate = (common.Environment.[(code "N").Value]).StartDate |> snd
            let startConditions = startValues startDate shrub
            let e = engine |> Bristlecone.withConditioning (Conditioning.Custom startConditions)

            // 2. Setup batches of dependent analyses
            for h in hypotheses do
                for _ in [ 1 .. Config.numberOfReplicates ] do
                    yield
                        async {
                            // A. Compute result
                            let result =
                                Bristlecone.tryFit e Config.endWhen (Bristlecone.fromDendro common) h.Model
                            // B. Save to file
                            match result with
                            | Ok r ->
                                Bristlecone.Data.EstimationResult.saveAll
                                    saveDirectory
                                    s.Identifier.Value
                                    h.ReferenceCode
                                    Config.thinTrace
                                    r

                                return r
                            | Error e -> return failwithf "Error in package: %s" e
                        }
    }

let work = workPackages dataset hypotheses engine Config.resultsDirectory

// Calling this function will run all of the analyses, but for this
// example script we won't run them here.
let run () =
    work
    |> Seq.iter (Orchestration.OrchestrationMessage.StartWorkPackage >> orchestrator.Post)
