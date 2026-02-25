(**
# Long-Term N-Shrub relations in Yamal, Russia

The following code demonstrates model-fitting 
and model-selection for the mechanisms controlling nutrient 
limitation to shrub individuals in the Arctic tundra.

This script completes model-selection on model fits
obtained by running `model-fitting.fsx`
for 24 shrub individuals.

First, we must load Bristlecone:
*)

#load "model-fitting.fsx"

open Bristlecone // Opens Bristlecone core library and estimation engine
open Bristlecone.Language // Open the language for writing Bristlecone models
open Bristlecone.Time

(**
### Model selection

After calling `run ()`, we should have a folder containing csv files, three
for each `EstimationResult`. We can load all of these in at once to
calculate model comparison statistics.

Functions for loading data are in the `Bristlecone.Data` namespace. You
may notice that in the work package definition above we used functions from
this namespace to save the results.
*)

module Settings =
    let resultsDirectory = "/Users/andrewmartin/Desktop/Bristlecone-3.0/Regional/"


type PlantAnnual = Dendro.PlantIndividual.PlantIndividual<Dendro.Units.millimetre,DatingMethods.Annual,int<year>,int<year>>
type Hypothesis = Hypotheses.Hypothesis<year>

let saveDiagnostics () =

    // 1. Get all results sliced by plant and hypothesis
    let results =
        let get (subject: PlantAnnual) (hypothesis: Hypothesis) =
            Data.EstimationResult.loadAll
                toSeries
                Settings.resultsDirectory
                subject.Identifier.Value
                hypothesis.Model
                hypothesis.ReferenceCode

        ModelSelection.ResultSet.arrangeResultSets ``Model-fitting``.dataset Model.hypotheses get

    // 2. Save convergence statistics to file
    results
    |> Diagnostics.Convergence.gelmanRubinAll
        10000
        (fun (s: PlantAnnual) -> s.Identifier.Value)
        (fun (h: Hypothesis) -> h.ReferenceCode)
    |> Data.Convergence.save Settings.resultsDirectory

    // 3. Save Akaike weights to file
    results
    |> ModelSelection.Akaike.akaikeWeightsForSet (fun (h: Hypothesis) -> h.ReferenceCode)
    |> Seq.map (fun (x, a, b, c) -> x.Identifier.Value, a, b, c)
    |> Data.ModelSelection.save Settings.resultsDirectory

    // 4. Save out logged components
    // --> TODO Update in Bristlecone to use here.

    // 5. One-step ahead predictions
    let oneStepPredictions =
        results
        |> Seq.filter(fun r -> r.BestResult.IsSome)
        |> Seq.map (fun resultSet ->

            let prepared = ``Model-fitting``.preTransform resultSet.Subject
            let oneStepResult = Bristlecone.oneStepAheadDendro prepared.Engine resultSet.Hypothesis.Model Bristlecone.FittingMethod.CumulativeGrowth Model.SR.Code resultSet.Subject id resultSet.BestResult.Value.Parameters

            // Save each n-step ahead result to a csv file
            Bristlecone.Data.NStepAhead.save
                Settings.resultsDirectory
                resultSet.Subject.Identifier.Value
                resultSet.Hypothesis.ReferenceCode
                resultSet.BestResult.Value.ResultId
                1
                oneStepResult

            resultSet.Subject.Identifier.Value, resultSet.Hypothesis.ReferenceCode, oneStepResult)

    Bristlecone.Data.NStepAhead.saveAllRMSE Settings.resultsDirectory oneStepPredictions

    ()
