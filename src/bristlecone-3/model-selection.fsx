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

#r "/Users/andrewmartin/Documents/GitHub Projects/bristlecone/src/Bristlecone/bin/Debug/net5.0/Bristlecone.dll"
#r "/Users/andrewmartin/Documents/GitHub Projects/bristlecone/src/Bristlecone.Dendro/bin/Debug/net5.0/Bristlecone.Dendro.dll"

open Bristlecone // Opens Bristlecone core library and estimation engine
open Bristlecone.Language // Open the language for writing Bristlecone models
open Bristlecone.Time

#load "model/model.fsx"

(**
### Model selection

After calling `run ()`, we should have a folder containing csv files, three
for each `EstimationResult`. We can load all of these in at once to
calculate model comparison statistics.

Functions for loading data are in the `Bristlecone.Data` namespace. You
may notice that in the work package definition above we used functions from
this namespace to save the results.
*)

let saveDiagnostics () =

    // 1. Get all results sliced by plant and hypothesis
    let results =
        let get (subject: Dendro.PlantIndividual.PlantIndividual) (hypothesis: Hypotheses.Hypothesis) =
            Bristlecone.Data.EstimationResult.loadAll
                Config.resultsDirectory
                subject.Identifier.Value
                hypothesis.Model
                hypothesis.ReferenceCode

        Bristlecone.ModelSelection.ResultSet.arrangeResultSets dataset hypotheses get

    // 2. Save convergence statistics to file
    results
    |> Diagnostics.Convergence.gelmanRubinAll
        10000
        (fun (s: PlantIndividual.PlantIndividual) -> s.Identifier.Value)
        (fun (h: Hypotheses.Hypothesis) -> h.ReferenceCode)
    |> Data.Convergence.save Config.resultsDirectory

    // 3. Save Akaike weights to file
    results
    |> ModelSelection.Akaike.akaikeWeightsForSet (fun (h: Hypotheses.Hypothesis) -> h.ReferenceCode)
    |> Seq.map (fun (x, a, b, c) -> x.Identifier.Value, a, b, c)
    |> Data.ModelSelection.save Config.resultsDirectory


    // // 4. Save out logged components
    // results
    // |> Seq.map(fun r ->
    //    Diagnostics.ModelComponents.calculateComponents fit engine r)

    // 5. One-step ahead predictions

    let bestFits =
        Seq.allPairs dataset hypotheses
        |> Seq.map (fun (s, h) ->
            s, h, Bristlecone.Data.MLE.loadBest Config.resultsDirectory s.Identifier.Value h.Model h.ReferenceCode)

    let oneStepPredictions =
        bestFits
        |> Seq.map (fun (s, h, mle) ->

            // 0. Convert x into biomass
            let preTransform (data: CodedMap<TimeSeries<float>>) =
                data
                |> Map.toList
                |> List.collect (fun (k, v) ->
                    if k.Value = "x" then
                        [ (k, v)
                          ((code "bs").Value,
                           v |> TimeSeries.map (fun (x, _) -> x * 1.<mm> |> Allometric.Proxies.toBiomassMM)) ]
                    else
                        [ (k, v) ])
                |> Map.ofList

            // 1. Arrange the subject and settings (same as in model-fitting)
            let shrub = s |> PlantIndividual.toCumulativeGrowth
            let common = shrub |> PlantIndividual.keepCommonYears
            let startDate = (common.Environment.[(code "N").Value]).StartDate |> snd
            let startConditions = startValues startDate shrub

            let e =
                engine
                |> Bristlecone.withConditioning (Bristlecone.Conditioning.Custom startConditions)

            let result =
                Bristlecone.oneStepAhead e h.Model preTransform (Bristlecone.fromDendro common) (mle |> snd |> snd)

            // Save each n-step ahead result to a csv file
            Bristlecone.Data.NStepAhead.save
                Config.resultsDirectory
                s.Identifier.Value
                h.ReferenceCode
                (mle |> fst)
                1
                result

            s.Identifier.Value, h.ReferenceCode, result)

    Bristlecone.Data.NStepAhead.saveAllRMSE Config.resultsDirectory oneStepPredictions

    ()
