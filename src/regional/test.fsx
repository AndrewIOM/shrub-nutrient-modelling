(**
### Testing the engine and model

Running a full test is strongly recommended. The test will demonstrate if the current
configuration can find known parameters for a model. If this step fails, there is an
issue with either your model, or the Bristlecone configuration.

First, load in the model (defined in the model folder) and open Bristlecone.
*)

#load "model/model.fsx"

open Bristlecone

(**
Next, run the test.
*)

let testSettings =
    Test.create
    |> Test.addNoise (Test.Noise.tryAddNormal "σ[y]" "N")
    |> Test.addNoise (Test.Noise.tryAddNormal "σ[x]" "bs")
    |> Test.addGenerationRules
        [ Test.GenerationRules.alwaysMoreThan -3. "N"
          Test.GenerationRules.alwaysLessThan 20. "N"
          Test.GenerationRules.alwaysMoreThan 0. "bs"
          Test.GenerationRules.monotonicallyIncreasing "x" ] // There must be at least 10mm of wood production
    |> Test.addStartValues [ "x", 5.0; "bs", 5.0<Dendro.mm> |> Allometric.Proxies.toBiomassMM; "N", 3.64 ]
    |> Test.withTimeSeriesLength 30
    |> Test.endWhen (Optimisation.EndConditions.afterIteration 1000)

(*** do-not-eval ***)
let testResult =
    hypotheses
    |> List.map ((fun h -> h.Model) >> Bristlecone.testModel engine testSettings)
