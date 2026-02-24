(**
### Testing the engine and model

Running a full test is strongly recommended. The test will demonstrate if the current
configuration can find known parameters for a model. If this step fails, there is an
issue with either your model, or the Bristlecone configuration.

First, load in the model (defined in the model folder) and open Bristlecone.
*)

#load "model/units.fsx"
#load "model/model.fsx"

open Bristlecone
open Bristlecone.Time
open Bristlecone.Language

(**
Next, we set some common configuration in the below Settings module, and repeat
the same engine definition used within the model-fitting step.
*)

module Settings =
    let output = Logging.Console.logger 100<iteration>
    let endCondition = Optimisation.EndConditions.Profiles.mcmc 100<iteration> output


let engine: EstimationEngine.EstimationEngine<DatingMethods.Annual,int<year>,year,1> =
    Bristlecone.mkContinuous ()
    |> Bristlecone.withContinuousTime Integration.RungeKutta.rk4
    |> Bristlecone.withTimeConversion DateMode.Conversion.Annual.toYears
    |> Bristlecone.withOutput Settings.output
    |> Bristlecone.withConditioning Conditioning.RepeatFirstDataPoint
    |> Bristlecone.withBristleconeOptimiser

(**
To configure the test, we use the built-in helpers within Bristlecone,
which exist as composable functions within the Test.* module (note:
you must open Bristlecone.Language to use these).
*)

let testStartBiomass = Constant 5.0<Dendro.Units.millimetre> |> ShrubModel.Allometry.Proxies.toBiomassMM |> ExpressionCompiler.compileSimple
let testStartN = Constant 3.64 |> Model.``δ15N -> N availability`` |> ExpressionCompiler.compileSimple

let mkTemperature =
    let meanT = 10.0<Dendro.Units.celsius>
    let phi, sigma = 0.5, 1.0<Dendro.Units.celsius>
    let series =
        Test.Synthetic.ar1 phi sigma engine.Random
        |> Seq.map (fun a -> meanT + a)
        |> Seq.cache
    fun (t: float<year>) ->
        let i = t |> int
        series |> Seq.item i |> Dendro.Units.celsiusToKelvin

let testSettings =
    Test.annualSettings
    |> Test.seriesLength 30
    |> Test.t1 (Require.measure Model.SR) 5.0<Units.mm>
    |> Test.t1 (Require.state Model.B) testStartBiomass
    |> Test.t1 (Require.state Model.N) testStartN
    |> Test.rule (Require.state Model.N) (Test.GenerationRules.between -3 20.)
    |> Test.rule (Require.state Model.B) (Test.GenerationRules.alwaysMoreThan (float testStartBiomass - 0.01)) // TODO Units fix in Bristlecone.
    |> Test.rule (Require.measure Model.SR) Test.GenerationRules.monotonicallyIncreasing
    |> Test.withObservationError (Require.state Model.N) (Test.Error.normal Model.σN)
    |> Test.withObservationError (Require.measure Model.SR) (Test.Error.normal Model.σSR)
    |> Test.withEnvironmentGenBySpan Model.TMax mkTemperature (fun ts -> ts |> float |> (*) 1.<year>)

(**
Here, we have set up the test such that it requires N state to stay within realistic bounds,
requires that wood is produced, and that the biomass of the plant is always more than just below the
starting biomass.

We set the initial hidden biomass as the biomass indicated by the allometric equations at 5mm
stem radius, set initial stem radius to 5mm, and initial N to d15N = 3.64.

We add observation error around the time-series in line with the configured starting bounds
for sigma N and sigma Bs respectively.

Finally, we may run a test for any particular hypothesis:
*)

let hypothesisToTest = Model.hypotheses.[4]

(*** do-not-eval ***)
let testResult =
    Bristlecone.testModel engine Settings.endCondition testSettings hypothesisToTest.Model

(**
Inspection of the test outputs will indicate if the model is identifiable given
the engine, test settings, and model setup.
*)
