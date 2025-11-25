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

(**
Next, run the test.
*)

let engine: EstimationEngine.EstimationEngine<int<Bristlecone.Time.year>,Time.year,1> =
    Bristlecone.mkContinuous ()
    |> Bristlecone.withContinuousTime Integration.RungeKutta.rk4
    |> Bristlecone.withTimeConversion (fun (i: int<Time.year>) -> float i * 1.<Time.year>)
    |> Bristlecone.withOutput (Logging.Console.logger 10<iteration>)
    |> Bristlecone.withConditioning Conditioning.RepeatFirstDataPoint
        |> Bristlecone.withCustomOptimisation (Optimisation.MonteCarlo.SimulatedAnnealing.fastSimulatedAnnealing 0.01<``optim-space``> true
            { HeatStepLength = Optimisation.EndConditions.atIteration 250<iteration>
              HeatRamp = fun t -> t * 1.10
              BoilingAcceptanceRate = 0.85
              TemperatureCeiling = Some 200.
              InitialTemperature = 1.00
              PreTuneLength = 5000<iteration>
              Tuning = { InitialScale = 0.001
                         TuneLength = 100000<iteration>
                         TuneN = 50<iteration> }
              AnnealStepLength = Optimisation.EndConditions.improvementCount 250 250<iteration> })

let testStartBiomass = Language.Constant 5.0<Dendro.Units.millimetre> |> ShrubModel.Allometry.Proxies.toBiomassMM |> Language.ExpressionCompiler.compileSimple

let testSettings =
    Test.annualSettings
    |> Test.addNoise (Test.Noise.tryAddNormal "σ[y]" "N")
    |> Test.addNoise (Test.Noise.tryAddNormal "σ[x]" "bs")
    |> Test.addGenerationRules
        [ Test.GenerationRules.alwaysMoreThan -3. "N"
          Test.GenerationRules.alwaysLessThan 20. "N"
          Test.GenerationRules.alwaysMoreThan 0. "bs"
          Test.GenerationRules.monotonicallyIncreasing "x" ] // There must be at least 10mm of wood production
    |> Test.addStartValues [ "x", 5.0; "bs", testStartBiomass * 1.</Units.g>; "N", 3.64 ]
    |> Test.withTimeSeriesLength 30
    |> Test.endWhen (Optimisation.EndConditions.atIteration 1000<iteration>)

(*** do-not-eval ***)
let testResult =
    Model.hypotheses
    |> List.map ((fun h -> h.Model) >> Bristlecone.testModel engine testSettings)
