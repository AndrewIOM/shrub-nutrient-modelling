#r "../packages/NETStandard.Library.NETFramework/build/net461/lib/netstandard.dll"
#load "../packages/Bristlecone/bristlecone.fsx"
#load "../packages/Bristlecone/charts.fsx"
#load "components/components.fsx"
#r "../packages/Bristlecone.Dendro/lib/netstandard2.0/bristlecone.Dendro.dll"

/////////////////////////////////////
/// Bristlecone Model Tests
/////////////////////////////////////

(* To ensure Bristlecone is effective at solving long-term ecological
   problems, we run a model-fitting-model-selection (MFMS) procedure
   for two ecological models. First, simple plant growth models are fit
   to plant growth data. Second, growth-resource models are fit.  *)

open Bristlecone
open Bristlecone.ModelSystem

// Test 1. Simple growth model (no noise)
// _____________________________

// A. Define a logistic growth model 
let vonBertalanffy' eta beta kappa mass =
    eta * mass ** beta - kappa * mass

let vonBertalanffy p t x environment =
    vonBertalanffy' 
      (p |> Pool.getEstimate "eta") 
      (p |> Pool.getEstimate "beta") 
      (p |> Pool.getEstimate "kappa") x

// B. Define model system
let hypothesis =
    { Equations  = [ code "x", vonBertalanffy ] |> Map.ofList
      Likelihood = ModelLibrary.Likelihood.sumOfSquares ["x"]
      Measures   = [] |> Map.ofList
      Parameters = [ code "eta",    parameter Unconstrained   0.001 1.00
                     code "beta",   parameter Unconstrained   0.001 1.00
                     code "kappa",  parameter Unconstrained   0.001 1.00 ] |> Map.ofList }

// C. Define fitting method
let logger = 
  let consolePost = Bristlecone.Logging.Console.logger()
  let graphLog = Bristlecone.Logging.RealTimeTrace.TraceGraph(Logging.Device.X11,20.,10000)
  (fun event -> consolePost event; graphLog.Log event)

open Optimisation.MonteCarlo.SimulatedAnnealing

let settings = {
      HeatStepLength = Optimisation.EndConditions.afterIteration 1000
      HeatRamp = fun t -> t * 1.10
      BoilingAcceptanceRate = 0.85
      InitialTemperature = 1.00
      AnnealStepLength = EndConditions.improvementCount 10000 100
      TemperatureCeiling = Some 500. }

let method =
      Bristlecone.mkContinuous
      |> Bristlecone.withCustomOptimisation (Optimisation.MonteCarlo.SimulatedAnnealing.fastSimulatedAnnealing 0.001 settings)
      |> Bristlecone.withOutput logger

// D. Define test methodology
let startValues = [code "x", 5.] |> Map.ofList
let generationRules = 
    [ code "x", fun (data:seq<float>) -> data |> Seq.pairwise |> Seq.sumBy (fun (a,b) -> b - a) < 50. ]
let endCondition = Optimisation.EndConditions.afterIteration 5000
let noNoise p x = x

// Test 2. Simple growth model (gaussian noise)
// _____________________________

let gaussianNegativeLogLikelihood key p data =
      let n = data |> Seq.length |> float
      let sos = ModelLibrary.Likelihood.sumOfSquares [key] p data
      let sigma = p |> Pool.getEstimate "sigma"
      -(n / 2.) * log ( 2. * System.Math.PI) - (n / 2.) * log (sigma**2.) - (1. / 2. * sigma ** 2.) * sos

let random = MathNet.Numerics.Random.MersenneTwister()

// Adds gaussian noise to all variables
let addNoise p data =
      let sigma = 0.43 //p |> Pool.getEstimate "sigma"
      let draw = Bristlecone.Statistics.Distributions.Normal.draw random 0. sigma
      data
      |> Map.map(fun key value -> value |> TimeSeries.map (fun (x,y) -> x + draw()))

let hypothesisNoisy =
    { Equations  = [ code "x", vonBertalanffy ] |> Map.ofList
      Measures   = [] |> Map.ofList
      Likelihood = gaussianNegativeLogLikelihood "x"
      Parameters = [ code "eta",    parameter PositiveOnly   0.001 1.00
                     code "beta",   parameter PositiveOnly   0.001 1.00
                     code "kappa",  parameter PositiveOnly   0.001 1.00
                     code "sigma",  parameter PositiveOnly   0.001 0.01 ] |> Map.ofList }


/// Predator-Prey Model
/// ________________________________

let ``population with single resource limitation`` =

    let dxdt' x y alpha beta =
        alpha * x - beta * x * y

    let dydt' y x beta delta gamma =
        delta * beta * x * y - gamma * y

    let dxdt p _ x (e:Environment) =
        dxdt' x (lookup e "y") (p |> Pool.getEstimate "alpha") (p |> Pool.getEstimate "beta")
    
    let dydt p _ y (e:Environment) =
        dydt' y (lookup e "x") (p |> Pool.getEstimate "beta") (p |> Pool.getEstimate "delta") (p |> Pool.getEstimate "gamma")

    { Equations =  [ ShortCode.create "x", dxdt
                     ShortCode.create "y", dydt ] |> Map.ofList
      Measures = [] |> Map.ofList
      Parameters = [ ShortCode.create "alpha",  Parameter.create Unconstrained 1.00 1.001
                     ShortCode.create "beta",   Parameter.create Unconstrained 0.20 0.201
                     ShortCode.create "delta",  Parameter.create Unconstrained 0.50 0.501
                     ShortCode.create "gamma",  Parameter.create Unconstrained 0.20 0.201 ] |> Map.ofList
                     //ShortCode.create "sigmax", Parameter.create Unconstrained -0.50 0.50
                     //ShortCode.create "sigmay", Parameter.create Unconstrained -0.50 0.50 
                     //ShortCode.create "rho",    Parameter.create Unconstrained -0.25 0.25 ] |> Map.ofList
      Likelihood = ModelLibrary.Likelihood.sumOfSquares ["x";"y"] } //ModelLibrary.Likelihood.bivariateGaussian "x" "y" }


// Run Tests
// ______________________________

// 1. Simple growth
let test1 =
    hypothesis
    |> Bristlecone.testModel method 30 startValues endCondition generationRules noNoise

// 2. Growth with noise
let test2 =
      hypothesis
      |> Bristlecone.testModel method 30 startValues endCondition generationRules addNoise

// 3. Population with noise
let test3 =
      let startValues = [ ShortCode.create "x", 1.00; ShortCode.create "y", 2.00 ] |> Map.ofList
      ``population with single resource limitation``
      |> Bristlecone.testModel method 30 startValues endCondition generationRules noNoise

