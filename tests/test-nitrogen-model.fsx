#r "../packages/NETStandard.Library.NETFramework/build/net461/lib/netstandard.dll"
#load "../packages/Bristlecone/bristlecone.fsx"
#load "../packages/Bristlecone/charts.fsx"
#load "components/components.fsx"
#r "../packages/Bristlecone.Dendro/lib/netstandard2.0/bristlecone.Dendro.dll"

////////////////////////////////////////////////////
/// Yuribei `Salix lanata` Shrub - Nitrogen Interactions
////////////////////////////////////////////////////

// Shrub ring width modelled with a single
// resource limitation.

open Bristlecone
open Bristlecone.ModelSystem
open Bristlecone.Workflow.Orchestration

// 1. Configure Options
// ----------------------------

module Options =
    let resultsDirectory = "/Users/andrewmartin/Desktop/Bristlecone Results/OptimisationTests/"
    let endWhen = Optimisation.EndConditions.afterIteration 500000 //Optimisation.MonteCarlo.SimulatedAnnealing.EndConditions.stoppedImproving 10
    let logger = Logging.RealTimeTrace.graphWithConsole 30. 10000 //Logging.Console.logger() 
    let engine =
        Bristlecone.mkContinuous 
        |> Bristlecone.withContinuousTime Integration.MathNet.integrate
        |> Bristlecone.withOutput logger
    let orchestrator = OrchestrationAgent(logger, 6)//System.Environment.ProcessorCount)


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
        if (newCumulativeMass - oldCumulativeMass) > 0.
        then newCumulativeMass |> ModelComponents.Proxies.toRadiusMM
        else lastRadius

    { Equations  = [ code "bs",        dbsdt
                     code "N",         dndt ] |> Map.ofList
      Measures   = [ code "x",         stemRadius ] |> Map.ofList
      Parameters = [ // for nitrogen dynamics
                     code "lambda",    parameter PositiveOnly   0.001 0.500   // Rate of nitrogen replenishment
                     code "gamma[n]",  parameter PositiveOnly   0.001 0.200   // Loss rate of nitrogen
                     // for shrub physiology
                     code "gamma[b]",  parameter PositiveOnly   0.001 0.200   // Loss rate of biomass
                     // for likelihood function
                     code "rho",       parameter Unconstrained  -0.50 0.500   // Covariance between growth and nitrogen
                     code "sigma[x]",  parameter PositiveOnly   0.001 0.001   // Standard deviation of x (biomass)
                     code "sigma[y]",  parameter PositiveOnly   0.001 0.001   // Standard deviation of y (nitrogen)
                    ] |> List.append additionalParameters |> Map.ofList
      Likelihood = ModelLibrary.Likelihood.bivariateGaussian "x" "N" }

let hypotheses =

    // [A] N may limited growth via combined N-limitations on (a) photosynthetic and (b) uptake rates
    let limitationModes =
        [ (fun p -> ModelComponents.GrowthLimitation.hollingDiscModel ((p |> Pool.getEstimate "a") / 1000.) ((p |> Pool.getEstimate "r") * 1000.) (p |> Pool.getEstimate "h")),
           [ code "a",      parameter PositiveOnly   0.100 4.000          // N-uptake efficiency
             code "h",      parameter PositiveOnly   0.100 4.000
             code "r",      parameter PositiveOnly   0.100 0.500 ]      // N-handling time (including uptake and incorporation)
          (fun p -> ModelComponents.GrowthLimitation.linear ((p |> Pool.getEstimate "a") / 1000.)), 
           [ code "a",      parameter PositiveOnly   0.001 0.010
             code "r",      parameter PositiveOnly   0.100 1.000 ]        // N-uptake efficiency
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

// 2. Test using fake data
// ----------------------------
let startValues = [ ShortCode.create "x", 5.; ShortCode.create "N", 3.64; ShortCode.create "bs", (5. |> ModelComponents.Proxies.toBiomassMM)] |> Map.ofList
let generationRules = 
    [ code "bs", fun data -> data |> Seq.min > 0.
      code "N", fun data -> data |> Seq.min > -3.0
      code "x", fun data -> data |> Seq.pairwise |> Seq.sumBy (fun (a,b) -> b - a) > 10.   // There must be at least 10mm of wood production
      code "N",  fun data -> data |> Seq.max < 20. ]                                        // N must not get to levels above 20 units

let noNoise p data =
    data

let addNoise p data =
    let random = System.Random()
    data |> Map.map(fun key value -> 
        let sigma =
            match key with
            | k when k = code "N" -> p |> Pool.getEstimate "sigma[y]"
            | k when k = code "bs" -> p |> Pool.getEstimate "sigma[x]"
            | _ -> invalidOp "No sigma for this!"
        let draw = Bristlecone.Statistics.Distributions.Normal.draw random 0. sigma
        value |> TimeSeries.map (fun (x,_) -> x + draw()))

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

// let vonBertalanffyRules = 
//     [ code "x", (fun data -> 
//         let cumulative = data |> Seq.pairwise |> Seq.sumBy (fun (a,b) -> b - a)
//         cumulative > 10. && cumulative < 10000.) ]

// B. Define model system
let simpleGrowthModel =
    { Equations  = [ code "x", vonBertalanffy ] |> Map.ofList
      Likelihood = ModelLibrary.Likelihood.sumOfSquares ["x"]
      Measures   = [] |> Map.ofList
      Parameters = [ code "eta",    parameter Unconstrained   0.001 0.50
                     code "beta",   parameter Unconstrained   0.001 0.50
                     code "kappa",  parameter Unconstrained   0.001 0.50 ] |> Map.ofList }



/// Setup test cases
/// __________________
/// 
/// We want to test the following:
/// A. Simulated Annealing (Cauchy, temp-independent)
/// B. Single Amoeba
/// C. Swarm Amoeba
/// D. Adaptive Metropolis
/// E. Tuned Random Walk?

let testEngines = [
    ("CauchyAnnealing", Options.engine |> Bristlecone.withCustomOptimisation (Optimisation.MonteCarlo.SimulatedAnnealing.fastSimulatedAnnealing 0.0005 { Optimisation.MonteCarlo.SimulatedAnnealing.AnnealSettings.Default with HeatRamp = (fun t -> t * 1.05); HeatStepLength = Optimisation.EndConditions.afterIteration 1000; AnnealStepLength = Optimisation.EndConditions.afterIteration 5000  }))
    //("AdaptiveMetropolis", Options.engine |> Bristlecone.withCustomOptimisation (Optimisation.MonteCarlo.adaptiveMetropolis 0.750 500))
    //("TunedMCMC", Options.engine |> Bristlecone.withTunedMCMC [ Optimisation.MonteCarlo.TuneMethod.CovarianceWithScaleTotalHistory 0.250, 1000, Optimisation.EndConditions.afterIteration 20000 ])
    //("SingleAmoeba", Options.engine |> Bristlecone.withGradientDescent)
]

let models = [
    //("NitrogenSaturatingWithFeedback", hypotheses.[0], noNoise, generationRules)
    //("BerfalanffyGrowthModel-Noise", simpleGrowthModel, addNoise, vonBertalanffyRules)
    //("BerfalanffyGrowthModel", simpleGrowthModel, noNoise, vonBertalanffyRules)
    ("NitrogenSaturatingWithFeedback-Noise", hypotheses.[0], addNoise, generationRules)
]

let workPackages models engines saveDirectory =
    seq {
        for (modelName,model,noise,rules) in models do
            for (engineName,engine) in engines do
                yield async {
                    // Returns result, original series, and original parameters
                    let result,originalSeries,originalTheta =
                        Bristlecone.testModel engine 50 startValues Options.endWhen rules noise model

                    // Save trace of optimisation
                    let jobId = System.Guid.NewGuid()
                    Bristlecone.Data.Cache.saveTrace saveDirectory engineName 1 jobId [ result ]

                    // Save expected and observed parameters of optimisation
                    let paramFileName = sprintf "%sparameters-%s-%s-%s.csv" saveDirectory engineName modelName (jobId.ToString())
                    System.IO.File.WriteAllLines(paramFileName, (originalTheta |> Seq.map(fun x -> sprintf "estimated,%s,%f" x.Key.Value (x.Value |> Parameter.getEstimate))))
                    System.IO.File.AppendAllLines(paramFileName, (result.Parameters |> Seq.map(fun x -> sprintf "actual,%s,%f" x.Key.Value (x.Value |> Parameter.getEstimate))))

                    // Save expected and observed parameters of optimisation
                    let seriesFileName = sprintf "%sseries-%s-%s-%s.csv" saveDirectory engineName modelName (jobId.ToString())
                    System.IO.File.WriteAllLines(seriesFileName, (result.Series |> Seq.collect(fun x -> x.Value.Expected |> Seq.mapi(fun i v -> sprintf "estimated,%s,%f,%i" x.Key.Value v i))))
                    System.IO.File.AppendAllLines(seriesFileName, (result.Series |> Seq.collect(fun x -> x.Value.Observed |> Seq.mapi(fun i v -> sprintf "actual,%s,%f,%i" x.Key.Value v i))))

                    // TODO Set up exporting of data here in Bristlecone, rather than a seperate script
                    return result
                }
    }

// Orchestrate the analyses
let work = workPackages models testEngines Options.resultsDirectory
work |> Seq.iter (OrchestrationMessage.StartWorkPackage >> Options.orchestrator.Post)

for w in work do w |> Async.RunSynchronously
