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
open Bristlecone.Dendro
open Bristlecone.Dendro.PlantIndividual
open Bristlecone.Data
open Bristlecone.Workflow.Orchestration

// 1. Configure Options
// ----------------------------

module Options =
    let resultsDirectory = "/Users/andrewmartin/Desktop/Bristlecone Results/YuribeiAnnual-Filzbach"
    let chains = 3
    let endWhen = Optimisation.EndConditions.afterIteration 100000
    let logger = Logging.RealTimeTrace.graphWithConsole 30. 10000 //Logging.Console.logger()
    let engine =
        Bristlecone.mkContinuous 
        |> Bristlecone.withContinuousTime Integration.MathNet.integrate
        |> Bristlecone.withOutput logger
        |> Bristlecone.withCustomOptimisation (Optimisation.MonteCarlo.SimulatedAnnealing.fastSimulatedAnnealing 0.0001 false
            { Optimisation.MonteCarlo.SimulatedAnnealing.AnnealSettings<float>.Default with 
                BoilingAcceptanceRate = 0.85
                HeatRamp = (fun t -> t + sqrt t); TemperatureCeiling = Some 500.
                HeatStepLength = Optimisation.EndConditions.afterIteration 1000
                AnnealStepLength = (fun x -> Optimisation.MonteCarlo.SimulatedAnnealing.EndConditions.improvementCount 5000 250 x || Optimisation.EndConditions.afterIteration 10000 x) }) //(Optimisation.EndConditions.afterIteration 10000) })

    let orchestrator = OrchestrationAgent(logger, System.Environment.ProcessorCount, false)


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
        dbsdt' bs ((e.[code "N"]) |> ModelComponents.Proxies.d15NtoAvailability)
            (p |> Pool.getEstimate "gamma[b]") ((p |> Pool.getEstimate "r")) (maxGrowthRate p) (nLimitation p)

    /// Bristlecone function for dN/dt
    let dndt p _ n (e:Environment) =
        dndt' (e.[code "bs"]) (n |> ModelComponents.Proxies.d15NtoAvailability) 
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
                     code "sigma[x]",  parameter PositiveOnly   0.001 0.100   // Standard deviation of x (biomass)
                     code "sigma[y]",  parameter PositiveOnly   0.001 0.100   // Standard deviation of y (nitrogen)
                    ] |> List.append additionalParameters |> Map.ofList
      Likelihood = ModelLibrary.Likelihood.bivariateGaussian "x" "N" }

let hypotheses =

    // [A] N may limited growth via combined N-limitations on (a) photosynthetic and (b) uptake rates
    let limitationModes =
        [ (fun p -> ModelComponents.GrowthLimitation.hollingDiscModelDual ((p |> Pool.getEstimate "a") / 1000.)  ((p |> Pool.getEstimate "r") * 1000.) (p |> Pool.getEstimate "h") 10.),
           [ code "a",      parameter PositiveOnly   0.100 0.400          // N-uptake efficiency
             code "h",      parameter PositiveOnly   0.100 0.400
             code "r",      parameter PositiveOnly   0.500 1.000 ]     // N-handling time (including uptake and incorporation)
          (fun p -> ModelComponents.GrowthLimitation.linear ((p |> Pool.getEstimate "a") / 1000.) 10.), 
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
// let startValues = [ code "x", 5.; code "N", 3.64; code "bs", (5. |> ModelComponents.Proxies.toBiomassMM)] |> Map.ofList
// let generationRules = 
//     [ code "bs", fun data -> data |> Seq.min > 0.
//       code "N", fun data -> data |> Seq.min > -3.0
//       code "x", fun data -> data |> Seq.pairwise |> Seq.sumBy (fun (a,b) -> b - a) > 10.   // There must be at least 10mm of wood production
//       code "N",  fun data -> data |> Seq.max < 20. ]                                        // N must not get to levels above 20 units

// let addNoise p data =
//     let random = System.Random()
//     data |> Map.map(fun key value -> 
//         let sigma =
//             match key with
//             | k when k = code "N" -> p |> Pool.getEstimate "sigma[y]"
//             | k when k = code "bs" -> p |> Pool.getEstimate "sigma[x]"
//             | _ -> invalidOp "No sigma for this!"
//         let draw = Bristlecone.Statistics.Distributions.Normal.draw random 0. sigma
//         value |> TimeSeries.map (fun (x,_) -> x + draw()))

// let testResult =
//     hypotheses
//     |> Seq.head
//     |> Bristlecone.testModel Options.engine 30 startValues Options.endWhen generationRules addNoise


// 3. Load Real Data and Estimate
// ----------------------------

let shrubs = 
    let yuribei = Data.PlantIndividual.loadRingWidths (__SOURCE_DIRECTORY__ + "/../data/yamal-rw.csv")
    let d15N = Data.PlantIndividual.loadLocalEnvironmentVariable (__SOURCE_DIRECTORY__ + "/../data/yuribei-d15N-imputed.csv")
    yuribei
    |> Seq.map (fun s -> (s.Identifier.Value, s))
    |> Seq.keyMatch d15N
    |> Seq.map (fun (_,plant,d15N) -> PlantIndividual.zipEnv (code "N") plant d15N)
    |> Seq.toList

let getStartValues (startDate:System.DateTime) (plant:PlantIndividual) =
    let initialRadius =
        match plant.Growth with
        | RingWidth s -> 
            match s with
            | Absolute c -> c.Head |> fst |> removeUnit
            | Cumulative c -> 
                let start = (c |> TimeSeries.trimStart (startDate - System.TimeSpan.FromDays(366.))).Values |> Seq.head |> removeUnit
                // printfn "Start cumulative growth = %f" start
                start
            | Relative _ -> invalidOp "Not implemented"
        | _ -> invalidOp "Not implemented 2"
    let initialMass = initialRadius |> removeUnit |> ModelComponents.Proxies.toBiomassMM
    let initialNitrogen = plant.Environment.[code "N"].Head |> fst
    [ (code "x", initialRadius)
      (code "N", initialNitrogen)
      (code "bs", initialMass) ] |> Map.ofList

let workPackages shrubs hypotheses engine saveDirectory =
    seq {
        for s in shrubs do

            // 1. Arrange the subject and settings
            let shrub = s |> PlantIndividual.toCumulativeGrowth
            let common = shrub |> PlantIndividual.keepCommonYears
            let startDate = (common.Environment.[code "N"]).StartDate |> snd
            let startConditions = getStartValues startDate shrub
            let e = engine |> Bristlecone.withConditioning (Custom startConditions)

            // 2. Setup batches of dependent analyses
            for h in [ 1 .. hypotheses |> List.length ] do
                for _ in [ 1 .. Options.chains ] do
                    yield async {
                            // A. Compute result
                            let result = Bristlecone.PlantIndividual.fit e Options.endWhen hypotheses.[h-1] common |> fst
                            // B. Save to file
                            Bristlecone.Data.EstimationResult.saveAll saveDirectory s.Identifier.Value h 1 result
                            return result }
    }

// Orchestrate the analyses
let work = workPackages shrubs hypotheses Options.engine Options.resultsDirectory
let run() = work |> Seq.iter (OrchestrationMessage.StartWorkPackage >> Options.orchestrator.Post)
