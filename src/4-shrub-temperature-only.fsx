#r "../packages/NETStandard.Library.NETFramework/build/net461/lib/netstandard.dll"
#load "../packages/Bristlecone/bristlecone.fsx"
#load "../packages/Bristlecone/charts.fsx"
#load "components/components.fsx"
#load "components/temperature.fsx"
#r "../packages/Bristlecone.Dendro/lib/netstandard2.0/bristlecone.Dendro.dll"

////////////////////////////////////////////////////
/// Yuribei `Salix lanata` Shrub - Nitrogen Interactions
////////////////////////////////////////////////////

// Shrub ring width modelled with a single
// resource limitation.

open Bristlecone
open Bristlecone.ModelSystem
open Bristlecone.PlantIndividual
open Bristlecone.Workflow.Orchestration

// 1. Configure Options
// ----------------------------

module Allometry =

    open ModelComponents

    // Empirically-derived parameters:
    let k5 = 19.98239 // Allometric fit to Yamal shrub BD-length data #1 (in centimetres)
    let k6 = 0.42092 // Allometric fit to Yamal shrub BD-length data #2 (in centimetres)

    // Constants from the literature:
    let a = 2. // the number of child branches added to previous branches (including the tops of the stems) for a shrub
    let p = 0.5 // the length of a child branch as a proportion of its parent branch/stem
    let lmin = 20. //cm. the length at which a stem or branch gets child branches
    let rtip = 0.1 //cm. the radius of the outermost tip of a stem or branch
    let b = 0.0075 // the ratio of the basal radius of a stem or branch and its length
    let salixWoodDensity = 0.5 // g / cm3 (from internet)
    let numberOfStems = 1.

    /// Radius in millimetres
    let toBiomassMM radiusMM = 
        radiusMM / 10. |> Allometrics.shrubBiomass b a rtip p lmin k5 k6 numberOfStems salixWoodDensity

    /// Biomass in grams
    let toRadiusMM biomassGrams = 
        let radiusCm = biomassGrams |> Allometrics.shrubRadius b a rtip p lmin k5 k6 numberOfStems salixWoodDensity
        radiusCm * 10.



module Options =
    let resultsDirectory = "/Users/andrewmartin/Desktop/Bristlecone Results/YuribeiAnnual-TemperatureDependent-Tuned/"
    let chains = 3
    let endWhen = Optimisation.EndConditions.afterIteration 100000
    let logger = Logging.RealTimeTrace.graphWithConsole 30. 10000
    let engine =
        Bristlecone.mkContinuous 
        |> Bristlecone.withContinuousTime Integration.MathNet.integrate
        |> Bristlecone.withOutput logger
        |> Bristlecone.withCustomOptimisation (Optimisation.MonteCarlo.adaptiveMetropolis 0.250 500)

    let orchestrator = OrchestrationAgent(logger, System.Environment.ProcessorCount, false)


// 2. Create Hypotheses
// ----------------------------

let ``base model`` maxGrowthRate temperatureDependency additionalParameters =

    let biomass b r gammab geom tempEffect : float =
        b * r * geom(b) * tempEffect - gammab * b

    /// Cumulative stem biomass [b].
    let dbsdt' (b:float) gammab r maxGrowthRate tempEffect =
        biomass b (r * 1000.) gammab maxGrowthRate tempEffect

    /// Bristlecone function for dBs/dt
    let dbsdt p _ b (e:Environment) =
        dbsdt' b (p |> Pool.getEstimate "gamma[b]") ((p |> Pool.getEstimate "r")) (maxGrowthRate p) (temperatureDependency p e)

    /// Measurement (Size) variable: stem radius
    let stemRadius lastRadius lastEnv env =
        let oldCumulativeMass = lookup lastEnv "bs"
        let newCumulativeMass = lookup env "bs"
        if (newCumulativeMass - oldCumulativeMass) > 0.
        then newCumulativeMass |> Allometry.toRadiusMM
        else lastRadius

    { Equations  = [ code "bs",        dbsdt ] |> Map.ofList
      Measures   = [ code "x",         stemRadius ] |> Map.ofList
      Parameters = [ code "gamma[b]",  parameter PositiveOnly   0.001 0.200   // Loss rate of biomass
                     // for likelihood function
                    //  code "rho",       parameter Unconstrained  -0.50 0.500   // Covariance between growth and nitrogen
                    //  code "sigma[x]",  parameter PositiveOnly   0.001 0.100   // Standard deviation of x (biomass)
                    //  code "sigma[y]",  parameter PositiveOnly   0.001 0.100   // Standard deviation of y (nitrogen)
                    ] |> List.append additionalParameters |> Map.ofList
      Likelihood = ModelLibrary.Likelihood.sumOfSquares ["x"] }

let hypotheses =

    // [A] A plant may be subject to mechanical constraints on its maximum size
    let growthModel = 
        [  (fun p -> ModelComponents.GeometricConstraint.none), []
           (fun p -> ModelComponents.GeometricConstraint.chapmanRichards ((p |> Pool.getEstimate "k") * 1000.)),
           [ code "k",  parameter PositiveOnly   3.00 5.00 ] ]     // Asymptotic biomass (grams)

    // [B] Net photosynthetic rate is temperature-dependent
    let temperature =
        [ (fun p e -> ModelComponents.Temperature.arrhenius (p |> Pool.getEstimate "A") (p |> Pool.getEstimate "Ea") (lookup e "T[max]")),
           [ code "A",              parameter PositiveOnly 1.00 2.00
             code "Ea",             parameter PositiveOnly 1.00 2.00  ]
          (fun _ -> ModelComponents.Temperature.none), [] ]

    List.allPairs growthModel temperature
    |> List.map (fun ((growth,gp),(temperature,tp)) -> 
        ``base model`` growth temperature (List.concat [gp; tp]))


// 3. Load Real Data and Estimate
// ----------------------------

let shrubs = 
    Data.PlantIndividual.loadRingWidths (__SOURCE_DIRECTORY__ + "/../data/yamal-rw.csv")

open FSharp.Data
type TemperatureData = CsvProvider<"/Users/andrewmartin/Projects/GitHub-Projects/shrub-nutrient-modelling/data/marre-sale-maxtemp.csv">

let temperatureData =
    let maxTemperatures = TemperatureData.Load "/Users/andrewmartin/Projects/GitHub-Projects/shrub-nutrient-modelling/data/marre-sale-maxtemp.csv"
    maxTemperatures.Rows
    |> Seq.map(fun r -> ( r.``T[max]``, r.Date))
    |> TimeSeries.fromObservations

temperatureData |> TimeSeries.resolution

let getStartValues (startDate:System.DateTime) (plant:PlantIndividual) =
    let initialRadius =
        match plant.Growth with
        | PlantIndividual.PlantGrowth.RingWidth s -> 
            match s with
            | GrowthSeries.Absolute c -> c.Head |> fst |> removeUnit
            | GrowthSeries.Cumulative c -> 
                let start = (c |> TimeSeries.trimStart (startDate - System.TimeSpan.FromDays(366.))).Values |> Seq.head |> removeUnit
                printfn "Start cumulative growth = %f" start
                start
            | GrowthSeries.Relative _ -> invalidOp "Not implemented"
        | _ -> invalidOp "Not implemented 2"
    let initialMass = initialRadius |> removeUnit |> Allometry.toBiomassMM
    let initialTmax = plant.Environment.[ShortCode.create "T[max]"].Head |> fst
    [ (code "x", initialRadius)
      (code "T[max]", initialTmax)
      (code "bs", initialMass) ] |> Map.ofList

let workPackages shrubs hypotheses engine saveDirectory =
    seq {
        for s in shrubs do

            // 1. Arrange the subject and settings
            let shrub = 
                s 
                |> PlantIndividual.toCumulativeGrowth 
                |> PlantIndividual.zipEnv (code "T[max]") temperatureData
            let common = shrub |> PlantIndividual.keepCommonYears
            let startDate = (common.Environment.[code "T[max]"]).StartDate |> snd
            let startConditions = getStartValues startDate shrub
            let e = engine |> Bristlecone.withConditioning (Custom startConditions)

            // 2. Setup batches of dependent analyses
            for h in [ 1 .. hypotheses |> List.length ] do
                for i in [ 1 .. Options.chains ] do
                    yield async {
                            // A. Compute result
                            let result = Bristlecone.PlantIndividual.fit e Options.endWhen hypotheses.[h-1] common
                            // B. Save to file
                            Bristlecone.Data.EstimationResult.saveAll saveDirectory s.Identifier.Value h result
                            return result }
    }

// Orchestrate the analyses
// let work = workPackages (shrubs |> Seq.where (fun x -> x.Identifier.Value = "BV05AB")) hypotheses Options.engine Options.resultsDirectory
let work = workPackages shrubs hypotheses Options.engine Options.resultsDirectory

// work |> Seq.iter (OrchestrationMessage.StartWorkPackage >> Options.orchestrator.Post)

(work |> Seq.head) |> Async.RunSynchronously
