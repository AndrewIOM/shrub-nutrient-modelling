#r "../packages/NETStandard.Library.NETFramework/build/net461/lib/netstandard.dll"
#load "../packages/Bristlecone/bristlecone.fsx"
#load "../packages/Bristlecone/charts.fsx"
#r "../packages/Bristlecone.Dendro/lib/netstandard2.0/bristlecone.Dendro.dll"
#load "components/components.fsx"
#load "components/temperature.fsx"
#load "components/snow.fsx"

////////////////////////////////////////////////////
/// Yamal Salix lanata Shrub - Nitrogen Interactions
////////////////////////////////////////////////////

// Shrub ring width modelled with a single
// resource limitation.

// TODO Specify index mode (e.g. interpolate / floor / replicate) in scripts

open Bristlecone
open Bristlecone.ModelSystem
open Bristlecone.PlantIndividual
open Bristlecone.Workflow.Orchestration

// 1. Configure Options
// ----------------------------

module Options =
    let resultsDirectory = "/Users/andrewmartin/Desktop/Bristlecone Results/Yuribei-Snow/"
    let endWhen = Optimisation.EndConditions.afterIteration 1000
    let chains = 8
    let logger = Logging.RealTimeTrace.graphWithConsole 30. 10000
    let engine =
        Bristlecone.mkContinuous 
        |> Bristlecone.withContinuousTime Integration.MathNet.integrate
        |> Bristlecone.withOutput logger
        |> Bristlecone.withCustomOptimisation (Optimisation.MonteCarlo.SimulatedAnnealing.fastSimulatedAnnealing 0.001 { Optimisation.MonteCarlo.SimulatedAnnealing.AnnealSettings.Default with HeatRamp = (fun t -> t + sqrt t); TemperatureCeiling = Some 5000.; HeatStepLength = Optimisation.EndConditions.afterIteration 1500; AnnealStepLength = Optimisation.EndConditions.afterIteration 5000  })
    let orchestrator = OrchestrationAgent(logger, System.Environment.ProcessorCount)


// 2. Create Hypotheses
// ----------------------------

module BaseEquations =

    /// Cumulative stem biomass [dBs/dt]
    let biomass b n r gammab geom f protectionEffect tempEffect : float =
        b * r * (f n) * geom(b) * tempEffect - gammab * protectionEffect(b)

    let soilNitrogen n b gamman y geom f feedback tempEffect : float =
        y - (geom(b) * b * (f n) * tempEffect) - gamman * n + feedback(b)

    let soilNitrogenNoUptake n b gamman y feedback : float =
        y - gamman * n + feedback(b)


let ``base model`` maxGrowthRate nLimitation nitrogenFeedback nReplenishment protectionEffect temperatureDependency additionalParameters =

    /// Cumulative stem biomass [b].
    let dbsdt' (b:float) n gammab r maxGrowthRate limit protectionEffect tempEffect =
        match limit with
        | Some l -> BaseEquations.biomass b n (r * 1000.) gammab maxGrowthRate l protectionEffect tempEffect
        | None -> BaseEquations.biomass b n r gammab maxGrowthRate (fun _ -> 1.) protectionEffect tempEffect

    /// Bioavailable soil nitrogen [N]
    let dndt' bs n gamman maxGrowthRate feedback limit nReplenishment tempEffect = 
        match limit with
        | Some l -> BaseEquations.soilNitrogen n bs gamman nReplenishment maxGrowthRate l feedback tempEffect
        | None -> BaseEquations.soilNitrogenNoUptake n bs gamman nReplenishment feedback

    /// Bristlecone function for dBs/dt
    let dbsdt p _ bs (e:Environment) =
        dbsdt' bs ((lookup e "N") |> ModelComponents.Proxies.d15NtoAvailability)
            (p |> Pool.getEstimate "gamma[b]") (p |> Pool.getEstimate "r") (maxGrowthRate p) (nLimitation p) (protectionEffect p e) (temperatureDependency p e)

    /// Bristlecone function for dN/dt
    let dndt p _ n (e:Environment) =
        dndt' (lookup e "bs") (n |> ModelComponents.Proxies.d15NtoAvailability) 
            (p |> Pool.getEstimate "gamma[n]") (maxGrowthRate p) (nitrogenFeedback p) (nLimitation p) (nReplenishment p e) (temperatureDependency p e)

    /// Measurement (Size) variable: stem radius
    let stemRadius lastRadius lastEnv env =
        let oldCumulativeMass = lookup lastEnv "bs"
        let newCumulativeMass = lookup env "bs"
        if (newCumulativeMass - oldCumulativeMass) > 0.
        then newCumulativeMass |> ModelComponents.Proxies.toRadiusMM
        else lastRadius


    /// Bristlecone function for dr/dt
    { Equations  = [ code "bs",        dbsdt
                     code "N",         dndt ] |> Map.ofList
      Measures   = [ code "x",         stemRadius ] |> Map.ofList
      Parameters = [ code "gamma[n]",  parameter PositiveOnly   0.001 0.200   // Loss rate of nitrogen
                     code "gamma[b]",  parameter PositiveOnly   0.001 0.200   // Loss rate of biomass
                     code "rho",       parameter Unconstrained  -0.50 0.500   // Covariance between growth and nitrogen
                     code "sigma[x]",  parameter PositiveOnly   0.001 0.100   // Standard deviation of x (biomass)
                     code "sigma[y]",  parameter PositiveOnly   0.001 0.100   // Standard deviation of y (nitrogen)
                    ] |> List.append additionalParameters |> Map.ofList
      Likelihood = ModelLibrary.Likelihood.bivariateGaussian "x" "N" }

let hypotheses =

    // [A] Increased snow levels protect shrub biomass from storm and other damage.
    let snowProtection =
        [ ((fun p e -> ModelComponents.SnowProtection.none), []) ]
          //(fun p e -> ModelComponents.SnowProtection.linearSnowProtection  (lookup e "bs" |> ModelComponents.Proxies.shrubHeightCm) (26. - (lookup e "d18O"))), [] ]

    // [B] Snow insulates soils, which increases the efficiency of N-mineralising microbes.
    let nitrogenReplenishment =
        [ //(fun p _ -> ModelComponents.Temperature.NitrogenReplenishment.linear (p |> Pool.getEstimate "lambda")),
          //  [ code "lambda",        parameter PositiveOnly   1.000 1.000 ]
          (fun p e -> ModelComponents.Temperature.NitrogenReplenishment.temperatureDependent (p |> Pool.getEstimate "lambda") (p |> Pool.getEstimate "soilEa") (p |> Pool.getEstimate "insulation") (lookup e "Tsummer") (lookup e "Twinter")),
            [ code "lambda",        parameter PositiveOnly   0.001 1.000
              code "soilEa",        parameter PositiveOnly   0.001 1.000
              code "insulation",    parameter PositiveOnly   0.001 1.000 ]
          (fun p e -> ModelComponents.Temperature.NitrogenReplenishment.temperatureDependentSnowInsulation (p |> Pool.getEstimate "lambda") (p |> Pool.getEstimate "soilEa") (p |> Pool.getEstimate "insulation") (26. - (lookup e "d18O")) (lookup e "Tsummer") (lookup e "Twinter")),
            [ code "lambda",        parameter PositiveOnly   0.001 1.000
              code "soilEa",        parameter PositiveOnly   0.001 1.000
              code "insulation",    parameter PositiveOnly   0.001 1.000 ] ]

    // [C] Net photosynthetic rate is temperature-dependent
    let temperature =
        [ //(fun _ -> ModelComponents.Temperature.none), []
          (fun p e -> ModelComponents.Temperature.arrhenius (p |> Pool.getEstimate "A") (p |> Pool.getEstimate "Ea") (lookup e "Tsummer")),
           [ code "A",              parameter PositiveOnly 0.01 1.00
             code "Ea",             parameter PositiveOnly 0.01 1.00  ] ]

    // [A] N may limited growth via combined N-limitations on (a) photosynthetic and (b) uptake rates
    let limitationModes =
        [ (fun p -> ModelComponents.GrowthLimitation.hollingDiscModel ((p |> Pool.getEstimate "a") / 100.) ((p |> Pool.getEstimate "r") * 1000.) ((p |> Pool.getEstimate "h") / 10.)),
           [ code "a",      parameter PositiveOnly   0.010 1.000          // N-uptake efficiency
             code "h",      parameter PositiveOnly   0.010 0.500
             code "r",      parameter PositiveOnly   0.100 0.500 ]      // N-handling time (including uptake and incorporation)
          (fun p -> ModelComponents.GrowthLimitation.linear ((p |> Pool.getEstimate "a") / 100.)), 
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

    // Create all combinations of H1-H5 
    List.combine6 geometricModes limitationModes feedbackModes nitrogenReplenishment snowProtection temperature
    |> List.map (fun ((growth,gp),(limit,lp),(feedback,fp), (replace,rp), (snow,sp), (temp,tp)) -> 
        ``base model`` growth limit feedback replace snow temp (List.concat [lp; fp; gp; rp; sp; tp]))


// 3. Load Real Data and Estimate
// ----------------------------

open FSharp.Data

[<Literal>]
let ClimateUrl = __SOURCE_DIRECTORY__ + "/../data/yamal-temperature-mean-quarters.csv"

type Climate = CsvProvider<ClimateUrl>
let regionalClimate = Climate.Load ClimateUrl
let summerTemperature = 
    let convert (x:decimal<_>) = decimal x
    regionalClimate.Rows
    |> Seq.map(fun r -> ((r.``JJA Mean Temperaure`` |> convert |> float), r.Year))
    |> TimeSeries.fromObservations

let winterTemperature = 
    let convert (x:decimal<_>) = decimal x
    regionalClimate.Rows
    |> Seq.map(fun r -> ((r.``DJF Mean Temperature`` |> convert |> float), r.Year))
    |> TimeSeries.fromObservations

// b) Load shrub individual data + environment
let shrubs = 
    let yuribei = Data.PlantIndividual.loadRingWidths (__SOURCE_DIRECTORY__ + "/../data/yamal-rw.csv")
    let d15N = Data.PlantIndividual.loadLocalEnvironmentVariable (__SOURCE_DIRECTORY__ + "/../data/yuribei-d15N-imputed.csv")
    let d18O = Data.PlantIndividual.loadLocalEnvironmentVariable (__SOURCE_DIRECTORY__ + "/../data/yuribei-d18O-imputed.csv")
    yuribei
    |> Seq.map (fun s -> s.Identifier.Value, s)
    |> Seq.keyMatch d15N
    |> Seq.map (fun (_,plant,d15N) -> PlantIndividual.zipEnv (code "N") plant d15N)
    |> Seq.map (fun s -> s.Identifier.Value, s)
    |> Seq.keyMatch d18O
    |> Seq.map (fun (_,plant,d18O) -> PlantIndividual.zipEnv (code "d18O") plant d18O)
    |> Seq.map (fun (plant) -> PlantIndividual.zipEnv (code "Tsummer") summerTemperature plant)
    |> Seq.map (fun (plant) -> PlantIndividual.zipEnv (code "Twinter") winterTemperature plant)
    |> Seq.toList


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
    let initialMass = initialRadius |> removeUnit |> ModelComponents.Proxies.toBiomassMM
    let initialNitrogen = plant.Environment.[ShortCode.create "N"].Head |> fst
    let initialWinterT = plant.Environment.[ShortCode.create "Twinter"].Head |> fst
    let initialSummerT = plant.Environment.[ShortCode.create "Tsummer"].Head |> fst
    let initialSnow = plant.Environment.[ShortCode.create "d18O"].Head |> fst
    [ ShortCode.create "x", initialRadius
      ShortCode.create "N", initialNitrogen 
      ShortCode.create "d18O", initialSnow
      ShortCode.create "Twinter", initialWinterT
      ShortCode.create "Tsummer", initialSummerT
      ShortCode.create "bs", initialMass ] |> Map.ofList

let workPackages shrubs hypotheses engine saveDirectory =
    seq {
        for s in shrubs do

            // 1. Arrange the subject and settings
            let shrub = s |> PlantIndividual.toCumulativeGrowth
            let common = shrub |> PlantIndividual.keepCommonYears
            let startDate = (common.Environment.[ShortCode.create "N"]).StartDate |> snd
            let startConditions = getStartValues startDate shrub
            let e = engine |> Bristlecone.withConditioning (Custom startConditions)
            
            // 2. Setup batches of dependent analyses
            for h in [ 1 .. hypotheses |> List.length ] do
                for i in [ 1 .. Options.chains ] do
                    let jobId = System.Guid.NewGuid()
                    yield async {
                            // A. Compute result
                            let result = Bristlecone.PlantIndividual.fit e Options.endWhen hypotheses.[h-1] common
                            // B. Save to file
                            Bristlecone.Data.Cache.saveTrace saveDirectory s.Identifier.Value h jobId [ result ]
                            return result }
    }

(shrubs |> List.map (fun s -> printfn "%s" s.Identifier.Value))

// Orchestrate the analyses
let work = workPackages (shrubs |> List.take 1) (hypotheses |> List.skip 0 |> List.take 1) Options.engine Options.resultsDirectory

let run () =
    work |> Seq.iter (OrchestrationMessage.StartWorkPackage >> Options.orchestrator.Post)
