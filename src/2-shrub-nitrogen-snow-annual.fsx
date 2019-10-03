#r "../packages/NETStandard.Library.NETFramework/build/net461/lib/netstandard.dll"
// #r "../packages/MathNet.Numerics.FSharp/lib/netstandard2.0/MathNet.Numerics.FSharp.dll"
#load "../packages/Bristlecone/bristlecone.fsx"
#r "../packages/Bristlecone.Dendro/lib/netstandard2.0/bristlecone.Dendro.dll"

#load "components/components.fsx"
#load "components/temperature.fsx"
#load "components/snow.fsx"

////////////////////////////////////////////////////
/// Yamal Salix lanata Shrub - Nitrogen Interactions
////////////////////////////////////////////////////

open Bristlecone
open Bristlecone.ModelSystem
open Bristlecone.Dendro
open Bristlecone.Dendro.PlantIndividual
open Bristlecone.Workflow.Orchestration
open Bristlecone.Diagnostics.ModelComponents

// .NET Core requires this workaround to bind to MathNet. Hopefully will be fixed in 3.0 RTM.
let x = MathNet.Numerics.Random.MersenneTwister() 

// 1. Configure Options
// ----------------------------

module Options =
    let resultsDirectory = "/Users/andrewmartin/Desktop/Bristlecone Results/Paper2-Snow-Short-Annual-SA-Cover/"
    let endWhen = Optimisation.EndConditions.afterIteration 25000
    let chains = 1
    let thin = 50
    let engine =
        Bristlecone.mkContinuous 
        |> Bristlecone.withContinuousTime Integration.MathNet.integrate
        // |> Bristlecone.withTunedMCMC [ Optimisation.MonteCarlo.TuneMethod.CovarianceWithScale 0.250, 2000, Optimisation.EndConditions.afterIteration 75000 ]
        // |> Bristlecone.withCustomOptimisation (Optimisation.MonteCarlo.Filzbach.filzbach 
        //     { TuneAfterChanges = 50
        //       MaxScaleChange = 100.00
        //       MinScaleChange = 0.0010
        //       BurnLength = Optimisation.EndConditions.afterIteration 250000 })
        |> Bristlecone.withCustomOptimisation (Optimisation.MonteCarlo.SimulatedAnnealing.fastSimulatedAnnealing 0.0001 false
            { Optimisation.MonteCarlo.SimulatedAnnealing.AnnealSettings<float>.Default with 
                BoilingAcceptanceRate = 0.85
                HeatRamp = (fun t -> t + sqrt t); TemperatureCeiling = Some 500.
                HeatStepLength = Optimisation.EndConditions.afterIteration 1000
                AnnealStepLength = (fun x -> Optimisation.MonteCarlo.SimulatedAnnealing.EndConditions.improvementCount 5000 250 x || Optimisation.EndConditions.afterIteration 5000 x) }) //(Optimisation.EndConditions.afterIteration 10000) })

    let orchestrator = OrchestrationAgent(engine.LogTo, System.Environment.ProcessorCount, false)


// 2. Create Hypotheses
// ----------------------------

let ``base model`` maxGrowthRate nLimitation nitrogenFeedback nReplenishment protectionEffect temperatureDependency additionalParameters (cLog:IComponentLogger<float>) =

    let biomass b n r gammab geom f (protectionEffect:float) tempEffect lightEffect cLog : float =
        if protectionEffect > (1. - 1e-03) then nan // Protection cannot be perfect
        else b * r * (f n) * geom(b) * tempEffect * lightEffect - gammab * (1. - protectionEffect) * b

    let soilNitrogen n b gamman y geom f feedback tempEffect lightEffect (protectionEffect:float) cLog : float =
        y - (geom(b) * b * (f n) * tempEffect * lightEffect) - gamman * n + feedback(b) * (1. - protectionEffect)

    let limit p =
        match nLimitation p with
        | Some l -> l
        | None -> invalidOp "N limitation is required for this model."

    /// Cumulative biomass [B].
    let dbsdt' (b:float) n gammab r maxGrowthRate limit protectionEffect tempEffect dayLength =
        biomass b n (r * 1000.) gammab maxGrowthRate limit protectionEffect tempEffect dayLength

    /// Bioavailable soil nitrogen [N]
    let dndt' bs n gamman maxGrowthRate feedback limit nReplenishment tempEffect dayLength protectionEffect = 
        soilNitrogen n bs gamman nReplenishment maxGrowthRate limit feedback tempEffect dayLength protectionEffect

    let seasonalLight t = 1.

    /// Bristlecone function for dBs/dt
    let dbsdt p t bs (e:Environment) =
        (dbsdt' bs ((lookup e "N") |> ModelComponents.Proxies.d15NtoAvailability)
            (p |> Pool.getEstimate "gamma[b]") (p |> Pool.getEstimate "r") (maxGrowthRate p) (limit p) (protectionEffect p e) (((temperatureDependency p e))) (seasonalLight t)) (fun s f -> cLog.StoreValue s t f |> ignore)

    /// Bristlecone function for dN/dt
    let dndt p t n (e:Environment) =
        dndt' (lookup e "bs") (n |> ModelComponents.Proxies.d15NtoAvailability) 
            (p |> Pool.getEstimate "gamma[n]") (maxGrowthRate p) (nitrogenFeedback p) (limit p) ((nReplenishment p e)) (temperatureDependency p e) (seasonalLight t) (protectionEffect p e) (fun s f -> cLog.StoreValue s t f |> ignore)

    /// Measurement (Size) variable: stem radius
    let stemRadius lastRadius lastEnv env =
        let oldCumulativeMass = lookup lastEnv "bs"
        let newCumulativeMass = lookup env "bs"
        if (newCumulativeMass - oldCumulativeMass) > 0.
        then 
            // This clause mandates that stem radius can never decrease.
            let newRadius = newCumulativeMass |> ModelComponents.Proxies.toRadiusMM
            if newRadius > lastRadius then newRadius else lastRadius
        else lastRadius

    /// Bristlecone function for dr/dt
    { Equations  = [ code "bs",        dbsdt
                     code "N",         dndt ] |> Map.ofList
      Measures   = [ code "x",         stemRadius ] |> Map.ofList
      Parameters = [ code "r",         parameter PositiveOnly   0.500 1.000
                     code "gamma[n]",  parameter PositiveOnly   0.001 0.200   // Loss rate of nitrogen
                     code "gamma[b]",  parameter PositiveOnly   0.001 0.200   // Loss rate of biomass
                     code "rho",       parameter Unconstrained  -0.50 0.500   // Covariance between growth and nitrogen
                     code "sigma[x]",  parameter PositiveOnly   0.001 0.100   // Standard deviation of x (biomass)
                     code "sigma[y]",  parameter PositiveOnly   0.001 0.100   // Standard deviation of y (nitrogen)
                    ] |> List.append additionalParameters |> Map.ofList
      Likelihood = ModelLibrary.Likelihood.bivariateGaussian "x" "N" }


module NReplenishment =

    let linear lambda = lambda

    /// The universal gas constant in J mol−1 K−1
    let gasConstant = 8.314

    /// An Arrhenius function to represent temperature limitation on growth.
    /// Form of equation from paper: https://pubag.nal.usda.gov/download/13565/PDF.
    /// Temperature limitation is between 0 and 1.
    /// When temperature increases above 25 degrees Celsius, temperature limitation = 1.
    let temperatureLimitation preExp activationEnergy temperature =
        preExp * (min (System.Math.E ** ((1000. * activationEnergy * (temperature - 298.)) / (298. * gasConstant * temperature))) 1.)

// The snow protection function requires height, but this computes very slowly.
// Here is a lookup table to make this faster (only approximate).

let biomassToHeight' =
    let x = [ 1. .. 5. .. 100000. ] |> List.map (fun x -> x, x |> ModelComponents.Proxies.toRadiusMM |> ModelComponents.Proxies.shrubHeightCm)
    List.append x [(System.Double.MaxValue, nan)]

let biomassToHeight biomass =
    let h = biomassToHeight' |> List.tryFind(fun (b,_) -> b > biomass)
    match h with
    | Some h -> h |> snd
    | None -> nan

let hypotheses =

    // [A] Snow insulates soils, which increases the efficiency of N-mineralising microbes.
    let nitrogenReplenishment =
        [ (fun p _ -> ModelComponents.Temperature.NitrogenReplenishment.linear (p |> Pool.getEstimate "lambda")),
            [ code "lambda",        parameter PositiveOnly   0.001 1.000 ]
          (fun p e -> ModelComponents.Temperature.NitrogenReplenishment.temperatureDependent (p |> Pool.getEstimate "soilExp") (p |> Pool.getEstimate "soilEa") (p |> Pool.getEstimate "conductivity") (lookup e "Tsummer") (lookup e "Twinter")),
            [ code "soilExp",       parameter PositiveOnly   0.001 1.000
              code "soilEa",        parameter PositiveOnly   10.00 20.00
              code "conductivity",    parameter PositiveOnly   0.001 1.000 ]
          (fun p e -> ModelComponents.Temperature.NitrogenReplenishment.temperatureDependentSnowInsulation (p |> Pool.getEstimate "soilExp") (p |> Pool.getEstimate "soilEa") (p |> Pool.getEstimate "insulation") (lookup e "snowDepth") (lookup e "Tsummer") (lookup e "Twinter")),
            [ code "soilExp",       parameter PositiveOnly   0.001 1.000
              code "soilEa",        parameter PositiveOnly   10.00 20.00
              code "insulation",    parameter PositiveOnly   0.001 1.000 ] ]

    // [B] Increased snow levels protect shrub biomass from storm and other damage.
    let snowProtection =
        [ (fun _ _ -> ModelComponents.SnowProtection.none), []
          (fun p e -> ModelComponents.SnowProtection.withShrubHeight (p |> Pool.getEstimate "spe") (lookup e "bs" |> biomassToHeight) (lookup e "snowDepth")),
            [ code "spe",        parameter PositiveOnly   0.0001 0.500 ] ]

    // [C] Net photosynthetic rate is temperature-dependent
    let temperature =
        [ (fun _ -> ModelComponents.Temperature.none), []
          (fun p e -> NReplenishment.temperatureLimitation 1. (*(p |> Pool.getEstimate "A")*) (p |> Pool.getEstimate "Ea") (lookup e "Tsummer")),
           [ //code "A",              parameter PositiveOnly 5.00 20.0
             code "Ea",             parameter PositiveOnly 05.00 40.00  ] ]

    // [D] N may limited growth via combined N-limitations on (a) photosynthetic and (b) uptake rates
    let limitationModes =
        [ (fun p -> ModelComponents.GrowthLimitation.hollingDiscModelDual ((p |> Pool.getEstimate "a") / 1000.) ((p |> Pool.getEstimate "r") * 1000.) ((p |> Pool.getEstimate "h") / 1.) 5.00),
           [ code "a",      parameter PositiveOnly   0.100 5.000          // N-uptake efficiency
             code "h",      parameter PositiveOnly   0.001 0.200 ]
          (fun p -> ModelComponents.GrowthLimitation.linear ((p |> Pool.getEstimate "a") / 1000.) 10.), 
           [ code "a",      parameter PositiveOnly   0.100 5.000 ] ]       // N-uptake efficiency

    // [E] Loss of plant material may feedback into the soil pool of available nitrogen (instant)
    let feedbackModes =
        [ (fun p -> ModelComponents.FeedbackToSoil.withBiomassLoss ((p |> Pool.getEstimate "alpha") / 100.) (p |> Pool.getEstimate "gamma[b]") ),
          [ code "alpha",  parameter PositiveOnly   0.001 0.400 ]
          (fun _ -> ModelComponents.FeedbackToSoil.none), [] ]    // N-recycling efficiency

    // [F] A plant may be subject to mechanical constraints on its maximum size
    let geometricModes = 
        [  (fun p -> ModelComponents.GeometricConstraint.chapmanRichards ((p |> Pool.getEstimate "k") * 1000.)),
           [ code "k",  parameter PositiveOnly   3.00 5.00 ]
           (fun p -> ModelComponents.GeometricConstraint.none), [] ]     // Asymptotic biomass (grams)

    // Create all combinations of H1-H5 
    List.combine6 geometricModes feedbackModes limitationModes nitrogenReplenishment snowProtection temperature
    |> List.map (fun ((growth,gp),(feedback,fp),(limit,lp),(replenish,rp), (snow,sp), (temp,tp)) -> 
        ``base model`` growth limit feedback replenish snow temp (List.concat [lp; fp; gp; rp; sp; tp]))


// 3. Load Real Data and Estimate
// ----------------------------
// A. Daily air temperature data from Marre Sale weather station.
// B. Daily snow depth data from earth observation.

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

[<Literal>] 
let SnowUrl = __SOURCE_DIRECTORY__ + "/../data/yuribei-snow-depth-annual.csv"

type SnowData = CsvProvider<SnowUrl>

let winterSnowDepthMean =
    let snowDepths = SnowData.Load SnowUrl
    snowDepths.Rows
    |> Seq.map(fun (r:SnowData.Row) -> (float r.``Winter (October - March) Mean``, r.Year))
    |> TimeSeries.fromObservations

// b) Load shrub individual data + environment
let shrubs = 
    let yuribei = Data.PlantIndividual.loadRingWidths (__SOURCE_DIRECTORY__ + "/../data/yamal-rw.csv")
    let d15N = Data.PlantIndividual.loadLocalEnvironmentVariable (__SOURCE_DIRECTORY__ + "/../data/yuribei-d15N-imputed.csv")
    let d18O = Data.PlantIndividual.loadLocalEnvironmentVariable (__SOURCE_DIRECTORY__ + "/../data/yuribei-d18O-imputed.csv")
    yuribei
    |> Seq.map (fun s -> s.Identifier.Value, s)
    |> Seq.keyMatch d15N
    |> Seq.map (fun (_,plant,d15N) -> 
        let s = PlantIndividual.zipEnv (code "N") plant d15N
        s.Identifier.Value, s)
    |> Seq.keyMatch d18O
    |> Seq.map (fun (_,plant,d18O) -> 
        PlantIndividual.zipEnv (code "d18O") plant d18O
        |> PlantIndividual.zipEnv (code "Tsummer") summerTemperature
        |> PlantIndividual.zipEnv (code "Twinter") winterTemperature
        |> PlantIndividual.zipEnv (code "snowDepth") winterSnowDepthMean)
    |> Seq.toList

/// How to generalise `getStartValues` for Bristlecone or
/// Bristlecone.Dendro?
/// 
/// - A. Determine if start time - 1 (what time-step?) is present. If so, use this. Otherwise, repeat start value.
/// - B. For PlantIndividual, must interpret the type of growth series.
/// - C. When doing allometric model, must set start mass as allometric transform.
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
    let initialNitrogen = plant.Environment.[code "N"].Head |> fst
    let initialWinterT = plant.Environment.[ShortCode.create "Twinter"].Head |> fst
    let initialSummerT = plant.Environment.[ShortCode.create "Tsummer"].Head |> fst
    let initialD18O = plant.Environment.[ShortCode.create "d18O"].Head |> fst
    let initialSnow = plant.Environment.[code "snowDepth"].Head |> fst
    [ code "x", initialRadius
      code "N", initialNitrogen 
      code "Twinter", initialWinterT
      code "Tsummer", initialSummerT
      code "d18O", initialD18O
      code "snowDepth", initialSnow
      code "bs", initialMass ] |> Map.ofList

/// How to generalise the `fit` function?
/// - Start Date = common 'variables' (as opposed to environmental data). Codes can be got from ModelSystem.
/// - End Date = common 'variables' (as opposed to environmental data).
/// - Ensure that start values are present (at t-1, including for env variables).
let fit s hypothesis engine = 
    let shrub = s |> PlantIndividual.toCumulativeGrowth
    let common = shrub |> PlantIndividual.keepCommonYears
    let startDate = (common.Environment.[code "N"]).StartDate |> snd
    let startConditions = getStartValues startDate shrub
    let e = engine |> Bristlecone.withConditioning (Custom startConditions)
    Bristlecone.PlantIndividual.fit e Options.endWhen hypothesis common

/// Use this type in Bristlecone itself.
type WorkPackage = Async<EstimationResult>

/// How to generalise this function?
/// - It is just a nested set of commands.
let workPackages shrubs hypotheses engine saveDirectory =
    seq {
        for s in shrubs do
            for h in [ 1 .. hypotheses |> List.length ] do
                for _ in [ 1 .. Options.chains ] do
                    if h > 60 && h < 73 then
                        yield async {
                                let cLog = PassThrough()
                                let result = fit s (hypotheses.[h-1] cLog) engine
                                Bristlecone.Data.EstimationResult.saveAll saveDirectory s.Identifier.Value h Options.thin (result |> fst)
                                return result |> fst }
    }

// Orchestrate the analyses
let shrubsWithIsotope = [ "YUSL03A"; "YUSL05A"; "YUSL26A"; "YUSL29A"; "YUSL39A" ]
let work = workPackages (shrubs |> Seq.where(fun s -> shrubsWithIsotope |> List.contains s.Identifier.Value)) hypotheses Options.engine Options.resultsDirectory |> Seq.toList
let run () = work |> Seq.rev |> Seq.iter (OrchestrationMessage.StartWorkPackage >> Options.orchestrator.Post)


// Temp
module Temp =

    open Bristlecone.Data

    /// Load an `EstimationResult` that has previously been saved as
    /// three seperate dataframes. Results will only be reconstructed
    /// when file names and formats are in original Bristlecone format.
    let loadAll directory subject (modelSystem:ModelSystem) modelId =
        let mles = MLE.load directory subject modelId |> Seq.map(fun (k,v) -> k.ToString(), v)
        let series = Series.load directory subject modelId |> Seq.map(fun (k,v) -> k.ToString(), v)
        let traces = Trace.load directory subject modelId |> Seq.map(fun (k,v) -> k.ToString(), v)

        printfn "%s %i" subject modelId

        let updateParameter (k:ShortCode) v newParameters =
            printfn "New key is %s" k.Value
            if k = code "conductivity" then 
                match newParameters |> Map.tryFind k with
                | Some i -> Parameter.setEstimate v i
                | None -> Parameter.setEstimate v (newParameters |> Map.find (code "insulation"))
            else Parameter.setEstimate v (newParameters |> Map.find k)

        mles
        |> Seq.keyMatch series
        |> Seq.map(fun (k,v1,v2) -> (k, (v1, v2)))
        |> Seq.keyMatch traces
        |> Seq.map(fun (k,t,(s,(l,p))) ->
            { ResultId = k |> System.Guid.Parse
              Likelihood = l
              Parameters = modelSystem.Parameters |> Map.map(fun k v -> updateParameter k v p)
              Series = s 
              Trace = t })


let saveDiagnostics () =

    let resultDir = "/Users/andrewmartin/Desktop/Bristlecone Results/Paper-2-Thesis-Version/Annual-RemoteSensed/"

    // 1. Get all results sliced by plant and hypothesis
    let results = 
        let hypotheses = hypotheses |> List.map(fun h -> h (PassThrough()))
        let get subject model modelId = Temp.loadAll resultDir subject.Identifier.Value model modelId
        Bristlecone.ModelSelection.ResultSet.arrangeResultSets shrubs hypotheses get

    // 2. Save convergence statistics to file
    results 
    |> Seq.map(fun (x,a,b,c) -> x.Identifier.Value,a,b,c)
    |> Seq.toList
    |> Diagnostics.Convergence.gelmanRubinAll 10000 3
    |> Data.Convergence.save Options.resultsDirectory

    // 3. Save Akaike weights to file
    results
    |> ModelSelection.Select.weights
    |> Seq.map(fun (x,a,b,c) -> x.Identifier.Value,a,b,c)
    |> Data.ModelSelection.save resultDir

    // 4. Save out logged components
    // results
    // |> Seq.map(fun r -> calculateComponents fit Options.engine r)
