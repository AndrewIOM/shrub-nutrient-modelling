#r "../packages/NETStandard.Library.NETFramework/build/net461/lib/netstandard.dll"
#load "../packages/Bristlecone/bristlecone.fsx"
// #load "../packages/Bristlecone/charts.fsx"
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

module CustomLog =

    open System.Threading
    open Bristlecone.Logging

    let print threadId (x:LogEvent) = 
        match x with
        | Bristlecone.Logging.OptimisationEvent e ->
            if e.Iteration % 750 = 0
            then printfn "##%i## At iteration %i (-logL = %f) %A" threadId e.Iteration e.Likelihood e.Theta
        | _ -> printfn "##%i## %A" threadId x

    let logger () =

        let agent = MailboxProcessor.Start(fun inbox -> 
            let rec messageLoop () = async {
                let! threadId,msg = inbox.Receive()
                print threadId msg
                return! messageLoop ()
                }
            messageLoop () )

        agent.Error.Add(fun e -> printfn "Error = %A" e)

        fun msg -> 
            let threadId = Thread.CurrentThread.ManagedThreadId
            agent.Post (threadId, msg)

// 1. Configure Options
// ----------------------------

module Options =
    let resultsDirectory = "/Users/andrewmartin/Desktop/Bristlecone Results/Yuribei-Snow/"
    let endWhen = Optimisation.EndConditions.afterIteration 25000
    let chains = 1
    let logger = CustomLog.logger() //Logging.Console.logger() //Logging.RealTimeTrace.graphWithConsole 30. 5000
    let filzbachOptions : Optimisation.MonteCarlo.Filzbach.FilzbachSettings<float> = {
        TuneAfterChanges = 50
        MaxScaleChange = 100.00
        MinScaleChange = 0.0010
        BurnLength = Optimisation.EndConditions.afterIteration 25000 }
    let engine =
        Bristlecone.mkContinuous 
        |> Bristlecone.withContinuousTime Integration.MathNet.integrate
        |> Bristlecone.withOutput logger
        |> Bristlecone.withCustomOptimisation (Optimisation.MonteCarlo.Filzbach.filzbach filzbachOptions)
        // |> Bristlecone.withTunedMCMC [ Optimisation.MonteCarlo.TuneMethod.CovarianceWithScale 0.200, 200, Optimisation.EndConditions.afterIteration 5000 ]
        // |> Bristlecone.withCustomOptimisation (Optimisation.MonteCarlo.SimulatedAnnealing.fastSimulatedAnnealing 0.0001 false
        //     { Optimisation.MonteCarlo.SimulatedAnnealing.AnnealSettings<float>.Default with 
        //         BoilingAcceptanceRate = 0.85
        //         TuneLength = 5000
        //         HeatRamp = (fun t -> t + sqrt t); TemperatureCeiling = Some 100.
        //         HeatStepLength = Optimisation.EndConditions.afterIteration 250
        //         AnnealStepLength = (fun x -> Optimisation.MonteCarlo.SimulatedAnnealing.EndConditions.improvementCount 250 250 x) }) //(Optimisation.EndConditions.afterIteration 10000) })
    let orchestrator = OrchestrationAgent(logger, System.Environment.ProcessorCount, true)
    let latitude = 68.91    // Yuribei North
    let longitude = -70.23  // Yuribei West
    let timezone = "Asia/Yekaterinburg"    // Yuribei is in Yekaterinburg Time


// 2. Create Hypotheses
// ----------------------------

module BaseEquations =

    let biomass b n r gammab geom f protectionEffect tempEffect lightEffect : float =
        b * r * (f n) * geom(b) * tempEffect * lightEffect - gammab * protectionEffect(b)

    let soilNitrogen n b gamman y geom f feedback tempEffect lightEffect : float =
        y - (geom(b) * b * (f n) * tempEffect * lightEffect) - gamman * n + feedback(b)

    // Newton's law of cooling / Fourier's law
    let soilTemperature soilT ambientT conductivity =
        if conductivity > 1. then nan
        else - conductivity * (soilT - ambientT)


type ComponentLogger() =

    let mutable (data:Map<string,Map<float,float>>) = [] |> Map.ofList
    with
        member __.StoreValue(componentId, t, v) =
            let existing = data |> Map.tryFind componentId
            match existing with
            | Some e -> data <- data |> Map.add componentId (e |> Map.add t v)
            | None -> data <- data |> Map.add componentId ([t,v] |> Map.ofList)
            v

        member __.GetAll() = data

let ``base model`` maxGrowthRate nLimitation nitrogenFeedback nReplenishment protectionEffect temperatureDependency soilConductivity additionalParameters (cLog:ComponentLogger) (startDate:System.DateTime) latitude longitude timeZone =

    let randomPrint t p e =
        ignore
        // if System.Random().Next(1, 10000) = 1
        // then 
        //     let date = startDate.AddMonths (int t)
        //     Options.logger <| Logging.LogEvent.GeneralEvent (sprintf "[%s] E = %A" (date.ToShortDateString()) e)

    /// Light limitation effect (linear between 0 and 1).
    let seasonalLight t =
        let date = startDate.AddMonths (int t)
        Sunrise.calculate date.Year date.Month date.Day latitude longitude timeZone
        |> Sunrise.dayFraction

    /// Cumulative biomass [B].
    let dbsdt' (b:float) n gammab r maxGrowthRate limit protectionEffect tempEffect dayLength =
        BaseEquations.biomass b n (r * 1000.) gammab maxGrowthRate limit protectionEffect tempEffect dayLength

    /// Bioavailable soil nitrogen [N]
    let dndt' bs n gamman maxGrowthRate feedback limit nReplenishment tempEffect dayLength = 
        BaseEquations.soilNitrogen n bs gamman nReplenishment maxGrowthRate limit feedback tempEffect dayLength

    /// Soil temperature [T_s]
    let dtsdt' ts ta conductivity =
        BaseEquations.soilTemperature ts ta conductivity

    /// Bristlecone function for dBs/dt
    let dbsdt p t bs (e:Environment) =
        (dbsdt' bs ((lookup e "N") |> ModelComponents.Proxies.d15NtoAvailability)
            (p |> Pool.getEstimate "gamma[b]") (p |> Pool.getEstimate "r") (maxGrowthRate p) (nLimitation p) (protectionEffect p e) (*(cLog.StoreValue("photosynthesis-temperature-limit",t,*)(((temperatureDependency p e))) (seasonalLight t))

    /// Bristlecone function for dN/dt
    let dndt p t n (e:Environment) =
        dndt' (lookup e "bs") (n |> ModelComponents.Proxies.d15NtoAvailability) 
            (p |> Pool.getEstimate "gamma[n]") (maxGrowthRate p) (nitrogenFeedback p) (nLimitation p) (*(cLog.StoreValue("nitrogen-replenishment",t,*)((nReplenishment p e)) (temperatureDependency p e) (seasonalLight t)

    /// Bristlecone function for dTs/dt
    let dtsdt p t ts (e:Environment) =
        randomPrint t p e 
        // cLog.StoreValue("day-fraction",t,(seasonalLight t)) |> ignore
        // cLog.StoreValue("soil-conductivity",t,(soilConductivity p e)) |> ignore
        // cLog.StoreValue("soil-temperature",t,ts + (dtsdt' ts (lookup e "Tair") (soilConductivity p e))) |> ignore
        dtsdt' ts (lookup e "Tair") (soilConductivity p e)

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
                     code "N",         dndt
                     code "Tsoil",     dtsdt ] |> Map.ofList
      Measures   = [ code "x",         stemRadius ] |> Map.ofList
      Parameters = [ code "lambda",    parameter PositiveOnly   0.050 0.500   // N-replenishment rate
                     code "gamma[b]",  parameter PositiveOnly   0.001 0.010   // Loss rate of biomass
                     code "gamma[n]",  parameter PositiveOnly   0.001 0.010   // Loss rate of nitrogen
                     code "rho",       parameter Unconstrained  -0.50 0.500   // Covariance between growth and nitrogen
                     code "sigma[x]",  parameter PositiveOnly   0.001 0.500   // Standard deviation of x (biomass)
                     code "sigma[y]",  parameter PositiveOnly   0.001 0.500   // Standard deviation of y (nitrogen)
                    ] |> List.append additionalParameters |> Map.ofList
      Likelihood = ModelLibrary.Likelihood.bivariateGaussian "x" "N" }

module NReplenishment =

    let linear lambda = lambda

    // let temperatureDependent a ea soilTemperature =
    //     ModelComponents.Temperature.arrhenius a ea soilTemperature

    /// The universal gas constant in J mol−1 K−1
    let gasConstant = 8.314

    /// An Arrhenius function to represent temperature limitation on growth.
    /// Form of equation from paper: https://pubag.nal.usda.gov/download/13565/PDF
    let temperatureLimitation preExp activationEnergy temperature =
        //System.Math.E ** (- ((activationEnergy * 1000.) / (8.314 * temperature)))
        preExp * System.Math.E ** ((1000. * activationEnergy * (temperature - 298.)) / (298. * gasConstant * temperature))

let hypotheses =

    // [A] Snow insulates soils, which increases the efficiency of N-mineralising microbes.
    let nitrogenReplenishment =
        [ // 1. Linear rate not affected by soil temperatures
          (fun p _ -> ModelComponents.Temperature.Conductivity.linear 1.),
          (fun p _ -> ModelComponents.Temperature.NitrogenReplenishment.linear (p |> Pool.getEstimate "lambda")), []
          // 2. Temperature-dependent microbial activity
          (fun p _ -> ModelComponents.Temperature.Conductivity.linear (p |> Pool.getEstimate "condct")),
          (fun p e -> NReplenishment.temperatureLimitation (p |> Pool.getEstimate "lambda") (p |> Pool.getEstimate "soilEa") (lookup e "Tsoil")),
            [ code "soilEa",        parameter PositiveOnly   5.000 10.00
              code "condct",    parameter PositiveOnly   0.001 0.200 ]
          // 3. Temperature-dependent microbial activity, with a snow-insulation effect on soil temperatures
          (fun p e -> ModelComponents.Temperature.Conductivity.snowConductivity (p |> Pool.getEstimate "condct") (lookup e "snowDepth") ),
          (fun p e -> NReplenishment.temperatureLimitation (p |> Pool.getEstimate "lambda") (p |> Pool.getEstimate "soilEa") (lookup e "Tsoil")),
            [ code "soilEa",    parameter PositiveOnly   20.00 30.00
              code "condct",    parameter PositiveOnly   0.001 0.200 ] ]

    // [B] Increased snow levels protect shrub biomass from storm and other damage.
    let snowProtection =
        [ (fun _ _ -> ModelComponents.SnowProtection.none), []
          (fun p e -> ModelComponents.SnowProtection.linearSnowProtection (p |> Pool.getEstimate "spe") (lookup e "bs" |> ModelComponents.Proxies.shrubHeightCm) (lookup e "snowDepth")),
            [ code "spe",        parameter PositiveOnly   0.001 1.000 ] ]

    // [C] Net photosynthetic rate is temperature-dependent
    let temperature =
        [ (fun _ -> ModelComponents.Temperature.none), []
          (fun p e -> NReplenishment.temperatureLimitation 1. (*(p |> Pool.getEstimate "A")*) (p |> Pool.getEstimate "Ea") (lookup e "Tair")),
           [ //code "A",              parameter PositiveOnly 5.00 20.0
             code "Ea",             parameter PositiveOnly 30.00 40.00  ] ]

    // [A] N may limited growth via combined N-limitations on (a) photosynthetic and (b) uptake rates
    let limitationModes =
        [ (fun p -> ModelComponents.GrowthLimitation.hollingDiscModelDual ((p |> Pool.getEstimate "a") / 1000.) ((p |> Pool.getEstimate "r") * 1000.) (p |> Pool.getEstimate "h") 5.),
           [ code "a",      parameter PositiveOnly   0.100 0.500          // N-uptake efficiency
             code "h",      parameter PositiveOnly   0.100 8.000
             code "r",      parameter PositiveOnly   2.000 8.000 ]      // N-handling time (including uptake and incorporation)
          (fun p -> ModelComponents.GrowthLimitation.linear ((p |> Pool.getEstimate "a") / 1000.) 5.), 
           [ code "a",      parameter PositiveOnly   0.010 0.100
             code "r",      parameter PositiveOnly   4.000 8.000 ] ]       // N-uptake efficiency

    // [B] Loss of plant material may feedback into the soil pool of available nitrogen (instant)
    let feedbackModes =
        [ (fun _ -> ModelComponents.FeedbackToSoil.none), []
          (fun p -> ModelComponents.FeedbackToSoil.withBiomassLoss ((p |> Pool.getEstimate "alpha") / 100.) (p |> Pool.getEstimate "gamma[b]") ),
          [ code "alpha",  parameter PositiveOnly   0.001 0.002 ] ]    // N-recycling efficiency

    // [C] A plant may be subject to mechanical constraints on its maximum size
    let geometricModes = 
        [  (fun p -> ModelComponents.GeometricConstraint.none), []
           (fun p -> ModelComponents.GeometricConstraint.chapmanRichards ((p |> Pool.getEstimate "k") * 1000.)),
           [ code "k",  parameter PositiveOnly   3.00 5.00 ] ]     // Asymptotic biomass (grams)

    // Create all combinations of H1-H5 
    List.combine6 geometricModes limitationModes feedbackModes nitrogenReplenishment snowProtection temperature
    |> List.map (fun ((growth,gp),(limit,lp),(feedback,fp), (conduct,replace,rp), (snow,sp), (temp,tp)) -> 
        ``base model`` growth limit feedback replace snow temp conduct (List.concat [lp; fp; gp; rp; sp; tp]))


// 3. Load Real Data and Estimate
// ----------------------------
// A. Daily air temperature data from Marre Sale weather station.
// B. Daily snow depth data from earth observation.

open FSharp.Data

[<Literal>] 
let ClimateUrl = __SOURCE_DIRECTORY__ + "/../data/marre-sale-meantemp.csv"
[<Literal>] 
let SnowUrl = __SOURCE_DIRECTORY__ + "/../data/yuribei-snow-depth.csv"

type DailyTemperature = CsvProvider<ClimateUrl>
type DailySnowDepth = CsvProvider<SnowUrl>

let monthlyTemperatures =
    let maxTemperatures = DailyTemperature.Load ClimateUrl
    maxTemperatures.Rows
    |> Seq.map(fun r -> ( (if System.Double.IsNaN r.``T[avg]`` then None else Some (r.``T[avg]`` + 273.15)), r.Date))
    |> TimeSeries.fromObservations
    |> TimeSeries.interpolate
    |> TimeSeries.generalise (Months 1) (fun x -> x |> Seq.averageBy fst)

let monthlySnow =
    let snowDepths = DailySnowDepth.Load SnowUrl
    snowDepths.Rows
    |> Seq.map(fun (r:DailySnowDepth.Row) -> (if  (float r.MeanSnowDepth) = -9999. then None else Some <| float r.MeanSnowDepth), r.DATE)
    |> TimeSeries.fromObservations
    |> TimeSeries.interpolate // Need to interpolate dates
    |> TimeSeries.generalise (Months 1) (fun x -> x |> Seq.averageBy fst)

let growthMap fn plant =
    match plant with
    | RingWidth rw ->
        match rw with
        | Absolute rws -> rws |> TimeSeries.toObservations |> Seq.map fn |> TimeSeries.fromObservations |> Absolute |> RingWidth
        | _ -> invalidOp "Not implemented"
    | _ -> invalidOp "Not implemented"

// b) Load shrub individual data + environment
let shrubs = 
    let yuribei = Data.PlantIndividual.loadRingWidths (__SOURCE_DIRECTORY__ + "/../data/yamal-rw.csv")
    let d15N = Data.PlantIndividual.loadLocalEnvironmentVariable (__SOURCE_DIRECTORY__ + "/../data/yuribei-d15N-imputed.csv")
    yuribei
    |> Seq.map (fun s -> s.Identifier.Value, s)
    |> Seq.keyMatch d15N
    |> Seq.map (fun (_,d15N,plant) -> 
        // Custom: move RW and d15N measurements to end of September (from end of December)
        let p = {plant with Growth = plant.Growth |> growthMap (fun (x,t) -> (x, t - System.TimeSpan.FromDays 92.)) }
        let n = d15N |> TimeSeries.toObservations |> Seq.map(fun (x,t) -> (x, t - System.TimeSpan.FromDays 92.)) |> TimeSeries.fromObservations
        PlantIndividual.zipEnv (code "N") n p
        |> PlantIndividual.zipEnv (code "Tair") monthlyTemperatures
        |> PlantIndividual.zipEnv (code "snowDepth") monthlySnow )
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
    let initialNitrogen = plant.Environment.[code "N"].Head |> fst
    let initialAirT = plant.Environment.[code "Tair"].Head |> fst
    let initialSnow = plant.Environment.[code "snowDepth"].Head |> fst
    [ code "x", initialRadius
      code "N", initialNitrogen 
      code "Tair", initialAirT
      code "Tsoil", initialAirT
      code "snowDepth", initialSnow
      code "bs", initialMass ] |> Map.ofList

let saveLog (cLog:ComponentLogger) (h:int) (guid:System.Guid) =
    let csv = System.IO.File.AppendText (sprintf "%slogger-%i-%s.csv" Options.resultsDirectory h (guid.ToString()))
    let x = cLog.GetAll()
    for m in x do
        for ts in m.Value do
            csv.WriteLine (sprintf "%s,%f,%f" m.Key ts.Key ts.Value)

let workPackages shrubs hypotheses engine saveDirectory =
    seq {
        for s in shrubs do

            // 1. Arrange the subject and settings
            let shrub = s |> PlantIndividual.toCumulativeGrowth
            let common = 
                shrub 
                |> PlantIndividual.keepCommonYears
            let startDate = (common.Environment.[code "N"]).StartDate |> snd
            let endDate = (common.Environment.[code "N"]) |> TimeSeries.endDate
            let shrubWithHighResEnvironment =
                common
                |> PlantIndividual.zipEnv (code "Tair") (monthlyTemperatures |> TimeSeries.bound (startDate - System.TimeSpan.FromDays 366.) endDate)
                |> PlantIndividual.zipEnv (code "snowDepth") (monthlySnow |> TimeSeries.bound (startDate - System.TimeSpan.FromDays 366.) endDate) // TODO Bristlecone should standardise environmental time-series of varying length
            let startConditions = getStartValues startDate shrub
            let e = engine |> Bristlecone.withConditioning (Custom startConditions)

            printfn "Shrub: %A" shrubWithHighResEnvironment

            // 2. Setup batches of dependent analyses
            for h in [ 1 .. hypotheses |> List.length ] do
                for _ in [ 1 .. Options.chains ] do
                    yield async {
                            let cLog = ComponentLogger()
                            let result = Bristlecone.PlantIndividual.fit e Options.endWhen (hypotheses.[h-1] cLog startDate Options.latitude Options.longitude Options.timezone) shrubWithHighResEnvironment
                            saveLog cLog h (result |> fst).ResultId

                            Bristlecone.Data.EstimationResult.saveAll saveDirectory s.Identifier.Value h 1 (result |> fst)
                            return result |> fst }
    }

// Orchestrate the analyses
let work = workPackages shrubs hypotheses Options.engine Options.resultsDirectory |> Seq.toList

let run () =
    work |> Seq.iter (OrchestrationMessage.StartWorkPackage >> Options.orchestrator.Post)
