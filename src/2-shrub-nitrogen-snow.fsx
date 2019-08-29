#r "../packages/NETStandard.Library.NETFramework/build/net461/lib/netstandard.dll"
#r "../packages/MathNet.Numerics.FSharp/lib/netstandard2.0/MathNet.Numerics.FSharp.dll"
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

let x = MathNet.Numerics.Random.MersenneTwister()  // For some reason, this is needed on .net core??

module CustomLog =

    open System.Threading
    open Bristlecone.Logging

    let print threadId (x:LogEvent) = 
        match x with
        | Bristlecone.Logging.OptimisationEvent e ->
            if e.Iteration % 1000 = 0
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
    let resultsDirectory = "/Users/andrewmartin/Desktop/Bristlecone Results/Paper2-Snow-Short-Reformulated4/"
    let endWhen = Optimisation.EndConditions.afterIteration 25000
    let chains = 1
    let thin = 50
    let logger = CustomLog.logger() //Logging.Console.logger() //Logging.RealTimeTrace.graphWithConsole 30. 5000
    let filzbachOptions : Optimisation.MonteCarlo.Filzbach.FilzbachSettings<float> = {
        TuneAfterChanges = 50
        MaxScaleChange = 100.00
        MinScaleChange = 0.0010
        BurnLength = Optimisation.EndConditions.afterIteration 75000 }
    let engine =
        Bristlecone.mkContinuous 
        |> Bristlecone.withContinuousTime Integration.MathNet.integrate
        |> Bristlecone.withOutput logger
        |> Bristlecone.withCustomOptimisation (Optimisation.MonteCarlo.Filzbach.filzbach filzbachOptions)
        // |> Bristlecone.withTunedMCMC [ Optimisation.MonteCarlo.TuneMethod.CovarianceWithScale 0.100, 250, Optimisation.EndConditions.afterIteration 20000 ]
        // |> Bristlecone.withCustomOptimisation (Optimisation.MonteCarlo.SimulatedAnnealing.fastSimulatedAnnealing 0.0001 false
        //     { Optimisation.MonteCarlo.SimulatedAnnealing.AnnealSettings<float>.Default with 
        //         BoilingAcceptanceRate = 0.85
        //         TuneLength = 2500
        //         HeatRamp = (fun t -> t + sqrt t); TemperatureCeiling = Some 100.
        //         HeatStepLength = Optimisation.EndConditions.afterIteration 250
        //         AnnealStepLength = (fun x -> Optimisation.MonteCarlo.SimulatedAnnealing.EndConditions.improvementCount 1000 250 x) }) //(Optimisation.EndConditions.afterIteration 10000) })
    let orchestrator = OrchestrationAgent(logger, System.Environment.ProcessorCount, false)
    let latitude = 68.91    // Yuribei North
    let longitude = -70.23  // Yuribei West
    let timezone = "Asia/Yekaterinburg"    // Yuribei is in Yekaterinburg Time


// 2. Create Hypotheses
// ----------------------------

type ComponentLogger(turnOff) =

    let mutable (data:Map<string,Map<float,float>>) = [] |> Map.ofList
    with
        member __.StoreValue(componentId, t, v) =
            if not turnOff then
                let existing = data |> Map.tryFind componentId
                match existing with
                | Some e -> data <- data |> Map.add componentId (e |> Map.add t v)
                | None -> data <- data |> Map.add componentId ([t,v] |> Map.ofList)
            v

        member __.GetAll() = data


type CachedSunlight(latitude, longitude, timeZone) =

    let mutable (lightData:Map<System.DateTime,Sunrise.DayLength>) = [] |> Map.ofList
    with
        member __.GetLight(date) =
            match lightData |> Map.tryFind date with
            | Some l -> l
            | None -> 
                let l = Sunrise.calculate date.Year date.Month date.Day latitude longitude timeZone
                lightData <- lightData |> Map.add date l
                l

// type LookupTable(fn,sigFigs) =

//     let roundToSigFigs sigFigs (x:float) = 
//         let absx = System.Math.Abs(x)
//         System.Math.Round(0.5 + System.Math.Log10(absx))
//         |> fun c -> float(c)-float(sigFigs)
//         |> fun d  -> System.Math.Round(absx * 10.0**(-d)) * 10.0**d
//         |> fun a -> a * float(System.Math.Sign(x))

//     let mutable (lookup:Map<float,float>) = [] |> Map.ofList
//     with
//         member __.Get(key) =
//             if key < 0. || System.Double.IsNaN key || System.Double.IsInfinity key || System.Double.IsNegativeInfinity key then nan
//             else
//                 let k = key |> roundToSigFigs sigFigs
//                 match lookup |> Map.tryFind k with
//                 | Some v -> v
//                 | None -> 
//                     let v = (key |> fn)
//                     lookup <- lookup |> Map.add k v
//                     v


module BaseEquations =

    let biomass b n r gammab geom f (protectionEffect:float) tempEffect lightEffect cLog : float =

        // cLog "stem" (b |> ModelComponents.Proxies.toRadiusMM |> ModelComponents.Proxies.shrubHeightCm)
        // cLog "biomassGrowth" (b * r * (f n) * geom(b) * tempEffect * lightEffect)
        // cLog "biomassLoss" (gammab * (1. - protectionEffect) * b)

        b * r * (f n) * geom(b) * tempEffect * lightEffect - gammab * (1. - protectionEffect) * b

    let soilNitrogen n b gamman y geom f feedback tempEffect lightEffect (protectionEffect:float) cLog : float =

        // cLog "nReplenishment" y
        // cLog "nUptake" (geom(b) * b * (f n) * tempEffect * lightEffect)
        // cLog "nLoss" (gamman * n + feedback(b) * (1. - protectionEffect))

        y - (geom(b) * b * (f n) * tempEffect * lightEffect) - gamman * n + feedback(b) * (1. - protectionEffect)

    // Newton's law of cooling / Fourier's law
    let soilTemperature soilT ambientT conductivity =
        if conductivity > 1. || conductivity < 0.05 then nan
        else - conductivity * (soilT - ambientT)


let ``base model`` maxGrowthRate nLimitation nitrogenFeedback nReplenishment protectionEffect temperatureDependency soilConductivity additionalParameters (cLog:ComponentLogger) (startDate:System.DateTime) latitude longitude timeZone =

    let limit p =
        match nLimitation p with
        | Some l -> l
        | None -> invalidOp "N limitation is required for this model."

    /// Light limitation effect (linear between 0 and 1).
    /// Light is cached in an object to avoid unnecessary computation.
    let lightFn = CachedSunlight(latitude, longitude, timeZone)
    let seasonalLight t =
        let date = startDate.AddMonths (int t)
        date |> lightFn.GetLight |> Sunrise.dayFraction

    /// Cumulative biomass [B].
    let dbsdt' (b:float) n gammab r maxGrowthRate limit protectionEffect tempEffect dayLength =
        BaseEquations.biomass b n (r * 1000.) gammab maxGrowthRate limit protectionEffect tempEffect dayLength

    /// Bioavailable soil nitrogen [N]
    let dndt' bs n gamman maxGrowthRate feedback limit nReplenishment tempEffect dayLength protectionEffect = 
        BaseEquations.soilNitrogen n bs gamman nReplenishment maxGrowthRate limit feedback tempEffect dayLength protectionEffect

    /// Soil temperature [T_s]
    let dtsdt' ts ta conductivity =
        BaseEquations.soilTemperature ts ta conductivity

    /// Bristlecone function for dBs/dt
    let dbsdt p t bs (e:Environment) =

        // cLog.StoreValue("protectionEffect",t,((protectionEffect p e))) |> ignore
        // cLog.StoreValue("biomass",t,((bs))) |> ignore
        // cLog.StoreValue("temperature-limit-on-photosynthesis",t,(temperatureDependency p e)) |> ignore

        (dbsdt' bs ((lookup e "N") |> ModelComponents.Proxies.d15NtoAvailability)
            (p |> Pool.getEstimate "gamma[b]") (p |> Pool.getEstimate "r") (maxGrowthRate p) (limit p) (protectionEffect p e) (((temperatureDependency p e))) (seasonalLight t)) (fun s f -> cLog.StoreValue(s,t,f) |> ignore)

    /// Bristlecone function for dN/dt
    let dndt p t n (e:Environment) =
        dndt' (lookup e "bs") (n |> ModelComponents.Proxies.d15NtoAvailability) 
            (p |> Pool.getEstimate "gamma[n]") (maxGrowthRate p) (nitrogenFeedback p) (limit p) ((nReplenishment p e)) (temperatureDependency p e) (seasonalLight t) (protectionEffect p e) (fun s f -> cLog.StoreValue(s,t,f) |> ignore)

    /// Bristlecone function for dTs/dt
    let dtsdt p t ts (e:Environment) =

        // cLog.StoreValue("daylight",t,(seasonalLight t)) |> ignore
        // cLog.StoreValue("conductivity",t,(soilConductivity p e)) |> ignore
        // cLog.StoreValue("soilT",t,ts + (dtsdt' ts (lookup e "Tair") (soilConductivity p e))) |> ignore

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
                     code "r",         parameter PositiveOnly   2.000 8.000   // Intrinsic growth rate
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
        [ // 1. Linear rate not affected by soil temperatures
          (fun p _ -> ModelComponents.Temperature.Conductivity.linear 1.),
          (fun p _ -> ModelComponents.Temperature.NitrogenReplenishment.linear (p |> Pool.getEstimate "lambda")), []
          // 2. Temperature-dependent microbial activity
          (fun p _ -> ModelComponents.Temperature.Conductivity.linear 1.), //ModelComponents.Temperature.Conductivity.linear (p |> Pool.getEstimate "condct")),
          (fun p e -> NReplenishment.temperatureLimitation (p |> Pool.getEstimate "lambda") (p |> Pool.getEstimate "soilEa") (lookup e "Tsoil")),
            [ code "soilEa",        parameter PositiveOnly   5.000 10.00 ]
          // 3. Temperature-dependent microbial activity, with a snow-insulation effect on soil temperatures
          (fun p e -> ModelComponents.Temperature.Conductivity.snowConductivity (p |> Pool.getEstimate "condct") (lookup e "snowDepth") ),
          (fun p e -> NReplenishment.temperatureLimitation (p |> Pool.getEstimate "lambda") (p |> Pool.getEstimate "soilEa") (lookup e "Tsoil")),
            [ code "soilEa",    parameter PositiveOnly   20.00 30.00
              code "condct",    parameter PositiveOnly   0.001 0.200 ] ]

    // [B] Increased snow levels protect shrub biomass from storm and other damage.
    let snowProtection =
        [ (fun _ _ -> ModelComponents.SnowProtection.none), []
          (fun p e -> ModelComponents.SnowProtection.withShrubHeight (p |> Pool.getEstimate "spe") (lookup e "bs" |> biomassToHeight) (lookup e "snowDepth")),
            [ code "spe",        parameter PositiveOnly   0.001 1.000 ] ]

    // [C] Net photosynthetic rate is temperature-dependent
    let temperature =
        [ (fun _ -> ModelComponents.Temperature.none), []
          (fun p e -> NReplenishment.temperatureLimitation 1. (*(p |> Pool.getEstimate "A")*) (p |> Pool.getEstimate "Ea") (lookup e "Tair")),
           [ //code "A",              parameter PositiveOnly 5.00 20.0
             code "Ea",             parameter PositiveOnly 30.00 40.00  ] ]

    // [D] N may limited growth via combined N-limitations on (a) photosynthetic and (b) uptake rates
    let limitationModes =
        [ (fun p -> ModelComponents.GrowthLimitation.hollingDiscModelDual ((p |> Pool.getEstimate "a") / 1000.) ((p |> Pool.getEstimate "r") * 1000.) ((p |> Pool.getEstimate "h") / 1.) 5.00),
           [ code "a",      parameter PositiveOnly   0.100 5.000          // N-uptake efficiency
             code "h",      parameter PositiveOnly   0.001 0.250 ]
          (fun p -> ModelComponents.GrowthLimitation.linear ((p |> Pool.getEstimate "a") / 1000.) 10.), 
           [ code "a",      parameter PositiveOnly   0.010 0.100] ]       // N-uptake efficiency

    // [E] Loss of plant material may feedback into the soil pool of available nitrogen (instant)
    let feedbackModes =
        [ (fun p -> ModelComponents.FeedbackToSoil.withBiomassLoss ((p |> Pool.getEstimate "alpha") / 100.) (p |> Pool.getEstimate "gamma[b]") ),
          [ code "alpha",  parameter PositiveOnly   0.001 0.002 ]
          (fun _ -> ModelComponents.FeedbackToSoil.none), [] ]    // N-recycling efficiency

    // [F] A plant may be subject to mechanical constraints on its maximum size
    let geometricModes = 
        [  (fun p -> ModelComponents.GeometricConstraint.chapmanRichards ((p |> Pool.getEstimate "k") * 1000.)),
           [ code "k",  parameter PositiveOnly   3.00 5.00 ]
           (fun p -> ModelComponents.GeometricConstraint.none), [] ]     // Asymptotic biomass (grams)

    // Create all combinations of H1-H5 
    List.combine6 geometricModes feedbackModes limitationModes nitrogenReplenishment snowProtection temperature
    |> List.map (fun ((growth,gp),(feedback,fp),(limit,lp),(conduct,replace,rp), (snow,sp), (temp,tp)) -> 
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
                    if h < 25 then //h < 13 then // Just catches the 'most complex' base models for this paper..
                        yield async {
                                let cLog = ComponentLogger(true)
                                let result = Bristlecone.PlantIndividual.fit e Options.endWhen (hypotheses.[h-1] cLog startDate Options.latitude Options.longitude Options.timezone) shrubWithHighResEnvironment
                                Bristlecone.Data.EstimationResult.saveAll saveDirectory s.Identifier.Value h Options.thin (result |> fst)
                                return result |> fst }
    }

// Orchestrate the analyses
let shrubsWithIsotope = [ "YUSL03A"; "YUSL05A"; "YUSL26A"; "YUSL29A"; "YUSL39A" ]
let work = workPackages (shrubs |> Seq.where(fun s -> shrubsWithIsotope |> List.contains s.Identifier.Value)) hypotheses Options.engine Options.resultsDirectory |> Seq.toList

let run () =
    work |> Seq.iter (OrchestrationMessage.StartWorkPackage >> Options.orchestrator.Post)


module Diagnostics =

    open Bristlecone.Data

    /// Load all EstimationResults from a directory.
    let loadAll' directory subject (modelSystem:ModelSystem) modelId =
        let mles = MLE.load directory subject modelId |> Seq.map(fun (k,v) -> k.ToString(), v)
        let series = Series.load directory subject modelId |> Seq.map(fun (k,v) -> k.ToString(), v)
        let trace = Trace.load directory subject modelId |> Seq.map(fun (k,v) -> k.ToString(), v)
        mles
        |> Seq.keyMatch series
        |> Seq.map(fun (k,v1,v2) -> (k, (v1, v2)))
        // |> Seq.keyMatch trace
        // |> Seq.map(fun (k,v1,(v2,v3)) -> (k, (v1, v2, v3)))
        |> Seq.choose(fun (k,(s,(l,p))) ->
            try
                let parameters = 
                    // Convert 'b' into 'h'
                    modelSystem.Parameters |> Map.map(fun k v -> 
                        if Option.isNone (p |> Map.tryFind k)
                        then
                            printfn "Found a 'b' model..."
                            // h = b / (r * a)
                            let h = (p |> Map.find (code "b")) / ((p |> Map.find (code "r")) * (p |> Map.find (code "a")))
                            Parameter.setEstimate v h
                        else Parameter.setEstimate v (p |> Map.find k))
                { ResultId = k |> System.Guid.Parse
                  Likelihood = l
                  Parameters = parameters
                  Series = s
                  Trace = [] } |> Some
            with _ -> None )

    type Result = (PlantIndividual * (ComponentLogger -> System.DateTime -> float -> float -> string -> ModelSystem) * int * (EstimationResult seq * EstimationResult) option)

    /// Load all EstimationResults for the analyses in this study.
    let loadAll () : Result list =
        List.allPairs shrubs (hypotheses |> List.take 24 |> List.mapi(fun i v -> (i+1,v)))
        |> List.map (fun (s,(hi,h)) ->
            let r = loadAll' Options.resultsDirectory s.Identifier.Value (h (ComponentLogger(false)) (System.DateTime.Now) Options.latitude Options.longitude Options.timezone) hi
            printfn "[%s %i] %A" s.Identifier.Value hi (r |> Seq.map(fun r -> r.Likelihood))
            if r |> Seq.isEmpty
            then (s, h, hi, None )
            else 
                let r' = r |> Seq.filter(fun x -> not (System.Double.IsNaN(x.Likelihood)))
                if Seq.isEmpty r'
                then (s, h, hi, None )
                else 
                    (s, h, hi, (r', r' |> Seq.minBy(fun x -> x.Likelihood)) |> Some))

    let saveLog (cLog:ComponentLogger) (h:int) plantcode (guid:System.Guid) =
        let csv = System.IO.File.AppendText (sprintf "%sbristlecone-%s-%i-components-%s.csv" Options.resultsDirectory plantcode h (guid.ToString()))
        csv.WriteLine "Component,Time,Value"
        let x = cLog.GetAll()
        for m in x do
            for ts in m.Value do
                csv.WriteLine (sprintf "%s,%f,%f" m.Key ts.Key ts.Value)

    // A1. Save out clog for each using its parameters.
    let logOutComponents (results:Result list) =
        let engine = Options.engine |> Bristlecone.withCustomOptimisation (Optimisation.None.passThrough)
        results
        |> List.iter(fun (s,h,hi,e) -> 
            match e with
            | None -> ()
            | Some (_,mle) ->
                let cLog = ComponentLogger(false)
                let p = (mle.Parameters |> Map.map(fun k v -> parameter Unconstrained (v |> Parameter.getEstimate) (v |> Parameter.getEstimate)))
                    
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
                // TODO Remove this whole repeat of 1. here.

                let hypothesis = { (h cLog startDate Options.latitude Options.longitude Options.timezone) with Parameters = p }
                let result = Bristlecone.PlantIndividual.fit e Options.endWhen hypothesis shrubWithHighResEnvironment
                saveLog cLog hi s.Identifier.Value mle.ResultId )


    // B. Calculate Akaike weights
    let weights results = 
        results
        |> Seq.groupBy(fun (s,_,_,_) -> s.Identifier)
        |> Seq.collect(fun (_,r) ->
            let weights =
                r 
                |> Seq.choose(fun (s,h,hi,r) -> r)
                |> Seq.map(fun (chains,mle) -> mle)
                |> fun x -> printfn "Seq is %A" x; x
                |> ModelSelection.Akaike.akaikeWeights 
            r
            |> Seq.zip weights
            |> Seq.map (fun (w,(s,h,hi,r)) -> s, h, hi, fst w, snd w)
            |> Seq.toList )
        |> Seq.toList

    let saveWeights results =
        ModelSelection.save Options.resultsDirectory (results |> weights |> List.map(fun (x,y,z,a,b) -> x.Identifier.Value,z,a,b))

    module Convergence =

        type GelmanRubinStats = CsvProvider<"Subject (string), Hypothesis (int), Parameter (string), RHat (float)">

        let toCsvRows result =
            result |> Seq.concat |> Seq.map GelmanRubinStats.Row

        let calculate nBurn n (results:Result list) =
            results
            |> Seq.choose(fun (s,b,hi,r) ->
                printfn "Calculating Rhat for %s %i" (s.Identifier.Value) hi
                match r with
                | Some (results, mle) ->
                    let chains = 
                        results
                        |> Seq.sortBy(fun x -> x.Likelihood)
                        |> (fun c -> if c |> Seq.length > n then c |> Seq.take n else c)
                        |> Seq.map(fun r -> r.Trace |> Seq.map snd)
                    let minChainLength = chains |> Seq.minBy(fun x -> x |> Seq.length) |> Seq.length
                    if (chains |> Seq.length) < 2 then None
                    else
                        let chains' = chains |> Seq.map(fun c -> c |> Seq.take minChainLength)
                        printfn "Using %i chains (with %i common iterations)" (chains |> Seq.length) minChainLength
                        mle.Parameters
                        |> Map.toList
                        |> List.mapi(fun i (code,p) ->
                            s.Identifier.Value,
                            hi,
                            code.Value, 
                            chains'
                            |> Seq.map(fun c -> c |> Seq.map(fun x -> x.[i]))
                            |> Statistics.Convergence.GelmanRubin.rHat
                        ) |> Some
                | None -> None )

        let filename directory =
            let path = System.IO.DirectoryInfo(directory)
            if path.Exists then sprintf "%sbristlecone-rhat.csv" path.FullName
            else invalidArg "directory" "The specified directory does not exist"

        let save nBurn takeBestChains result =
            let csv = new GelmanRubinStats (result |> calculate nBurn takeBestChains |> toCsvRows)
            let filePath = filename Options.resultsDirectory
            csv.Save(filePath)


let saveDiagnostics () =
    let results = Diagnostics.loadAll()
    results |> Diagnostics.saveWeights
    results |> Diagnostics.logOutComponents
    results |> Diagnostics.Convergence.save 10000 3
