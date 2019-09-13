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


// TEMP

// Filzbach with Cauchy distribution?

/// An adaptation of the Filzbach method (originally by Drew Purves)
module Filzbach =

    open Bristlecone.Optimisation
    open Bristlecone.Optimisation.MonteCarlo
    open Bristlecone.Logging
    open Bristlecone.Statistics.Distributions
    
    type FilzbachSettings<'a> = {
        TuneAfterChanges: int
        MaxScaleChange: float
        MinScaleChange: float
        BurnLength: EndCondition<'a>
    }

    let filzbach' settings (theta:float[]) random writeOut (sampleEnd:EndCondition<float>) (domain:Domain) (f:Objective<float>) =
        writeOut <| GeneralEvent (sprintf "[Optimisation] Starting Filzbach-style MCMC optimisation")
        let sample sd = Normal.draw random 0. sd
        // let sample sd = fun () -> MathNet.Numerics.Distributions.Cauchy.Sample(random, 0., sd)
        let scaleRnd = Bristlecone.Statistics.Distributions.ContinuousUniform.draw random 0. 1.
        let paramRnd = MathNet.Numerics.Distributions.DiscreteUniform(0, theta.Length - 1, random)
        let initialScale = domain |> Array.map(fun (l,u,_) -> (u - l) / 6.)
        let l1 = f theta

        let rec step burning (p:(float*float[])[]) endWhen (l1, theta1) d =
            // Change one to many parameter at once (Filzbach-style)
            let scalesToChange = 
                if scaleRnd() < 0.670
                then
                    // Choose one parameter to change
                    let rnd = paramRnd.Sample()
                    // Change also nearby parameters with probability 1/2:
                    p |> Array.mapi(fun i x -> 
                        (x, if i = rnd then true
                            else if i = (rnd - 1) then scaleRnd() < 0.5
                            else if i = (rnd + 1) then scaleRnd() < 0.5
                            else false))
                else
                    // Probability of change:
                    let pChange = exp(4.0 * (scaleRnd() - 0.50))
                    // Try and allocate random changes to array
                    let rec changeRandom p =
                        let r = p |> Array.mapi(fun i x -> (x, scaleRnd() < pChange)) 
                        if r |> Array.where(fun (_,b) -> b) |> Array.isEmpty
                        then changeRandom p
                        else r
                    changeRandom p

            let propose theta =
                Array.zip3 theta scalesToChange domain
                |> Array.map(fun (x,((ti,n),shouldChange),(_,_,con)) -> 
                    if shouldChange
                    then constrainJump x (sample ti ()) 1. con
                    else x )

            // Metropolis step here
            let result = SimulatedAnnealing.tryMove propose (SimulatedAnnealing.Machines.boltzmann 1.) random f (l1, theta1)
            // End metropolis step
            
            // Tune Scales (burnin only)
            let newScaleInfo = 
                if not burning then p
                else 
                    scalesToChange 
                    |> Array.zip (result |> snd)
                    |> Array.mapi(fun parameteri (v, ((ti,previous),changed)) ->
                        let ti, previous = 
                            if changed then (ti, (previous |> Array.append [|v|]))  // Append new parameter values to previous ones
                            else (ti, previous)
                        if previous |> Array.length = settings.TuneAfterChanges
                        then
                            let changes = previous |> Array.pairwise |> Array.where(fun (a,b) -> a <> b) |> Array.length
                            match (float changes) / (float settings.TuneAfterChanges) with
                            | ar when ar < 0.25 -> 
                                if (ti * 0.80) < (initialScale.[parameteri] * settings.MinScaleChange)
                                then (initialScale.[parameteri] * settings.MinScaleChange, Array.empty)
                                else (ti * 0.80, Array.empty)
                            | ar when ar > 0.25 -> 
                                if (ti * 1.20) > (initialScale.[parameteri] * settings.MaxScaleChange)
                                then (initialScale.[parameteri] * settings.MaxScaleChange, Array.empty)
                                else (ti * 1.20, Array.empty)
                            | _ -> (ti, Array.empty)
                        else (ti, previous) )
            // End Tune Scales (burnin only)

            let newResult = result::d
            if endWhen d
            then (newResult, newScaleInfo |> Array.map fst)
            else 
                writeOut <| OptimisationEvent { Iteration = newResult |> List.length; Likelihood = result |> fst; Theta = result |> snd }
                step burning newScaleInfo endWhen result newResult

        writeOut <| GeneralEvent (sprintf "[Filzbach] Starting burn-in at point %A (L = %f)" theta l1)
        let burnResults,burnScales = step true (initialScale |> Array.map(fun x -> x, Array.empty)) settings.BurnLength (l1, theta) []
        writeOut <| GeneralEvent (sprintf "[Filzbach] Burn-in complete. Starting sampling at point %A" (burnResults |> Seq.head))
        let results,_ = step false (burnScales |> Array.map(fun t -> (t, Array.empty))) sampleEnd (burnResults |> Seq.head) []
        [ results; burnResults ] |> List.concat

    let filzbachCauchy settings writeOut n domain (f:Objective<float>) : (float * float[]) list =
        let random = MathNet.Numerics.Random.MersenneTwister(true)
        match tryGenerateTheta f domain random 10000 with
        | Ok theta ->
            writeOut <| GeneralEvent (sprintf "[Optimisation] Initial theta is %A" theta)
            filzbach' settings theta random writeOut n domain f
        | Error _ -> invalidOp "Could not generate theta"

// END TEMP






// 1. Configure Options
// ----------------------------

module Options =
    let resultsDirectory = "/Users/andrewmartin/Desktop/Bristlecone Results/Paper2-Snow-Short-Reformulated4/"
    let endWhen = Optimisation.EndConditions.afterIteration 25000
    let chains = 1
    let thin = 50
    let engine =
        Bristlecone.mkContinuous 
        |> Bristlecone.withContinuousTime Integration.MathNet.integrate
        // |> Bristlecone.withTunedMCMC [ Optimisation.MonteCarlo.TuneMethod.CovarianceWithScale 0.250, 500, Optimisation.EndConditions.afterIteration 75000 ]
        |> Bristlecone.withCustomOptimisation (Filzbach.filzbachCauchy 
            { TuneAfterChanges = 50
              MaxScaleChange = 100.00
              MinScaleChange = 0.0010
              BurnLength = Optimisation.EndConditions.afterIteration 25000 })

    let orchestrator = OrchestrationAgent(engine.LogTo, System.Environment.ProcessorCount, false)
    let latitude = 68.91    // Yuribei North
    let longitude = -70.23  // Yuribei West
    let timezone = "Asia/Yekaterinburg"    // Yuribei is in Yekaterinburg Time


// 2. Create Hypotheses
// ----------------------------

module BaseEquations =

    let biomass b n r gammab geom f (protectionEffect:float) tempEffect lightEffect cLog : float =

        cLog "stem" (b |> ModelComponents.Proxies.toRadiusMM |> ModelComponents.Proxies.shrubHeightCm)
        cLog "biomassGrowth" (b * r * (f n) * geom(b) * tempEffect * lightEffect)
        cLog "biomassLoss" (gammab * (1. - protectionEffect) * b)

        b * r * (f n) * geom(b) * tempEffect * lightEffect - gammab * (1. - protectionEffect) * b

    let soilNitrogen n b gamman y geom f feedback tempEffect lightEffect (protectionEffect:float) cLog : float =

        cLog "nReplenishment" y
        cLog "nUptake" (geom(b) * b * (f n) * tempEffect * lightEffect)
        cLog "nLoss" (gamman * n + feedback(b) * (1. - protectionEffect))

        y - (geom(b) * b * (f n) * tempEffect * lightEffect) - gamman * n + feedback(b) * (1. - protectionEffect)

    // Newton's law of cooling / Fourier's law
    let soilTemperature soilT ambientT conductivity =
        if conductivity > 1. || conductivity < 0.05 then nan
        else - conductivity * (soilT - ambientT)


let ``base model`` maxGrowthRate nLimitation nitrogenFeedback nReplenishment protectionEffect temperatureDependency soilConductivity additionalParameters (cLog:IComponentLogger<float>) ((startDate:System.DateTime),latitude,longitude,timeZone) =

    let limit p =
        match nLimitation p with
        | Some l -> l
        | None -> invalidOp "N limitation is required for this model."

    /// Light limitation effect (linear between 0 and 1).
    /// Light is cached in an object to avoid unnecessary computation.
    let lightFn = Sunrise.DayLengthCache(latitude, longitude, timeZone)
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

        cLog.StoreValue "protectionEffect" t ((protectionEffect p e)) |> ignore
        cLog.StoreValue "biomass" t bs |> ignore
        cLog.StoreValue "temperature-limit-on-photosynthesis" t (temperatureDependency p e) |> ignore

        (dbsdt' bs ((lookup e "N") |> ModelComponents.Proxies.d15NtoAvailability)
            (p |> Pool.getEstimate "gamma[b]") (p |> Pool.getEstimate "r") (maxGrowthRate p) (limit p) (protectionEffect p e) (((temperatureDependency p e))) (seasonalLight t)) (fun s f -> cLog.StoreValue s t f |> ignore)

    /// Bristlecone function for dN/dt
    let dndt p t n (e:Environment) =
        dndt' (lookup e "bs") (n |> ModelComponents.Proxies.d15NtoAvailability) 
            (p |> Pool.getEstimate "gamma[n]") (maxGrowthRate p) (nitrogenFeedback p) (limit p) ((nReplenishment p e)) (temperatureDependency p e) (seasonalLight t) (protectionEffect p e) (fun s f -> cLog.StoreValue s t f |> ignore)

    /// Bristlecone function for dTs/dt
    let dtsdt p t ts (e:Environment) =

        cLog.StoreValue "daylight" t (seasonalLight t) |> ignore
        cLog.StoreValue "conductivity" t (soilConductivity p e) |> ignore
        cLog.StoreValue "soilT" t (ts + (dtsdt' ts (lookup e "Tair") (soilConductivity p e))) |> ignore

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

// Plant process...
// ====================

// 1. Allometric fit
// 2. Measurement variable: ring width


module TreeRing =

    [<Measure>] type g

    // Convert plant radius into biomass and vice versa.
    // Validate allometrics so that they can never be negative.
    // Validate allometrics so that they do the symmetrical transform.
    type Allometrics = {
        StemRadiusToBiomass: float<mm> -> float<g>
        BiomassToStemRadius: float<g> -> float<mm>
    }

    /// A base `ModelSystem` that, given a biomass equation, outputs
    /// ring widths as a measurement series. 
    let treeRing biomassCode biomassFn parameters (allometric:Allometrics) l =

        /// Measurement (Size) variable: stem radius
        let stemRadius lastRadius lastEnv env =
            let oldCumulativeMass = lookup lastEnv biomassCode
            let newCumulativeMass = lookup env biomassCode
            if (newCumulativeMass - oldCumulativeMass) > 0.
            then 
                let newRadius = newCumulativeMass * 1.<g> |> allometric.BiomassToStemRadius
                if (newRadius |> removeUnit) > lastRadius then (newRadius |> removeUnit) else lastRadius
            else lastRadius

        // A. There should be no repeat parameters.
        let parameters = parameters |> List.concat

        { Equations  = [ code biomassCode, biomassFn ] |> Map.ofList
          Measures   = [ code "x", stemRadius ] |> Map.ofList
          Parameters = parameters |> Map.ofList
          Likelihood = l }

    let requiresVariable c eq h =
        { h with Equations = h.Equations |> Map.add (code c) eq }



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

type NestedComponent<'value> = {
    Parameters: Parameter list
    Function: ParameterPool -> Environment -> 'value
}


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
    let meanTemperatures = DailyTemperature.Load ClimateUrl
    meanTemperatures.Rows
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
    let initialAirT = plant.Environment.[code "Tair"].Head |> fst
    let initialSnow = plant.Environment.[code "snowDepth"].Head |> fst
    [ code "x", initialRadius
      code "N", initialNitrogen 
      code "Tair", initialAirT
      code "Tsoil", initialAirT
      code "snowDepth", initialSnow
      code "bs", initialMass ] |> Map.ofList

/// How to generalise the `fit` function?
/// - Start Date = common 'variables' (as opposed to environmental data). Codes can be got from ModelSystem.
/// - End Date = common 'variables' (as opposed to environmental data).
/// - Ensure that start values are present (at t-1, including for env variables).
let fit s hypothesis engine = 
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
    Bristlecone.PlantIndividual.fit e Options.endWhen (hypothesis (startDate,Options.latitude,Options.longitude,Options.timezone)) shrubWithHighResEnvironment

let setupHypothesis s hypothesis =
    let shrub = s |> PlantIndividual.toCumulativeGrowth
    let common = shrub |> PlantIndividual.keepCommonYears
    let startDate = (common.Environment.[code "N"]).StartDate |> snd
    hypothesis (startDate,Options.latitude,Options.longitude,Options.timezone)


/// Use this type in Bristlecone itself.
type WorkPackage = Async<EstimationResult>

/// How to generalise this function?
/// - It is just a nested set of commands.
let workPackages shrubs hypotheses engine saveDirectory =
    seq {
        for s in shrubs do
            for h in [ 1 .. hypotheses |> List.length ] do
                for _ in [ 1 .. Options.chains ] do
                    if h = 12 || h = 24 then
                        yield async {
                                let cLog = PassThrough()
                                let result = fit s (hypotheses.[h-1] cLog) engine
                                Bristlecone.Data.EstimationResult.saveAll saveDirectory s.Identifier.Value h Options.thin (result |> fst)
                                return result |> fst }
    }

// Orchestrate the analyses
let shrubsWithIsotope = [ "YUSL03A"; "YUSL05A"; "YUSL26A"; "YUSL29A"; "YUSL39A" ]
let work = workPackages (shrubs |> Seq.where(fun s -> shrubsWithIsotope |> List.contains s.Identifier.Value)) hypotheses Options.engine Options.resultsDirectory |> Seq.toList

let run () =
    work |> Seq.map(fun x -> x,System.Random().Next()) |> Seq.sortBy snd |> Seq.map fst |> Seq.iter (OrchestrationMessage.StartWorkPackage >> Options.orchestrator.Post)


let saveDiagnostics () =

    // 1. Get all results sliced by plant and hypothesis
    let results = 
        let hypotheses = hypotheses |> List.map(fun h -> h (PassThrough()) ((System.DateTime.Now,Options.latitude,Options.longitude,Options.timezone)))
        let get subject model modelId = Data.EstimationResult.loadAll Options.resultsDirectory subject.Identifier.Value model modelId
        ModelSelection.ResultSet.arrangeResultSets shrubs hypotheses get

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
    |> Data.ModelSelection.save Options.resultsDirectory

    // 4. Save out logged components
    // results
    // |> Seq.map(fun r -> calculateComponents fit Options.engine r)


// Take the base case and run an optimisation routine from there.
// - AKA hypothesis 13 for every individual
// - Introduce new factor in with effective zero effect


// Start chains that have already finished and run again.
// Uses MLE for each shrub and hypothesis.


// 1. Get all results sliced by plant and hypothesis
let results = 
    let get subject model modelId = Data.EstimationResult.loadAll Options.resultsDirectory subject.Identifier.Value (model (PassThrough()) ((System.DateTime.Now,Options.latitude,Options.longitude,Options.timezone))) modelId
    Bristlecone.ModelSelection.ResultSet.arrangeResultSets shrubs hypotheses get

let workPackage (shrub:PlantIndividual) hi hypothesis engine saveDirectory =
    async {
        let cLog = PassThrough()
        let result = fit shrub (hypothesis cLog) engine |> fst
        Bristlecone.Data.EstimationResult.saveAll saveDirectory shrub.Identifier.Value hi 100 result
        return result
    }

// Make all analyses start with best H13 parameter set.
let basicModel = 13



let work2 =
    results
    |> Seq.collect(fun (plant,h,hi,r) ->
        match r with
        | Some r ->
            let sndBestr = if r |> fst |> Seq.length > 1 then r |> fst |> Seq.sortBy(fun x -> x.Likelihood) |> Seq.skip 1 |> Seq.head else snd r
            let mleToBounds mlePool = mlePool |> Map.map(fun k v -> Parameter.create (Parameter.detatchConstraint v |> snd) (v |> Parameter.getEstimate) (v |> Parameter.getEstimate))
            let hypothesisMle a b = { h a b with Parameters = mleToBounds sndBestr.Parameters }
            [ workPackage plant hi hypothesisMle Options.engine Options.resultsDirectory ]
        | None -> [] )

let runAgain () =
    work2 |> Seq.iter (OrchestrationMessage.StartWorkPackage >> Options.orchestrator.Post)
