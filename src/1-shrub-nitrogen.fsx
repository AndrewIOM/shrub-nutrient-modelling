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
open Bristlecone.PlantIndividual
open Bristlecone.Workflow.Orchestration

// 1. Configure Options
// ----------------------------

module Options =
    let resultsDirectory = "/Users/andrewmartin/Desktop/Bristlecone Results/YuribeiAnnual-TemperatureDependent-Tuned/"
    let chains = 3
    let endWhen = Optimisation.EndConditions.afterIteration 1000000
    let logger = Logging.RealTimeTrace.graphWithConsole 30. 10000
    let engine =
        Bristlecone.mkContinuous 
        |> Bristlecone.withContinuousTime Integration.MathNet.integrate
        |> Bristlecone.withOutput logger
        |> Bristlecone.withCustomOptimisation (Optimisation.MonteCarlo.SimulatedAnnealing.fastSimulatedAnnealing 0.0001 
            { Optimisation.MonteCarlo.SimulatedAnnealing.AnnealSettings<float>.Default with 
                BoilingAcceptanceRate = 0.60
                HeatRamp = (fun t -> t + sqrt t); TemperatureCeiling = Some 500.
                HeatStepLength = Optimisation.EndConditions.afterIteration 1000
                AnnealStepLength = (fun x -> Optimisation.MonteCarlo.SimulatedAnnealing.EndConditions.improvementCount 5000 250 x || Optimisation.EndConditions.afterIteration 15000 x) }) //(Optimisation.EndConditions.afterIteration 10000) })

    let orchestrator = OrchestrationAgent(logger, System.Environment.ProcessorCount)


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
                     code "sigma[x]",  parameter PositiveOnly   0.001 0.100   // Standard deviation of x (biomass)
                     code "sigma[y]",  parameter PositiveOnly   0.001 0.100   // Standard deviation of y (nitrogen)
                    ] |> List.append additionalParameters |> Map.ofList
      Likelihood = ModelLibrary.Likelihood.bivariateGaussian "x" "N" }

let hypotheses =

    // [A] N may limited growth via combined N-limitations on (a) photosynthetic and (b) uptake rates
    let limitationModes =
        [ (fun p -> ModelComponents.GrowthLimitation.lizzy ((p |> Pool.getEstimate "a") / 1000.) ((p |> Pool.getEstimate "h") / 1.)), //.hollingDiscModel ((p |> Pool.getEstimate "a") / 1000.) 1. ((p |> Pool.getEstimate "h") / 1.)), //Removed r: ((p |> Pool.getEstimate "r") * 1000.)
           [ code "a",      parameter PositiveOnly   0.100 0.400          // N-uptake efficiency
             code "h",      parameter PositiveOnly   0.100 0.400
             code "r",      parameter PositiveOnly   0.500 1.000 ]     // N-handling time (including uptake and incorporation)
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
// let startValues = [ ShortCode.create "x", 5.; ShortCode.create "N", 3.64; ShortCode.create "bs", (5. |> ModelComponents.Proxies.toBiomassMM)] |> Map.ofList
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
    |> Seq.map (fun (_,plant,d15N) -> PlantIndividual.zipEnv (ShortCode.create "N") plant d15N)
    |> Seq.toList

let getStartValues (startDate:System.DateTime) (plant:PlantIndividual) =
    let initialRadius =
        match plant.Growth with
        | PlantIndividual.PlantGrowth.RingWidth s -> 
            match s with
            | GrowthSeries.Absolute c -> c.Head |> fst |> removeUnit
            | GrowthSeries.Cumulative c -> 
                let start = (c |> TimeSeries.trimStart (startDate - System.TimeSpan.FromDays(366.))).Values |> Seq.head |> removeUnit
                // printfn "Start cumulative growth = %f" start
                start
            | GrowthSeries.Relative _ -> invalidOp "Not implemented"
        | _ -> invalidOp "Not implemented 2"
    let initialMass = initialRadius |> removeUnit |> ModelComponents.Proxies.toBiomassMM
    let initialNitrogen = plant.Environment.[ShortCode.create "N"].Head |> fst
    [ (ShortCode.create "x", initialRadius)
      (ShortCode.create "N", initialNitrogen)
      (ShortCode.create "bs", initialMass) ] |> Map.ofList

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

// Orchestrate the analyses
let work = workPackages (shrubs |> Seq.skip 0) hypotheses Options.engine Options.resultsDirectory
// work |> Seq.take 3 |> Seq.iter (OrchestrationMessage.StartWorkPackage >> Options.orchestrator.Post)

// // work 
// // |> Seq.toArray
// // |> Array.chunkBySize 1
// // |> Array.take 1
// // |> Array.map(fun p -> p |> Array.Parallel.map(fun w -> w |> Async.RunSynchronously))

// open RProvider
// open RProvider.graphics

// let pi = System.Math.PI

// /// Total shrub volume given height and number of stems
// let vol b a rtip p lmin k5 k6 n h =

//     let radius = ModelComponents.NiklasAndSpatz_Allometry.basalRadius k5 k6
//     let mainStemVolume =
//         match radius h with
//         | r when r > rtip -> n * pi * h * ((radius h) ** 2. + (radius h) * rtip + rtip ** 2.) / 3.
//         | r -> n * pi * h * rtip ** 2.

//     let mutable volume = mainStemVolume
//     let mutable k = 0.

//     while (p ** k * h > lmin * 2./3.) do
//         let volToAdd =
//             printfn "[1] - %f" (n * a * ((a + 1.) ** (float k)) * pi * 3. * p * ((p ** k) * h - 2. * lmin / 3.) * ((radius (3. * p * ((p ** k) * h - 2. * lmin / 3.))) * (radius (3. * p * ((p ** k) * h - 2. * lmin / 3.)) * rtip + (rtip ** 2.))) / 3.)
//             printfn "[2] - %f" (n * a * ((a + 1.) ** (float k)) * 3. * p * ((p ** k) * h - 2. * lmin / 3.) * pi * (rtip ** 2.))
//             printfn "[3] - %f" (n * a * ((a + 1.) ** (float k)) * pi * (p ** (k + 1.)) * h * (((radius ((p ** (k+1.)) * h)) ** 2.) + (radius ((p ** (k + 1.)) * h)) * rtip + (rtip ** 2.)) / 3.)
//             printfn "[4] - %f" (n * a * ((a + 1.) ** (float k)) * (p ** (k + 1.)) * h * pi * (rtip ** 2.))
            
//             // If the branch to add is less than min branch length (lmin),
//             // then scale the amount to add continuously so that the length
//             // of the child branch is 3*p*(p^k*h-2*l_min/3).
//             printfn "[%f] Branch length to add: %f" mainStemVolume (p ** k * h)
//             match (p ** k * h < lmin) with
//             | true ->
//                 // Adding a branch with length less than l_min
//                 // Add the volume of the child branches to every
//                 // parent branch / stem.
//                 printfn "[Test 1] - %f" ((b * 3. * p * (p ** k * h - 2. * lmin / 3.)))
//                 match (b * 3. * p * (p ** k * h - 2. * lmin / 3.) > rtip) with
//                 | true ->
//                     // Length of child branch is more than r_tip
//                     // Model child branch as a cone
//                     printfn "1 - Radius = %f" (radius h)
//                     n * a * (a + 1.) ** (float k) * pi * 3. * p * (p ** k * h - 2. * lmin / 3.) * ((radius (3. * p * (p ** k * h - 2. * lmin / 3.))) * (radius (3. * p * (p ** k * h - 2. * lmin / 3.)) * rtip + rtip ** 2.)) / 3.
//                 | false ->
//                     // Length of child branch is less than r_tip
//                     // Model child branch as a cylinder
//                     printfn "2 - Radius = %f" (radius h)
//                     //V_s = V_s + n*a_s*(a_s+1)^k * 3*p*(p^k*h_s-2*l_min/3)*pi*r_tip^2;
//                     n * a * (a + 1.) ** (float k) * 3. * p * (p ** k * h - 2. * lmin / 3.) * pi * rtip ** 2.
//             | false ->
//                 printfn "[Test 2] - %f" (radius (p ** (k + 1.) * h))
//                 match (radius (p ** (k + 1.) * h) > rtip) with
//                 | true ->
//                     printfn "Ratio = %f" (radius (p ** (k + 1.) * h))
//                     // b * branch length > r_tip (max radius of branch tip)
//                     // Model child branch as a truncated cone
//                     printfn "3 - Radius = %f" (radius h)
//                     // V_s=V_s + n*a_s*(a_s+1)^k*pi*p^(k+1)*h_s*((rad(p^(k+1)*h_s))^2+rad(p^(k+1)*h_s)*r_tip+r_tip^2)/3;
//                     printfn "Component1 = %f" (n * a * (a + 1.) ** (float k) * pi * p ** (k + 1.) * h)
//                     printfn "Component2 = %f" (((radius (p ** (k+1.) * h)) ** 2. + (radius (p ** (k + 1.) * h)) * rtip + rtip ** 2.))
//                     n * a * (a + 1.) ** (float k) * pi * p ** (k + 1.) * h * ((radius (p ** (k+1.) * h)) ** 2. + (radius (p ** (k + 1.) * h)) * rtip + rtip ** 2.) / 3.
//                 | false ->
//                     // Model child branch as a cylinder
//                     printfn "4 - Radius = %f" (radius h)
//                     n * a * (a + 1.) ** (float k) * p ** (k + 1.) * h * pi * rtip ** 2.
//         printfn "Adding %f" volToAdd
//         volume <- volume + volToAdd
//         k <- k + 1.

//     k, volume
 
// let radius = [ 1. .. 1. .. 100. ] // in mm
// let height = radius |> List.map (fun r -> ModelComponents.Allometrics.shrubHeight Constants.Allometrics.k5 Constants.Allometrics.k6 (r / 10.))
// let volume = height |> List.map (fun r -> vol Constants.Allometrics.b Constants.Allometrics.a Constants.Allometrics.rtip Constants.Allometrics.p Constants.Allometrics.lmin Constants.Allometrics.k5 Constants.Allometrics.k6 2. r)
// let biomass = radius |> List.map (fun r -> ModelComponents.Proxies.toBiomassMM r)
// volume |> List.zip radius
// let rgr = biomass |> List.pairwise |> List.map(fun (b1,b2) -> (b2 - b1) / b1)

// R.plot(radius, height)
// R.plot(volume |> List.map snd, height)
// R.lines(radius, volume |> List.map fst)
// R.plot(radius, biomass)

// // A. Volume vs Height
// R.plot(volume |> List.map snd, height)

// // D. Volume vs Twig Count
// let twigCount = volume |> List.map (fun (t,_) -> 3. * 2. ** t)
// R.plot(volume |> List.map snd, twigCount)

// R.plot(radius |> List.tail, rgr)
