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
open FSharp.Data

type BristleconeResult = CsvProvider<Sample = "PlantCode (string), Hypothesis (int), Iteration (int), Chain (int), Parameter (string), Likelihood (float), value (float)", CacheRows = true>


#r "../lib/DiffSharp.dll"
module HMC =

    open DiffSharp.AD.Float64
    open Bristlecone.Logging
    open Bristlecone.Optimisation.MonteCarlo
    open Bristlecone.Optimisation

    /// Leapfrog integrator
    /// `u`: potential energy function
    /// `k`: kinetic energy function
    /// `d`: integration step size
    /// `steps`: number of integration steps
    /// `(x0, p0)`: initial position and momentum vectors
    let leapFrog (u:DV->D) (k:DV->D) (d:D) steps (x0, p0) =
        let hd = d / 2.
        [1 .. steps] 
        |> List.fold (fun (x, p) _ ->
            let p' = p - hd * grad u x
            let x' = x + d * grad k p'
            x', p' - hd * grad u x') (x0, p0)

    // let Rnd = new System.Random()

    // // Uniform random number ~U(0, 1)
    // let rnd() = Rnd.NextDouble()

    // // Standard normal random number ~N(0, 1), via Box-Muller transform
    // let rec rndn() =
    //     let x, y = rnd() * 2.0 - 1.0, rnd() * 2.0 - 1.0
    //     let s = x * x + y *y
    //     if s > 1.0 then rndn() else x * sqrt (-2.0 * (log s) / s)

    // /// Hamiltonian Monte Carlo
    // /// n: number of samples wanted
    // /// hdelta: step size for Hamiltonian dynamics
    // /// hsteps: number of steps for Hamiltonian dynamics
    // /// x0: initial state
    // /// f: target distribution function
    // let hmc' n hdelta hsteps (x0:DV) (f:DV->D) =
    //     let u x = -log (f x) // potential energy
    //     let k p = (p * p) / D 2. // kinetic energy
    //     let hamilton x p = u x + k p
    //     let x = ref x0
    //     [|for _ in 1..n do
    //         let p = DV.init x0.Length (fun _ -> rndn()) // PROPOSAL: Lift onto random level set
    //         //printfn "P: %A" p
    //         let x', p' = leapFrog u k hdelta hsteps (!x, p) // PRPOPSAL: Explore with hamiltonian trajectory, and project back down to target parameter space
    //         //printfn "X+P: %A %A" x' p'
    //         if rnd() < float (exp ((hamilton !x p) - (hamilton x' p'))) then x := x' // TEST: Metropolis test as per other algorithms
    //         yield !x|] // yields theta, but not likelihood

    // // Just get prposal using Hamiltonian dynamics (DiffSharp)
    // let propose hdelta hsteps (x0:DV) f =
    //     let u x = 
    //         printfn "Potential energy calc on %A" x
    //         try f x |> ignore with
    //         | e -> printfn "AAH %A" e
    //         -log (f x) // potential energy
    //     let k p = 
    //         printfn "Kinetic energy calc on %A" p
    //         (p * p) / D 2. // kinetic energy
    //     let hamilton x p = u x + k p
    //     let p = DV.init x0.Length (fun _ -> rndn()) // PROPOSAL: Lift onto random level set
    //     printfn "Proposal: %A" p
    //     let x', p' = leapFrog u k hdelta hsteps (x0, p) // PRPOPSAL: Explore with hamiltonian trajectory, and project back down to target parameter space
    //     x'

    let hmc (hdelta:float) hsteps writeOut n domain (f:Point->float) : (float * float[]) list =
        writeOut <| GeneralEvent (sprintf "[Optimisation] Started Hamiltonian Monte Carlo")

        // Get initial theta
        let random = MathNet.Numerics.Random.MersenneTwister(true)
        let theta = initialise domain random
        let l = f theta
        writeOut <| GeneralEvent (sprintf "[Optimisation] Initial theta is %A" theta)

        // HMC
        let hmc' n hdelta hsteps (x0:DV) (f:DV->D) =
            let u x = (*-log*) (f x) // potential energy
            let k p = (p * p) / D 2. // kinetic energy
            let hamilton x p = u x + k p // energy function
            let x = ref x0
            [|for i in 1..n do
                let p = DV.init x0.Length (fun _ -> random.NextDouble()) // PROPOSAL: Lift onto random level set

                let randSteps = MathNet.Numerics.Distributions.DiscreteUniform(1,hsteps).Sample()
                let x', p' = leapFrog u k hdelta randSteps (!x, p) // PRPOPSAL: Explore with hamiltonian trajectory, and project back down to target parameter space
                
                let rand = MathNet.Numerics.Distributions.ContinuousUniform(0.,1.).Sample()
                let l1 = hamilton !x p
                let l2 = hamilton x' p'
                if rand < float (exp ((hamilton !x p) - (hamilton x' p'))) 
                then 
                    printfn "ACCEPT. %A -> %A" l1 l2
                    x := x' // TEST: Metropolis test as per other algorithms
                else
                    printfn "REJECT. %A -> %A" l1 l2
                yield !x

                let l = u !x |> float
                let thetaWithLikelihood =
                    x.Value 
                    |> DV.toArray 
                    |> Array.map float
                    |> Array.append [|l|]
                writeOut <| OptimisationEvent { Iteration = i; Likelihood = l; Theta = thetaWithLikelihood }                
                |]
            // TODO It's going up, not down (aka for -logL)

        let fd = fun (t:DV) -> 
            t |> DV.toArray |> Array.map float |> f |> D

        let x = hmc' n (D hdelta) hsteps (toDV theta) fd

        printfn "X is %A" x
        failwith "Not finished"
        //let likelihoods = x |> Array.unzip |> fst |> Array.map float
        //let thetas = x |> Array.unzip |> snd |> Array.map DV.toArray
        //y

        // // Algorithm for generating proposals
        // let proposeJump _ (theta:float[]) =
        //     try
        //         let fd = fun (t:DV) -> 
        //             printfn "Starting f..."
        //             let r = t |> DV.toArray |> Array.map float |> f |> D
        //             printfn "R = %A" r
        //             r
        //         printfn "Got to here"
        //         let theta1 = propose (D hdelta) hsteps (toDV theta) fd
        //         printfn "Got even further"
        //         theta1 |> DV.toArray |> Array.map float
        //     with
        //     | e -> printfn "%A" e; invalidOp "Crashed"

        // // Run an MCMC with the HMC proposal mechanism
        // let result = metropolisHastings' writeOut proposeJump TuningMode.none f theta l n [] ()
        // result |> fst

// open DiffSharp.AD.Float64

// let multiNormal (mu:DV) (sigma:DM) (x:DV) =
//     let s = sigma |> DM.inverse
//     exp (-((x - mu) * s * (x - mu)) / D 2.)

// let iterations = 10000
// let hStepSize = 0.1
// let hStepCount = 10
// let targetDistributionFuction = 
//     multiNormal (toDV [0.; 0.]) (toDM [[1.; 0.8]; [0.8; 1.]])

// let theta = toDV [0.; 0.]
// let samples = 
//     HMC.hmc' iterations (D hStepSize) hStepCount theta targetDistributionFuction

// // Plot HMC sample
// open RProvider
// open RProvider.graphics
// R.plot(samples |> Array.map (fun v -> float v.[0]), samples |> Array.map (fun v -> float v.[1]))


// 1. Configure Options
// ----------------------------

module Options =
    let resultsDirectory = "/Users/andrewmartin/Desktop/"
    let iterations = 100000
    let chains = 6
    let logger = 
        //let consolePost = Bristlecone.Logging.Console.logger(30.)
        let graphLog = Bristlecone.Logging.RealTimeTrace.TraceGraph(Logging.Device.X11,60.,5000)
        (fun event -> (*consolePost event;*) graphLog.Log event)

    let engine =
        Bristlecone.mkContinuous 
        |> Bristlecone.withContinuousTime Integration.MathNet.integrate
        |> Bristlecone.withOutput logger
        |> Bristlecone.withCustomOptimisation (HMC.hmc 0.000001 10)


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

    let randomLog environment parameters =
        if System.Random().Next(0,500000) = 1 then
            printfn "RANDOM LOG: Environment = %A; Parameters = %A" environment parameters

    /// Cumulative stem biomass [b].
    let dbsdt' (b:float) n gammab r maxGrowthRate limit =
        match limit with
        | Some l -> BaseEquations.biomass b n r gammab maxGrowthRate l
        | None -> BaseEquations.biomass b n r gammab maxGrowthRate (fun _ -> 1.)

    /// Bioavailable soil nitrogen [N]
    let dndt' bs n lambda gamman maxGrowthRate feedback limit = 
        match limit with
        | Some l -> BaseEquations.soilNitrogen n bs gamman lambda maxGrowthRate l feedback
        | None -> BaseEquations.soilNitrogenNoUptake n bs gamman lambda feedback

    /// Measurement variable: stem radius [rw].
    let drwdt' bs n r gammab maxGrowthRate limit =
        let biomassStemChange = dbsdt' bs n gammab r maxGrowthRate limit
        if biomassStemChange > 0.
            then
                let oldRadius = bs |> ModelComponents.Proxies.toRadiusMM
                let newRadius = (bs + biomassStemChange) |> ModelComponents.Proxies.toRadiusMM
                newRadius - oldRadius
            else 0.

    /// Bristlecone function for dBs/dt
    let dbsdt p _ bs (e:Environment) =
        randomLog e p
        dbsdt' bs ((e.[ShortCode.create "N"]) |> ModelComponents.Proxies.d15NtoAvailability)
            (p |> Pool.getEstimate "gamma[b]") (p |> Pool.getEstimate "r") (maxGrowthRate p) (nLimitation p)

    /// Bristlecone function for dN/dt
    let dndt p _ n (e:Environment) =
        dndt' (e.[ShortCode.create "bs"]) (n |> ModelComponents.Proxies.d15NtoAvailability) 
            (p |> Pool.getEstimate "lambda") (p |> Pool.getEstimate "gamma[n]") (maxGrowthRate p) (nitrogenFeedback p) (nLimitation p)

    /// Bristlecone function for dr/dt
    let drwdt p _ _ (e:Environment) =
        drwdt' (e.[ShortCode.create "bs"]) ((e.[ShortCode.create "N"]) |> ModelComponents.Proxies.d15NtoAvailability) 
            (p |> Pool.getEstimate "r") (p |> Pool.getEstimate "gamma[b]") (maxGrowthRate p) (nLimitation p)

    { Equations  = [ code "x",         drwdt
                     code "bs",        dbsdt
                     code "N",         dndt ] |> Map.ofList
      Parameters = [ // for nitrogen dynamics
                     code "lambda",    parameter PositiveOnly   0.001 0.500   // Rate of nitrogen replenishment
                     code "gamma[n]",  parameter PositiveOnly   0.001 0.200   // Loss rate of nitrogen
                     // for shrub physiology
                     code "gamma[b]",  parameter PositiveOnly   0.001 0.200   // Loss rate of biomass
                     //code "r",         parameter PositiveOnly   0.500 1.000   // Photosynthetic efficiency (either N-limited or N-unlimited)
                     // for likelihood function
                     code "rho",       parameter Unconstrained  -0.50 0.500   // Covariance between growth and nitrogen
                     code "sigma[x]",  parameter PositiveOnly   0.100 1.200   // Standard deviation of x (biomass)
                     code "sigma[y]",  parameter PositiveOnly   0.250 0.750   // Standard deviation of y (nitrogen)
                    ] |> List.append additionalParameters |> Map.ofList
      Likelihood = ModelLibrary.Likelihood.bivariateGaussian "x" "N" }

let hypotheses =

    // [A] N may limited growth via combined N-limitations on (a) photosynthetic and (b) uptake rates
    let limitationModes =
        [ (fun p -> ModelComponents.GrowthLimitation.hollingDiscModel (p |> Pool.getEstimate "a") (p |> Pool.getEstimate "r") (p |> Pool.getEstimate "h")),
           [ code "a",      parameter PositiveOnly   0.001 0.010          // N-uptake efficiency
             code "h",      parameter PositiveOnly   0.001 4.000
             code "r",      parameter PositiveOnly   300.0 500.0 ]      // N-handling time (including uptake and incorporation)
          (fun p -> ModelComponents.GrowthLimitation.linear (p |> Pool.getEstimate "a")), 
           [ code "a",      parameter PositiveOnly   0.001 0.010
             code "r",      parameter PositiveOnly   100.0 500.0 ]        // N-uptake efficiency
          (fun _ -> ModelComponents.GrowthLimitation.none), 
          [ code "r",       parameter PositiveOnly   0.001 1.000 ] ]

    // [B] Loss of plant material may feedback into the soil pool of available nitrogen (instant)
    let feedbackModes =
        [ (fun _ -> ModelComponents.FeedbackToSoil.none), []
          (fun p -> ModelComponents.FeedbackToSoil.withBiomassLoss (p |> Pool.getEstimate "alpha") (p |> Pool.getEstimate "gamma[b]") ),
          [ code "alpha",  parameter PositiveOnly   0.0001 0.0010 ] ]    // N-recycling efficiency

    // [C] A plant may be subject to mechanical constraints on its maximum size
    let geometricModes = 
        [  (fun p -> ModelComponents.GeometricConstraint.none), []
           (fun p -> ModelComponents.GeometricConstraint.chapmanRichards (p |> Pool.getEstimate "k")),
           [ code "k",  parameter PositiveOnly   3000.00 5000.00 ] ]     // Asymptotic biomass (grams)

    List.combine3 geometricModes limitationModes feedbackModes
    |> List.map (fun ((growth,gp),(limit,lp),(feedback,fp)) -> 
        ``base model`` growth limit feedback (List.concat [lp; fp; gp]))


// 3. Load Real Data and Estimate
// ----------------------------

let shrubs = 
    let yuribei = Data.PlantIndividual.loadRingWidths (__SOURCE_DIRECTORY__ + "/../data/yuribei-rw.csv")
    let d15N = Data.PlantIndividual.loadLocalEnvironmentVariable (__SOURCE_DIRECTORY__ + "/../data/yuribei-d15N-imputed.csv")
    yuribei
    |> Seq.map (fun s -> s.Identifier.Value, s)
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
                let start = (c |> TimeSeries.trimStart (startDate - System.TimeSpan.FromDays(366.))).Values |> Array.head |> removeUnit
                printfn "Start cumulative growth = %f" start
                start
            | GrowthSeries.Relative _ -> invalidOp "Not implemented"
        | _ -> invalidOp "Not implemented 2"
    let initialMass = initialRadius |> removeUnit |> ModelComponents.Proxies.toBiomassMM
    let initialNitrogen = plant.Environment.[ShortCode.create "N"].Head |> fst
    [ ShortCode.create "x", initialRadius
      ShortCode.create "N", initialNitrogen 
      ShortCode.create "bs", initialMass ] |> Map.ofList

let everyNth n seq = 
    seq |> Seq.mapi (fun i el -> el, i)              // Add index to element
        |> Seq.filter (fun (el, i) -> i % n = n - 1) // Take every nth element
        |> Seq.map fst                               // Drop index from the result

let workPackages shrubs (hypotheses:ModelSystem list) engine saveDirectory =
    seq {
        for s in shrubs do
            for hi in [ 1.. hypotheses.Length ] do
                    let estimate() =
                        printfn "Starting estimate for shrub %s (H%i)" s.Identifier.Value hi
                        let estimates =
                            [hypotheses.[hi-1]]
                            |> List.toArray
                            |> Array.map(fun h ->
                                [| 1 .. Options.chains |]
                                |> Array.Parallel.mapi(fun i _ ->
                                    [s] |> List.map (fun s ->
                                        let shrub = s |> PlantIndividual.toCumulativeGrowth
                                        let common = shrub |> PlantIndividual.keepCommonYears
                                        let startDate = common.Environment.[ShortCode.create "N"] |> TimeSeries.start
                                        let startConditions = getStartValues startDate shrub
                                        try 
                                            let e = 
                                                engine 
                                                |> Bristlecone.withConditioning (Custom startConditions)         
                                            let result = common |> Bristlecone.PlantIndividual.fit e Options.iterations h
                                            printfn "Result %A" result
                                            Some (s.Identifier, h, result)
                                        with
                                        | e ->
                                            printfn "A chain failed with expection: %A" e
                                            None
                                    )))
                        printfn "Done with estimates for %s. Saving..." s.Identifier.Value

                        // PlantCode, Hypothesis, Iteration, Chain, Parameter, Likelihood
                        let chainsDataFrame =
                            estimates
                            |> Array.map (fun hypothesis -> 
                                hypothesis
                                |> Array.mapi (fun chainNumber chain ->
                                    chain
                                    |> List.choose id
                                    |> List.map (fun (plantCode,model,result) ->
                                        result.Trace
                                        |> List.rev
                                        |> Seq.mapi (fun iterationNumber (likelihood,values) ->
                                            result.Parameters
                                            |> Map.toList
                                            |> List.mapi(fun i (name,p) -> 
                                                plantCode.Value,
                                                hi,
                                                iterationNumber + 1,
                                                chainNumber + 1,
                                                name.Value,
                                                likelihood,
                                                values.[i] )
                                        )
                                        |> everyNth 5 // Thinning by this amount
                                        |> List.concat
                                        |> Seq.toList
                                    )
                                    |> List.concat )
                                |> Array.toList
                                |> List.concat )
                            |> Array.toList
                            |> List.concat

                        // Save as CSV file
                        let buildRowFromObject = fun (a,b,c,d,e,f,g) -> BristleconeResult.Row(a,b,c,d,e,f,g)
                        let buildTableFromObjects = (Seq.map buildRowFromObject) >> Seq.toList >> BristleconeResult
                        let myCsv = chainsDataFrame |> buildTableFromObjects
                        myCsv.Save(sprintf "%sdphil-shrub-output-%s-H%i.csv" saveDirectory s.Identifier.Value hi)
                    yield estimate
            }

(shrubs |> List.map (fun s -> printfn "%s" s.Identifier.Value))

let run() =
    workPackages (shrubs |> List.take 1) (hypotheses |> List.skip 2 |> List.take 1) Options.engine Options.resultsDirectory
    |> Seq.chunkBySize 1
    |> Seq.toArray
    |> Array.map (fun f -> f |> Array.Parallel.map(fun g -> g()))
