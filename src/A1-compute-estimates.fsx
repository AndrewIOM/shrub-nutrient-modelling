// #load "1-shrub-nitrogen.fsx"
#load "3-shrub-nitrogen-longterm.fsx"

open Bristlecone
open FSharp.Data
// open ``1-shrub-nitrogen``
open ``3-shrub-nitrogen-longterm``

let loadFrom = "/Users/andrewmartin/Desktop/Bristlecone Results/LongTermThreeYearly/"

// A. Load in pre-computed results
// ________________________________

type BristleconeResult = CsvProvider<"/Users/andrewmartin/Projects/GitHub-Projects/bristlecone/src/Bristlecone/templates/saved-data.csv">

let data =
    System.IO.Directory.GetFiles(loadFrom, "*.csv")

let bestParameterPool step (constraints:Constraint[]) (r: BristleconeResult) =
    let lowestLikelihood = (r.Rows |> Seq.minBy (fun x -> x.NegativeLogLikelihood)).NegativeLogLikelihood
    printfn "Lowest likelihood: %f" lowestLikelihood
    r.Rows
    |> Seq.sortBy(fun x -> x.NegativeLogLikelihood)
    |> Seq.takeWhile(fun x -> x.NegativeLogLikelihood = lowestLikelihood)
    |> Seq.take constraints.Length
    |> Seq.mapi(fun i row -> (ShortCode.create row.ParameterCode, Parameter.create (constraints.[i]) (row.ParameterValue) (row.ParameterValue + step)))
    |> Map.ofSeq

let hypothesisNumberFromFilename (name:string) =
    name.Split('H').[1].Split('-').[0] |> int

let fit engineAdditions endCondition s h =
    let shrub = s |> PlantIndividual.toCumulativeGrowth
    let common = shrub |> PlantIndividual.keepCommonYears
    let startDate = (common.Environment.[ShortCode.create "N"]).StartDate |> snd
    let startConditions = getStartValues startDate shrub
    let e = Bristlecone.mkContinuous |> Bristlecone.withOutput Options.logger |> Bristlecone.withContinuousTime Integration.MathNet.integrate |> Bristlecone.withTunedMCMC [] |> Bristlecone.withConditioning (Custom startConditions) |> engineAdditions
    common |> Bristlecone.PlantIndividual.fit e endCondition h

let fitData =
    data
    |> Array.map(fun fileName ->
        let data = BristleconeResult.Load fileName
        let headRow = data.Rows |> Seq.head
        let hypothesisNumber = headRow.ModelId
        let realShrub = shrubs |> List.find (fun s -> s.Identifier.Value = headRow.Subject)
        let constraints = hypotheses.[hypothesisNumber-1].Parameters |> Seq.map(fun x -> x.Value |> Parameter.detatchConstraint |> snd) |> Seq.toArray
        let hypothesis =  { hypotheses.[hypothesisNumber-1] with Parameters = bestParameterPool 0. constraints data }
        printfn "Hypothesis (%s) %i" headRow.Subject hypothesisNumber
        let a = fit id (Optimisation.EndConditions.afterIteration 0) realShrub hypothesis
        (realShrub, hypothesisNumber, a))
    |> Array.groupBy(fun (plant,h,_) -> (plant, h))
    |> Array.map(fun ((plant,h),r) -> r |> Seq.minBy(fun (_,_,e) -> e.Likelihood))
    |> Array.groupBy(fun (plant,_,_) -> plant)
    |> Array.collect(fun (plant,data) -> 
        let weights = data |> Array.map (fun (x,y,z) -> z) |> ModelSelection.Akaike.akaikeWeights |> Seq.toArray
        printfn "Weights are %A for data %A" (weights |> Array.map snd) data
        Array.zip weights data )
    |> Array.collect(fun ((e,w),(plant, h, estimate)) ->
        // Recompute the time index to assign to observed and predicted series
        let shrub = plant |> PlantIndividual.toCumulativeGrowth
        let common = shrub |> PlantIndividual.keepCommonYears
        let startDate = (common.Environment.[ShortCode.create "N"]).StartDate |> snd
        let timeIndex = common.Environment.[ShortCode.create "N"] |> TimeSeries.toObservations |> Seq.map snd |> Seq.toArray
        let startConditions = getStartValues startDate shrub
        estimate.Series 
        |> Map.toArray
        |> Array.collect (fun (s,v) ->
            // Convert values from Cumulative radius to RGR (mm)
            let expected,observed = 
                if s.Value = "x" 
                then
                    let e = 
                        v.Expected 
                        |> Array.scan(fun (mass,_) newMass -> newMass, ((newMass - mass) / mass)) (startConditions.[code "x"],0.)
                        |> Array.tail // Remove t0
                        |> Array.map snd
                    let o = 
                        v.Observed 
                        |> Array.scan(fun (mass,_) newMass -> newMass, ((newMass - mass) / mass)) (startConditions.[code "x"],0.)
                        |> Array.tail // Remove t0
                        |> Array.map snd
                    (e, o)
                else (v.Expected, v.Observed)
            expected
            |> Array.zip observed
            |> Array.mapi (fun i (o,e) -> plant.Identifier.Value, h, s.Value, timeIndex.[i], e, o, estimate.Likelihood, w ) ) )

// TODO Make time a date rather than an integer

type Prediction = CsvProvider<Sample = "PlantCode (string), Hypothesis (int), Variable (string), Time (date), Expected (float), Observed (float), Likelihood (float), AkaikeWeight (float)">

// Save as CSV file
let buildRowFromObject = fun (a,b,c,d,e,f,g,h) -> Prediction.Row(a,b,c,d,e,f,g,h)
let buildTableFromObjects = (Seq.map buildRowFromObject) >> Seq.toList >> Prediction
let myCsv = fitData |> buildTableFromObjects
myCsv.Save(sprintf "%sdphil-shrub-predictions-paper3-3yearly.csv" "/Users/andrewmartin/Desktop/")


// Profile Likelihood Method
// _______________________________

// Given the MLE (i.e. from Simulated Annealing), 
// run an MCMC algorithm that samples around the MLE
// Collect any results where the difference in likelihoods is less than 2

type Interval = {
    Lower: float
    Upper: float
}

and ConfidenceInterval = {
    ``95%``: Interval
    ``68%``: Interval
    Estimate: float
}

/// The difference in likelihood at 68% confidence
let lowerBound = 0.49447324 // qchisq(0.68,1)/2
/// The difference in likelihood at 95% confidence
let upperBound = 1.92072941 // qchisq(0.95,1)/2

// let randomWalk theta writeOut n domain f : (float * float[]) list =
//     let random = MathNet.Numerics.Random.MersenneTwister(true)
//     let initialCovariance = Bristlecone.Optimisation.MonteCarlo.TuningMode.covarianceFromBounds 10000 domain random
//     let tune = (Optimisation.MonteCarlo.TuneMethod.Covariance 0.250, 500, Optimisation.EndConditions.afterIteration 20000)
//     Bristlecone.Optimisation.MonteCarlo.randomWalk' initialCovariance 1. theta [ tune ] random writeOut n domain f |> fst

open Bristlecone.Logging

// Given a candidate distribution + machine, run base SA algorithm
let tunedSearch scale settings annealEnd machine (jump:System.Random->float->float->unit->float) cool writeOut domain f =

    // 1. Initial conditions
    let random = MathNet.Numerics.Random.MersenneTwister(true)
    let draw' = jump random
    let theta1 = Bristlecone.Optimisation.MonteCarlo.initialise domain random
    let l1 = f theta1

    // 2. Tune individual parameter scales (for many univariate distributions)
    let initialScale = 
        [| 1 .. theta1.Length |] |> Array.map (fun _ -> scale)

    let kMax = 60000//20000 * domain.Length
    let rec tune (p:(float*float[])[]) k results (l1, theta1) =
        let chance = (float k) / (float kMax)
        let parameterToChange = random.Next(0, (p |> Array.length) - 1)
        let scalesToChange = p |> Array.mapi (fun i x -> (x, random.NextDouble() < chance || i = parameterToChange))
        let propose theta =
            Array.zip3 theta scalesToChange domain
            |> Array.map(fun (x,((ti,n),shouldChange),(_,_,con)) -> 
                if shouldChange
                then Bristlecone.Optimisation.MonteCarlo.constrainJump x (draw' ti 1. ()) 1. con
                else x )
        let result = Bristlecone.Optimisation.MonteCarlo.SimulatedAnnealing.tryMove propose (machine 1.) random f (l1, theta1)
        let tuneN = 50
        let newScaleInfo = 
            scalesToChange 
            |> Array.zip (result |> snd)
            |> Array.map(fun (v, ((ti,previous),changed)) ->
                if changed then (ti, (previous |> Array.append [|v|]))  // Append new parameter values to previous ones
                else (ti, previous) )
            |> Array.map(fun (ti,previous) ->
                if previous |> Array.length = tuneN
                then
                    let changes = previous |> Array.pairwise |> Array.where(fun (a,b) -> a <> b) |> Array.length
                    match (float changes) / (float tuneN) with
                    | ar when ar < 0.35 -> (ti * 0.80, Array.empty)
                    | ar when ar > 0.50 -> (ti * 1.20, Array.empty)
                    | _ -> (ti, Array.empty)
                else (ti, previous) )

        if k % 1000 = 0 then
            writeOut <| GeneralEvent (sprintf "Tuning is at %A (k=%i/%i)" (newScaleInfo |> Array.map fst) k kMax)

        if k < kMax
        then
            writeOut <| OptimisationEvent { Iteration = k; Likelihood = result |> fst; Theta = result |> snd } 
            tune newScaleInfo (k + 1) (results |> List.append [result]) result
        else newScaleInfo |> Array.map fst

    let result = tune (initialScale |> Array.map(fun t -> (t, Array.empty))) 1 [] (l1, theta1)
    writeOut <| GeneralEvent (sprintf "Tuned = %A" result)
    [(l1, theta1)]

let classic scale settings writeOut n domain f =
    let gaussian rnd scale t = Bristlecone.Statistics.Distributions.Normal.draw rnd 0. (scale * sqrt (1.))
    tunedSearch scale settings n Optimisation.MonteCarlo.SimulatedAnnealing.Machines.boltzmann gaussian (Optimisation.MonteCarlo.SimulatedAnnealing.CoolingSchemes.exponential 0.05) writeOut domain f

let optimise =
    classic 0.0001
        { Optimisation.MonteCarlo.SimulatedAnnealing.AnnealSettings<float>.Default with 
            BoilingAcceptanceRate = 0.10
            HeatRamp = (fun t -> t + sqrt t); TemperatureCeiling = Some 0.99
            HeatStepLength = Optimisation.EndConditions.afterIteration 1
            AnnealStepLength = (fun x -> Optimisation.MonteCarlo.SimulatedAnnealing.EndConditions.improvementCount 5000 250 x || Optimisation.EndConditions.afterIteration 50000 x) }


type OptimOutput() =
    let mutable (events:Bristlecone.Logging.ModelFitState list) = []
    with
        member __.SaveEvent(e) =
            match e with
            | Bristlecone.Logging.LogEvent.OptimisationEvent o -> events <- (events |> List.append [o])
            | _ -> ()
            Options.logger e

        member __.GetAll() = events


/// `n` is the number of samples to generate within the likelihood profile.
let profileLikelihood (fileName:string) =

    // 1. Load in saved Bristlecone data using filename
    let data = BristleconeResult.Load fileName
    let headRow = data.Rows |> Seq.head
    let mle = (data.Rows |> Seq.minBy(fun r -> r.NegativeLogLikelihood)).NegativeLogLikelihood
    printfn "MLE is %f" mle
    let hypothesisNumber = headRow.ModelId
    let realShrub = shrubs |> List.find (fun s -> s.Identifier.Value = headRow.Subject)
    let constraints = hypotheses.[hypothesisNumber-1].Parameters |> Seq.map(fun x -> x.Value |> Parameter.detatchConstraint |> snd) |> Seq.toArray
    let hypothesis =  { hypotheses.[hypothesisNumber-1] with Parameters = bestParameterPool 0. constraints data}
    printfn "Hypothesis (%s) %i" headRow.Subject hypothesisNumber    
    let theta = hypothesis.Parameters |> Seq.map(fun k -> k.Value |> Parameter.detatchConstraint |> fst |> Parameter.bounds |> fst) |> Seq.toArray
    printfn "Parameter pool is %A" hypothesis.Parameters
    printfn "Theta is %A" theta

    // 2. Generate a trace of at least n samples that deviate in L less than 2.0
    let results = OptimOutput()
    let rec fit' currentTrace =
        let a = fit (Bristlecone.withCustomOptimisation optimise >> Bristlecone.withOutput results.SaveEvent) (Optimisation.EndConditions.afterIteration 20000) realShrub hypothesis
        let validTrace = a.Trace |> List.filter(fun (l,_) -> (l - mle) < 2.00 && (l - mle) > 0.00) |> List.distinct
        printfn "Valid trace was %i elements long..." (validTrace |> List.length)
        // if currentTrace |> List.append validTrace |> List.length < n
        // then currentTrace |> List.append validTrace |> fit'
        // else 
        currentTrace |> List.append validTrace
    fit' [] |> ignore

    let trace = 
        results.GetAll()
        |> List.map(fun s -> (s.Likelihood, s.Theta |> Seq.toArray))
    printfn "Actual trace was %A" trace

    // 3. Calculate min and max at the specified limit for each parameter
    let interval limit =
        trace
        |> List.filter(fun (l,p) -> (l - mle) < limit)
        |> List.fold(fun best (l,p) -> best |> Array.zip p |> Array.map(fun (v,(mn,mx)) -> ((min v mn), (max v mx)))) 
            (Array.init hypothesis.Parameters.Count (fun _ -> (System.Double.MaxValue, System.Double.MinValue)))

    let lowerInterval = interval lowerBound
    let upperInterval = interval upperBound

    (headRow.Subject, headRow.ModelId, (theta
    |> Seq.zip3 lowerInterval upperInterval
    |> Seq.mapi(fun i ((l1,l2),(u1,u2),e) -> 
        (hypothesis.Parameters |> Seq.toArray |> Seq.item i).Key, {
            Estimate = e
            ``68%`` = { Lower = l1; Upper = l2 }
            ``95%`` = { Lower = u1; Upper = u2 }})
    |> Map.ofSeq))

// Save as CSV file
type ConfidenceData = CsvProvider<Sample = "PlantCode (string), Hypothesis (int), Parameter (string), Value (float), 68% Lower (float), 68% Upper (float), 95% Lower (float), 95% Upper (float)">

let bestAnswers =
    data
    |> Array.map(fun fileName ->
        let data = BristleconeResult.Load fileName
        let headRow = data.Rows |> Seq.head
        let hypothesisNumber = headRow.ModelId
        let mle = (data.Rows |> Seq.minBy (fun x -> x.NegativeLogLikelihood)).NegativeLogLikelihood
        (headRow.Subject, hypothesisNumber, fileName, mle))
    |> Array.groupBy(fun (plant,h,_,_) -> (plant, h))
    |> Array.chunkBySize 8
    |> Array.collect(fun g ->
        g |> Array.Parallel.map(fun (_,r) -> 
            r 
            |> Seq.minBy(fun (_,_,_,mle) -> mle)
            |> fun (_,_,d,_) -> profileLikelihood d
            |> fun x -> printfn "Bounds were: %A" x; x ))
    |> Array.collect(fun (s,h,ci) -> ci |> Seq.map(fun p -> ConfidenceData.Row(s, h, p.Key.Value, p.Value.Estimate, p.Value.``68%``.Lower, p.Value.``68%``.Upper, p.Value.``95%``.Lower, p.Value.``95%``.Upper)) |> Seq.toArray )

let d = new ConfidenceData(bestAnswers)
d.Save(sprintf "%sdphil-shrub-predictions-paper3-tuned-intervals.csv" "/Users/andrewmartin/Desktop/")



// One-step-ahead predictions
// __________________________

// Predict the next point from the current observation.
// Calculate the Root Mean Square Error (RMSE) between expected and next observed

let rootMeanSquareDeviation pairs =
    List.averageBy (fun (a,b) -> (a - b)**2.0) pairs |> System.Math.Sqrt

/// Given a time series, finds the model prediction from each point to the next point.
let oneStepAhead fit t0 observations =
    observations
    |> List.append t0
    |> List.map(fun obs-> fit 1. obs)

// Output:
// - One-step-ahead predicted time-series
// - Overall RMSD

let oneStepPredictions =
    data
    |> Array.map(fun fileName ->
        let data = BristleconeResult.Load fileName
        let headRow = data.Rows |> Seq.head
        let hypothesisNumber = headRow.ModelId
        let mle = (data.Rows |> Seq.minBy (fun x -> x.NegativeLogLikelihood)).NegativeLogLikelihood
        (headRow.Subject, hypothesisNumber, fileName, mle))
    |> Array.groupBy(fun (plant,h,_,_) -> (plant, h))
    |> Array.map(fun (_,r) -> 
        r 
        |> Seq.minBy(fun (_,_,_,mle) -> mle)
        |> fun (_,_,d,_) ->
            Bristlecone.generateFixedSeries  )