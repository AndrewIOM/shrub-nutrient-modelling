#load "1-shrub-nitrogen.fsx"

open Bristlecone
open FSharp.Data
open ``1-shrub-nitrogen``

// A. Load in pre-computed results
// ________________________________

type BristleconeResult = CsvProvider<"/Users/andrewmartin/Projects/GitHub-Projects/bristlecone/src/Bristlecone/templates/saved-data.csv">

let data =
    System.IO.Directory.GetFiles(Options.resultsDirectory, "*.csv")

let bestParameterPool (r: BristleconeResult) =
    let lowestLikelihood = (r.Rows |> Seq.minBy (fun x -> x.NegativeLogLikelihood)).NegativeLogLikelihood
    printfn "Lowest likelihood: %f" lowestLikelihood
    r.Rows
    |> Seq.sortBy(fun x -> x.NegativeLogLikelihood)
    |> Seq.takeWhile(fun x -> x.NegativeLogLikelihood = lowestLikelihood)
    |> Seq.map(fun row -> (ShortCode.create row.ParameterCode, Parameter.create Unconstrained row.ParameterValue row.ParameterValue))
    |> Map.ofSeq

let hypothesisNumberFromFilename (name:string) =
    name.Split('H').[1].Split('-').[0] |> int

let fit s h =
    let shrub = s |> PlantIndividual.toCumulativeGrowth
    let common = shrub |> PlantIndividual.keepCommonYears
    let startDate = common.Environment.[ShortCode.create "N"] |> TimeSeries.start
    let startConditions = getStartValues startDate shrub
    let e = Bristlecone.mkContinuous |> Bristlecone.withContinuousTime Integration.MathNet.integrate |> Bristlecone.withTunedMCMC [] |> Bristlecone.withConditioning (Custom startConditions)
    common |> Bristlecone.PlantIndividual.fit e (Optimisation.EndConditions.afterIteration 0) h

let fitData =
    data
    |> Array.map(fun fileName ->
        let data = BristleconeResult.Load fileName
        let headRow = data.Rows |> Seq.head
        let hypothesisNumber = headRow.ModelId
        let realShrub = shrubs |> List.find (fun s -> s.Identifier.Value = headRow.Subject)
        let hypothesis =  { hypotheses.[hypothesisNumber-1] with Parameters = bestParameterPool data}
        printfn "Hypothesis (%s) %i" headRow.Subject hypothesisNumber
        let a = fit realShrub hypothesis
        (headRow.Subject, hypothesisNumber, a))
    |> Array.groupBy(fun (plant,h,_) -> (plant, h))
    |> Array.map(fun ((plant,h),r) -> r |> Seq.minBy(fun (_,_,e) -> e.Likelihood))
    |> Array.groupBy(fun (plant,_,_) -> plant)
    |> Array.collect(fun (plant,data) -> 
        let weights = data |> Array.map (fun (x,y,z) -> z) |> ModelSelection.Akaike.akaikeWeights |> Seq.toArray
        printfn "Weights are %A for data %A" (weights |> Array.map snd) data
        Array.zip weights data )
    |> Array.collect(fun ((e,w),(shrubId, h, estimate)) ->
        estimate.Series 
        |> Map.toArray
        |> Array.collect (fun (s,v) ->
            v.Expected |> Array.mapi (fun i x -> shrubId, h, s.Value, i, x, estimate.Likelihood, w ) ) )

type Prediction = CsvProvider<Sample = "PlantCode (string), Hypothesis (int), Variable (string), Time (int), Expected (float), Likelihood (float), AkaikeWeight (float)">

// Save as CSV file
let buildRowFromObject = fun (a,b,c,d,e,f,g) -> Prediction.Row(a,b,c,d,e,f,g)
let buildTableFromObjects = (Seq.map buildRowFromObject) >> Seq.toList >> Prediction
let myCsv = fitData |> buildTableFromObjects
myCsv.Save(sprintf "%sdphil-shrub-predictions.csv" "/Users/andrewmartin/Desktop/")
