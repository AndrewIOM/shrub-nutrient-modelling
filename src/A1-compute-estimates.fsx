#load "1-shrub-nitrogen.fsx"

open Bristlecone
open FSharp.Data
open ``1-shrub-nitrogen``

/////////////////////////////
/// PLOTTING BELOW
/// //////////////////////////

// A. Load in pre-computed results
// ________________________________

let data =
    System.IO.Directory.GetFiles(Options.resultsDirectory, "*.csv")

let bestParameterPool (r: BristleconeResult) =
    let lowestLikelihood = (r.Rows |> Seq.minBy (fun x -> x.Likelihood)).Likelihood
    printfn "Lowest likelihood: %f" lowestLikelihood
    r.Rows
    |> Seq.sortBy(fun x -> x.Likelihood)
    |> Seq.takeWhile(fun x -> x.Likelihood = lowestLikelihood)
    |> Seq.map(fun row -> ShortCode.create row.Parameter, Parameter.create Unconstrained row.Value row.Value)
    |> Map.ofSeq

let hypothesisNumberFromFilename (name:string) =
    name.Split('H').[1].Split('.').[0] |> int

let fit s h =
    let shrub = s |> PlantIndividual.toCumulativeGrowth
    let common = shrub |> PlantIndividual.keepCommonYears
    let startDate = common.Environment.[ShortCode.create "N"] |> TimeSeries.start
    let startConditions = getStartValues startDate shrub
    let e = Bristlecone.mkContinuous |> Bristlecone.withContinuousTime Integration.MathNet.integrate |> Bristlecone.withConditioning (Custom startConditions)
    common |> Bristlecone.PlantIndividual.fit e 1 h

let fitData =
    data
    |> Array.rev
    |> Array.choose(fun fileName ->
        let data = BristleconeResult.Load fileName
        let headRow = data.Rows |> Seq.head
        let hypothesisNumber = hypothesisNumberFromFilename fileName
        let realShrub = shrubs |> List.find (fun s -> s.Identifier.Value = headRow.PlantCode)
        let hypothesis =  { hypotheses.[hypothesisNumber-1] with Parameters = bestParameterPool data}
        printfn "Hypothesis = %A" hypothesis
        try
            let a = fit realShrub hypothesis
            (headRow.PlantCode, hypothesisNumber, a) |> Some
        with
        | _ -> None )
    |> Seq.groupBy(fun (plant,_,_) -> plant)
    |> Seq.collect(fun (plant,data) -> 
        let weights = data |> Seq.map (fun (x,y,z) -> z) |> ModelSelection.Akaike.akaikeWeights
        Seq.zip weights data )
    |> Seq.collect(fun ((e,w),(shrubId, h, estimate)) ->
        estimate.Series 
        |> Map.toArray
        |> Array.collect (fun (s,v) ->
            v.Expected |> Array.mapi (fun i x -> shrubId, h, s.Value, i, x, estimate.Likelihood, w ) ) )
    |> Seq.toArray

type Prediction = CsvProvider<Sample = "PlantCode (string), Hypothesis (int), Variable (string), Time (int), Expected (float), Likelihood (float), AkaikeWeight (float)">

// Save as CSV file
let buildRowFromObject = fun (a,b,c,d,e,f,g) -> Prediction.Row(a,b,c,d,e,f,g)
let buildTableFromObjects = (Seq.map buildRowFromObject) >> Seq.toList >> Prediction
let myCsv = fitData |> buildTableFromObjects
myCsv.Save(sprintf "%sdphil-shrub-predictions.csv" "/Users/andrewmartin/Desktop/")
