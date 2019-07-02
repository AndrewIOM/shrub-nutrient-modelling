#r "../packages/NETStandard.Library.NETFramework/build/net461/lib/netstandard.dll"
#load "../packages/Bristlecone/bristlecone.fsx"

open FSharp.Data

type SnowData = CsvProvider<"/Users/andrewmartin/Documents/DPhil Zoology/Data/Snow Data (By Peter L)/snow-with-dates.csv">

let rawData = SnowData.Load "/Users/andrewmartin/Documents/DPhil Zoology/Data/Snow Data (By Peter L)/snow-with-dates.csv"


// Snow depth (cm) mean
// 1. ..in summer
let summer =
    rawData.Rows
    |> Seq.groupBy(fun r -> r.Year)
    |> Seq.map(fun (y,v) ->
        let values = v |> Seq.where(fun r -> r.Month >= 4 && r.Month <= 9)
        if values |> Seq.isEmpty then y, None
        else y, values |> Seq.averageBy(fun r -> r.``Chang.1987.sd..cm.``) |> Some)
    |> Seq.toList

// 2. .. in previous winter
let previousWinter =
    rawData.Rows
    |> Seq.groupBy(fun r -> r.Year)
    |> Seq.pairwise
    |> Seq.map(fun ((y1,v1),(y2,v2)) ->
        let winterStart = v1 |> Seq.where(fun r -> r.Month > 12)
        let winterEnd = v2 |> Seq.where(fun r -> r.Month > 3 && r.Month < 5)
        let values = Seq.concat [ winterStart; winterEnd ]
        if values |> Seq.isEmpty then y2, None
        else y2, (values |> Seq.averageBy(fun r -> r.``Aschbacher.1989.sd..cm.``)) |> Some)
    |> Seq.toList

for p in previousWinter do
    printfn "%i,%f" (p |> fst) ((p |> snd).Value |> float)

type MonthPrecip = CsvProvider<"/Users/andrewmartin/Desktop/cru-yurbei-precip-monthly.csv">

let m = MonthPrecip.Load "/Users/andrewmartin/Desktop/cru-yurbei-precip-monthly.csv"

let previousWinterMonth =
    m.Rows
    |> Seq.sortBy(fun r -> r.Year)
    |> Seq.pairwise
    |> Seq.map(fun (y1,y2) ->
        printfn "Y1 = %A Y2 = %A" y1.Year y2.Year
        let winter = [ y1.Sep; y1.Oct; y1.Nov; y1.Dec ]
        y2.Year, (winter |> Seq.average))
    |> Seq.toList

for p in previousWinterMonth do
    printfn "%i,%f" (fst p) (snd p |> float)