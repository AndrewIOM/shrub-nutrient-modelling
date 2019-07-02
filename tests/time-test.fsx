#r "../packages/NETStandard.Library.NETFramework/build/net461/lib/netstandard.dll"
#load "../packages/Bristlecone/bristlecone.fsx"
#load "components/components.fsx"
#r "../packages/Bristlecone.Dendro/lib/netstandard2.0/bristlecone.Dendro.dll"

/////////////////////////////////////
/// Bristlecone Model Tests
/////////////////////////////////////

(* To ensure Bristlecone is effective at solving long-term ecological
   problems, we run a model-fitting-model-selection (MFMS) procedure
   for two ecological models. First, simple plant growth models are fit
   to plant growth data. Second, growth-resource models are fit.  *)

open Bristlecone

// A. Time series are ordered

let data =
    let random = System.Random()
    seq {
        for i in [1. .. 21. ] do
            yield (random.NextDouble(), (System.DateTime.Create 01 01 (1939 + (int i))))
    } |> Seq.toList

let ts1 = TimeSeries.fromObservations data
let ts2 = TimeSeries.fromSeq (System.DateTime.Create 01 01 1940) (Years 1) (data |> Seq.map fst)

// Test both methods return same series
ts1 = ts2
ts1.Resolution
ts2.Resolution

// Trim start
ts1 |> TimeSeries.trimStart (System.DateTime.Create 31 12 1946)

// Trim end
ts1 |> TimeSeries.trimEnd (System.DateTime.Create 31 12 1946)

// Bound
ts1 |> TimeSeries.bound (System.DateTime.Create 01 01 1943) (System.DateTime.Create 31 12 1946)

let ts3 = TimeSeries.fromSeq (System.DateTime.Create 01 01 1940) (3 |> Years) (data |> Seq.map fst)
ts3 |> TimeSeries.findExact (System.DateTime.Create 01 01 1943)
ts3.Resolution

// Bouding - appears to work
ts1 |> TimeSeries.bound (System.DateTime.Create 01 01 1943) (System.DateTime.Create 31 12 1946)

let hello = TimeSeries.fromSeq (System.DateTime.Create 01 01 1946) (24 |> Months) (data |> Seq.map fst)
hello.Resolution
// Time modes of the time series
// A. Discrete - the readings occur at the exact moments given in time
// B. Period - the readings occur within the following time interval
    // e.g. the start (t1) occurs before t2 (when the first TimeSpan has elapsed)
    // This means that the last point does not have a span in which it occurs.

// Unify two time series of differing resolutions
// A. Daily and Yearly

// B. Yearly and Yearly

let averageBin obs = obs |> Seq.averageBy(fun (v,_) -> v)
let upscaled = ts1 |> TimeSeries.generalise (3 |> Years) averageBin

upscaled |> TimeSeries.dates |> Seq.map(fun x -> x.Year) |> Seq.toList
upscaled.Resolution


// Testing out time indexes

let baseline = System.DateTime.Create 01 01 1900
let index = TimeIndex.TimeIndex(baseline, Days 1, ts1)
ts1 |> TimeSeries.toObservations |> Seq.map(fun (v,d) -> d.ToShortDateString()) |> Seq.toList
index.Values |> Seq.toList
//index.Item 0.

// Compare two time series

let x = TimeSeries.commonTimeline [ ts1; ts2 ]


// Mismatched resolution data
