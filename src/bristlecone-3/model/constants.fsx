module Constants

#r "/Users/andrewmartin/Documents/GitHub Projects/bristlecone/src/Bristlecone/bin/Debug/net5.0/Bristlecone.dll"
#load "units.fsx"

open Bristlecone
open Bristlecone.Language
open Units

module Allometrics =

    // Empirically-derived parameters:
    let k5 = Constant 19.98239 // Allometric fit to Yamal shrub BD-length data #1 (in centimetres)
    let k6 = Constant 0.42092<cm> // Allometric fit to Yamal shrub BD-length data #2 (in centimetres)

    // Constants from the literature:
    let a = Constant 2. // the number of child branches added to previous branches (including the tops of the stems) for a shrub
    let p = Constant 0.5 // the length of a child branch as a proportion of its parent branch/stem
    let lmin = Constant 20.<cm> // the length at which a stem or branch gets child branches
    let rtip = Constant 0.1<cm> // the radius of the outermost tip of a stem or branch
    let b = Constant 0.0075 // the ratio of the basal radius of a stem or branch and its length
    let salixWoodDensity = Constant 0.5<g/cm^3>
    let numberOfStems = Constant 2.2
