module Units

// #r "nuget: Bristlecone.Dendro, 3.0.0-beta1"
#r "/Users/andrewmartin/Documents/GitHub Projects/bristlecone/src/Bristlecone.Dendro/bin/Debug/net10.0/Bristlecone.dll"
#r "/Users/andrewmartin/Documents/GitHub Projects/bristlecone/src/Bristlecone.Dendro/bin/Debug/net10.0/Bristlecone.Dendro.dll"
// #r "nuget: DiffSharp-cpu, v=1.0.7"
#r "nuget: FSharp.Data"
#r "nuget: MathNet.Numerics, v=5.0.0"

// Use a local fork of DiffSharp
#I "../../../lib"
#r "DiffSharp.Core.dll"
#r "DiffSharp.Backends.Reference.dll"

open DiffSharp

// dsharp.config(backend=Backend.Torch, device=Device.CPU)

open Bristlecone.Language
open FSharp.Data.UnitSystems.SI

[<Measure>] type g
[<Measure>] type kg = UnitNames.kilogram
[<Measure>] type cm

[<Measure>] type δ15N
[<Measure>] type nutrient

[<Measure>] type km
[<Measure>] type area = km^2
[<Measure>] type year = Bristlecone.Time.year
[<Measure>] type mm = Bristlecone.Dendro.Units.millimetre
[<Measure>] type m = UnitNames.metre

[<Measure>] type J = UnitNames.joule
[<Measure>] type kJ
[<Measure>] type mol = UnitNames.mole
[<Measure>] type K = UnitNames.kelvin
[<Measure>] type DegC


// Conversion factors
let gramsPerKg = Constant 1000.0<g/kg>
let cmPerM = Constant 100.0<cm/m>
let mmPerCm = Constant 10.<mm/cm>
let joulesPerKilojoule = Constant 1000.<J/kJ>

// Explicit conversions
let toMetres (x: ModelExpression<cm>) = x / cmPerM
let toGrams (x: ModelExpression<kg>) = x * gramsPerKg
let toCentimetres (x: ModelExpression<m>) = x * cmPerM
let mmToCm (x:ModelExpression<mm>) = x / mmPerCm
let cmToMm (x:ModelExpression<cm>) = x * mmPerCm
let toJoule (x:ModelExpression<kJ * 'u>) : ModelExpression<J * 'u> = x * joulesPerKilojoule