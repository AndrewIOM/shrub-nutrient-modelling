module Units

#r "/Users/andrewmartin/Documents/GitHub Projects/bristlecone/src/Bristlecone/bin/Debug/net5.0/Bristlecone.dll"
#r "/Users/andrewmartin/Documents/GitHub Projects/bristlecone/src/Bristlecone.Dendro/bin/Debug/net5.0/Bristlecone.Dendro.dll"

open Bristlecone
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