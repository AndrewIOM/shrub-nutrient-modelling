#load "../src/bristlecone.fsx"

////////////////////////////////////////////////////
/// Yamal Salix lanata Shrub - Nitrogen Interactions
////////////////////////////////////////////////////

// Shrub ring width modelled with a single
// resource limitation. The models presented are based on
// Tilman (1990 / 1988).

open Bristlecone
open Bristlecone.ModelSystem
open Bristlecone.PlantIndividual
open FSharp.Data

type BristleconeResult = CsvProvider<Sample = "PlantCode (string), Hypothesis (int), Iteration (int), Chain (int), Parameter (string), Likelihood (float), value (float)", CacheRows = true>

// 0. Configure Options
// ----------------------------

module Options =
    let iterations = 10
    let chains = 3
    let testSeriesLength = 30
    let engine =
        Bristlecone.mkContinuous 
        |> Bristlecone.withContinuousTime Integration.MathNet.integrate
        |> Bristlecone.withTunedMCMC [ Optimisation.MonteCarlo.TuneMethod.CovarianceWithScale 0.750, 500, 200000 ]
    let resultsDirectory = "/Volumes/Macbook Time Machine/Testing-Temperature-Other-Hypotheses"



module Constants =

    // Empirically-derived parameters:
    let k5 = 19.98239 // Allometric fit to Yamal shrub BD-length data #1 (in centimetres)
    let k6 = 0.42092 // Allometric fit to Yamal shrub BD-length data #2 (in centimetres)

    // Constants from the literature:
    let a = 2. // the number of child branches added to previous branches (including the tops of the stems) for a shrub
    let p = 0.5 // the length of a child branch as a proportion of its parent branch/stem
    let lmin = 20. //cm. the length at which a stem or branch gets child branches
    let rtip = 0.1 //cm. the radius of the outermost tip of a stem or branch
    let b = 0.0075 // the ratio of the basal radius of a stem or branch and its length
    let salixWoodDensity = 0.5 // g / cm3 (from internet)
    let numberOfStems = 2.2


// 1. Setup model components
// ----------------------------

module ModelComponents =

    let pi = System.Math.PI

    module NiklasAndSpatz_Allometry =

        let nthroot n A =
            let rec f x =
                let m = n - 1.
                let x' = (m * x + A/x**m) / n
                match abs(x' - x) with
                | t when t < abs(x * 1e-9) -> x'
                | _ -> f x'
            f (A / double n)

        /// Gives the basal radius in centimeters of a stem/branch given its length in centimeters. Function from Niklas and Spatz (2004). 
        let basalRadius k5 k6 stemLength =
            100. * (( 0.01 * stemLength + k6) / k5) ** (3. / 2.) / 2.

        /// Inverse equation of basalRadius, rearranged using Wolfram Alpha
        /// http://www.wolframalpha.com/input/?i=solve+r+%3D+100*((0.01*h%2Bk_6)%2Fk_5)%5E(3%2F2)%2F2+for+h
        let stemLength k5 k6 radius =
            2. * ((nthroot 3. 2.) * 5. ** (2./3.) * k5 * radius ** (2./3.) - 50. * k6)

    module Götmark2016_ShrubModel =

        /// Total shrub volume given height and number of stems
        let shrubVolume b a rtip p lmin k5 k6 n h =

            let radius = NiklasAndSpatz_Allometry.basalRadius k5 k6
            let mainStemVolume =
                match radius h with
                | r when r > rtip -> n * pi * h * ((radius h) ** 2. + (radius h) * rtip + rtip ** 2.) / 3.
                | _ -> n * pi * h * rtip ** 2.

            let mutable volume = mainStemVolume
            let mutable k = 0.

            while (p ** k * h > lmin * 2./3.) do
                let volToAdd =
                    match (p ** k * h < lmin) with
                    | true ->
                        match (b * 3. * p * (p ** k * h - 2. * lmin / 3.) > rtip) with
                        | true ->
                            n * a * (a + 1.) ** (float k) * pi * 3. * p * (p ** k * h - 2. * lmin / 3.) * ((radius (3. * p * (p ** k * h - 2. * lmin / 3.))) * (radius (3. * p * (p ** k * h - 2. * lmin / 3.)) * rtip + rtip ** 2.)) / 3.
                        | false ->
                            n * a * (a + 1.) ** (float k) * 3. * p * (p ** k * h - 2. * lmin / 3.) * pi * rtip ** 2.
                    | false ->
                        match (radius (p ** (k + 1.) * h) > rtip) with
                        | true ->
                            n * a * (a + 1.) ** (float k) * pi * p ** (k + 1.) * h * ((radius (p ** (k+1.) * h)) ** 2. + (radius (p ** (k + 1.) * h)) * rtip + rtip ** 2.) / 3.
                        | false ->
                            n * a * (a + 1.) ** (float k) * p ** (k + 1.) * h * pi * rtip ** 2.
                
                volume <- volume + volToAdd
                k <- k + 1.

            k, volume


    module Allometrics =

        open Götmark2016_ShrubModel

        let mass woodDensity volume =
            volume * woodDensity

        let massToVolume woodDensity mass =
            mass / woodDensity

        let shrubBiomass b a rtip p lmin k5 k6 n woodDensity radius =
            radius
            |> NiklasAndSpatz_Allometry.stemLength k5 k6
            |> shrubVolume b a rtip p lmin k5 k6 n |> snd
            |> mass woodDensity

        let shrubRadius b a rtip p lmin k5 k6 n woodDensity mass =
            let findRadius volume =
                let v x = x |> NiklasAndSpatz_Allometry.stemLength k5 k6 |> shrubVolume b a rtip p lmin k5 k6 n |> snd
                let f = (fun x -> (v x) - volume )
                Optimisation.RootFinding.bisect 0 200 f 0.01 1000000.00 1e-8 // Assumption that shrub radius is between 0.01 and 100.0cm.
            mass
            |> massToVolume woodDensity
            |> findRadius

        let shrubHeight k5 k6 radius =
            radius |> NiklasAndSpatz_Allometry.stemLength k5 k6


    module GrowthLimitation =

        /// A rearranged version of a Monod model
        let hollingDiscModel a h (r:float) =
            (a * r) / (1. + (a * h * r))

        /// TEST: An integrated supply and use model
        let supplyAndUse a h b i (r:float) =
            (a * b * r) / (1. + (a * b * h * r) + (i * b * r))

        /// Is saturation ever reached?
        let linear a (r:float) =
            a * r

        /// The resource enforces no limitation on growth, and is negated
        let none (r:float) =
            1.

        /// **Description**
        /// A monotonically increasing function of a resource `r`.
        /// **Parameters**
        ///   * `h` - soil resource concentration required for growth at half the maximum rate
        ///   * `r` - the current resource concentration
        let michaelisMenten h (r:float) =
            r / (h + r)

        /// From Jabot and Pottier 2012
        let monod k (r:float) =
            r / (k + r)


    module AbioticResource =

        /// **Description**
        /// A standard chemostat-type model for replenshment of an abiotic resource.
        /// **Parameters**
        ///   * `d` - a rate constant
        ///   * `s` - resource concentration of the inflow
        ///   * `n` - the current resource concentration
        let chemostat d s (n:float) =
            d * (s - n)


// 2. Create Hypotheses
// ----------------------------

module Proxies =

    /// Radius in millimetres
    let toBiomassMM radiusMM = 
        radiusMM / 10. |> ModelComponents.Allometrics.shrubBiomass Constants.b Constants.a Constants.rtip Constants.p Constants.lmin Constants.k5 Constants.k6 Constants.numberOfStems Constants.salixWoodDensity

    /// Biomass in grams
    let toRadiusMM biomassGrams = 
        let radiusCm = biomassGrams |> ModelComponents.Allometrics.shrubRadius Constants.b Constants.a Constants.rtip Constants.p Constants.lmin Constants.k5 Constants.k6 Constants.numberOfStems Constants.salixWoodDensity
        radiusCm * 10.

    let shrubHeightCm radiusMM =
        radiusMM / 10. |> ModelComponents.Allometrics.shrubHeight Constants.k5 Constants.k6

    /// d15N to N availability. From Craine 2009, as shown in Craine 2015 (Plant and Soil).
    /// Assuming d15N is a linear index of N availability, the minimum supported value of d15N is -3.09, as 0 N availability.
    let d15NtoAvailability d15N =
        (100. * d15N + 309.) / 359.


module GeometricConstraint = 

    /// Linear growth rate in dM/dt form.
    let none _ = 1.

    /// A dM/dt form of the von Bertalanffy monomollecular growth function, where M = mass.
    let vonBertalanffy k m = 
        (k / m)

    /// The dM/dt form of the Chapman-Richards growth function, where M = mass.
    let chapmanRichards k m =
        (1. - (m / k))

module TemperatureDependence =

    // NB The parameter A is our parameter 'r' or 'Q', hence its absence from the below equations.
    let none _ = 1.
    let linear a t = a * t

    /// Activation energy is in KJ rather than J
    let arrhenius a ea t = a * System.Math.E ** (- ((ea * 1000.) / (8.314 * t)))


module ProtectionEffect =

    let none b = b

    let linearSnowProtection protectionEffect shrubHeightCentimetres snowMass b = 
        // Returns protected mass, which is then multiplied by gammaB
        if shrubHeightCentimetres < (snowMass * protectionEffect)
        then 1.
        else 1.


module SoilTemperature =

    // Soil temperature is a function of a flux between summer and winter mean temperatures.
    // The flux is moderated by a snow regulation factor. 
    // Snow depth = snow mass / snow density.
    // BUT, here the density effect is collapsed down into the insulation factor, leaving one term to be estimated.

    /// See: 10.1002/2016JG003725
    let snowMassRegulation insulationFactor snow =
        System.Math.E ** (-insulationFactor * snow)

    /// Fouriers Law for heat flux
    let localHeatFlux conductivity outsideTemp insideTemp =
        - conductivity * (outsideTemp - insideTemp)


module NitrogenReplenishment =

    let linear lambda = lambda

    // There is a heat flux between the air and soil.
    // Here, it is modelled simply as a gradient between summer and winter mean air temperatures.
    // It is modified by the mass of snow on the soil, which provides an insulation layer.
    let temperatureDependentSnowInsulation a ea insulationFactor snowMass summerTemp winterTemp = 
        let conductivityViaSnow = snowMass |> SoilTemperature.snowMassRegulation insulationFactor
        let soilTemperature = SoilTemperature.localHeatFlux conductivityViaSnow summerTemp winterTemp
        TemperatureDependence.arrhenius a ea soilTemperature


module FeedbackToSoil =

    let none b : float = 0. * b
    let withBiomassLoss alpha gammab b : float = alpha * b * gammab


let ``base model`` maxGrowthRate nLimitation nitrogenFeedback nReplenishment protectionEffect temperatureDependency additionalParameters =

    /// Cumulative stem biomass [bs].
    let dbsdt' (b:float) n gammab maxGrowthRate limit protectionEffect tempEffect =
        limit(n) * b * maxGrowthRate(b) * tempEffect - gammab * protectionEffect(b)

    /// Bioavailable soil nitrogen [N]
    let dndt' b n gamman q maxGrowthRate feedback limit nReplenishment tempEffect : float = 
        nReplenishment - q * (maxGrowthRate(b) * b * limit(n) * tempEffect) - gamman * n + feedback(b)

    /// Measurement variable: stem radius [rw].
    let drwdt' bs n gammab maxGrowthRate limit protectionEffect tempEffect =
        let biomassStemChange = dbsdt' bs n gammab maxGrowthRate limit protectionEffect tempEffect
        if biomassStemChange > 0.
            then
                let oldRadius = bs |> Proxies.toRadiusMM
                let newRadius = (bs + biomassStemChange) |> Proxies.toRadiusMM
                newRadius - oldRadius
            else 0.

    /// Bristlecone function for dBs/dt
    let dbsdt p _ bs (e:Environment) =
        dbsdt' bs ((e.[ShortCode.create "N"]) |> Proxies.d15NtoAvailability)
            (p |> Pool.getEstimate "gamma[b]") (maxGrowthRate p) (nLimitation p) (protectionEffect p e) (temperatureDependency p e)

    /// Bristlecone function for dN/dt
    let dndt p _ n (e:Environment) =
        dndt' (e.[ShortCode.create "bs"]) (n |> Proxies.d15NtoAvailability) (p |> Pool.getEstimate "gamma[n]") 
            (p |> Pool.getEstimate "q") (maxGrowthRate p) (nitrogenFeedback p) (nLimitation p) (nReplenishment p e) (temperatureDependency p e)

    /// Bristlecone function for dr/dt
    let drwdt p _ _ (e:Environment) =
        drwdt' (e.[ShortCode.create "bs"]) ((e.[ShortCode.create "N"]) |> Proxies.d15NtoAvailability)
            (p |> Pool.getEstimate "gamma[b]") (maxGrowthRate p) (nLimitation p) (protectionEffect p e) (temperatureDependency p e)

    { Equations  = [ code "x",      drwdt
                     code "bs",     dbsdt
                     code "N",      dndt ] |> Map.ofList
      Parameters = [ code "gamma[n]", parameter PositiveOnly   0.001 2.500   // Loss rate of nitrogen
                     code "gamma[b]", parameter PositiveOnly   0.001 0.250  // Loss rate of biomass
                     code "q",        parameter PositiveOnly   0.001 0.002
                     code "rho",      parameter Unconstrained  -0.50 0.500   // Covariance between growth and nitrogen
                     code "sigma[x]", parameter PositiveOnly   0.100 1.200   // Standard deviation of x (biomass)
                     code "sigma[y]", parameter PositiveOnly   0.250 0.750   // Standard deviation of y (nitrogen)
                    ] |> List.append additionalParameters |> Map.ofList
      Likelihood = ModelLibrary.Likelihood.bivariateGaussian "x" "N" }


let hypotheses =

    // H1. Growth is N-dependent, due to (a) uptake and (b) photosynthetic constraints.
    let limitationModes =
        [ (fun p -> ModelComponents.GrowthLimitation.hollingDiscModel (p |> Pool.getEstimate "a") (p |> Pool.getEstimate "h")),
          [ code "a",      parameter PositiveOnly   0.500 1.500
            code "h",      parameter PositiveOnly   0.100 3.000 ] ]
          //(fun p -> ModelComponents.GrowthLimitation.linear (p |> Pool.getEstimate "a")), 
          //[ code "a",      parameter PositiveOnly   1.000 1.05 ]
          //(fun _ -> ModelComponents.GrowthLimitation.none), [] ]

    // H2. There is a positive feedback of nitrogen to soils, as plant environmental losses increase.
    let feedbackModes =
        [ //(fun _ -> FeedbackToSoil.none), []
          (fun p -> FeedbackToSoil.withBiomassLoss (p |> Pool.getEstimate "alpha") (p |> Pool.getEstimate "gamma[b]") ),
          [ code "alpha",  parameter PositiveOnly   0.0001 0.0002 ] ]

    // H3. The maximum size of a shrub is constrained by geometry.
    let geometricModes = 
        [ (fun p -> GeometricConstraint.none), [] ]
          //(fun p -> GeometricConstraint.chapmanRichards (p |> Pool.getEstimate "k")),
          // [ code "k",  parameter PositiveOnly   3000.00 5000.00 ] ]

    // H4. Increased snow levels protect shrub biomass from storm and other damage.
    let snowProtection =
        [ (fun p e -> ProtectionEffect.none), [] ]
          //(fun p e -> ProtectionEffect.linearSnowProtection (lookup e "bs" |> Proxies.shrubHeightCm) (26. - (lookup e "d18O"))), [] ]

    // H5. Snow insulates soils, which increases the efficiency of N-mineralising microbes.
    let nitrogenReplenishment =
        [ (fun p _ -> NitrogenReplenishment.linear (p |> Pool.getEstimate "lambda")),
            [ code "lambda", parameter PositiveOnly   1.000 10.000 ]
          (fun p e -> NitrogenReplenishment.temperatureDependentSnowInsulation (p |> Pool.getEstimate "lambda") (p |> Pool.getEstimate "soilEa") (p |> Pool.getEstimate "insulation") (26. - (lookup e "d18O")) (lookup e "Tsummer") (lookup e "Twinter")),
            [ code "lambda",        parameter PositiveOnly   0.001 0.100
              code "soilEa",        parameter PositiveOnly   0.001 0.100
              code "insulation",    parameter PositiveOnly   0.001 0.100 ] ]

    // H7. Net photosynthetic rate is temperature-dependent
    let temperature =
        [ (fun _ -> TemperatureDependence.none), []
          //(fun p e -> TemperatureDependence.linear (p |> Pool.getEstimate "A") (lookup e "Tsummer")),
          // [ code "A",  parameter PositiveOnly 0.001 0.100 ]
          (fun p e -> TemperatureDependence.arrhenius (p |> Pool.getEstimate "A") (p |> Pool.getEstimate "Ea") (lookup e "Tsummer")),
           [ code "A",              parameter PositiveOnly 1.00 2.00
             code "Ea",             parameter PositiveOnly 1.00 2.00  ] ]

    // Create all combinations of H1-H5 
    List.combine6 geometricModes limitationModes feedbackModes nitrogenReplenishment snowProtection temperature
    |> List.map (fun ((growth,gp),(limit,lp),(feedback,fp), (replace,rp), (snow,sp), (temp,tp)) -> 
        ``base model`` growth limit feedback replace snow temp (List.concat [lp; fp; gp; rp; sp; tp]))


// 3. Load Real Data and Estimate
// ----------------------------

// a) Load regional (shared) environment
type Climate = CsvProvider<"../samples/data/yamal-temperature-means.csv">
let regionalClimate = Climate.Load("../samples/data/yamal-temperature-means.csv")
let summerTemperature = 
    let convert (x:decimal<_>) = decimal x
    regionalClimate.Rows
    |> Seq.map(fun r -> r.Year, r.``JJA Mean Temperaure`` |> convert |> float) 
    |> TimeSeries.createVarying

let winterTemperature = 
    let convert (x:decimal<_>) = decimal x
    regionalClimate.Rows
    |> Seq.map(fun r -> r.Year, r.``DJF Mean Temperature`` |> convert |> float) 
    |> TimeSeries.createVarying

// b) Load shrub individual data + environment
let shrubs = 
    let yuribei = DataAccess.Shrub.loadRingWidths (__SOURCE_DIRECTORY__ + "/data/yuribei-rw.csv")
    let d15N = DataAccess.Shrub.loadLocalEnvironmentVariable (__SOURCE_DIRECTORY__ + "/data/yuribei-d15N-imputed.csv")
    let d18O = DataAccess.Shrub.loadLocalEnvironmentVariable (__SOURCE_DIRECTORY__ + "/data/yuribei-d18O-imputed.csv")
    yuribei
    |> Seq.map (fun s -> s.Identifier.Value, s)
    |> Seq.keyMatch d15N
    |> Seq.map (fun (_,plant,d15N) -> PlantIndividual.zipEnv (code "N") plant d15N)
    |> Seq.map (fun s -> s.Identifier.Value, s)
    |> Seq.keyMatch d18O
    |> Seq.map (fun (_,plant,d18O) -> PlantIndividual.zipEnv (code "d18O") plant d18O)
    |> Seq.map (fun (plant) -> PlantIndividual.zipEnv (code "Tsummer") summerTemperature plant)
    |> Seq.map (fun (plant) -> PlantIndividual.zipEnv (code "Twinter") winterTemperature plant)
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
    let initialMass = initialRadius |> removeUnit |> Proxies.toBiomassMM
    let initialNitrogen = plant.Environment.[ShortCode.create "N"].Head |> fst
    let initialWinterT = plant.Environment.[ShortCode.create "Twinter"].Head |> fst
    let initialSummerT = plant.Environment.[ShortCode.create "Tsummer"].Head |> fst
    let initialSnow = plant.Environment.[ShortCode.create "d18O"].Head |> fst
    [ ShortCode.create "x", initialRadius
      ShortCode.create "N", initialNitrogen 
      ShortCode.create "d18O", initialSnow
      ShortCode.create "Twinter", initialWinterT
      ShortCode.create "Tsummer", initialSummerT
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
                                |> Array.Parallel.map(fun _ ->
                                    [s] |> List.map (fun s ->
                                        let shrub = s |> PlantIndividual.toCumulativeGrowth
                                        let common = shrub |> PlantIndividual.keepCommonYears
                                        let startDate = common.Environment.[ShortCode.create "N"] |> TimeSeries.start
                                        let startConditions = getStartValues startDate shrub
                                        try 
                                            let e = engine |> Bristlecone.withConditioning (Custom startConditions)
                                            let result = common |> Bristlecone.PlantIndividual.fit e Options.iterations h
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

workPackages shrubs hypotheses Options.engine Options.resultsDirectory
|> Seq.chunkBySize 3
|> Seq.toArray
|> Array.map (fun f -> f |> Array.Parallel.map(fun g -> g()))




/////////////////////////////
/// PLOTTING BELOW
/// //////////////////////////

// A. Load in pre-computed results
// ________________________________

let data =
    System.IO.Directory.GetFiles(Options.resultsDirectory, "*.csv")
    |> Array.map BristleconeResult.Load

let bestParameterPool (r: BristleconeResult) =
    let lowestLikelihood = (r.Rows |> Seq.minBy (fun x -> x.Likelihood)).Likelihood
    printfn "Lowest likelihood: %f" lowestLikelihood
    r.Rows
    |> Seq.sortBy(fun x -> x.Likelihood)
    |> Seq.takeWhile(fun x -> x.Likelihood = lowestLikelihood)
    |> Seq.map(fun row -> ShortCode.create row.Parameter, Parameter.create Unconstrained row.Value row.Value)
    |> Map.ofSeq


let bestFits =
  
    // TODO remove duplication of this function
    let fit s h =
        let shrub = s |> PlantIndividual.toCumulativeGrowth
        let common = shrub |> PlantIndividual.keepCommonYears
        let startDate = common.Environment.[ShortCode.create "N"] |> TimeSeries.start
        let startConditions = getStartValues startDate shrub
        let e = Bristlecone.mkContinuous |> Bristlecone.withContinuousTime Integration.MathNet.integrate |> Bristlecone.withGradientDescent |> Bristlecone.withConditioning (Custom startConditions)
        common |> Bristlecone.PlantIndividual.fit e 1 h

    data
    |> Array.groupBy(fun i -> (i.Rows |> Seq.head).PlantCode, (i.Rows |> Seq.head).Hypothesis)
    |> Array.map(fun ((shrub,hi),d) -> shrub, hi, { hypotheses.[hi-1] with Parameters = bestParameterPool (d |> Seq.head) } )
    |> Array.map(fun (shrubId,hi,h) -> 
        let realShrub = shrubs |> List.find (fun s -> s.Identifier.Value = shrubId)
        shrubId, hi, fit realShrub h )

let weights =
    bestFits
    |> Seq.map(fun (s,i,r) -> r)
    |> ModelSelection.Akaike.akaikeWeights
    |> Seq.toArray

// ShrubId, Hypothesis, Variable, Time, Expected Value, Weight
let fitData =
    bestFits
    |> Array.zip weights
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
myCsv.Save(sprintf "%sdphil-shrub-predictions-snow.csv" "/Users/andrewmartin/Desktop/")
