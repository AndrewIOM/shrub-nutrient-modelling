#r "../packages/NETStandard.Library.NETFramework/build/net461/lib/netstandard.dll"
#r "../packages/FSharp.Data/lib/net45/FSharp.Data.dll"
#r "../packages/MathNet.Numerics/lib/net40/MathNet.Numerics.dll"
#r "../packages/Bristlecone/lib/netstandard2.0/bristlecone.dll"
#r "../packages/Bristlecone.Dendro/lib/netstandard2.0/bristlecone.Dendro.dll"
#load "components/components.fsx"
#load "components/temperature.fsx"
#load "components/snow.fsx"

////////////////////////////////////////////////////
/// Yamal Salix lanata Shrub - Nitrogen Interactions
////////////////////////////////////////////////////

// Shrub ring width modelled with a single
// resource limitation.

open Bristlecone
open Bristlecone.ModelSystem
open Bristlecone.PlantIndividual
open FSharp.Data

type BristleconeResult = CsvProvider<Sample = "PlantCode (string), Hypothesis (int), Iteration (int), Chain (int), Parameter (string), Likelihood (float), value (float)", CacheRows = true>


// 1. Configure Options
// ----------------------------

module Options =
    let resultsDirectory = "/Users/andrewmartin/Desktop/Snow-temperature models/"
    let iterations = 10
    let chains = 3
    let engine =
        Bristlecone.mkContinuous 
        |> Bristlecone.withContinuousTime Integration.MathNet.integrate
        |> Bristlecone.withOutput (Bristlecone.Logging.Console.logger())
        |> Bristlecone.withTunedMCMC [ Optimisation.MonteCarlo.TuneMethod.CovarianceWithScale 0.750, 500, 100000 ]


// 2. Create Hypotheses
// ----------------------------

module BaseEquations =

    /// Cumulative stem biomass [dBs/dt]
    let biomass b n r gammab geom f protectionEffect tempEffect : float =
        b * r * (f n) * geom(b) * tempEffect - gammab * protectionEffect(b)

    let soilNitrogen n b gamman y geom f feedback tempEffect : float =
        y - (geom(b) * b * (f n) * tempEffect) - gamman * n + feedback(b)

    let soilNitrogenNoUptake n b gamman y feedback : float =
        y - gamman * n + feedback(b)


let ``base model`` maxGrowthRate nLimitation nitrogenFeedback nReplenishment protectionEffect temperatureDependency additionalParameters =

    let randomLog environment parameters =
        if System.Random().Next(0,50000) = 1 then
            printfn "RANDOM LOG: Environment = %A; Parameters = %A" environment parameters

    /// Cumulative stem biomass [b].
    let dbsdt' (b:float) n gammab r maxGrowthRate limit protectionEffect tempEffect =
        match limit with
        | Some l -> BaseEquations.biomass b n r gammab maxGrowthRate l protectionEffect tempEffect
        | None -> BaseEquations.biomass b n r gammab maxGrowthRate (fun _ -> 1.) protectionEffect tempEffect

    /// Bioavailable soil nitrogen [N]
    let dndt' bs n gamman maxGrowthRate feedback limit nReplenishment tempEffect = 
        match limit with
        | Some l -> BaseEquations.soilNitrogen n bs gamman nReplenishment maxGrowthRate l feedback tempEffect
        | None -> BaseEquations.soilNitrogenNoUptake n bs gamman nReplenishment feedback

    /// Measurement variable: stem radius [rw].
    let drwdt' bs n r gammab maxGrowthRate limit protectionEffect tempEffect =
        let biomassStemChange = dbsdt' bs n gammab r maxGrowthRate limit protectionEffect tempEffect
        if biomassStemChange > 0.
            then
                let oldRadius = bs |> ModelComponents.Proxies.toRadiusMM
                let newRadius = (bs + biomassStemChange) |> ModelComponents.Proxies.toRadiusMM
                newRadius - oldRadius
            else 0.

    /// Bristlecone function for dBs/dt
    let dbsdt p _ bs (e:Environment) =
        randomLog e p
        dbsdt' bs ((lookup e "N") |> ModelComponents.Proxies.d15NtoAvailability)
            (p |> Pool.getEstimate "gamma[b]") (p |> Pool.getEstimate "r") (maxGrowthRate p) (nLimitation p) (protectionEffect p e) (temperatureDependency p e)

    /// Bristlecone function for dN/dt
    let dndt p _ n (e:Environment) =
        dndt' (lookup e "bs") (n |> ModelComponents.Proxies.d15NtoAvailability) 
            (p |> Pool.getEstimate "gamma[n]") (maxGrowthRate p) (nitrogenFeedback p) (nLimitation p) (nReplenishment p e) (temperatureDependency p e)

    /// Bristlecone function for dr/dt
    let drwdt p _ _ (e:Environment) =
        drwdt' (lookup e "bs") ((lookup e "N") |> ModelComponents.Proxies.d15NtoAvailability) 
            (p |> Pool.getEstimate "r") (p |> Pool.getEstimate "gamma[b]") (maxGrowthRate p) (nLimitation p) (protectionEffect p e) (temperatureDependency p e)

    { Equations  = [ code "x",         drwdt
                     code "bs",        dbsdt
                     code "N",         dndt ] |> Map.ofList
      Parameters = [ code "gamma[n]",  parameter PositiveOnly   0.001 0.100   // Loss rate of nitrogen
                     code "gamma[b]",  parameter PositiveOnly   0.001 0.010   // Loss rate of biomass
                     code "r",         parameter PositiveOnly   0.010 0.100   // Photosynthetic efficiency (either N-limited or N-unlimited)
                     code "rho",       parameter Unconstrained  -0.50 0.500   // Covariance between growth and nitrogen
                     code "sigma[x]",  parameter PositiveOnly   0.100 1.200   // Standard deviation of x (biomass)
                     code "sigma[y]",  parameter PositiveOnly   0.250 0.750   // Standard deviation of y (nitrogen)
                    ] |> List.append additionalParameters |> Map.ofList
      Likelihood = ModelLibrary.Likelihood.bivariateGaussian "x" "N" }


let hypotheses =

    // [A] Increased snow levels protect shrub biomass from storm and other damage.
    let snowProtection =
        [ (fun p e -> ModelComponents.SnowProtection.none), [] ]
          //(fun p e -> ModelComponents.SnowProtection.linearSnowProtection  (lookup e "bs" |> ModelComponents.Proxies.shrubHeightCm) (26. - (lookup e "d18O"))), [] ]

    // [B] Snow insulates soils, which increases the efficiency of N-mineralising microbes.
    let nitrogenReplenishment =
        [ //(fun p _ -> ModelComponents.Temperature.NitrogenReplenishment.linear (p |> Pool.getEstimate "lambda")),
          //  [ code "lambda",        parameter PositiveOnly   1.000 10.000 ]
          (fun p e -> ModelComponents.Temperature.NitrogenReplenishment.temperatureDependent (p |> Pool.getEstimate "lambda") (p |> Pool.getEstimate "soilEa") (p |> Pool.getEstimate "insulation") (lookup e "Tsummer") (lookup e "Twinter")),
            [ code "lambda",        parameter PositiveOnly   0.001 10.00
              code "soilEa",        parameter PositiveOnly   0.001 10.00
              code "insulation",    parameter PositiveOnly   0.001 1.000 ]
          (fun p e -> ModelComponents.Temperature.NitrogenReplenishment.temperatureDependentSnowInsulation (p |> Pool.getEstimate "lambda") (p |> Pool.getEstimate "soilEa") (p |> Pool.getEstimate "insulation") (26. - (lookup e "d18O")) (lookup e "Tsummer") (lookup e "Twinter")),
            [ code "lambda",        parameter PositiveOnly   0.001 10.00
              code "soilEa",        parameter PositiveOnly   0.001 10.00
              code "insulation",    parameter PositiveOnly   0.001 1.000 ] ]

    // [C] Net photosynthetic rate is temperature-dependent
    let temperature =
        [ //(fun _ -> ModelComponents.Temperature.none), []
          (fun p e -> ModelComponents.Temperature.arrhenius (p |> Pool.getEstimate "A") (p |> Pool.getEstimate "Ea") (lookup e "Tsummer")),
           [ code "A",              parameter PositiveOnly 1.00 2.00
             code "Ea",             parameter PositiveOnly 1.00 2.00  ] ]

    // H1. Growth is N-dependent, due to (a) uptake and (b) photosynthetic constraints.
    let limitationModes =
        [ //(fun p -> ModelComponents.GrowthLimitation.hollingDiscModel (p |> Pool.getEstimate "a") (p |> Pool.getEstimate "h")),
          //[ code "a",      parameter PositiveOnly   0.500 1.500
          //  code "h",      parameter PositiveOnly   0.100 3.000 ] ]
          (fun p -> ModelComponents.GrowthLimitation.linear (p |> Pool.getEstimate "a")), 
          [ code "a",      parameter PositiveOnly   1.000 1.05 ] ]
          //(fun _ -> ModelComponents.GrowthLimitation.none), [] ]

    // H2. There is a positive feedback of nitrogen to soils, as plant environmental losses increase.
    let feedbackModes =
        [ (fun _ -> ModelComponents.FeedbackToSoil.none), [] ]
          //(fun p -> ModelComponents.FeedbackToSoil.withBiomassLoss (p |> Pool.getEstimate "alpha") (p |> Pool.getEstimate "gamma[b]") ),
          //[ code "alpha",  parameter PositiveOnly   0.0001 0.0002 ] ]

    // H3. The maximum size of a shrub is constrained by geometry.
    let geometricModes = 
        [ (fun p -> ModelComponents.GeometricConstraint.none), [] ]
          //(fun p -> GeometricConstraint.chapmanRichards (p |> Pool.getEstimate "k")),
          // [ code "k",  parameter PositiveOnly   3000.00 5000.00 ] ]

    // Create all combinations of H1-H5 
    List.combine6 geometricModes limitationModes feedbackModes nitrogenReplenishment snowProtection temperature
    |> List.map (fun ((growth,gp),(limit,lp),(feedback,fp), (replace,rp), (snow,sp), (temp,tp)) -> 
        ``base model`` growth limit feedback replace snow temp (List.concat [lp; fp; gp; rp; sp; tp]))


// 3. Load Real Data and Estimate
// ----------------------------

type Climate = CsvProvider<"../data/yamal-temperature-mean-quarters.csv">
let regionalClimate = Climate.Load("../data/yamal-temperature-mean-quarters.csv")
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
    let yuribei = Data.PlantIndividual.loadRingWidths (__SOURCE_DIRECTORY__ + "/../data/yuribei-rw.csv")
    let d15N = Data.PlantIndividual.loadLocalEnvironmentVariable (__SOURCE_DIRECTORY__ + "/../data/yuribei-d15N-imputed.csv")
    let d18O = Data.PlantIndividual.loadLocalEnvironmentVariable (__SOURCE_DIRECTORY__ + "/../data/yuribei-d18O-imputed.csv")
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
    let initialMass = initialRadius |> removeUnit |> ModelComponents.Proxies.toBiomassMM
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
    workPackages shrubs hypotheses Options.engine Options.resultsDirectory
    |> Seq.chunkBySize 2
    |> Seq.toArray
    |> Array.map (fun f -> f |> Array.Parallel.map(fun g -> g()))
