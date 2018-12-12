#load "../src/bristlecone.fsx"

open Bristlecone
open Bristlecone.ModelSystem
open Bristlecone.PlantIndividual
open FSharp.Data

type BristleconeResult = CsvProvider<Sample = "PlantCode (string), Hypothesis (int), Iteration (int), Chain (int), Parameter (string), Likelihood (float), value (float)">

// 0. Configure Options
// ----------------------------

module Options =

    let resolution = Annual
    let iterations = 10
    let burn = 10
    let chains = 3
    let testSeriesLength = 40


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

    module GrowthLimitation =

        /// **Description**
        /// A monotonically increasing function of a resource `r`.
        /// **Parameters**
        ///   * `h` - soil resource concentration required for growth at half the maximum rate
        ///   * `r` - the current resource concentration
        let michaelisMenten h (r:float) =
            (r / (h + r))

        /// From Jabot and Pottier 2012
        let monod k (r:float) =
            r / (k + r)

        /// Is saturation ever reached?
        let linear a (r:float) =
            a * r

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

    /// d15N to N availability. From Craine 2009, as shown in Craine 2015 (Plant and Soil).
    /// Assuming d15N is a linear index of N availability, the minimum supported value of d15N is -3.09, as 0 N availability.
    let d15NtoAvailability d15N =
        (100. * d15N + 309.) / 359.


module Allocation =

    type PlantMass =
    | TotalMass of float
    | Compartments of PlantBiomass

    and PlantBiomass = {
        Leaf: float
        Root: float
        Stem: float
        LeafToStemRatio: float
        RootToStemRatio: float }

    /// All biomass is assumed to be photosynthetic, with no allocation based on tissue function. 
    let none b = TotalMass b

    /// Roots, leaves, and stem are in equal 1:1:1 ratio. 
    let equal bs = Compartments { Stem = bs; Root = bs; Leaf = bs; LeafToStemRatio = 1.; RootToStemRatio = 1. }

    /// Floating allocation. 
    let floating al ar bs = Compartments { Stem = bs; Leaf = (bs * al); Root = (bs * ar); LeafToStemRatio = al; RootToStemRatio = ar }


module GrowthRate = 

    /// Linear growth rate in dM/dt form.
    let linear r _ = r

    /// A dM/dt form of the von Bertalanffy monomollecular growth function, where M = mass.
    let vonBertalanffy r k m = 
        r * (k / m)

    /// The dM/dt form of the Chapman-Richards growth function, where M = mass.
    let chapmanRichards r k m =
        r * m * (1. - (m / k))


module FeedbackToSoil =

    let none b : float = 0. * b
    let withBiomassLoss alpha gammab b : float = alpha * b * gammab


// module BaseEquationsOld =
    
//     open Allocation

//     /// Uptake rate (modelled as 'effective nutrient availability')
//     let effectiveNutrient n rootMass leafMass : float = 
//         n * rootMass / leafMass

//     /// **Description**
//     /// Cumulative stem biomass [db/dt].
//     /// If bs>0 and bl and br are set at zero, allocation is disabled and bs represents total plant biomass.
//     /// If bl and br are >0, allocation enabled and equation returns stem biomass production.
//     /// 
//     /// **Parameters**
//     ///   * `b` - Current `PlantBiomass`
//     ///   * `n` - Current nitrogen availability
//     ///   * `r` - Maximum rate of growth
//     ///   * `limit` - Limiting effect of an environmental factor on growth
//     ///   * `m` - Loss rate of biomass due to disturbance
//     ///   * `resp` - Respiration rate per unit biomass
//     ///
//     /// **Output Type**
//     ///   * `float` - Change in [stem] biomass during this time point
//     let plantBiomass (b:PlantBiomass) n r limit m resp = 
//         // printfn "Stem fraction: %f" ((b.Stem / (b.Leaf + b.Root + b.Stem)))
//         let netPhotosynthesis = b.Leaf * (r b.Leaf) * (limit (effectiveNutrient n b.Root b.Leaf)) - resp * (b.Stem + b.Leaf + b.Root)
//         (b.Stem / (b.Leaf + b.Root + b.Stem)) * (netPhotosynthesis - m * (b.Stem + b.Leaf + b.Root))
//         // Stem fraction * [ (leaf capacity - respiration) - environmental loss ]
//         //b.Stem * (f - resp * (b.Stem + b.Leaf + b.Root)) - m * (b.Stem + b.Leaf + b.Root)

//     let soilResource (b:PlantBiomass) (n:float) lambda gamman q r limit resp feedback : float =
//         let netPhotosynthesis = b.Leaf * (r b.Leaf) * (limit (effectiveNutrient n b.Root b.Leaf)) - resp * (b.Stem + b.Leaf + b.Root)
//         lambda - q * netPhotosynthesis - gamman * n + feedback (b.Stem + b.Leaf + b.Root)


module BaseEquations =

    open Allocation

    module ProportionalAllocation =
    
        /// In this model, nutrient aquisition depends on root mass, while photosynthesis depends on leaf mass.
        /// Stem mass serves no function.

        /// Effective nutrient availability, depending on the plant's allocation to nutrient-seeking tissues
        let effectiveNutrient n rootMass leafMass : float = 
            n * rootMass / leafMass

        /// Cumulative stem biomass [dBs/dt]
        let stemBiomass b n al ar gammab rr r f =
            (1. / (1. + al + ar)) * ((b.Leaf * r(b.Leaf) * f(n) - rr * (b.Leaf + b.Stem + b.Root)) - gammab * (b.Leaf + b.Stem + b.Root))

        /// Soil nitrogen availability
        let soilNitrogen n b q rr gamman y r f =
            y - q * (b.Leaf * r(b.Leaf) * f(n) - rr * (b.Leaf + b.Stem + b.Root)) - gamman * n

    module NoAllocation =

        /// Cumulative stem biomass [dBs/dt]
        let biomass b n rr gammab r f : float =
            b * r(b) * f(n) - rr * b - gammab * b //NB Can this solve with both respiration and other losses?
        
        let soilNitrogen n b q rr gamman y r f : float =
            y - q * (b * r(b) * f(n) - rr * b) - gamman * n

let ``base model`` maxGrowthRate nLimitation allocationMode nitrogenFeedback additionalParameters =

    /// Cumulative stem biomass [bs].
    let dbsdt' (bs:float) n gammab rr maxGrowthRate allocationMode limit = 
        match bs |> allocationMode with
        | Allocation.TotalMass b -> BaseEquations.NoAllocation.biomass b n rr gammab maxGrowthRate limit
        | Allocation.Compartments b -> BaseEquations.ProportionalAllocation.stemBiomass b n b.LeafToStemRatio b.RootToStemRatio gammab rr maxGrowthRate limit

    /// Bioavailable soil nitrogen [N]
    let dndt' bs n lambda gamman q rr maxGrowthRate allocationMode feedback limit = 
        match bs |> allocationMode with
        | Allocation.TotalMass b -> BaseEquations.NoAllocation.soilNitrogen n b q rr gamman lambda maxGrowthRate limit
        | Allocation.Compartments b -> BaseEquations.ProportionalAllocation.soilNitrogen n b q rr gamman lambda maxGrowthRate limit

    /// Measurement variable: stem radius [rw].
    let drwdt' bs n m r maxGrowthRate allocationMode limit = 
        let biomassStemChange = dbsdt' bs n m r maxGrowthRate allocationMode limit
        if biomassStemChange > 0.
            then
                let oldRadius = bs |> Proxies.toRadiusMM
                let newRadius = (bs + biomassStemChange) |> Proxies.toRadiusMM
                newRadius - oldRadius
            else 0.

    /// Bristlecone function for dBs/dt
    let dbsdt p _ bs (e:Environment) =
        dbsdt' bs ((e.[ShortCode.create "N"]) |> Proxies.d15NtoAvailability)
            (p |> Pool.getEstimate "gammab") (p |> Pool.getEstimate "resp") (maxGrowthRate p) (allocationMode p) (nLimitation p)

    /// Bristlecone function for dN/dt
    let dndt p _ n (e:Environment) =
        dndt' (e.[ShortCode.create "x"]) (n |> Proxies.d15NtoAvailability) (p |> Pool.getEstimate "lambda") (p |> Pool.getEstimate "gamman") 
            (p |> Pool.getEstimate "q") (p |> Pool.getEstimate "resp") (maxGrowthRate p) (allocationMode p) (nitrogenFeedback p) (nLimitation p)

    /// Bristlecone function for dr/dt
    let drwdt p _ _ (e:Environment) =
        drwdt' (e.[ShortCode.create "bs"]) ((e.[ShortCode.create "N"]) |> Proxies.d15NtoAvailability)
            (p |> Pool.getEstimate "gammab") (p |> Pool.getEstimate "resp") (maxGrowthRate p) (allocationMode p) (nLimitation p)

    { Equations  = [ ShortCode.create "x",      drwdt
                     ShortCode.create "bs",     dbsdt
                     ShortCode.create "N",      dndt ] |> Map.ofList
      Parameters = [ // for nitrogen replenishment
                     ShortCode.create "lambda", Parameter.create PositiveOnly   0.400 0.800   // Rate of nitrogen replenishment
                     ShortCode.create "gamman", Parameter.create PositiveOnly   0.300 1.000   // Loss rate of nitrogen
                     // for shrub allocation and physiology
                     ShortCode.create "resp",   Parameter.create PositiveOnly   0.001 0.100   // Respiration cost per unit biomass
                     ShortCode.create "q",      Parameter.create PositiveOnly   0.001 0.100   // Nutrient requirement per unit biomass
                     ShortCode.create "gammab", Parameter.create PositiveOnly   0.001 0.750   // Loss rate of biomass
                     // for likelihood function
                     ShortCode.create "rho",    Parameter.create Unconstrained  -0.50 0.50     // Covariance between growth and nitrogen
                     ShortCode.create "sigmax", Parameter.create PositiveOnly   0.200 1.000    // Standard deviation of x (biomass)
                     ShortCode.create "sigmay", Parameter.create PositiveOnly   0.200 0.600    // Standard deviation of y (nitrogen)
                    ] |> List.append additionalParameters |> Map.ofList
      Likelihood = ModelLibrary.Likelihood.bivariateGaussian "x" "N" }

let hypotheses =

    let growthRateModes = 
        [ (fun p -> GrowthRate.linear (p |> Pool.getEstimate "r")), 
            [ ShortCode.create "r",  Parameter.create PositiveOnly   0.500 1.000 ] ]
        //   (fun p -> GrowthRate.vonBertalanffy (p |> Pool.getEstimate "r") (p |> Pool.getEstimate "k")),
        //     [ ShortCode.create "r",  Parameter.create PositiveOnly   0.001 1.000
        //       ShortCode.create "k",  Parameter.create PositiveOnly   0.001 1.000 ]
        //   (fun p -> GrowthRate.chapmanRichards (p |> Pool.getEstimate "r") (p |> Pool.getEstimate "k")),
        //     [ ShortCode.create "r",  Parameter.create PositiveOnly   0.001 1.000
        //       ShortCode.create "k",  Parameter.create PositiveOnly   0.001 1.000 ] ]

    let limitationModes =
        [ //(fun p -> ModelComponents.GrowthLimitation.linear (p |> Pool.getEstimate "a")), 
          //  [ ShortCode.create "a",      Parameter.create PositiveOnly   0.200 1.600 ]
          (fun p -> ModelComponents.GrowthLimitation.michaelisMenten (p |> Pool.getEstimate "h")),
            [ ShortCode.create "h",      Parameter.create PositiveOnly   0.001 0.750 ] ]   // Soil resource availability for growth at half of maximum rate 

    let feedbackModes =
        [ (fun p -> FeedbackToSoil.none), [] ]
          //(fun p -> FeedbackToSoil.withBiomassLoss (p |> Pool.getEstimate "alpha") (p |> Pool.getEstimate "gammab") ),
          // [ ShortCode.create "alpha",  Parameter.create PositiveOnly   0.001 0.005 ] ]
          //Additional hypothesis: lagged input to soils?

    let allocationModes =
        [ (fun _ -> Allocation.none), []
        //   (fun _ -> Allocation.equal), [] ]
          (fun p -> Allocation.floating (p |> Pool.getEstimate "al") (p |> Pool.getEstimate "ar")),
          [ ShortCode.create "al",     Parameter.create PositiveOnly   0.900 1.000 
            ShortCode.create "ar",     Parameter.create PositiveOnly   0.900 1.000 ] ]

    List.combine4 growthRateModes limitationModes feedbackModes allocationModes
    |> List.map (fun ((growth,gp),(limit,lp),(feedback,fp),(allocation,ap)) -> 
        ``base model`` growth limit allocation feedback (List.concat [lp; ap; fp; gp]))


// 3. Load Real Data and Estimate
// ----------------------------

let shrubs = 
    let yuribei = DataAccess.Shrub.loadRingWidths (__SOURCE_DIRECTORY__ + "/data/yuribei-rw.csv")
    let d15N = DataAccess.Shrub.loadLocalEnvironmentVariable (__SOURCE_DIRECTORY__ + "/data/yuribei-d15N-imputed.csv")
    yuribei
    |> Seq.map (fun s -> s.Identifier.Value, s)
    |> Seq.keyMatch d15N
    |> Seq.map (fun (_,plant,d15N) -> PlantIndividual.zipEnv (ShortCode.create "N") plant d15N)
    |> Seq.toList

// TODO implement getting start values properly

let getStartValues startDate (plant:PlantIndividual) =
    let initialRadius =
        match plant.Growth with
        | PlantIndividual.PlantGrowth.RingWidth s -> 
            match s with
            | GrowthSeries.Absolute c -> c.Head |> fst |> removeUnit
            | GrowthSeries.Cumulative c -> 6.3 //c.Head |> fst |> removeUnit
            | _ -> invalidOp "Not implemented"
        | _ -> invalidOp "Not implemented 2"
    let initialMass = initialRadius |> removeUnit |> Proxies.toBiomassMM
    let initialNitrogen = plant.Environment.[ShortCode.create "N"].Head |> fst
    [ ShortCode.create "x", initialRadius
      ShortCode.create "N", initialNitrogen 
      ShortCode.create "bs", initialMass ] |> Map.ofList


// Generate results for all models, chains, and plants
let estimates =
    hypotheses
    |> List.toArray
    |> Array.Parallel.map(fun h ->
        [| 1 .. Options.chains |]
        |> Array.Parallel.map(fun _ ->
            [shrubs.[0]] |> List.map (fun s ->
                let shrub = s |> PlantIndividual.toCumulativeGrowth
                let common = shrub |> PlantIndividual.keepCommonYears
                let startDate = common.Environment.[ShortCode.create "N"] |> TimeSeries.start
                let startConditions = getStartValues startDate shrub
                s.Identifier, h, common |> Bristlecone.PlantIndividual.fit (Bristlecone.mkContinuous |> Bristlecone.withConditioning (Custom startConditions)) Options.iterations Options.burn h )))

printfn "Done with estimates. Saving..."


// 4. Save Results to R data file
// ----------------------------
let burnin = 0

let everyNth n seq = 
    seq |> Seq.mapi (fun i el -> el, i)              // Add index to element
        |> Seq.filter (fun (el, i) -> i % n = n - 1) // Take every nth element
        |> Seq.map fst                               // Drop index from the result

// PlantCode, Hypothesis, Iteration, Chain, Parameter, Likelihood
let chainsDataFrame =
    estimates
    |> Array.mapi (fun hypothesisNumber hypothesis -> 
        hypothesis
        |> Array.mapi (fun chainNumber chain ->
            chain
            |> List.map (fun (plantCode,model,result) ->
                result.Trace
                |> List.rev
                |> List.skip burnin // Skip burn in period
                |> Seq.mapi (fun iterationNumber (likelihood,values) ->
                    result.Parameters
                    |> Map.toList
                    |> List.mapi(fun i (name,p) -> 
                        plantCode.Value,
                        hypothesisNumber,
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
myCsv.Save("dphil-shrub-output.csv")