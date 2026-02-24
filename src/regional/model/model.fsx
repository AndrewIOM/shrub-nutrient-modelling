(**
# Long-Term N-Shrub relations in Yamal, Russia

The following code defines models to compete to determine
the mechanisms controlling shrub growth in Arctic tundra.

The hypotheses are tested on a per-shrub basis. For each shrub,
there are four components that may vary in the model, which
when multiplied together results in 12 possible combinations
of mechanisms for each shrub.

First, we must load Bristlecone and the seperate script that
defines the allometric equations that we will use later.
*)

#load "units.fsx"
#load "components.fsx"

open Bristlecone.Language // Open the language for writing Bristlecone models
open Bristlecone.Time
open Units

(**
### Defining a 'base model'

Then, we define a `base model` system, into which we can nest
the components that vary (which represent different hypotheses).

Our base model consists of the following parts:

* **An empirical transform.** As this model uses long-term ecological data, 
   there is a function that translates from an ecological proxy (stable nitrogen isotope)
   into the desired variable - soil nitrogen availability.
* **Two ordinary differential equations (ODEs).** These represent the soil N availability
   and plant biomass through time.
* **A *measurement* variable.** Bristlecone 'measurements' are effectively
   derived time-series that are calulated from the integrated ODE output. Here,
   we use a measurement variable to transform biomass into wood ring production,
   where rings are only produced where plant biomass has increased during the year.

The code for the parts of the base model is shown below.
*)

// Empirical transform from δ15N to N availability.
let ``δ15N -> N availability`` n =
    (Constant 100. * n + Constant 309.) / Constant 359.

// Parameters (rates):
let r = parameter "r" Positive 1e3<g/1> 1e4<g/1> // gram per unit N (N is dimensionless)
let lambda = parameter "λ" (Bounded 1e-12</year> 50.</year>) 0.100</year> 2.000</year>
let gammaN = parameter "γ[N]" (Bounded 1e-12</year> 10.</year>) 0.100</year> 10.00</year>
let gammaB = parameter "γ[b]" Positive 0.001</year> 1.000</year>

// States:
let N = state     "N"
let B = state<g> "bs"
let SR = measure<mm> "x"

// Environmental forcings
let TMax = environment<K> "T[max]"

// f(N) represents the combined effect of nitrogen uptake and conversion into plant biomass
type NLimitation =
  { UptakeRatePerBiomass : ModelExpression<1> -> ModelExpression<1/g/year>
    UptakeMultiplier: ModelExpression<1> }

/// ODE1. Cumulative stem biomass
/// NB. A parameter constraint is added using a `Conditional` expression.
/// NB. The parameter r is multiplied by 1000 when applying an N-limitation purely to aid model-fitting.
let ``db/dt`` (geomLimit:ModelExpression<g> -> ModelExpression<1>) (nLimitation:NLimitation) (envLimit:ModelExpression<1>) : ModelExpression<g/year> =
    This<g> * P r
        * nLimitation.UptakeRatePerBiomass (``δ15N -> N availability`` (State N))
        * geomLimit This<g>
        * envLimit
    - P gammaB * This

/// ODE2. A dimensionless state variable representing "soil nitrogen availability"
let ``dN/dt`` (geomLimit:ModelExpression<g> -> ModelExpression<1>) nLimitation (envLimit:ModelExpression<1>) (feedback:ModelExpression<g> -> ModelExpression</year>) : ModelExpression</year> =
    P lambda
    - P gammaN * ``δ15N -> N availability`` This<1>
    + feedback (State B)
    - nLimitation.UptakeMultiplier * ((geomLimit (State B)) * (State B) * (nLimitation.UptakeRatePerBiomass (``δ15N -> N availability`` This<1>)) * (envLimit))

/// Measurement variable: stem radius
let stemRadius : ModelExpression<mm> =
    let oldCumulativeMass = StateAt (-1<``time index``>, B)
    let newCumulativeMass = StateAt (0<``time index``>, B)
    Conditional (newCumulativeMass - oldCumulativeMass .> Constant 0.<g>)
        (newCumulativeMass |> ShrubModel.Allometry.Proxies.toRadiusMM)
        This

/// Define a function to set the initial biomass, based on the initial
/// radius.
let initialMass = ShrubModel.Allometry.Proxies.toBiomassMM (Measure SR)

(**
Once we have defined the components, we can scaffold them into a model system.
We can plug in the nestable hypotheses (defined further below) by defining
the base model as a function that takes parameters representing the
alternative hypotheses.
*)

let ρ = parameter "ρ" NoConstraints -0.500 0.500
let σSR = parameter "σ_SR" (Bounded 1e-12<mm> 5.<mm>) 0.001<mm> 0.100<mm>
let σN = parameter "σ_N" (Bounded 1e-12 3.) 0.001 0.100

let NLL = Bristlecone.ModelLibrary.NegLogLikelihood.BivariateNormal (Require.measure SR) (Require.state N) σSR σN ρ

let ``base model``
    (geometricConstraint: ModelExpression<g> -> ModelExpression<1>)
    (plantSoilFeedback: ModelExpression<g> -> ModelExpression</year>)
    (envLimit: ModelExpression<1>)
    (nLimitation: NLimitation) =
    Model.empty
    |> Model.addRateEquation B (``db/dt`` geometricConstraint nLimitation envLimit)
    |> Model.addRateEquation N (``dN/dt`` geometricConstraint nLimitation envLimit plantSoilFeedback)
    |> Model.addMeasure SR stemRadius
    |> Model.initialiseHiddenStateWith B initialMass
    |> Model.estimateParameter lambda
    |> Model.estimateParameter gammaN
    |> Model.estimateParameter gammaB
    |> Model.useLikelihoodFunction NLL

(**
### Defining the competing hypotheses

Here, we define 24 alternative hypotheses by defining four interchangeable 
components:

 - Asymptotic plant size (2 types);
 - Plant-soil feedback presence / absence (2 types);
 - Nitrogen limitation form (3 types); and
 - Temperature-limiting effect on photosynthesis (2 types).

Once we have defined each of the three components, we can take the
product of them with the base model which forms 24
alternative hypotheses, each represented as a `ModelSystem`.

#### Geometric constraint

Plants do not grow indefinitely, but eventually reach an asymptotic mass
owing either to geometric or resource constraints. Here, we define two
competing hypotheses: that a shrub does not show evidence of nearing its
asymptote, or that it does (based on a Chapman-Richards growth function).

In Bristlecone, we use `modelComponent` to construct a pluggable component
into a model system. We pass `modelComponent` a list of `subComponent`s, which
each have a name and an equation. In its equation, a model component can 
take one or many parameters, but these must match the signiture required 
by the hole in the base model. For example, the geometric constraint here 
takes the current `mass` only.

In addition, a `subComponent` may require additional estimatable parameters
over the base model. In this case, the Chapman-Richards model requires an
extra parameter *K*, which represents the asymptotic biomass. These may be
added to a `subComponent` by using `|> estimateParameter` afterwards, as below.
*)

let ``geometric constraint``: Components.ModelComponent<(ModelExpression<g> -> ModelExpression<1>)> =
    let K = parameter "K" (Bounded 1e-12<kg> 100.<kg>) 3.00<kg> 5.00<kg> // Asymptote biomass
    Components.modelComponent
        "Geometric constraint"
        [ Components.subComponent "None" (fun _ -> Constant 1.)
          Components.subComponent "Chapman-Richards" (fun (mass: ModelExpression<g>) ->
            Constant 1. - (mass / Units.toGrams (P K)))
          |> Components.estimateParameter K ]

(**
#### Plant-soil feedback

The plant-soil feedback is the flow of nutrients from plant biomass into the
soil available nitrogen pool. Here, we effectively turn on or off N input into
the soil pool on biomass loss. In the base model, density-dependent biomass
loss occurs. Turning on the feedback here creates an N input into soil based
on the biomass lost multiplied by a conversion factor `ɑ`.

The plant-soil feedback functions should require a biomass in grams, and return
a dimensionless contribution to the N‑availability index.
*)

let ``plant-soil feedback`` =

    let alpha = parameter "ɑ" Positive 0.01</g> 1.00</g>

    let biomassLoss (biomass:ModelExpression<g>) : ModelExpression</year> =
        (P alpha / Constant 100.) * biomass * P gammaB

    Components.modelComponent
        "Plant-Soil Feedback"
        [ Components.subComponent "None" (fun _ -> Constant 0.</year>)
          Components.subComponent "Biomass Loss" biomassLoss
          |> Components.estimateParameter alpha ] // N-recycling efficiency

(**
#### Nitrogen limitation

We specify three plausable mechanisms for nutrient limitation to shrub
growth: (1) that growth is independent on soil N availability; (2) that
growth is dependent on soil N availability in a linear way; or (3) that
a mechanistic model of root-foraging (saturating) best represents
N-limitation of shrub growth.

A new concept here is the `Conditional` element in an equation. This
term exposes a `ModelExpression -> float` (compute), allowing a calculation
to be conditional on the state of parameters or values. In this example,
we use it to restrict the models such that the N-limiting effect cannot
be zero.

Note: in the saturating (Holling) model, r appears both in the uptake
function and in biomass growth, so the model is more sensitive to r
than in the linear or N-independent cases.
*)

let ``N-limitation to growth`` =

    // a = how effectively roots capture available N from the soil (1/g/year).
    // b / r = efficiency of incorporating captured N into biomass (dimensionless or g/unit N).
    // h = integrated handling time — the bottleneck that emerges when both processes (foraging and incorporation) are operating together (years).
    // note: ecologically, h is the saturation effect: even if roots forage efficiently and incorporation is efficient, there’s still a finite rate at which N can be processed.
    let a = parameter "a" Positive 1e-4</g/year> 5e-3</g/year>
    let h = parameter "h" Positive 1e-3<year> 0.250<year>
    let rIndependent = parameter "r" Positive 0.10<g/1> 10.00<g/1> // gram per unit N (N is dimensionless)

    let holling =

        let hollingDiscModelDual =
            ShrubModel.GrowthLimitation.hollingDiscModelDual
                (P a)
                (P r)
                (P h)

        Components.subComponent "Saturating (Holling disc)"
            { UptakeMultiplier = Constant 1. // Turn uptake on.
              UptakeRatePerBiomass      = hollingDiscModelDual }
        |> Components.estimateParameter a
        |> Components.estimateParameter h
        |> Components.estimateParameter r

    let linear =
        Components.subComponent "Linear"
            { UptakeMultiplier = Constant 1. // Turn uptake on.
              UptakeRatePerBiomass      = ShrubModel.GrowthLimitation.linear (P a) }
        |> Components.estimateParameter a
        |> Components.estimateParameter r

    let none =
        Components.subComponent "None" { UptakeRatePerBiomass = (fun _ -> Constant 1.</g/year>); UptakeMultiplier = Constant 0. }
        |> Components.estimateParameter rIndependent

    Components.modelComponent "N-limitation" [ holling; linear; none ]

(**
#### Environmental limitation: temperature

Two model options are stated to determine the role of air temperature in controlling
the growth rate of shrubs. 
*)

let ``temperature limitation to growth`` =

    /// The universal gas constant in J mol−1 K−1
    let gasConstant = Constant 8.314<J mol^1 K^1>

    /// An Arrhenius function to represent temperature limitation on growth.
    /// Form of equation from paper: https://pubag.nal.usda.gov/download/13565/PDF
    let arrhenius (activationEnergy: ModelExpression<J K mol>) (temperature: ModelExpression<K>) =
        Constant System.Math.E ** ((activationEnergy * (temperature - Constant 298.<K>)) / (Constant 298. * gasConstant * temperature))

    let Ea = parameter "Ea" (Bounded 1e-12<kJ mol^1 K^1> 150.<kJ mol^1 K^1>) 10.<kJ mol^1 K^1> 30.<kJ mol^1 K^1>

    Components.modelComponent
        "Temperature limiting effect on photosynthetic rate"
        [
            Components.subComponent "None" (Constant 1.)
            Components.subComponent "Arrhenius" (arrhenius (P Ea |> toJoule) (Environment TMax))
            |> Components.estimateParameter Ea
        ]


(**
#### Putting the hypotheses together

We use functions from the `Hypothesis` module of Bristlecone to scaffold together 
all of the varying model components. Calling `Hypothesis.compile` creates a list
of hypotheses, which reflect the product of all of the components.
*)

let hypotheses =
    ``base model``
    |> Hypotheses.createFromModel
    |> Hypotheses.apply ``geometric constraint``
    |> Hypotheses.apply ``plant-soil feedback``
    |> Hypotheses.apply ``temperature limitation to growth``
    |> Hypotheses.apply ``N-limitation to growth``
    |> Hypotheses.compile
