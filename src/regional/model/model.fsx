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

#r "nuget: Bristlecone.Dendro,2.0.0"

open Bristlecone // Opens Bristlecone core library and estimation engine
open Bristlecone.Language // Open the language for writing Bristlecone models
open Bristlecone.Time

#load "components.fsx"

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

/// ODE1. Cumulative stem biomass
/// NB. A parameter constraint is added using a `Conditional` expression.
/// NB. The parameter r is multiplied by 1000 when applying an N-limitation purely to aid model-fitting.
let ``db/dt`` geomLimit nLimitName nLimitation environmentLimitation =
    Conditional <| fun compute ->
        if compute (Parameter "r") > 100000
        then Invalid
        else
            This * (if nLimitName = "None" then Parameter "r" else Parameter "r" * Constant 1000.) * 
                (nLimitation (``δ15N -> N availability`` (Environment "N"))) * (geomLimit This) * (environmentLimitation (Constant 1.))
                - Parameter "γ[b]" * This

/// ODE2. Soil nitrogen availability
let ``dN/dt`` geomLimit feedback limitationName nLimitation environmentLimitation =
    if limitationName = "None" then
        Parameter "λ" - Parameter "γ[N]" * ``δ15N -> N availability`` This + feedback (Environment "bs")
    else
        Parameter "λ" - Parameter "γ[N]" * ``δ15N -> N availability`` This + feedback (Environment "bs")
        - ((geomLimit (Environment "bs")) * (Environment "bs") * (nLimitation (``δ15N -> N availability`` This))
            * (environmentLimitation (Constant 1.)))

/// Measurement variable: stem radius
let stemRadius lastRadius lastEnv env =
    let oldCumulativeMass = lastEnv |> lookup "bs"
    let newCumulativeMass = env |> lookup "bs"

    if (newCumulativeMass - oldCumulativeMass) > 0. then
        newCumulativeMass |> ModelComponents.Proxies.toRadiusMM // from components.fsx file
    else
        lastRadius

(**
Once we have defined the components, we can scaffold them into a model system.
We can plug in the nestable hypotheses (defined further below) by defining
the base model as a function that takes parameters representing the
alternative hypotheses.
*)

let ``base model`` geometricConstraint plantSoilFeedback environmentLimit (nLimitMode, nLimitation) =
    Model.empty
    |> Model.addEquation "bs" (``db/dt`` geometricConstraint nLimitMode nLimitation environmentLimit)
    |> Model.addEquation "N" (``dN/dt`` geometricConstraint plantSoilFeedback nLimitMode nLimitation environmentLimit)
    |> Model.includeMeasure "x" stemRadius
    |> Model.estimateParameter "λ" notNegative 0.001 0.500
    |> Model.estimateParameter "γ[N]" notNegative 0.001 0.200
    |> Model.estimateParameter "γ[b]" notNegative 0.001 0.200
    |> Model.useLikelihoodFunction (ModelLibrary.Likelihood.bivariateGaussian "x" "N")
    |> Model.estimateParameter "ρ" noConstraints -0.500 0.500
    |> Model.estimateParameter "σ[x]" notNegative 0.001 0.100
    |> Model.estimateParameter "σ[y]" notNegative 0.001 0.100

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

let ``geometric constraint`` =
    modelComponent
        "Geometric constraint"
        [ subComponent "None" (Constant 1. |> (*))
          subComponent "Chapman-Richards" (fun mass -> Constant 1. - (mass / (Parameter "k" * Constant 1000.)))
          |> estimateParameter "K" notNegative 3.00 5.00 ] // Asymptotic biomass (in kilograms)

(**
#### Plant-soil feedback

The plant-soil feedback is the flow of nutrients from plant biomass into the
soil available nitrogen pool. Here, we effectively turn on or off N input into
the soil pool on biomass loss. In the base model, density-dependent biomass
loss occurs. Turning on the feedback here creates an N input into soil based
on the biomass lost multiplied by a conversion factor `ɑ`.
*)

let ``plant-soil feedback`` =

    let biomassLoss biomass =
        (Parameter "ɑ" / Constant 100.) * biomass * Parameter "γ[b]"

    modelComponent
        "Plant-Soil Feedback"
        [ subComponent "None" (Constant 1. |> (*))
          subComponent "Biomass Loss" biomassLoss
          |> estimateParameter "ɑ" notNegative 0.01 1.00 ] // N-recycling efficiency

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

Note: when using the saturating (holling disc) model, the scale of the
parameter r (intrinsic growth rate) becomes different.
*)

let ``N-limitation to growth`` =

    let hollingDiscModelDual =
        ModelComponents.GrowthLimitation.hollingDiscModelDual
            (Parameter "a" / Constant 1000.) // To aid model-fitting
            (Parameter "r")
            (Parameter "h")
            (Constant 5.00)

    let linear =
        ModelComponents.GrowthLimitation.linear
            (Parameter "a" / Constant 1000.) // To aid model-fitting
            (Constant 5.00)

    modelComponent
        "N-limitation"
        [ subComponent "Saturating (Holling disc)" hollingDiscModelDual
          |> estimateParameter "a" notNegative 0.100 0.400
          |> estimateParameter "h" notNegative 0.100 0.400
          |> estimateParameter "r" notNegative 0.500 1.000
          subComponent "Linear" linear
          |> estimateParameter "a" notNegative 0.100 0.400
          |> estimateParameter "r" notNegative 0.500 1.000
          subComponent "None" (Constant 1. |> (*))
          |> estimateParameter "r" notNegative 0.500 1.000 ]

(**
#### Environmental limitation: temperature

Two model options are stated to determine the role of air temperature in controlling
the growth rate of shrubs. 
*)

let ``temperature limitation to growth`` =

    /// The universal gas constant in J mol−1 K−1
    let gasConstant = Constant 8.314

    /// An Arrhenius function to represent temperature limitation on growth.
    /// Form of equation from paper: https://pubag.nal.usda.gov/download/13565/PDF
    let arrhenius activationEnergy temperature =
        Constant System.Math.E ** ((Constant 1000. * activationEnergy * (temperature - Constant 298.)) / (Constant 298. * gasConstant * temperature))

    modelComponent
        "Temperature limiting effect on photosynthetic rate"
        [
            subComponent "None" (Constant 1. |> (*))
            subComponent "Arrhenius" (arrhenius (Parameter "Ea") (Environment "T[max]") |> (*))
            |> estimateParameter "Ea" notNegative 10. 30.
        ]


(**
#### Putting the hypotheses together

We use functions from the `Hypothesis` module of Bristlecone to scaffold together 
all of the varying model components. Calling `Hypothesis.compile` creates a list
of hypotheses, which reflect the product of all of the components.
*)

let hypotheses =
    ``base model``
    |> Hypotheses.createFromComponent ``geometric constraint``
    |> Hypotheses.useAnother ``plant-soil feedback``
    |> Hypotheses.useAnother ``temperature limitation to growth``
    |> Hypotheses.useAnotherWithName ``N-limitation to growth``
    |> Hypotheses.compile
