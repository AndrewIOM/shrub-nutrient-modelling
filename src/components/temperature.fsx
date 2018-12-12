module ModelComponents

module TemperatureDependence =

    // NB The parameter A is our parameter 'r' or 'Q', hence its absence from the below equations.
    let none _ = 1.
    let linear a t = a * t

    /// Activation energy is in KJ rather than J
    let arrhenius a ea t = a * System.Math.E ** (- ((ea * 1000.) / (8.314 * t)))


module SoilTemperature =

    // Soil temperature is a function of a flux between summer and winter mean temperatures.
    // The flux is moderated by a snow regulation factor. 
    // Snow depth = snow mass / snow density.
    // BUT, here the density effect is collapsed down into the insulation factor, leaving one term to be estimated.

    /// See: 10.1002/2016JG003725
    let snowMassEffect insulationFactor snow =
        System.Math.E ** (-insulationFactor * snow)

    /// Fouriers Law for heat flux
    let localHeatFlux outsideTemp insideTemp conductivity : float =
        - conductivity * (outsideTemp - insideTemp)

module NitrogenReplenishment =

    let linear lambda = lambda

    let temperatureDependent a ea conductivity summerTemp winterTemp =
        conductivity
        |> SoilTemperature.localHeatFlux summerTemp winterTemp // NB are these the right way around?
        |> TemperatureDependence.arrhenius a ea

    let temperatureDependentSnowInsulation a ea insulationFactor snowMass summerTemp winterTemp = 
        snowMass 
        |> SoilTemperature.snowMassEffect insulationFactor
        |> fun con -> temperatureDependent a ea con summerTemp winterTemp

