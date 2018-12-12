// Shrub model: transposed from MATLAB

module Götmark2016_ShrubModel =

    let pi = System.Math.PI


    /// Gives the basal radius in centimeters of a stem/branch given its length in centimeters. Function from Niklas and Spatz (2004). 
    let basalRadius k5 k6 stemLength =
        100. * (( 0.01 * stemLength + k6) / k5) ** (3. / 2.) / 2.


    /// Total shrub volume given height and number of stems
    let shrubVolume b a rtip p lmin k5 k6 n h =

        let radius = basalRadius k5 k6

        // Starting value for shrub volume
        let mainStemVolume =
            match radius h with
            | r when r > rtip -> 
                // Modelling stem as a tuncated cone
                n * pi * h * ((radius h) ** 2. + (radius h) * rtip + rtip ** 2.) / 3.
            | _ ->
                // Modelling stem as a cylinder
                n * pi * h * rtip ** 2.

        // Keep adding branches while the length of the parent stem / branch is longer than (2/3)*lmin
        let mutable volume = mainStemVolume
        let mutable k = 0. // Branch generations

        while (p ** k * h > lmin * 2./3.) do

            let volToAdd =
                match (p ** k * h < lmin) with
                | true ->
                    // When a branch is to be added to a stem/branch whose length is less than 
                    // l_min, scale the length of the child branch so that it is zero when the length of parent branch 
                    // (p^k*h) is (2/3)*l_min and p*l_min when the length of the parent branch (p^k*h) is l_min; that
                    // is, the length of the child branch will be 3*p*(p^k*h-2*l_min/3).
                    // This is to make the function continuous. 

                    match (b * 3. * p * (p ** k * h - 2. * lmin / 3.) > rtip) with
                    | true ->
                        // Add the volume of the child branches to every parent branch/stem. 
                        // If b*(length of child branch) is larger than r_tip, model the child branch as a
                        // truncated cone. 
                        n * a * (a + 1.) ** (float k) * pi * 3. * p * (p ** k * h - 2. * lmin / 3.) * ((radius (3. * p * (p ** k * h - 2. * lmin / 3.))) * (radius (3. * p * (p ** k * h - 2. * lmin / 3.)) * rtip + rtip ** 2.)) / 3.

                    | false ->
                        // Add the volume of the child branches to every parent
                        // branch/stem. If b*(length of child branch) is smaller than r_tip, 
                        // model the child branch as a cylinder. 
                        n * a * (a + 1.) ** (float k) * 3. * p * (p ** k * h - 2. * lmin / 3.) * pi * rtip ** 2.

                | false ->

                    match (radius (p ** (k + 1.) * h) > rtip) with
                    | true ->
                        // Add the volume of the child branches to every parent branch/stem. 
                        // If b*(length of child branch) is larger than r_tip, model the child branch as a
                        // truncated cone. 
                        // V_s=V_s + n*a*(a+1)^k*pi*p^(k+1)*h_s*((rad(p^(k+1)*h_s))^2+rad(p^(k+1)*h_s)*r_tip+r_tip^2)/3;
                        n * a * (a + 1.) ** (float k) * pi * p ** (k + 1.) * h * ((radius (p ** (k+1.) * h)) ** 2. + (radius (p ** (k + 1.) * h)) * rtip + rtip ** 2.) / 3.
                    | false ->
                        // V_s=V_s + n*a*(a+1)^k * p^(k+1)*h_s*pi*r_tip^2;
                        // Add the volume of the child branches to every parent
                        // branch/stem. If b*(length of child branch) is smaller than r_tip, 
                        // model the child branch as a cylinder. 
                        n * a * (a + 1.) ** (float k) * p ** (k + 1.) * h * pi * rtip ** 2.
            
            volume <- volume + volToAdd
            k <- k + 1.

        k, volume


    let shrubCrossSectionArea b rtip height numberOfStems =
        if (b * height > rtip) 
        then numberOfStems * pi * (b * height) ** 2.
        else numberOfStems * pi * rtip ** 2.

    let basalAreaToRadius a =
        sqrt (a / System.Math.PI)


let toBasalArea r =
    System.Math.PI * r ** 2.

let toStemLength k5 k6 radius =
    k5 * (2. * radius) ** (2. / 3.) - k6

let toMass woodDensity volume =
    volume * woodDensity

let shrubBiomass b a rtip p lmin k5 k6 n woodDensity radius =
    radius
    |> toStemLength k5 k6
    |> Götmark2016_ShrubModel.shrubVolume b a rtip p lmin k5 k6 n
    |> snd
    |> toMass woodDensity


// Parameters from paper
let a = 2. // the number of child branches added to previous branches (including the tops of the stems) for a shrub
let p = 0.5 // the length of a child branch as a proportion of its parent branch/stem, for parent branches/stems longer than l_min 
let lmin = 20. //cm. the length at which a stem or branch gets child branches
let rtip = 0.1 //cm. the radius of the outermost tip of a stem or branch
let b = 0.0075 // the ratio of the basal radius of a stem or branch and its length

// K5 and K6 are derived empirically from the Yamal stem length - basal radius data
let k5 = 80.200334 // Obtained using least squares fitting
let k6 = 9.918091 // Obtained using least squares fitting
let salixWoodDensity = 0.5 // g / cm3 (from internet)
let numberOfStems = 2.2

let biomass = shrubBiomass b a rtip p lmin k5 k6 numberOfStems salixWoodDensity


