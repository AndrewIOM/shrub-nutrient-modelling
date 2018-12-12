module ModelComponents

module ProtectionEffect =

    let none b = b

    let linearSnowProtection (protectionEffect:float) shrubHeightCentimetres snowMass b = 
        // Returns protected mass, which is then multiplied by gammaB
        if shrubHeightCentimetres < (snowMass * protectionEffect)
        then 1.
        else 1.
