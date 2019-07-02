module ModelComponents.SnowProtection

let none b = b

let linearSnowProtection (protectionEffect:float) shrubHeightCentimetres snowMass b = 
    // Returns protected mass, which is then multiplied by gammaB
    // if shrubHeightCentimetres < (snowMass * protectionEffect)
    // then protectionEffect * snowMass
    // else protectionEffect * snowMass
    (snowMass * protectionEffect * shrubHeightCentimetres) / 1. + (snowMass * protectionEffect * shrubHeightCentimetres)
