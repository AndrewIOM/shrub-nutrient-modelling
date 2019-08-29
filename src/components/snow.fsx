module ModelComponents.SnowProtection

let none = 0.

let withShrubHeight (protectionEffect:float) shrubHeightCentimetres snowMass = 
    (snowMass * protectionEffect) / (shrubHeightCentimetres + (snowMass * protectionEffect))
