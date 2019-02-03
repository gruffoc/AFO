const λ  = 0.002 #in mm => 150GHz
const θᵦ = (π / 180) * 0.05#   0.17 deg for STRIP W-band
const ω₀ = λ / ( π * θᵦ ) # θᵦ is the beam opening angle in radians
const z_atm = 3E4 #depends on the observation site
const Tᵧ = 240 # ground temperature in kelvin
const χ₁_0 = 1 # in m^{-2}
const χ₂_0 = 1 # in m^{-2}
const z_0  = 1E3 # the quote where water vapor vanish
const L₀   = 674 # in m, the turbulance correlation length
