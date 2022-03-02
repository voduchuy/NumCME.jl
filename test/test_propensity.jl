using EtFsp 
using Test 

p1 = propensity((G, RNA) -> RNA)
@test !istimevarying(p1)

p2 = propensity_timevarying((t, G, RNA) -> (1.0+sin(t))*G*RNA)
@test istimevarying(p2)
@test !istimeseparable(p2)

p3 = propensity_timevarying(t -> (1.0+sin(t)), (G, RNA)->G*RNA)
@test istimevarying(p3)
@test istimeseparable(p3)

