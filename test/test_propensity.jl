using Julifsp
using Test 

p1 = propensity(x -> x[2])
@test !istimevarying(p1)

p2 = propensity_timevarying((t, x) -> (1.0+sin(t))*x[1]*x[2])
@test istimevarying(p2)
@test !istimeseparable(p2)

p3 = propensity_timevarying(t -> (1.0+sin(t)), x->x[1]*x[2])
@test istimevarying(p3)
@test istimeseparable(p3)

