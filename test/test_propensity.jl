using Julifsp
using Test 

p1 = propensity((x,p) -> x[2])
@test !istimevarying(p1)

p2 = propensity((t, x,p) -> (1.0+sin(t))*x[1]*x[2])
@test istimevarying(p2)
@test !istimeseparable(p2)

p3 = propensity((x,p)->x[1]*x[2], (t,p) -> (1.0+sin(t)))
@test istimevarying(p3)
@test istimeseparable(p3)

