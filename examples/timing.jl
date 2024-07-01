include("goldstein_price.jl")

@time using IntervalArithmetic
@time using ForwardDiff

∇(f) = X -> ForwardDiff.gradient(f, X.v)

@time X = IntervalBox(3..4, 5..6)

@time f(X)

@time ∇(f)(X)
