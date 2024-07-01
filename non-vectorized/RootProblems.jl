# traits:

# IsPolynomial()

# HasPointSolutions()

# HasMultipleRoots()

# Do we actually want a type?
# Could be a callable object

# Use a function wrapper?

# Automatically categorise based on the expression

# See https://github.com/JuliaHomotopyContinuation/PolynomialTestSystems.jl/blob/master/src/systems.jl

# Maybe add a field `is_polynomial`

# Or make RootsProblem abstract

# abstract type AbstractRootProblem
# struct PolynomialRootProblem
# struct RootProblem

# If don't want to specify roots, make problems parametrised by n

# Maybe have a `solutions` function instead of storing solutions in the problem

using StaticArrays
using IntervalArithmetic
"""
Expresses a root-finding problem in dimension N.
Equations should be a function accepting an SVector and returning an SVector?
"""
struct RootProblem{N,F}
    equation::F
    solutions::Vector{Interval{Float64}}
end

RootProblem(N, equation::F) where {F} = RootProblem{N,F}(equation, Interval{Float64}[])


intersection_circle_line = RootProblem(2,
    ( (x, y), ) -> SVector(x^2 + y^2 - 1, x - y)
    # [SVector(1 / √2, 1 / √2),
    #  SVector(-1 / √2, -1 / √2)
    # ]
)


intersection_circle_ellipse = RootProblem(2,
    ( (x, y), ) -> SVector(x^2 + y^2 - 1,
                           ((x - 0.2) / 2)^2 + ((y - 0.1) / 0.5)^2 - 1
                           )
)


# https://www.juliahomotopycontinuation.org/guides/introduction/
homotopy_continuation_example = RootProblem(2,
    ( (x, y), ) -> SVector(
                    (x^4 + y^4 - 1) * (x^2 + y^2 - 2) + x^5 * y,
                    x^2+2x*y^2 - 2y^2 - 1/2
                    )
    )




using ModelingToolkit

@variables x[1:20]

symbolic(p::RootProblem{N}) where {N} = p.equation(x[1:N])

symbolic(intersection_circle_line)


include("interval_newton_methods.jl")

solve(p::RootProblem, X) = branch_and_prune(p.equation, X, K, 1e-8)

X = IntervalBox(-10..10, 2)
solve(intersection_circle_line, X)

solve(intersection_circle_ellipse, X)

X = IntervalBox(-Inf..Inf, 2)
@time solve(homotopy_continuation_example, X)
