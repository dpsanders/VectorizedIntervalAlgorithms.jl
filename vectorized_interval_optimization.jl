# Vectorized interval optimization:

# using CUDA,
using StaticArrays
using IntervalArithmetic # using IntervalOptimisation
using IntervalConstraintProgramming, ModelingToolkit
using ForwardDiff
using LinearAlgebra
using GPUifyLoops
using CuArrays

IntervalArithmetic.configure!(directed_rounding=:fast, powers=:fast)

# cd("/Users/dpsanders/Dropbox/code/interval_root_finding")

# include("forwarddiff_intervalbox.jl")

# ForwardDiff.gradient(f, X::Interval) = ForwardDiff.derivative(f, X)

function bisect_array(X)
    if isempty(X)
        return X
    end

    bisected = bisect.(X)

    return vcat(first.(bisected), last.(bisected))
end


wrap(x::Float64) = Interval(x)

wrap(X::SVector{N,T}) where {N,T} = IntervalBox(Interval.(X))

interval_mid(X::Interval{T}) where {T} = Interval(mid(X))
interval_mid(X::IntervalBox{N,T}) where {N,T} = wrap(mid(X))
interval_mid(X::SVector{N,T}) where {N,T} = wrap(mid(X))

import IntervalArithmetic: mid
mid(X::SVector{N,T}) where {N,T} = mid.(X)


vecminimize(f, X, numsteps=10) = vecminimize(f, [X], numsteps)

"Standard mean-value form for a scalar-valued function f"
# function mean_value_form(f, X)
#     m = interval_mid(X)
#     return f(m) + gradient(f, X) ⋅ (X - m)
# end
#
# mean_value_form(f) = X -> mean_value_form(f, X)
#
# f(x) = x^2 - 2x
# X = -1..1
#
# f(X)
#
#
# mean_value_form(f, X)


X = IntervalBox(1..2, 3..4)
inf.(X)

using Random

function take_random_sample(v::AbstractVector{Interval{T}}) where {T}

    ds = diam.(v)
    rs = similar(ds)

    return Interval.(inf.(v) .+ ds .* rs)
end


"Sample one random point in every box"
function take_random_sample(v::AbstractVector{IntervalBox{N,T}}) where {N,T}
    intervals = reinterpret(Interval{T}, v)  # split up intervalboxes into intervals

    return reinterpret(IntervalBox{N,T}, take_random_sample(intervals))
end

X = IntervalBox(-10..10, 2)

intervals = mince(X, 4)

# using CUDA
take_random_sample(CuArray(intervals))

"Minimize the function f over the boxes in Xs. A known upper bound
`m` may be specified. A contractor C may be given"
function vecminimize(f, Xs::AbstractVector, numsteps = 10; m = +∞, C=nothing)
    intervals = Xs

    for i in 1:numsteps

        # find possible new upper bound for global minimum
        # by evaluating f at the midpoint of each box:

        midpoint_images = sup.(f.(interval_mid.(intervals)))
        midpoint_images = midpoint_images[midpoint_images .> -Inf]
        # > -Inf ignores boxes whose images are empty

        if !isempty(midpoint_images)
            m = min( minimum(midpoint_images), m )
        end


        # evaluate at one random point in each box
        random_points = sup.(f.(take_random_sample(intervals)))
        random_points = random_points[random_points .> -Inf]

        if !isempty(random_points)
            m = min( minimum(random_points), m )
        end

        if C != nothing  # use contractor
            Z = IntervalBox(Interval(-∞, m))
            intervals = C.(Ref(Z), intervals)
            intervals = filter(!isempty, intervals)

        else

            # find which intervals possibly contain the global minimum:
            interval_images = inf.(f.(intervals))
            intervals = intervals[interval_images .≤ m]

        end

        # bisect to produce new interval list:

        intervals = bisect_array(intervals)


        # @show i, length(intervals), m
        # @show diam(reduce(∪, to_bisect, init=emptyinterval(X))[1])

    end

    images = f.(intervals)
    bound = Interval(minimum(inf.(images)), maximum(sup.(images)))
    # bound = reduce(∪, f.(intervals), init=emptyinterval(X))
    bound = bound ∩ Interval(-∞, m)

    return m, bound, intervals

end
#
# function vecminimize(C::IntervalConstraintProgramming.AbstractContractor, Xs::AbstractVector{IntervalBox{N,T}}, numsteps = 10, m = ∞) where {N,T}
#     intervals = Xs
#
#     for i in 1:numsteps
#         # obtain possible new upper bound for global minimum:
#         midpoint_images = sup.(C.(interval_mid.(intervals)))
#         m = min( minimum(midpoint_images[midpoint_images .> -Inf]), m )
#
#         # find which intervals possibly contain the global minimum:
#         contracted = similar(intervals)
#         contracted .= C.(Interval{T}(-Inf, m), intervals)
#         to_bisect = filter(!isempty, contracted)
#
#         # bisect to produce new interval list:
#         #bisected = bisect.(to_bisect)
#         intervals = bisect_array(to_bisect)
#         # intervals = vcat(first.(bisected), last.(bisected))
#
#         @show i, length(intervals), m
#         # @show to_bisect
#         # @show diam.(reduce(∪, to_bisect, init=emptyinterval(X)))
#     end
#
#     return m, intervals
#
# end

vecminimize(x->(x^2 - 2)^2 - 2, [-10..10])

# vecminimize(x->(x^2 - 2)^2 - 2, CuArray([-10..10]))

# combined_mean_value(f, X) = f(X) ∩ mean_value_form(f, X)
# combined_mean_value(f) = X -> combined_mean_value(f, X)

include("lennard_jones.jl")


#
# X3 = IntervalBox(0..6, 3)
# vecminimize(V3, X3, 20)
#
# @time vecminimize(V3, X3, 60)
#
# @time vecminimize(combined_mean_value(V3), X3, 60)
#
#
# X4 = IntervalBox(0..6, 6)
# @time vecminimize(V4, X4, 40)
# @time vecminimize(combined_mean_value(V4), X4, 40)
#
# X4 = IntervalBox([0.1*i..6 for i in 1:6])
# @time vecminimize(V4, X4, 40)
#
#
# X5 = IntervalBox([0.01*i..3 for i in 1:9])
# # @time vecminimize(mean_value_form(V5), X5, 30)

#
# vars = @variables x2, y2, y3
#
# C = Contractor(vars, V3)
#
# X = IntervalBox(0..6, 3)
#
# v = CuArray([X])
#
# C.(-Inf..10, v)
#
# vecminimize(C, v, 10)
#
#
# f(x) = (x^2 - 2)^2
#
# X = -10..10
#
# minimizers = vecminimize(f, X, 30)
#
#
#
# diam.(minimizers)
#
# g(x, y) = (x^2 - 2)^2 + 2*x*y + 3y^2
# g(xx) = g(xx[1], xx[2])
#
# X = IntervalBox(-10..10, 2)
# @time vecminimize(g, X, 30)
#
# @time vecminimize(mean_value_form_scalar(g), IntervalBox(-10..10, 2), 40)
#
# using IntervalOptimisation
#
# # @time IntervalOptimisation.minimise(g, X)
#
# # @time IntervalOptimisation.minimise(mean_value_form_scalar(g), X)
#
#     @time vecminimize(mean_value_form_scalar(g), X, 30)



vars = (@variables x[1:9])[1]

C = Contractor(vars, V5)

X = IntervalBox(-6..6, 9)

@time vecminimize(V5, [X], C=C, 16)

using CuArrays


make_basic_contractor(C) = BasicContractor(contextualize(C.forward.f), contextualize(C.backward.f))

C2 = make_basic_contractor(C)

inf(3..4)

2^16

@time vecminimize(V5, CuArray([X]), 10, C=C2)


v = CuArray([X])

Z = IntervalBox(1..2)

C2.(Ref(Z), v)

C = Contractor(vars[1:3], V3)

C2 = make_basic_contractor(C)

X3 = IntervalBox(-6..6, 3)

v = CuArray([X3])

C2.(Ref(Z), v)
