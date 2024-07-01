


f(x) = (x^2 - 2)^2

X = -10..10

rts = vecroots_bisection(f, [X], 20)


g( (x, y) ) = SVector(x^2 + y^2 - 1, x - y)
X = IntervalBox(-5..5, 2)

rts, unknown = vecroots(g, [X], 20)

Y = rts[1]
Y = K(g, Y)

rts, unknown = vecroots(g, CuArray([X]), 20)




# function h(X::IntervalBox{N,T}) where {N,T}
#     x1,x2,x3,x4,x5,x6,x7,x8,x9,x10 = X
#
#     p = SVector{23,T}(3600, 18, 18, 3600, 3600, 3600, 1800, 3600, 18, 18, 18, 1800, 18, 36, 11, 180, 0.7, 0.4, 30, 0.2, 4, 4.5, 0.4)
#
#     return IntervalBox(
#
#         -p[17] * x1 + -2 * p[1] * (x1 ^ 2 / 2) + 2 * p[2] * x2 + p[21] * p[18] * ((1 + p[19] * x7) / (p[20] + x7)),
#         -p[17] * x2 + p[1] * (x1 ^ 2 / 2) + -1 * p[2] * x2 + -1 * p[4] * x2 * x4 + p[9] * x3 + p[14] * x3 + -1 * p[6] * x2 * x7 + p[11] * x8,
#         -p[17] * x3 + p[4] * x2 * x4 + -1 * p[9] * x3 + -1 * p[5] * x3 * x4 + p[10] * x5 + -1 * p[14] * x3 + p[15] * x5 + p[7] * x8 * x4 + -1 * p[12] * x3 * x7,
#         -p[17] * x4 + -1 * p[4] * x2 * x4 + p[9] * x3 + -1 * p[5] * x3 * x4 + p[10] * x5 + -1 * p[7] * x8 * x4 + p[12] * x3 * x7 + p[16] * x9 + p[22] * p[18] * ((1 + p[19] * x7) / (p[20] + x7)),
#         -p[17] * x5 + p[5] * x3 * x4 + -1 * p[10] * x5 + -1 * p[15] * x5,
#         -p[17] * x6 + p[14] * x3 + p[15] * x5 + -1 * p[8] * x6 * x10 + p[13] * x9,
#         -p[17] * x7 + -1 * p[6] * x2 * x7 + p[11] * x8 + p[7] * x8 * x4 + -1 * p[12] * x3 * x7 + p[18] * ((1 + p[19] * x7) / (p[20] + x7)),
#         -p[17] * x8 + p[6] * x2 * x7 + -1 * p[11] * x8 + -1 * p[7] * x8 * x4 + p[12] * x3 * x7,
#         -p[17] * x9 + p[8] * x6 * x10 + -1 * p[13] * x9 + -1 * p[16] * x9,
#         x10 + x9 - p[23]
#         )
# end

function h(X::StaticVector)
    x1,x2,x3,x4,x5,x6,x7,x8,x9,x10 = X

    p = SVector(3600, 18, 18, 3600, 3600, 3600, 1800, 3600, 18, 18, 18, 1800, 18, 36, 11, 180, 0.7, 0.4, 30, 0.2, 4, 4.5, 0.4)

    return IntervalBox(

        -p[17] * x1 + -2 * p[1] * (x1 ^ 2 / 2) + 2 * p[2] * x2 + p[21] * p[18] * ((1 + p[19] * x7) / (p[20] + x7)),
        -p[17] * x2 + p[1] * (x1 ^ 2 / 2) + -1 * p[2] * x2 + -1 * p[4] * x2 * x4 + p[9] * x3 + p[14] * x3 + -1 * p[6] * x2 * x7 + p[11] * x8,
        -p[17] * x3 + p[4] * x2 * x4 + -1 * p[9] * x3 + -1 * p[5] * x3 * x4 + p[10] * x5 + -1 * p[14] * x3 + p[15] * x5 + p[7] * x8 * x4 + -1 * p[12] * x3 * x7,
        -p[17] * x4 + -1 * p[4] * x2 * x4 + p[9] * x3 + -1 * p[5] * x3 * x4 + p[10] * x5 + -1 * p[7] * x8 * x4 + p[12] * x3 * x7 + p[16] * x9 + p[22] * p[18] * ((1 + p[19] * x7) / (p[20] + x7)),
        -p[17] * x5 + p[5] * x3 * x4 + -1 * p[10] * x5 + -1 * p[15] * x5,
        -p[17] * x6 + p[14] * x3 + p[15] * x5 + -1 * p[8] * x6 * x10 + p[13] * x9,
        -p[17] * x7 + -1 * p[6] * x2 * x7 + p[11] * x8 + p[7] * x8 * x4 + -1 * p[12] * x3 * x7 + p[18] * ((1 + p[19] * x7) / (p[20] + x7)),
        -p[17] * x8 + p[6] * x2 * x7 + -1 * p[11] * x8 + -1 * p[7] * x8 * x4 + p[12] * x3 * x7,
        -p[17] * x9 + p[8] * x6 * x10 + -1 * p[13] * x9 + -1 * p[16] * x9,
        x10 + x9 - p[23]
        )
end



function f(X::SVector{N,T}, p) where {N,T}
    x1,x2,x3,x4,x5,x6,x7,x8,x9,x10 = X

    return SVector(

        -p[17] * x1 + -2 * p[1] * (x1 ^ 2 / 2) + 2 * p[2] * x2 + p[21] * p[18] * ((1 + p[19] * x7) / (p[20] + x7)),
        -p[17] * x2 + p[1] * (x1 ^ 2 / 2) + -1 * p[2] * x2 + -1 * p[4] * x2 * x4 + p[9] * x3 + p[14] * x3 + -1 * p[6] * x2 * x7 + p[11] * x8,
        -p[17] * x3 + p[4] * x2 * x4 + -1 * p[9] * x3 + -1 * p[5] * x3 * x4 + p[10] * x5 + -1 * p[14] * x3 + p[15] * x5 + p[7] * x8 * x4 + -1 * p[12] * x3 * x7,
        -p[17] * x4 + -1 * p[4] * x2 * x4 + p[9] * x3 + -1 * p[5] * x3 * x4 + p[10] * x5 + -1 * p[7] * x8 * x4 + p[12] * x3 * x7 + p[16] * x9 + p[22] * p[18] * ((1 + p[19] * x7) / (p[20] + x7)),
        -p[17] * x5 + p[5] * x3 * x4 + -1 * p[10] * x5 + -1 * p[15] * x5,
        -p[17] * x6 + p[14] * x3 + p[15] * x5 + -1 * p[8] * x6 * x10 + p[13] * x9,
        -p[17] * x7 + -1 * p[6] * x2 * x7 + p[11] * x8 + p[7] * x8 * x4 + -1 * p[12] * x3 * x7 + p[18] * ((1 + p[19] * x7) / (p[20] + x7)),
        -p[17] * x8 + p[6] * x2 * x7 + -1 * p[11] * x8 + -1 * p[7] * x8 * x4 + p[12] * x3 * x7,
        -p[17] * x9 + p[8] * x6 * x10 + -1 * p[13] * x9 + -1 * p[16] * x9,
        x10 + x9 - p[23]
        )
end

function f(X::IntervalBox{N,T}, p) where {N,T}
    x1,x2,x3,x4,x5,x6,x7,x8,x9,x10 = X

    return IntervalBox(

        -p[17] * x1 + -2 * p[1] * (x1 ^ 2 / 2) + 2 * p[2] * x2 + p[21] * p[18] * ((1 + p[19] * x7) / (p[20] + x7)),
        -p[17] * x2 + p[1] * (x1 ^ 2 / 2) + -1 * p[2] * x2 + -1 * p[4] * x2 * x4 + p[9] * x3 + p[14] * x3 + -1 * p[6] * x2 * x7 + p[11] * x8,
        -p[17] * x3 + p[4] * x2 * x4 + -1 * p[9] * x3 + -1 * p[5] * x3 * x4 + p[10] * x5 + -1 * p[14] * x3 + p[15] * x5 + p[7] * x8 * x4 + -1 * p[12] * x3 * x7,
        -p[17] * x4 + -1 * p[4] * x2 * x4 + p[9] * x3 + -1 * p[5] * x3 * x4 + p[10] * x5 + -1 * p[7] * x8 * x4 + p[12] * x3 * x7 + p[16] * x9 + p[22] * p[18] * ((1 + p[19] * x7) / (p[20] + x7)),
        -p[17] * x5 + p[5] * x3 * x4 + -1 * p[10] * x5 + -1 * p[15] * x5,
        -p[17] * x6 + p[14] * x3 + p[15] * x5 + -1 * p[8] * x6 * x10 + p[13] * x9,
        -p[17] * x7 + -1 * p[6] * x2 * x7 + p[11] * x8 + p[7] * x8 * x4 + -1 * p[12] * x3 * x7 + p[18] * ((1 + p[19] * x7) / (p[20] + x7)),
        -p[17] * x8 + p[6] * x2 * x7 + -1 * p[11] * x8 + -1 * p[7] * x8 * x4 + p[12] * x3 * x7,
        -p[17] * x9 + p[8] * x6 * x10 + -1 * p[13] * x9 + -1 * p[16] * x9,
        x10 + x9 - p[23]
        )
end

g(X::SVector{N,T}) where {N,T} = f(X, SVector{23,T}(3600, 18, 18, 3600, 3600, 3600, 1800, 3600, 18, 18, 18, 1800, 18, 36, 11, 180, 0.7, 0.4, 30, 0.2, 4, 4.5, 0.4))

g(X::IntervalBox{N,T}) where {N,T} = f(X, SVector{23,T}(3600, 18, 18, 3600, 3600, 3600, 1800, 3600, 18, 18, 18, 1800, 18, 36, 11, 180, 0.7, 0.4, 30, 0.2, 4, 4.5, 0.4))


#g(X::IntervalBox) = g(X.v)
#
# h(X::IntervalBox) = h(X.v)
#
# let λ=3.5, α=@SMatrix rand(7, 7)
#     global gg = X -> f(X, λ, α)
# end
#
# using IntervalArithmetic, StaticArrays
#
# X = IntervalBox(0..10, 7)

import Base: *
*(M::AbstractMatrix, X::IntervalBox) = IntervalBox(M * X.v)

using LinearAlgebra


isnotempty(X) = !(isempty(X))


function apply_krawczyk(K, X)
    KX = K.(X)

    is_root = KX .⪽ X

    found_roots = X[is_root]

    @show length(found_roots)

    non_roots = is_root .== false

    Xnew = KX[non_roots] .∩ X[non_roots]

    return found_roots, Xnew[isnotempty.(Xnew)]
end


# From Yingbo Ma:

"""
    jacvec(f, x, v) -> u
``J(f(x))*v``
"""
function jacvec(f, x, v)
    du = ForwardDiff.Dual{:___jacvec_tag}.(x, v)
    ForwardDiff.partials.(f(du), 1)
end

function cycle(K, found_roots, rts, n=2)
           rts = vecroots(gg, rts, 1)

           for i in 1:n

               new_found, rts = apply_krawczyk(K, rts)

               @show new_found

               if !isempty(new_found)
                   push!(found_roots, new_found)
               end

           end
           @show length(rts)



           return found_roots, rts
       end

f(N, λ, α) = N .- ( (λ .* N) ./ (1.0 .+ α * N) )

function make_functions(n)
    λ = 3.5
    α = @SMatrix rand(n, n)
    g = X -> f(X, λ, α)
    J(X) = jacobian(g, X)
    K(X) = 𝒦(g, J, X, 0.499754)

    return g, J, K
end

XX(n) = CuArray([IntervalBox(-0.01..1e4, n)])

for i in 1:10
    global Y = K(g, Y)
    @show Y
end




g( (x, y) ) = (x - 1)^3 + y^3 - 2y

g.(v)

# K(∇(g), v[1])

K.(∇(g), v)

vv = cu(v)

K.(∇(g), vv)





CUDA.@time rts, unknown = vecroots(g, CuArray([X]), 20)

@time vecroots(g, ([X]), 20)

include("nishi_transistor_example.jl")

X = IntervalBox(-10..10, 4)

@time rts, unknown = vecroots(X -> nishi(X, dT, G, J, factor), [X], 100)

maximum(diam.(unknown))

vv = CuArray([X])

CUDA.allowscalar(true)

nishi.(vv)

@time rts, unknown = vecroots(X -> nishi(X, T, G, J, factor), CuArray([X]), 100)
