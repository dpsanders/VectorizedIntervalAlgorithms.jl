
using Revise
using IntervalArithmetic
using StaticArrays

γ(a, b, x) = (b / a) ∩ x   #  only when 0 ∉ a

function gauss_seidel(A, b, X::IntervalBox{N,T}) where {N,T}

    for i in 1:N

        # @inbounds new_val = γ(A[i, i],
        #             b[i] - sum(A[i, j] * X[j] for j in 1:N if j ≠ i),
        #             X[i]
        #             )

        new_val = zero(X[i])

        # new_val = b[i] - sum(ifelse(i==j, zero(X[j]), A[i, j] * X[j]) for j in 1:N)

        @inbounds for j in 1:N
            j == i && continue
            new_val += A[i, j] * X[j]
        end

        # new_val = dot(A[i, :], X) - A[i, i] * X[i]

        @inbounds new_val = ( (b[i] - new_val) / A[i, i] ) ∩ X[i]


        @inbounds X = setindex(X, new_val, i)
    end

    return X
end

A = interval.(SA[1 2; 3 4])
b = interval.(SA[5, 6])

X = IntervalBox(-5..5, 2)

@code_warntype gauss_seidel(A, b, X)

gauss_seidel(A, b, X)

using BenchmarkTools

@btime gauss_seidel($A, $b, $X)




        # new_val = zero(X[i])

        # new_val = b[i] - sum(ifelse(i==j, zero(X[j]), A[i, j] * X[j]) for j in 1:N)

        # @inbounds for j in 1:N
        #     j == i && continue
        #     new_val += A[i, j] * X[j]
        # end
        #
        # @inbounds new_val = ( (b[i] - new_val) / A[i, i] ) ∩ X[i]


#
# function standard_gauss_seidel!(A, b, X)
#     N = length(X)
#
#     for i in 1:N
#
#         X[i] = (b[i] - sum(A[i, j] * X[j] for j in 1:N if j ≠ i)) / A[i, i]
#
#     end
#
#     return X
# end
#
# X = [3.0, 4.0]
# standard_gauss_seidel!(A, b, X)
#
#
# @code_warntype standard_gauss_seidel!(A, b, X)
#
#
#
# struct GaussSeidelDiff{U, V}
#     A::U
#     X::V
# end
#
# (gsd::GaussSeidelDiff)(i, j) = gsd.A[i, j] * gsd.X[j]
#
# function gauss_seidel(A, b, X::IntervalBox{N,T}) where {N,T}
#     gsd = GaussSeidelDiff(A, X)
#     for i in 1:N
#         @inbounds new_val = γ(A[i, i],
#                     b[i] - sum(gsd(i, j) for j in 1:N if j ≠ i),
#                     X[i]
#                     )
#         X = setindex(X, new_val, i)
#     end
#     return X
# end

gauss_seidel(A, b, X)

@code_warntype gauss_seidel(A, b, X)

@btime gauss_seidel($A, $b, $X)
#
#
# mutable struct GaussSeidelDiff2{U, V}
#     A::U
#     X::V
# end
# (gsd::GaussSeidelDiff2)(i, j) = gsd.A[i, j] * gsd.X[j]
# function gauss_seidel(A, b, X::IntervalBox{N,T}) where {N,T}
#     gsd = GaussSeidelDiff2(A, X)
#     for i in 1:N
#         @inbounds new_val = γ(A[i, i],
#             b[i] - sum(gsd(i, j) for j in 1:N if j ≠ i),
#             X[i]
#         )
#         gsd.X = setindex(gsd.X, new_val, i)
#     end
#     return gsd.X
# end

@btime gauss_seidel($A, $b, $X)

@code_warntype gauss_seidel(A, b, X)

@btime setindex($X, 1..2, 2)


## Example from Moore-Kearfott-Cloud pg. 92

A = [3 1; 3 2]
b = [1, 0]

X = IntervalBox(-1.12..1.12, 2)

Y = A^-1

A2 = Y * A
b2 = Y * b
X = gauss_seidel(A, b, X)
