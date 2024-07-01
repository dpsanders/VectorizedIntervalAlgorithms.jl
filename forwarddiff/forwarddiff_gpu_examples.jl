include("forwarddiff_gpu.jl")


v = [IntervalBox(1..2, 3..4)]
vv = CuArray(v)

g( (x, y) ) = x^2 + y

g.(vv)

gradient.(g, vv)

hessian.(g, vv)


## Krawczyk

using LinearAlgebra

IntervalArithmetic.mid(J::AbstractMatrix, α) = mid.(J, α)

"Krawczyk operator"
function K(f, X::IntervalBox{N,T}) where {N,T}

    α = T(0.4792215678)
    m = IntervalBox(mid(X, α))

    J = jacobian(f, X)

    M = mid(J, α)

    if abs(det(M)) < 1e-10
        return X
    end

    # LU = lu(M)  # NB: Not inv(M)

    # if abs(det(LU)) < 1e-5  # singular
    #     return X
    # end

    Y = inv(M)

    return m - Y*f(m) + (I - Y*J) * (X - m)   # if Y == inv(M)

    # return m - LU \ f(m) + (I - LU \ J) * (X - m)
end

# diagonal preconditioner:
# Y = inv(Diagonal(J))
# d = X - m
# return m - Y*f(m) + (I - Y*J) * d

v = [IntervalBox(-5..5, 2)]


g( (x, y) ) = (x - 1)^3 + y^3 - 2y

g.(v)

# K(∇(g), v[1])

K.(∇(g), v)

vv = cu(v)

K.(∇(g), vv)
