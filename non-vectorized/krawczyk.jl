
using LinearAlgebra

using ForwardDiff

jacobian(f, X::IntervalBox) = ForwardDiff.jacobian(f, X.v)

IntervalArithmetic.mid(J::AbstractMatrix, α) = mid.(J, α)

"Krawczyk operator"
function K(f::F, X::IntervalBox{N,T}) where {F,N,T}

    α = 0.4792215678   #  so works with Float32?
    m = IntervalBox(Interval.(mid(X, α)))

    J = jacobian(f, X)
    M = mid(J, α)

    if abs(det(M)) < 1e-10
        return X
    end

    # LU = lu(M)  # NB: Not inv(M)

    # if abs(det(LU)) < 1e-5  # singular
    #     return X
    # end

    Y = inv(M)   # M is a matrix of floats

    return m - Y * f(m) + (I - Y * J) * (X - m)   # if Y == inv(M)

    # return m - LU \ f(m) + (I - LU \ J) * (X - m)
end

# diagonal preconditioner:
# Y = inv(Diagonal(J))
# d = X - m
# return m - Y*f(m) + (I - Y*J) * d
