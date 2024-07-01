using IntervalArithmetic
using ForwardDiff
using StaticArrays
using LinearAlgebra

# setrounding(Interval, :accurate)

# ForwardDiff.jacobian(f, X::IntervalBox) = ForwardDiff.jacobian(f, X.v)
# ForwardDiff.gradient(f, X::IntervalBox) = ForwardDiff.gradient(f, X.v)

# include("forwarddiff_gpu.jl")
include("forwarddiff_intervalbox.jl")

"Branch and prune with only bisection"
function simple_branch_and_prune(f, X, ϵ=0.01)

    working = [X]
    possible_roots = typeof(X)[]

    while !isempty(working)
        X = popfirst!(working)

        # check containment
        if !contains_zero(f(X))
            continue
        end

        if diam(X) < ϵ
            push!(possible_roots, X)
            continue
        end

        push!(working, bisect(X)...)
    end

    return possible_roots
end


simple_branch_and_prune(x->x^2 - 2x, -10..10)

mid(0..1, 0.6)

"Interval Newton"
function N(f, fp, X::Interval)
    m = Interval(mid(X, 0.4992212234))

    return m - f(m) / fp(X)
end

N(f, X::Interval) = N(f, x->ForwardDiff.derivative(f, x), X)

N(x->x^2 - 2x, 2..5)

function N(f, fp, X::IntervalBox{M,T}) where {M,T}
    m = mid(X, convert(T, 0.489221))

    J = fp(X)

    return IntervalBox(m .- (J \ f(m)))
end


function KK(f, f′, X::Interval, α=0.489221)
    m = Interval(mid(X, α))
    m = mid(X, α)
    mm = Interval(m)

    if isempty(f(mm))  # outside domain of f
        return emptyinterval(X)
    end

    Y = 1 / f′(m)

    return IntervalBox(m - Y*f(m) + (1 - Y*f′(X)) * (X - m))
end


KK(f, X::Interval) = KK(f, x->ForwardDiff.derivative(f, x), X)


"""
Multi-variable Krawczyk operator.
"""
function KK(f, jacobian, X::IntervalBox{N,T}) where {N,T}

    α = 0.4792215678

    m = IntervalBox(mid(X, T(α)))

    J = jacobian(X)

    Y = inv(mid.(J, α))
    return m - Y*f(m) + (I - Y*J) * (X - m)

    # Y = inv(Diagonal(J))
    # d = X - m
    # return m - Y*f(m) + (I - Y*J) * d
end


N(f, X::IntervalBox) = N(f, x->ForwardDiff.jacobian(f, x), X)
KK(f, X::IntervalBox) = KK(f, x->ForwardDiff.jacobian(f, x), X)

f( (x, y) ) = SVector(x^2 + y^2 - 1, y - x)

X = (0.5..1) × (0.5..1)

X = X ∩ N(f, X)
X = X ∩ N(f, X)

# @code_warntype N(f, X)

"Branch and prune with a uniqueness contractor"
function branch_and_prune(f, X, C, ϵ=0.001)

    working = [X]
    possible_roots = typeof(X)[]
    unique_roots = typeof(X)[]

    while !isempty(working)
        X = pop!(working)

        # check containment
        if !contains_zero(f(X))
            continue
        end

        Y = C(f, X)  # apply contractor e.g. Newton

        if isempty(Y ∩ X)   # no root
            continue


        elseif Y ⊆ X  # exists unique root, so refine it

            push!(unique_roots, refine(f, C, Y))
            continue

        else    # roots are in intersection
            X = Y ∩ X
        end


        if diam(X) < ϵ
            push!(possible_roots, X)
            continue
        end

        push!(working, bisect(X, 0.4899123)...)
    end

    return unique_roots, possible_roots
end

"Refine an interval where have already proved there's a unique solution"
function refine(f, C, X, tol=1e-20)
    while diam(X) > tol  # avoid problem with tiny floating-point numbers if 0 is a root
        NX = C(f, X) ∩ X
        NX == X && break  # reached limit of precision
        X = NX
    end

    return X
end
