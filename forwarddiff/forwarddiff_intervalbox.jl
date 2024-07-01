using ForwardDiff


import ForwardDiff: gradient, jacobian, hessian

∇(f::F) where {F} = x -> gradient(f, x)

@inline gradient(f::F, X::IntervalBox{N,T}) where {F,N,T} = gradient(f, X.v)

@inline function jacobian(f::F, X::IntervalBox{N,S}) where {F,N,S}
    x = X.v
    T = typeof(ForwardDiff.Tag(f, eltype(x)))
    return ForwardDiff.extract_jacobian(T, ForwardDiff.static_dual_eval(T, f, x), x)
end

@inline hessian(f::F, X::IntervalBox{N,T}) where {F,N,T} = jacobian(x -> gradient(f, x), X)
