using CUDA, IntervalArithmetic

setrounding(Interval, :accurate)

X = -5..5

"Split interval `X` into `n` pieces"
function mince(X, n)
    nodes = range(X.lo, X.hi, length=n+1)
    return [interval(nodes[i], nodes[i+1]) for i in 1:length(nodes)-1]
end


Xs = mince(X, 10)
f(x) = x * sin(x)

# CPU:
f.(Xs)  # broadcast: apply f element-wise

# GPU:
v = CuArray(Xs)  # send to GPU
f.(v)

# Find roots:
0 .∈  f.(v)
