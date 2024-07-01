using CUDA, IntervalArithmetic   # dps/configure  branch of IntervalArithmetic

make_box(x) = IntervalBox(Interval.(x))   # wrap a vector in a box

function minimize(f, v, numsteps = 10)

    m = +∞       # upper bound for global minimum

    for i in 1:numsteps
        fs = f.(make_box.(mid.(v)))    # interval evaluation at midpoints
        m = min(minimum(sup.(fs)), m)  # update upper bound

        v = v[inf.(f.(v)) .≤ m]        # eliminate if cannot contain global minimum 

        bisected = bisect.(v)             
        v = vcat(first.(bisected), last.(bisected))
    end

    return m, v
end

f( (x, y) ) = (x^2 - 2)^2 + (y^2 - 3)^2

v = [IntervalBox(-5..5, 2)]
minval, minimizers = minimize(f, v, 10);

vv = CuArray(v)
IntervalArithmetic.configure!(directed_rounding=:fast, powers=:fast)

minval, minimizers = minimize(f, vv, 50)
julia> minval, minimizers = minimize(f, vv, 50)
(2.5227748965770077e-14, 
IntervalBox{2,Float64}[[-1.41422, -1.41421] × [1.73205, 1.73206], 
[-1.41422, -1.41421] × [-1.73206, -1.73205], 
...)