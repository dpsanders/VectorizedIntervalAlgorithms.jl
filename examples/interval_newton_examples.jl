
include("interval_newton_methods.jl")

# 1D
branch_and_prune(x->x^2 - 2x, -5..5, N)

# 2D
f( (x, y) ) = SVector(x^2 + y^2 - 1, y - x)}
@time branch_and_prune(f, IntervalBox(-5..5, 2), N)
@time branch_and_prune(f, IntervalBox(-5..5, 2), K, 1e-5)


"Griewank function"
# G(X) = 1 + sum(abs2, X) / 4000 -
#  prod( cos(X[i] / √(Interval(Float32(i)))) for i in 1:length(X) )

G(X) = 1 + sum(abs2, X) / 4000 -
  prod( cos(X[i] / √(Interval(i))) for i in 1:length(X) )


∇(G) = x -> ForwardDiff.gradient(G, x)

#x = Interval{Float32}(-100..100)
#typeof(x)

#IntervalBox(x, x)

# @time branch_and_prune(∇(G), IntervalBox(-600..600, 1), N, 0.001)

# @time branch_and_prune(∇(G), IntervalBox(Interval{Float32}(-20..20), 2), N, 0.001)
@time branch_and_prune(∇(G), IntervalBox(-10..10, 2), N, 1e-8)

@time branch_and_prune(∇(G), IntervalBox(-100..100, 2), N, 1e-8)
@time branch_and_prune(∇(G), IntervalBox(-100..100, 2), K, 1e-8)


@time branch_and_prune(∇(G), IntervalBox(-600..600, 2), N, 1e-8)


@time branch_and_prune(∇(G), IntervalBox(-600..600, 2), N, 1e-8)
@time branch_and_prune(∇(G), IntervalBox(Interval{Float32}(-100..100), 2), K, 1e-4)
@time branch_and_prune(∇(G), IntervalBox(-600..600, 2), N, 1e-8)
@time branch_and_prune(∇(G), IntervalBox(big(-10..10), 2), K, 1e-8)


@time branch_and_prune(∇(G), IntervalBox(-600..600, 2), N, 1e-8)

@time simple_branch_and_prune(∇(G), IntervalBox(-600..600, 2), 1.0)


using IntervalRootFinding


@time roots(∇(G), IntervalBox(-600..600, 1))
@time roots(∇(G), IntervalBox(-20..20, 2))

@time roots(∇(G), IntervalBox(-100..100, 2), Newton, 1e-5)

X = IntervalBox(Interval{Float32}.(X))

@code_warntype N(f, X)


X

typeof(cos(X[1]))

typeof(ans)



function kin1(X)
    x1, x2, x3, x4, x5, x6 = X

    return SVector(
    -0.4077 + cos(x2)*cos(x6) + cos(x3)*cos(x6) + cos(x4)*cos(x6) + cos(x5)*sin(x2)*sin(x6) - cos(x5)*sin(x3)*sin(x6) - cos(x5)*sin(x4)*sin(x6), 
    -1.9115 + cos(x5)*sin(x1) + cos(x1)*cos(x2)*sin(x5) + cos(x1)*cos(x3)*sin(x5) + cos(x1)*cos(x4)*sin(x5),
    -1.9791 + sin(x2)*sin(x5) + sin(x3)*sin(x5) + sin(x4)*sin(x5),
    -4.0616 + 3*cos(x1)*cos(x2) + 2*cos(x1)*cos(x3) + cos(x1)*cos(x4),
    -1.7172 + 3*cos(x2)*sin(x1) + 2*cos(x3)*sin(x1) + cos(x4)*sin(x1),
    -3.9701 + 3*sin(x2) + 2*sin(x3) + sin(x4)
    )
end

X = IntervalBox(0..2π, 6)


@time basic_branch_prune(kin1, X, K, 1e-6)


