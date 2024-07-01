include("vectorized_interval_root_finding.jl")

using CUDA
CUDA.allowscalar(false)

@time vecroots(x->SVector(x[1]^2 - 2), [IntervalBox(-100..100)])

@time vecroots(x->SVector(x[1]^2 - 2), [IntervalBox(-100f0..100f0)])


@time meta_branch_and_bound(x->SVector(x[1]^2 - 2), [IntervalBox(-10..10)], 0.000001, 1)

@time meta_branch_and_bound(x->SVector(x[1]^2 - 2),
    [IntervalBox(-10f0..10f0)], 0.000001, 1)



G(X) = 1 + sum(abs2, X) / 4000 -
  prod( cos(X[i] / √(Float32(i))) for i in 1:length(X) )

@time rts = meta_branch_and_bound(∇(G), [IntervalBox(-100..100, 2)])

@time rts = meta_branch_and_bound(∇(G),
    [IntervalBox(-100f0..100f0, 2)])

using CUDA
CUDA.allowscalar(false)

G2( (x, y) ) = 1 + (x^2 + y^2) / 4000 - cos(x) * cos(y / √2)

G3( (x, y, z) ) = 1 + (x^2 + y^2 + z^2) / 4000 - cos(x) * cos(y / √2) * cos(z / √3)

p( (x, y) ) = sin(x) * cos(y)

X = IntervalBox(-60..60, 2)
@time rts = meta_branch_and_bound(∇(G2), [X])

CUDA.@time frue)

CUDA.@time rts_gpu, unknown = meta_branch_and_bound(∇(G2), [IntervalBox(-3000..3000, 2)], 1e-8, 10^7, cuda=true)


CUDA.@time rts_gpu, unknown = meta_branch_and_bound(∇(p), [IntervalBox(-1200..1200, 2)], 1e-8, 10^7, cuda=true)

CUDA.@time rts_gpu, unknown = meta_branch_and_bound(∇(p), [IntervalBox(-1200..1200, 2)], 1e-8, 10^7, cuda=false)

length(rts_gpu), length(unknown)

X = IntervalBox(-3000..3000, 2)

G2(X)


rts[1][rts_gpu[1] .!= rts[1]]
# they differ only with very small values

X = IntervalBox(-30..30, 3)
# IntervalBox(-600..600, 3) gives 45375381 roots
# in 345 seconds on GPU

@time rts = meta_branch_and_bound(∇(G3), [X], 1e-8)


CUDA.@time rts_gpu_3d, unknowns = meta_branch_and_bound(∇(G3), [X], 1e-8, 10^7, cuda=true)


rts_gpu_3d

rts_gpu_3d

unknowns

@time rts_gpu_3d_2, unknowns = meta_branch_and_bound(∇(G3), unknowns, 1e-8, cuda=true)

vecroots(∇(G3), unknowns, 1e-8)
