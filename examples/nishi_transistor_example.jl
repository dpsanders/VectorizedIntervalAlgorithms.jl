# Reference:
# https://projecteuclid.org/euclid.jjiam/1265033785

# Japan J. Indust. Appl. Math.
# Volume 26, Number 2-3 (2009), 327-336.
# Numerical Existence Proof of Five Solutions for Certain Two-Transistor Circuit Equations
# Yusuke Nakaya, Tetsuo Nishi, Shin'ichi Oishi, and Martin Claus

# Compare KV example/test-nishi.cc
using IntervalArithmetic
using IntervalRootFinding
using StaticArrays

include("interval_newton_methods.jl")


const Rb = 10e3
const Rc = 5e3
const Vcc = -5

const αf = 0.99
const αr = 0.5

# case a
const Is = 1e-9
const VT = 0.053
const Vs = -0.64

# case b
# const Is = 1e-6
# const VT = 0.102
# const Vs = -0.44

const Gb = 1/Rb
const Gc = 1/Rc

const T = SA[
			 1   -αr  0   0
			-αf   1   0   0
			 0    0   1  -αr
			 0    0  -αf  1
			 ]

const G = SA[
			2Gb + Gc   -(Gb + Gc)  -2Gb   Gb
			-(Gb + Gc)  Gb + Gc   Gb   0
			-2Gb   Gb   2Gb + Gc   -(Gb + Gc)
			Gb   0   -(Gb + Gc)   Gb + Gc
			]

const J = SA[Gc*Vcc, Gb*Vs - Gc*Vcc, Gc*Vcc, Gb*Vs - Gc*Vcc]

const factor = Is ./ SA[αf, αr, αf, αr]

ffff(V) = factor .* (exp.(V ./ VT) .- 1)

function nishi(V)
	return T * ffff(V) + G*V + J
end



X = IntervalBox(-10..10, 4)

setrounding(Interval, :tight)

@show @time branch_and_prune(nishi, X, K, 1e-8)
