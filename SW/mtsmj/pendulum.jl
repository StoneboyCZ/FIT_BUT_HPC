using LinearAlgebra 
using Plots
using DifferentialEquations
using BenchmarkTools
using StaticArrays 
using Printf


d = 1;
m = 1;
g = 9.81; 
L = 1;

k = d/m;
a = g/L;

A = @SMatrix BigFloat[    
    -k 0 -a 0
    1 0 0 0
    0 0 0 0
    0 0 0 0
    ]


ij = @SMatrix Int[
    4 1
    3 1
]

B2 = @SMatrix BigFloat[
    0 0
    0 0
    1 0
    0 -1
    ]

B3 = nothing
B4 = nothing
B5 = nothing

ijk = nothing
ijkl = nothing
ijklm = nothing

b = @SVector BigFloat[0,0,0,0]

y0 = BigFloat[0,pi/3,sin(pi/3),cos(pi/3)]

t0 = Float64(0)
tmax = Float64(2*pi)
dt = Float64(tmax/10)
eps = Float64(1e-9)
maxORD = Int8(64)


include("solvers/mtsm_nonlinear.jl")

T, Y, ORD = mtsm_nonlinear(A,b,B2,B3,B4,B5,ij,ijk,ijkl,ijklm,y0,t0,tmax,dt,eps, maxORD)

display(T)
display(Y)
display(ORD)


#=
A = zeros(4,4);
A(1,1) = -data.k;
A(1,3) = -data.a;
A(2,1) = 1;

A2 = zeros(4,2);
A2(3,1) = 1;
A2(4,2) = -1;

ij=[4,1;
    3,1]; =#