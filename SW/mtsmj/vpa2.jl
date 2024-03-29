using LinearAlgebra 
using Plots
using DifferentialEquations
using BenchmarkTools
using StaticArrays 
using Printf

include("solvers/mtsm_linear.jl")

PREC = 128
setprecision(BigFloat, PREC)

omega = BigInt(1); 

A = @SMatrix BigFloat[
    0 omega
    -omega 0
    ]

b = @SVector BigFloat[0,0]

y0 = BigFloat[0,1]

t0 = Float64(0)

#dt = 0.06283185307179587
#tmax = 2*dt
tmax = Float64(6)
dt = Float64(0.1)

eps = Float64(1e-20)

T,Y,ORD = mtsm_linear(A,b,y0,t0,tmax,dt,eps) 

ANAL = sin.(omega*T)

ERR = Y[1,:]-ANAL

ORD