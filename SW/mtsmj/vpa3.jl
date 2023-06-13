using LinearAlgebra 
using Plots
using DifferentialEquations
using BenchmarkTools
using StaticArrays 
using Printf

include("solvers/mtsm_linear.jl")

PREC = 128
setprecision(BigFloat, PREC)

omega = BigInt(100); 

A = @SMatrix BigFloat[
        1
    ]

b = @SVector BigFloat[0]

y0 = BigFloat[1]

t0 = Float64(0)

#dt = 0.06283185307179587
#tmax = 2*dt
tmax = Float64(2*pi)
dt = Float64(pi/100)

eps = Float64(1e-20)

T,Y,ORD = mtsm_linear(A,b,y0,t0,tmax,dt,eps) 

ANAL = sin.(omega*T)


ERR = Y[1,:]-ANAL

ORD