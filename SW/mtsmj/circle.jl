using LinearAlgebra
using Plots
using DifferentialEquations
using BenchmarkTools
using StaticArrays 

include("solvers/mtsm_linear.jl")

function sincos!(du,u,p,t)
    omega = p;
    du[1] = omega*u[2]
    du[2] = -omega*u[1]
end

omega = 1; 

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

eps = Float64(1e-9)

T,Y,ORD = mtsm_linear(A,b,y0,t0,tmax,dt,eps) 
#display(T)
#display(Y[1,:])

p = omega
prob = ODEProblem(sincos!, BigFloat[0;1], (t0, tmax),p)

sol = solve(prob,DP5())

plot(sol)

plot(T,Y[1,:])
plot(T,Y[2,:])
scatter(T,ORD)
