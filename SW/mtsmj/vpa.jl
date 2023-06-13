using LinearAlgebra 
using Plots
using DifferentialEquations
using BenchmarkTools
using StaticArrays 
using Printf

include("solvers/mtsm_linear.jl")

PREC = 512
setprecision(BigFloat, PREC)

left = BigInt(10)
right = BigInt(10)
h = BigFloat(0.01)

x = range(start=-left*h,step=h,stop=right*h)
I = sin.(x)

indexes = vcat(collect(range(start=-left,step=1,stop=-1)), collect(range(start=1,step=1,stop=right)))

n = size(indexes)[1]

A = zeros(BigFloat,n,n)

A[:,1] = indexes
for i in range(start=2,stop=20,step=1)
    A[:,i] = indexes.^i    
end

b = zeros(BigFloat,n)

for i in range(start=1,stop=20,step=1)
    if left - i >= 0
        b[i] = I[i] - I[left+1]
    else
        b[i] = I[i+1] - I[left+1]
    end
end

DY = A\b

derivOrd = left + right


dy = zeros(BigFloat,derivOrd)
z = zeros(BigFloat,derivOrd)

for i in range(start=1,step=1,stop=derivOrd)
    dy[i] = DY[i] / ((h^i)/factorial(i))

    r = mod(i,4)
    if r == 1
        z[i] = cos(BigFloat(0))
    elseif r == 2
        z[i] = -sin(BigFloat(0))
    elseif r == 3
        z[i] = -cos(BigFloat(0))
    elseif r == 0
        z[i] = sin(BigFloat(0))
    end
end 
    
errorVec = broadcast(abs,z - dy)

@printf("Derivations: Numerical/Analytical(l=%d, k=%d, h=%f, prec=%d)",left,right,h,PREC)
for i in range(start=1,step=1,stop=derivOrd)
    @printf("%d\t%f\t%f\n", i, dy[i], z[i])
end

println("error:",errorVec)