"""
    mtsm_linear(A,b,y0,t0,tmax,dt,eps)    

MTSM linear solver. 

...
# Arguments
- `n::Matrix`: Jacobian matrix.
- `b::Vector`: RHS vector.
- `y0::Vector`: initial conditions.
- `t0::Float`: initial time.
- `tmax::Float`: maximum time.
- `dt::Float`: size of the integration step.
- `eps::Float`: accuracy of the calculation.
...

# Examples
```julia-repl
julia> 
1
```
"""
function mtsm_linear(A,b,y0,t0,tmax,dt,eps)
    T = Float64[t0]
    Y = y0
    ORD = [0]
    t = t0+dt
    nStopping = 3 

    while true
        y = y0
        Ay = A*y+b
        DY = dt*Ay
        y+= DY

        stopping = ones(BigFloat,nStopping)*1e10
        stopping[1] = abs(reduce((x,y) -> max.(x,y), DY))
        
        k = 2
        while norm(abs.(DY)) > eps
            DY = BigFloat(dt/k)*A*DY
            y+=DY
           
            stopping[mod(k-1,nStopping)+1] = abs(reduce((x,y) -> max.(x,y), DY))
            k+=1
        end # calculation in the current time ends

        Y = [Y y] # save the result
        push!(ORD,k) # save the value of ORD 
        push!(T,t) # save the current time

        if t == tmax # if the final time was achieved, end the calculation
            break
        end
        
        if (t+dt > tmax) # if the following step would overstep
            dt = tmax - t # correct the step size to the appropriate size
        end

        t = t+dt
        y0 = y


        #=
        if (t+dt < tmax)
            Y = [Y y]
            push!(ORD,k)
            push!(T,t)
            
            t = t+dt
            y0 = y
        else
            break
        end
        =#
    end
    
    return (T,Y,ORD)
end
