function mtsm_nonlinear(A,b,B2,B3,B4,B5,ij,ijk,ijkl,ijklm,y0,t0,tmax,dt,eps,maxORD)
    T = Float64[t0]
    Y = y0
    ORD = Int8[0]
    t = t0+dt
    nStopping = 3
    ne = length(y0) 

    # prepare the indexes of multiplication terms
    if ij === nothing
        B2y = zeros(BigFloat,ne,1)
    else
        i2 = ij[:, 1]
        j2 = ij[:, 2]
    end

    if ijk === nothing
        B3y = zeros(BigFloat,ne,1)
    else
        i3 = ijk[:, 1]
        j3 = ijk[:, 2]
        k3 = ijk[:, 3]
    end

    if ijkl === nothing
        B4y = zeros(BigFloat,ne,1) 
    else
        i4 = ijkl[:, 1]
        j4 = ijkl[:, 2]
        k4 = ijkl[:, 3]
        l4 = ijkl[:, 4]
    end
    
    if ijklm === nothing
        B5y = zeros(BigFloat,ne,1) 
    else
        i5 = ijklm[:, 1]
        j5 = ijklm[:, 2]
        k5 = ijklm[:, 3]
        l5 = ijklm[:, 4]
        m5 = ijklm[:, 5]
    end

    while true
        y = y0        
        
        DY = zeros(BigFloat, ne, 5)
        display(DY)
        DY[:,1] = y0
        display(DY)

        stopping = ones(BigFloat,nStopping)*1e10
        stopping[1] = abs(reduce((x,y) -> max.(x,y), DY[:,1]))
        
        Ay = (A*DY[:,1])+b
        display(Ay)

        #= Non-linear indexes =#
        ij1 = [2,1]
        ij2 = [1,2]

        ijk2 = ij1
        ijk3 = ij2

        ijkl2 = ij1
        ijkl3 = ij2

        ijklm2 = [2,1,1]
        ijklm3 = [1,2,1]
        ijklm4 = [1,1,2]

        if ij !== nothing
            B2y=B2*(DY[i2,1] .* DY[j2,1])            
        end

        DY[:, 2] = dt*(Ay+B2y)
        display(DY)
        y = y + DY[:,2]
        display(y)
        
        stopping[2] = abs(reduce((x,y) -> max.(x,y), DY[:,2]))

        k = Int8(3)
        s = Int8(3)
        while norm(abs.(stopping)) > eps
            # linear term
            Ay=A*DY[:,k-1] 
            
            if ij !== nothing
                B2y=B2*sum(DY[i2,ij1].*DY[j2,ij2],dims=2)
                B2y = reshape(B2y, length(B2y))
                display(B2y)
                display(length(B2y))
            end

            DY[:,k] = (dt/k)*(Ay+B2y)
            display(k)
            display(dt/k)
            display(DY[:,k])
            y = y + DY[:,k]
            stopping[s] = abs(reduce((x,y) -> max.(x,y), DY[:,k]))

            if s == nStopping
                s = 1
            else
                s = s+1    
            end
           
            ij1 = cat(k,ij1;dims=1);
            ij2 = cat(ij2,k; dims=1);
            display(ij1) 
            display(ij2)
        
            k = k+1
        end

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
    end

    return (T,Y,ORD)
end