maxORD = 3
multORD = 2
coefs = []

# init
init=[0] * multORD

# get derivatives
DY = []
for o in range(0,maxORD+1,1):
    if o == 0:
        #coefs.append(init)
        DY.append(init)
    else:
        for i,n in enumerate(DY[o-1]):
            coefs = []

            # current digit copies the previous result
            coefs.append(DY[o-1].copy())

            # change the appropriate value by incrementing 1 (simulate derivative)
            coefs[i][i] = coefs[i][i]+1 

        DY.append(coefs)

print(DY)    