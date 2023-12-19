''' PROBLEM
y' = (sin(t)*cos(t))/exp(t); y(0) = 1; TMAX = 1; h = 0.1


ODE
y1' = y6 y1(0) = 1
y2' = -y3 y2(0) = 1 cos(t)
y3' = y2 y3(0) = 0 sin(t)
y4' = y4 y4(0) = 1 e^t

DAE
y5 = y2*y4;
y6 = y5/y4;
'''

h = 0.1
tmax = 1
maxord = 64
eps = 1e-12

# results
y1 = []
y2 = []
y3 = []
y4 = []
y5 = []
y6 = []

# initial conditions
y1.append(1)
y2.append(1)
y3.append(0)
y4.append(1)
y5.append(y2[-1]*y3[-1])
y6.append(y5[-1]/y4[-1])

steps = round(tmax/h)

#print(f'i {0}: y1 {y1[-1]} \t y2 {y2[-1]} \t y3 {y3[-1]} \t y4 {y4[-1]} \t y5 {y5[-1]}')
t = 0+h
for i in range(1,steps+1):
    k = 0

    y1.append(y1[-1])
    y2.append(y2[-1])
    y3.append(y3[-1])
    y4.append(y4[-1])
    y5.append(y2[-1]*y3[-1])
    y6.append(y5[-1]/y4[-1])

    DY1 = []
    DY2 = []
    DY3 = []
    DY4 = []
    DY5 = []
    DY6 = []

    # set initial conditions
    DY1.append(y1[-1])
    DY2.append(y2[-1])
    DY3.append(y3[-1])
    DY4.append(y4[-1])
    DY5.append(y5[-1])
    DY6.append(y6[-1])
    

    #print(f'k {k}: DY1 {DY1[-1]} \t DY2 {DY2[-1]} \t DY3 {DY3[-1]} \t DY4 {DY4[-1]} ')

    #Y0 = DY2[0]/DY4[0]

    k = 1
    DY1.append(h*(1*DY6[0]))
    DY2.append(h*(-1*DY3[0]))
    DY3.append(h*(1*DY2[0]))
    DY4.append(h*(1*DY4[0]))

    DY5.append(DY2[k]*DY3[0] + DY2[0]*DY3[k])
    
    DY6.append((1/DY4[0])*(DY5[k] - DY6[0]*DY4[k]))

    y1[-1] = y1[-1] + DY1[-1]
    y2[-1] = y2[-1] + DY2[-1]
    y3[-1] = y3[-1] + DY3[-1]
    y4[-1] = y4[-1] + DY4[-1]
    y5[-1] = y5[-1] + DY5[-1]
    y6[-1] = y6[-1] + DY6[-1]


    #print(f'k {k}: DY1 {DY1[-1]} \t DY2 {DY2[-1]} \t DY3 {DY3[-1]} \t DY4 {DY4[-1]} ')

    #while (abs(DY1[-1]) > eps) and (abs(DY2[-1]) > eps) and (abs(DY3[-1]) > eps) and (abs(DY4[-1]) > eps):
    while k < 40:
        k = k+1
        print(f'k: {k}')
        DY1.append((h/k)*(1*DY6[k-1]))
        DY2.append((h/k)*(-1*DY3[k-1]))
        DY3.append((h/k)*(1*DY2[k-1]))
        DY4.append((h/k)*(1*DY4[k-1]))
        
        l = k
        m = 0
        sum1 = 0
        for r in range(0,k+1):
            #print(f'{r}: m {m} l {l}')
            sum1 = sum1 + (DY2[l]*DY3[m])
            #print(f'{mixed}')
            l = l-1
            m = m+1
        
        DY5.append(sum1)
        
        l = 1
        m = k-1
        mixed = 0
        for r in range(1,k+1):
            #print(f'{r}: m {m} l {l}')
            mixed = mixed + (DY6[m]*DY4[l])
            #print(f'{mixed}')
            l = l+1
            m = m-1
        
        DY6.append((1/DY4[0])*(DY5[k]-mixed))


        y1[-1] = y1[-1] + DY1[-1]
        y2[-1] = y2[-1] + DY2[-1]
        y3[-1] = y3[-1] + DY3[-1]
        y4[-1] = y4[-1] + DY4[-1]   
        y5[-1] = y5[-1] + DY5[-1]
        y6[-1] = y6[-1] + DY6[-1]
        
    
    print(f'i {i} t{t}: y1 {y1[-1]} \t y2 {y2[-1]} \t y3 {y3[-1]} \t y4 {y4[-1]} \t {y5[-1]}  ')    
    t = t+h    
    
print(y1)