''' PROBLEM
y' = cos(t)/exp(t); y(0) = 1; TMAX = 1; h = 0.1

AUX
y1' = y2/y4 y1(0) = 1
y2' = -y3 y2(0) = 1 cos(t)
y3' = y2 y3(0) = 0 sin(t)
y4' = y4 y4(0) = 1 e^t
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

# initial conditions
y1.append(1)
y2.append(1)
y3.append(0)
y4.append(1)

steps = round(tmax/h)

print(f'i {0}: y1 {y1[-1]} \t y2 {y2[-1]} \t y3 {y3[-1]} \t y4 {y4[-1]} ')

for i in range(1,steps+1):
    k = 0
    y1.append(y1[-1])
    y2.append(y2[-1])
    y3.append(y3[-1])
    y4.append(y4[-1])


    DY1 = []
    DY2 = []
    DY3 = []
    DY4 = []

    # set initial conditions
    DY1.append(y1[-1])
    DY2.append(y2[-1])
    DY3.append(y3[-1])
    DY4.append(y4[-1])
    
    #print(f'k {k}: DY1 {DY1[-1]} \t DY2 {DY2[-1]} \t DY3 {DY3[-1]} \t DY4 {DY4[-1]} ')

    Y0 = DY2[0]/DY4[0]

    k = 1
    DY1.append(h*Y0)
    DY2.append(h*(-1*DY3[0]))
    DY3.append(h*(1*DY2[0]))
    DY4.append(h*(1*DY4[0]))

    y1[-1] = y1[-1] + DY1[-1]
    y2[-1] = y2[-1] + DY2[-1]
    y3[-1] = y3[-1] + DY3[-1]
    y4[-1] = y4[-1] + DY4[-1]

    #print(f'k {k}: DY1 {DY1[-1]} \t DY2 {DY2[-1]} \t DY3 {DY3[-1]} \t DY4 {DY4[-1]} ')

    #while (abs(DY1[-1]) > eps) and (abs(DY2[-1]) > eps) and (abs(DY3[-1]) > eps) and (abs(DY4[-1]) > eps):
    while k < 15:
        k = k+1
        print(f'k: {k}')
        DY2.append((h/k)*(-1*DY3[k-1]))
        DY3.append((h/k)*(1*DY2[k-1]))
        DY4.append((h/k)*(1*DY4[k-1]))

        l = 0
        m = k-1
        mixed = 0
        for r in range(1,k-1):
            #print(f'{r}: y1 {l} y4 {m}')
            mixed = mixed + (DY1[l]*DY4[m])
            #print(f'{mixed}')
            l = l+1
            m = m-1
        
        DY1.append((h/(k*DY4[0]))*(DY2[k-1]-mixed-Y0*DY4[k-1]))

        y1[-1] = y1[-1] + DY1[-1]
        y2[-1] = y2[-1] + DY2[-1]
        y3[-1] = y3[-1] + DY3[-1]
        y4[-1] = y4[-1] + DY4[-1]
    
        
    print(f'i {i}: y1 {y1[-1]} \t y2 {y2[-1]} \t y3 {y3[-1]} \t y4 {y4[-1]} ')    


print(y1)