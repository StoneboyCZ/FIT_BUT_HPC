'''
dy(1) = y(3);
dy(2) = y(4);
dy(3) = -y(11);
dy(4) = -y(12);    
dy(5) = 3*y(9) + 3*y(10);
dy(6) = y(13) + y(14);

y(7) = y(1)*y(3)
y(8) = y(2)*y(4)
y(9) = y(6)*y(7)
y(10) = y(6)*y(8)
y(11) = y(1)/y(5)
y(12) = y(2)/y(5)
y(13) = y(7)/y(6)
y(14) = y(8)/y(6)
'''

from math import sqrt

h = 0.1
tmax = 1
maxord = 64
eps = 1e-12
e = 0.75

# results
y1 = []
y2 = []
y3 = []
y4 = []
y5 = []
y6 = []
y7 = []
y8 = []
y9 = []
y10 = []
y11 = []
y12 = []
y13 = []
y14 = []

y01 = 1-e
y02 = 0
y03 = 0
y04 = sqrt((1+e)/(1-e))
y05 = sqrt((y01*y01 + y02*y02)*(y01*y01 + y02*y02)*(y01*y01 + y02*y02))
y06 = sqrt((y01*y01 + y02*y02))
y07 = y01*y03
y08 = y02*y04
y09 = y06*y07
y010 = y06*y08
y011 = y01/y05
y012 = y02/y05
y013 = y07/y06
y014 = y08/y06

# initial conditions
y1.append(y01)
y2.append(y02)
y3.append(y03)
y4.append(y04)
y5.append(y05)
y6.append(y06)
y7.append(y07)
y8.append(y08)
y9.append(y09)
y10.append(y010)
y11.append(y011)
y12.append(y012)
y13.append(y013)
y14.append(y014)

steps = round(tmax/h)

#print(f'i {0}: y1 {y1[-1]} \t y2 {y2[-1]} \t y3 {y3[-1]} \t y4 {y4[-1]} \t y5 {y5[-1]}')
t = 0+h
for i in range(1,steps+1):
    k = 0

    y1.append(y1[-1])
    y2.append(y2[-1])
    y3.append(y3[-1])
    y4.append(y4[-1])
    y5.append(y5[-1])
    y6.append(y6[-1])

    y7.append(y1[-1]*y3[-1])
    y8.append(y2[-1]*y4[-1])
    y9.append(y6[-1]*y7[-1])
    y10.append(y6[-1]*y8[-1])
    
    y11.append(y1[-1]/y5[-1])
    y12.append(y2[-1]/y5[-1])
    y13.append(y7[-1]/y6[-1])
    y14.append(y8[-1]/y6[-1])
    
    DY1 = []
    DY2 = []
    DY3 = []
    DY4 = []
    DY5 = []
    DY6 = []
    DY7 = []
    DY8 = []
    DY9 = []
    DY10 = []
    DY11 = []
    DY12 = []
    DY13 = []
    DY14 = []
    
    # set initial conditions
    DY1.append(y1[-1])
    DY2.append(y2[-1])
    DY3.append(y3[-1])
    DY4.append(y4[-1])
    DY5.append(y5[-1])
    DY6.append(y6[-1])
    DY7.append(y7[-1])
    DY8.append(y8[-1])
    DY9.append(y9[-1])
    DY10.append(y10[-1])
    DY11.append(y11[-1])
    DY12.append(y12[-1])
    DY13.append(y13[-1])
    DY14.append(y14[-1])
    

    #print(f'k {k}: DY1 {DY1[-1]} \t DY2 {DY2[-1]} \t DY3 {DY3[-1]} \t DY4 {DY4[-1]} ')

    #Y0 = DY2[0]/DY4[0]

    k = 1
    
    DY1.append(h*(1*DY3[0]))
    DY2.append(h*(1*DY4[0]))
    DY3.append(h*(-1*DY11[0]))
    DY4.append(h*(-1*DY12[0]))
    DY5.append(h*(3*DY9[0] + 3*DY10[0]))
    DY6.append(h*(1*DY13[0] + 1*DY14[0]))

    DY7.append(DY1[k]*DY3[0] + DY1[0]*DY3[k])
    DY8.append(DY2[k]*DY4[0] + DY2[0]*DY4[k])
    DY9.append(DY6[k]*DY7[0] + DY6[0]*DY7[k])
    DY10.append(DY6[k]*DY8[0] + DY6[0]*DY8[k])
    
    DY11.append((1/DY5[0])*(DY1[k] - DY11[0]*DY5[k]))
    DY12.append((1/DY5[0])*(DY2[k] - DY12[0]*DY5[k]))
    DY13.append((1/DY6[0])*(DY7[k] - DY13[0]*DY6[k]))
    DY14.append((1/DY6[0])*(DY8[k] - DY14[0]*DY6[k]))
    
    y1[-1] = y1[-1] + DY1[-1]
    y2[-1] = y2[-1] + DY2[-1]
    y3[-1] = y3[-1] + DY3[-1]
    y4[-1] = y4[-1] + DY4[-1]
    y5[-1] = y5[-1] + DY5[-1]
    y6[-1] = y6[-1] + DY6[-1]
    y7[-1] = y7[-1] + DY7[-1]
    y8[-1] = y8[-1] + DY8[-1]
    y9[-1] = y9[-1] + DY9[-1]
    y10[-1] = y10[-1] + DY10[-1]
    y11[-1] = y11[-1] + DY11[-1]
    y12[-1] = y12[-1] + DY12[-1]
    y13[-1] = y13[-1] + DY13[-1]
    y14[-1] = y14[-1] + DY14[-1]
    
    #print(f'k {k}: DY1 {DY1[-1]} \t DY2 {DY2[-1]} \t DY3 {DY3[-1]} \t DY4 {DY4[-1]} ')

    #while (abs(DY1[-1]) > eps) and (abs(DY2[-1]) > eps) and (abs(DY3[-1]) > eps) and (abs(DY4[-1]) > eps):
    while k < 15:
        k = k+1
        print(f'k: {k}')
        
        DY1.append((h/k)*(1*DY3[k-1]))
        DY2.append((h/k)*(1*DY4[k-1]))
        DY3.append((h/k)*(-1*DY11[k-1]))
        DY4.append((h/k)*(-1*DY12[k-1]))
        DY5.append((h/k)*(3*DY9[k-1] + 3*DY10[k-1]))
        DY6.append((h/k)*(1*DY13[k-1] + 1*DY14[k-1]))
            
        # multiplications       
        l = k
        m = 0
        sum = 0
        for r in range(0,k+1):
            #print(f'{r}: m {m} l {l}')
            sum = sum + (DY1[l]*DY3[m])
            #print(f'{mixed}')
            l = l-1
            m = m+1
        
        DY7.append(sum)


        l = k
        m = 0
        sum = 0
        for r in range(0,k+1):
            #print(f'{r}: m {m} l {l}')
            sum = sum + (DY2[l]*DY4[m])
            #print(f'{mixed}')
            l = l-1
            m = m+1
        
        DY8.append(sum)

        l = k
        m = 0
        sum = 0
        for r in range(0,k+1):
            #print(f'{r}: m {m} l {l}')
            sum = sum + (DY6[l]*DY7[m])
            #print(f'{mixed}')
            l = l-1
            m = m+1
        
        DY9.append(sum)


        l = k
        m = 0
        sum = 0
        for r in range(0,k+1):
            #print(f'{r}: m {m} l {l}')
            sum = sum + (DY6[l]*DY8[m])
            #print(f'{mixed}')
            l = l-1
            m = m+1
        
        DY10.append(sum)

        # divisions
        l = 1
        m = k-1
        mixed = 0
        for r in range(1,k+1):
            #print(f'{r}: m {m} l {l}')
            mixed = mixed + (DY11[m]*DY5[l])
            #print(f'{mixed}')
            l = l+1
            m = m-1
        
        DY11.append((1/DY5[0])*(DY1[k]-mixed))

        l = 1
        m = k-1
        mixed = 0
        for r in range(1,k+1):
            #print(f'{r}: m {m} l {l}')
            mixed = mixed + (DY12[m]*DY5[l])
            #print(f'{mixed}')
            l = l+1
            m = m-1
        
        DY12.append((1/DY5[0])*(DY2[k]-mixed))


        l = 1
        m = k-1
        mixed = 0
        for r in range(1,k+1):
            #print(f'{r}: m {m} l {l}')
            mixed = mixed + (DY13[m]*DY6[l])
            #print(f'{mixed}')
            l = l+1
            m = m-1
        
        DY13.append((1/DY6[0])*(DY7[k]-mixed))

        l = 1
        m = k-1
        mixed = 0
        for r in range(1,k+1):
            #print(f'{r}: m {m} l {l}')
            mixed = mixed + (DY14[m]*DY6[l])
            #print(f'{mixed}')
            l = l+1
            m = m-1
        
        DY14.append((1/DY6[0])*(DY8[k]-mixed))

    y1[-1] = y1[-1] + DY1[-1]
    y2[-1] = y2[-1] + DY2[-1]
    y3[-1] = y3[-1] + DY3[-1]
    y4[-1] = y4[-1] + DY4[-1]
    y5[-1] = y5[-1] + DY5[-1]
    y6[-1] = y6[-1] + DY6[-1]
    y7[-1] = y7[-1] + DY7[-1]
    y8[-1] = y8[-1] + DY8[-1]
    y9[-1] = y9[-1] + DY9[-1]
    y10[-1] = y10[-1] + DY10[-1]
    y11[-1] = y11[-1] + DY11[-1]
    y12[-1] = y12[-1] + DY12[-1]
    y13[-1] = y13[-1] + DY13[-1]
    y14[-1] = y14[-1] + DY14[-1]
    
    print(f'i {i} t{t}: y1 {y1[-1]} \t y2 {y2[-1]} ')    
    t = t+h    
    
print(y1)