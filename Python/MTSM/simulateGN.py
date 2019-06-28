# MTSM implementation -- Python
# Solves linear systems defined as y' = Ay + b

### Imports ### 
import argparse
import numpy as np
import scipy.sparse as sp

#from matplotlib.pyplot import figure, show
# https://matplotlib.org/examples/pylab_examples/spy_demos.html
# fig = figure()
# fig = figure()
# ax1 = fig.add_subplot(221)
# ax2 = fig.add_subplot(222)
# ax3 = fig.add_subplot(223)
# ax4 = fig.add_subplot(224)

# x = numpy.random.randn(20, 20)
# x[5] = 0.
# x[:, 12] = 0.

# ax1.spy(x, markersize=5)
# ax2.spy(x, precision=0.1, markersize=5)

# ax3.spy(x)
# ax4.spy(x, precision=0.1)

# show()

# ----------------------------------
# Gets Parameters
def GetParams():
    parser = argparse.ArgumentParser()
    
    parser.add_argument('-N', action='store', dest='N', default=10, type=int,
                        help='Number of segments')
    
    parser.add_argument('-R1', action='store', dest='R1', default=100, type=float,
                        help='Value of R1 [Ohm]')
    
    parser.add_argument('-R2', action='store', dest='R2', default=100, type=float,
                        help='Value of R2 [Ohm]')

    parser.add_argument('-L', action='store', dest='L', default=1e-8, type=float,
                        help='Inductance [H]')

    parser.add_argument('-C', action='store', dest='C', default=1e-12, type=float,
                        help='Capacitance [F]')

    parser.add_argument('-om', action='store', dest='om', default=3e9, type=float,
                        help='Omega [rad/s]')

    parser.add_argument('-maxord', action='store', dest='maxord', default=64, type=int,
                        help='Maximum order of MTSM')
    
    parser.add_argument('-eps', action='store', dest='eps', default=1e-8, type=float,
                        help='Required accuracy')
    
    
    params = parser.parse_args()
    #print 'N     =', params.N
    #print 'R1   =', params.R1
    #print 'R2   =', params.R2
    return params

# ----------------------------------
# Generate inputs for the Telegraph Line problem
def GenerateInputsTelegraphLine(N, L, C, om, R1, R2):
	# matrix A11
	data_A11=np.array([[N-1,N-1]])
	rows = data_A11[:,0]
	cols = data_A11[:,1] 

	data_A11 = np.ones(rows.shape, dtype=float) * (-1/(R2*C))
	A_11 = sp.csr_matrix((data_A11, (rows, cols)), shape=(N,N))#.todense()
	print "Matrix A11"
	print A_11
	
	# -------------------
	# matrix A12
	diagonals = np.zeros((2, N))   # 2 diagonals
	diagonals[0,:] = np.ones(N, dtype=float) * (1/C)
	diagonals[1,:] = np.ones(N, dtype=float) * (-1/C)

	A_12 = sp.spdiags(diagonals, [0,1], N, N+2, format='csr')#.todense()
	print "Matrix A12 (1/C, -1/C)"
	print A_12
	#A12(1:nRLC, 1:nRLC) = spdiags([e*1/C e*-1/C], [0 1], nRLC, nRLC);
	
	# -------------------
	# matrix A21
	diagonals = np.zeros((2, N))   # 2 diagonals
	diagonals[0,:] = np.ones(N, dtype=float) * (1/L)
	diagonals[1,:] = np.ones(N, dtype=float) * (-1/L)

	A_21 = sp.spdiags(diagonals, [-1,0], N+2, N, format='csr')#.todense()
	print "Matrix A21 (1/L, -1/L)"
	print A_21
	#A21(1:nRLC, 1:nRLC) = spdiags([e*1/L e*-1/L], [-1 0], nRLC, nRLC);
	
	# -------------------
	# matrix A22
	#data_A22=np.array([[0,0], [1,N+1], [N+1,N+2], [N+2,N+1]])
	data_A22=np.array([[0,0], [0,N], [N,N+1], [N+1,N]])
	rows = data_A22[:,0]
	cols = data_A22[:,1]
	data_A22 = np.ones(rows.shape, dtype=float) * ([-R1/L, 1/L, om, -om])


	A_22 = sp.csr_matrix((data_A22, (rows, cols)), shape=(N+2,N+2))#.todense()
	print "Matrix A22"
	print A_22

	# -------------------
	# assembly matrix A
	A = sp.bmat([[A_11, A_12], [A_21, A_22]])

	print "Matrix A"
	print A

	# -------------------
	# vector ic
	ic = np.array(np.zeros(N*2+2, dtype=float))
	ic[N*2+1]=1.0
	print "Vector ic"
	print ic

	# -------------------
	# vector b
	b = np.array(np.zeros(N*2+2, dtype=float))
	#b[N*2+1]=1.0
	print "Vector b"
	print b


	return A, ic, b

# ----------------------------------
# Explicit MTSM
def explicitTaylor(dt, tmax, eps, maxord, N, A, ic, b):
	# stoping rule: maximum number DY terms are smaller then eps
	stopping = 3
	# last step tolerance
	ls_tol = dt/10

	i = 0;
	t = []
	#Y = []
	Y = np.array(np.zeros((maxord,2*N+2),dtype=float))
	ORD = []
	ADy = []

	#t[i] = 0.0
	t.append(0.0)
	Y[i,:] = ic
	print "Y[i,:]"
	print Y[i,:]
	#Y.append(ic)
	#Y(i, mslice[:]).lvalue = ic

	i = 1;
	t.append(t[i-1] + dt)
	#t[i] = t[i-1] + dt

	AT = np.transpose(A)

	# main cycle
	while t[i-1]+ls_tol < tmax:
		DY = np.array(np.zeros((maxord,2*N+2),dtype=float))
		print DY
		Y[i,:] = Y[i-1,:]; 
		AY = A.dot(Y[i,:])
		print AY
		DY[1,:] = dt*(AY+b)
		print "DY"
		print DY
		Y[i,:] = Y[i,:] + DY[1,:]
		k = 1

		one_arr = np.ones((3), dtype=int)*1e10 #+ np.zeros(maxord+stopping)
		zero_arr = np.zeros(maxord+stopping)
		maxDY = np.concatenate((one_arr, zero_arr), axis=0)
		print "maxDY"
		print maxDY
		k = 2

		# compute taylor series terms
		while np.linalg.norm(maxDY[k-1:stopping+k-2]) > eps:
		#while sp.spatial.distance.euclidean(maxDY[k-1:stopping+k-2]) > eps:
			ADy = DY[k-1,:]*AT      
			DY[k,:] = (dt/k)*ADy

			Y[i,:] = Y[i,:]+DY[k,:]
			print 'Y[i,:] = Y[i,:]+DY[k,:]'
			print Y[i,:]

			maxDY[k+stopping-1] = np.absolute(np.max(DY[k,:]))
			k = k+1;

			if k > maxord:
				pass
				#break
			# end while inner

			# maxord reached - halving dt
			if k > maxord:
				dt = dt/2; 
			else:
				#DY_all{i-1}=DY;
				#ORD[i] = k;
				ORD.append(k)
				i = i+1;
				#t[i] = t[i-1]+dt;
				t.append(t[i-1]+dt)
        
    		# end while outer    
	return Y, t

# -------------------------------
# plot graphs
# def plotResults(T, Y, N):
# 	fig = figure()
# 	# MTSM - UC1
# 	plot(T,Y[:,1], '-r', 'LineWidth',1.5)
# 	#grid on;
# 	legend('UC1')
# 	#TITLE=sprintf('MATLAB Explicit Taylor solution');
# 	title('MATLAB Explicit Taylor solution')

# 	hold(True) 

# 	# MTSM - UC100
# 	plot(T,Y[:,N], '--b','LineWidth',1.5)
# 	#grid on;
# 	legend('UC1' ,['UC',num2str(nRLC,'%u')])
# 	#TITLE=sprintf('MATLAB Explicit Taylor solution');
# 	title('MATLAB Explicit Taylor solution')

# -------------------------------
# Main program
def main():
	par = GetParams()
	print 'main: GetParams()'
	print 'N     =', par.N
	print 'R1   =', par.R1
	print 'R2   =', par.R2
	print 'L     =', par.L
	print 'C   =', par.C
	print 'om   =', par.om
	print 'maxord   =', par.maxord
	print 'eps   =', par.eps
	

	#calculate the parameters - dt, tmax
	dt = np.sqrt(par.L*par.C);  		# size of the integration step

	tdelay = par.N*np.sqrt(par.L*par.C); 	# total propagation delay of the line 

	# set the maximum simulation time
	if par.R1 == 100 and par.R2 == 100:
		#adjusted line
		tmax = 2*tdelay;#2e-8;
	elif par.R1 == 100 and par.R2 == 1e12:
		#amplified bounced
		tmax = 2*tdelay;#2e-8; # for voltageType = sinus, for puslse tmax = 4e-8;

	print 'dt   =', dt
	print 'tdelay   =', tdelay
	print 'tmax   =', tmax

	# -----------------------
	# generate inputs
	print 'main: GenerateInputsTelegraphLine()'
	A, ic, b = GenerateInputsTelegraphLine(par.N, par.L, par.C, par.om, par.R1, par.R2)
	print 'main: GenerateInputsTelegraphLine() - final vals'
	print 'A     =', A
	print 'ic   =', ic
	print 'b   =', b

	# -----------------------
	# solve the problem using MTSM
	print 'main: explicitTaylor()'
	Y, t = explicitTaylor(dt, tmax, par.eps, par.maxord, par.N, A, ic, b)
	print 'main: explicitTaylor() - results'
	#for y,t  in zip(Y,t):
	#	print f'{t}:{y}'
	print Y[:,par.N*2+1]
	print t


# -------------------------------
# MAIN
# -------------------------------
main()