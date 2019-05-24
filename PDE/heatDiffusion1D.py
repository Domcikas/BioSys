import numpy as np
import matplotlib.pyplot as plt
import math
import random
from math import sin, cos, tan, pi
length = 1
numberOfNodes = 50
leftBoundary = 0
rightBoundary = 0
finalTime = 60

D0 = np.ones(numberOfNodes)*0.1
energyS = np.zeros(numberOfNodes)
def forX(x):
	#return 2-1.5*x+sin(math.pi*x)
	if x >= 0.45 and x <= 0.55:
		return 100	
	return 0
	#return 2*math.exp(-(x-0.5)**2/4)
	#return 4*x-4*x*x

def calcTimeStep(dx, D):
	return 0.5*dx**2/D

def calcTimeStep2D(dx,dy,D):
	return ((dx**2)*(dy**2))/(2*D*((dx**2)+(dy**2)))

def calcTimeStepPolar(r, nS, n, D):
	dr = r/n
	#trapezoid points
	#v1
	x1 = r
	y1 = 0
	#v2
	x2 = r-dr
	y2 = 0
	#v3
	x3 = (r-dr)*cos(2*np.pi*1/nS)
	y3 = (r-dr)*sin(2*np.pi*1/nS)
	#v4
	x4 = r*cos(2*np.pi*1/nS)
	y4 = r*sin(2*np.pi*1/nS)
	#distances
	v2v3 = np.sqrt((x2-x3)**2+(y2-y3)**2)
	v1v4 = np.sqrt((x1-x4)**2+(y1-y4)**2)
	v1v2 = np.sqrt((x1-x2)**2+(y1-y2)**2)
	h = np.sqrt(v1v2**2-((v1v4-v2v3)/2)**2)
	diameter = np.sqrt(h**2+(v2v3+(v1v4-v2v3)/2)**2)
	return 0.001*diameter/D

		

def internalHeatSource(x, sourceExists):
	return 2*sourceExists

def distributeD(x, l):
	return 0.1
	#if x < (l*0.9)/2:
	#	return 0.001
	#elif x > (l*1.1)/2:
	#	return 0.01
	#else:
	#	return 0.001+x-(l*0.9)/2

def integrateCartesian(matrix, dx, dy):
	summ = 0
	for i in matrix:
		for j in i:
			summ += j*(dx*dy)
	return summ

def integratePolar(matrix, dr, nphi):
	summ = 0
	summ += matrix[0][0]*np.pi*(dr/2)**2
	for i in matrix:
		for j in range(1,len(i)-1):
			S = (np.pi*(j*dr+dr/2)**2-np.pi*(j*dr-dr/2)**2)/nphi
			summ += i[j]*S
		S = (np.pi*(j*dr)**2-np.pi*(j*dr-dr/2)**2)/nphi
		summ += i[-1]*S
	return summ

def calculateS(index, dr, nphi, nr):
	if index == 0:
		return np.pi*(dr/2)**2
	if index == nr:
		return (np.pi*(index*dr)**2-np.pi*(index*dr-dr/2)**2)/nphi
	else:
		return (np.pi*(index*dr+dr/2)**2-np.pi*(index*dr-dr/2)**2)/nphi

def finiteVolumeMethod(leftBoundary, rightBoundary, finalTime, length = 1, numberOfNodes = 10, thermalCoeffs = [0.1] * numberOfNodes, energySource = [0] * numberOfNodes):
	#Getting the parameters
	L = length
	n = numberOfNodes
	T = np.ones(numberOfNodes)
	T1S = leftBoundary
	T2S = rightBoundary
	dx = float(L)/float(n)
	x = np.linspace(dx/2, L-dx/2, n)
	heatSource = energySource
	
	#Applying energy source function	
	for i in range(0, len(x)):
		heatSource[i] = internalHeatSource(x[i], energySource[i])

	#Applying initial condition to X
	for i in range(0, len(x)):
		T[i] = forX(x[i])

	t_final = finalTime
	D = thermalCoeffs
	timeSteps = []

	#Applying thermal constant D to each x node
	for i in range(0,len(x)):
		D[i] = distributeD(x[i], L)

	#Calculating suitable dt
	for d in D:
		timeSteps.append(calcTimeStep(dx, float(d)))	
	dt = min(timeSteps)
	
	#Finding minimum and maximum values of temperature
	#They're either in initVector or one boundaries
	tMin = min(T)
	tMin = tMin if tMin < leftBoundary else leftBoundary
	tMin = tMin if tMin < rightBoundary else rightBoundary

	tMax = max(T)
	tMax = tMax if tMax > leftBoundary else leftBoundary
	tMax = tMax if tMax > rightBoundary else rightBoundary
	t = np.arange(0, t_final, dt)
	dTdt = np.empty(n)
	temperature = []
	s = 0
	for tem in T:
		s+=tem*dx
	print s
	#Running the simulation
	for j in range(1,len(t)):
		#T1S = 2*sin(t[j])
		#T2S = 2*cos(t[j])
		for i in range(1,n-1):
			dTdt[i] = D[i]*(-(T[i]-T[i-1])/dx**2+(T[i+1]-T[i])/dx**2)+heatSource[i]
		#dTdt[0] = D[i]*(-(T[0]-T1S)/dx**2+(T[1]-T[0])/dx**2)+heatSource[0]
		#dTdt[n-1] = D[i]*(-(T[n-1]-T[n-2])/dx**2+(T2S-T[n-1])/dx**2)+heatSource[n-1]
		#Boundaries are insulated
		dTdt[0] = D[i]*(-(T[0]-T[1])/dx**2+(T[1]-T[0])/dx**2)+heatSource[0]
		dTdt[n-1] = D[i]*(-(T[n-1]-T[n-2])/dx**2+(T[n-2]-T[n-1])/dx**2)+heatSource[n-1]
		T += dTdt*dt
		temperature.append(T.tolist())
		a = 0
		for tem in T:
			a+=tem*dx
		if j == len(t)-1:
			print j
			print a
			print T
			return 0
		
	
	#Plotting the simulation to heatmap
	fig, ax = plt.subplots()
	c = ax.pcolormesh(x, t, temperature, cmap='RdBu_r', vmin=tMin, vmax=tMax)
	ax.set_title('1D heat diffusion')
	ax.axis([dx/2, L-dx/2, 0, t_final])
	plt.xlabel("x")
	plt.ylabel("t")
	fig.colorbar(c, ax=ax)
	
	plt.show()

def FTCSexample(leftBoundary, rightBoundary, finalTime, length = 1, numberOfNodes = 10, thermalCoeffs = [0.1] * numberOfNodes, energySource = [0] * numberOfNodes):
	#Getting the parameters
	L = length
	n = numberOfNodes
	T = np.ones(numberOfNodes)
	T1S = leftBoundary
	T2S = rightBoundary
	dx = float(L)/float(n)
	x = np.linspace(dx/2, L-dx/2, n)
	heatSource = energySource
	
	#Applying energy source function	
	for i in range(0, len(x)):
		heatSource[i] = internalHeatSource(x[i], energySource[i])

	#Applying initial condition to X
	for i in range(0, len(x)):
		T[i] = forX(x[i])

	t_final = finalTime
	D = thermalCoeffs
	timeSteps = []

	#Applying thermal constant D to each x node
	for i in range(0,len(x)):
		D[i] = distributeD(x[i], L)

	#Calculating suitable dt
	for d in D:
		timeSteps.append(calcTimeStep(dx, float(d)))	
	dt = min(timeSteps)
	
	#Finding minimum and maximum values of temperature
	#They're either in initVector or one boundaries
	tMin = min(T)

	tMax = max(T)
	t = np.arange(0, t_final, dt)
	dTdt = np.empty(n)
	temperature = []
	for ts in t:
		tMin = tMin if tMin < 2*sin(ts) else 2*sin(ts)
		tMin = tMin if tMin < 2*cos(ts) else 2*cos(ts)
		tMax = tMax if tMax > 2*sin(ts) else 2*sin(ts)
		tMax = tMax if tMax > 2*cos(ts) else 2*cos(ts)	
	#Running the simulation
	for j in range(1,len(t)):
		#Dynamic boundaries
		T1S = 2*sin(t[j])
		T2S = 2*cos(t[j])

		for i in range(1,n-1):
			dTdt[i] = D[i]*(-(T[i]-T[i-1])/dx**2+(T[i+1]-T[i])/dx**2)+heatSource[i]
		dTdt[0] = D[i]*(-(T[0]-T1S)/dx**2+(T[1]-T[0])/dx**2)+heatSource[0]
		dTdt[n-1] = D[i]*(-(T[n-1]-T[n-2])/dx**2+(T2S-T[n-1])/dx**2)+heatSource[n-1]
		T += dTdt*dt
		temperature.append(T.tolist())
	
	#Plotting the simulation to heatmap
	fig, ax = plt.subplots()
	c = ax.pcolormesh(x, t, temperature, cmap='RdBu_r', vmin=tMin, vmax=tMax)
	ax.set_title('1D heat diffusion')
	ax.axis([dx/2, L-dx/2, 0, t_final])
	plt.xlabel("x")
	plt.ylabel("t")
	fig.colorbar(c, ax=ax)
	
	plt.show()


def FTCSmethod2D(lengthX, lengthY, numberOfNodesX, numberOfNodesY, initMatrix, leftBoundaryX, leftBoundaryY, rightBoundaryX, rightBoundaryY, thermalConst, finalTime):
	#Getting the parameters
	Lx = lengthX
	Ly = lengthY
	nx = numberOfNodesX
	ny = numberOfNodesY
	T = initMatrix
	T[5][6] = 1200
	T1Sx = leftBoundaryX
	T1Sy = leftBoundaryY
	T2Sx = rightBoundaryX
	T2Sy = rightBoundaryY
	dx = float(Lx)/float(nx)
	dy = float(Ly)/float(ny)
	x = np.linspace(dx/2, Lx-dx/2, nx)
	y = np.linspace(dy/2, Ly-dy/2, ny)
	X,Y = np.meshgrid(x,y)
	t_final = finalTime
	alpha = thermalConst

	#Calculating suitable dt
	dt = calcTimeStep2D(dx, dy, float(alpha))
	
	
	#Finding minimum and maximum values of temperature
	#They're either in initMatrix or one boundaries
	tMin = initMatrix[0][0]
	for i in initMatrix:
		tMin = tMin if tMin < min(i) else min(i)
	tMin = tMin if tMin < leftBoundaryX else leftBoundaryX
	tMin = tMin if tMin < rightBoundaryX else rightBoundaryX
	tMin = tMin if tMin < leftBoundaryY else leftBoundaryY
	tMin = tMin if tMin < rightBoundaryY else rightBoundaryY

	tMax = initMatrix[0][0]
	for i in initMatrix:
		tMax = tMax if tMax > max(i) else max(i)
	tMax = tMax if tMax > leftBoundaryX else leftBoundaryX
	tMax = tMax if tMax > rightBoundaryX else rightBoundaryX
	tMax = tMax if tMax > leftBoundaryY else leftBoundaryY
	tMax = tMax if tMax > rightBoundaryY else rightBoundaryY

	t = np.arange(0, t_final, dt)
	dTdt = np.empty((nx, ny))
	dAdt = np.empty((nx, ny))
	integrals = [[1200],[0]]
	#Reaction parameters
	k_cat_A_T = 4.2
	KM_A_T = 0.052
	EM = 0.012
	A = np.zeros((numberOfNodesX, numberOfNodesY))
	#Running the simulation
	for j in range(1,len(t)):
		for i in range(1,nx-1):
			for k in range(1, ny-1):
				#Inner square
				#T
				dTdt[i][k] = alpha*((-(T[i][k]-T[i-1][k])/dx**2+(T[i+1][k]-T[i][k])/dx**2)+(-(T[i][k]-T[i][k-1])/dy**2+(T[i][k+1]-T[i][k])/dy**2))
				dTdt[i][k] += (-k_cat_A_T*EM*T[i][k])#/(KM_A_T+T[i][k])
				#A
				dAdt[i][k] = alpha*((-(A[i][k]-A[i-1][k])/dx**2+(A[i+1][k]-A[i][k])/dx**2)+(-(A[i][k]-A[i][k-1])/dy**2+(A[i][k+1]-A[i][k])/dy**2))
				dAdt[i][k] += (k_cat_A_T*EM*T[i][k])#/(KM_A_T+T[i][k])
			#Left and right borders
			#T
			#dTdt[i][0] = alpha*((-(T[i][0]-T[i-1][0])/dx**2+(T[i+1][0]-T[i][0])/dx**2)+(-(T[i][0]-leftBoundaryY)/dy**2+(T[i][1]-T[i][0])/dy**2))
			#insulated boundaries
			dTdt[i][0] = alpha*((-(T[i][0]-T[i-1][0])/dx**2+(T[i+1][0]-T[i][0])/dx**2)+(-(T[i][0]-T[i][1])/dy**2+(T[i][1]-T[i][0])/dy**2))
			dTdt[i][0] += (-k_cat_A_T*EM*T[i][0])#/(KM_A_T+T[i][0])
			#dTdt[i][ny-1] = alpha*((-(T[i][ny-1]-T[i-1][ny-1])/dx**2+(T[i+1][ny-1]-T[i][ny-1])/dx**2)+(-(T[i][ny-1]-T[i][ny-2])/dy**2+(rightBoundaryY-T[i][ny-1])/dy**2))
			#insulated boundaries
			dTdt[i][ny-1] = alpha*((-(T[i][ny-1]-T[i-1][ny-1])/dx**2+(T[i+1][ny-1]-T[i][ny-1])/dx**2)+(-(T[i][ny-1]-T[i][ny-2])/dy**2+(T[i][ny-2]-T[i][ny-1])/dy**2))
			dTdt[i][ny-1] += (-k_cat_A_T*EM*T[i][ny-1])#/(KM_A_T+T[i][ny-1])
			#A
			#dAdt[i][0] = alpha*((-(A[i][0]-A[i-1][0])/dx**2+(A[i+1][0]-A[i][0])/dx**2)+(-(A[i][0]-leftBoundaryY)/dy**2+(A[i][1]-A[i][0])/dy**2))
			#insulated boundaries
			dAdt[i][0] = alpha*((-(A[i][0]-A[i-1][0])/dx**2+(A[i+1][0]-A[i][0])/dx**2)+(-(A[i][0]-A[i][1])/dy**2+(A[i][1]-A[i][0])/dy**2))
			dAdt[i][0] += (k_cat_A_T*EM*T[i][0])#/(KM_A_T+T[i][0])
			#dAdt[i][ny-1] = alpha*((-(A[i][ny-1]-A[i-1][ny-1])/dx**2+(A[i+1][ny-1]-A[i][ny-1])/dx**2)+(-(A[i][ny-1]-A[i][ny-2])/dy**2+(rightBoundaryY-A[i][ny-1])/dy**2))
			#insulated boundaries
			dAdt[i][ny-1] = alpha*((-(A[i][ny-1]-A[i-1][ny-1])/dx**2+(A[i+1][ny-1]-A[i][ny-1])/dx**2)+(-(A[i][ny-1]-A[i][ny-2])/dy**2+(A[i][ny-2]-A[i][ny-1])/dy**2))
			dAdt[i][ny-1] += (k_cat_A_T*EM*T[i][ny-1])#/(KM_A_T+T[i][ny-1])

		for y in range(1, ny-1):
			#Top and bottom borders
			#T
			#dTdt[0][y] = alpha*((-(T[0][y]-leftBoundaryX)/dx**2+(T[1][y]-T[0][y])/dx**2)+(-(T[0][y]-T[0][y-1])/dy**2+(T[0][y+1]-T[0][y])/dy**2))
			#dTdt[nx-1][y] = alpha*((-(T[nx-1][y]-T[nx-2][y])/dx**2+(rightBoundaryX-T[nx-1][y])/dx**2)+(-(T[nx-1][y]-T[nx-1][y-1])/dy**2+(T[nx-1][y+1]-T[nx-1][y])/dy**2))
			#insulated boundaries
			dTdt[0][y] = alpha*((-(T[0][y]-T[1][y])/dx**2+(T[1][y]-T[0][y])/dx**2)+(-(T[0][y]-T[0][y-1])/dy**2+(T[0][y+1]-T[0][y])/dy**2))
			dTdt[nx-1][y] = alpha*((-(T[nx-1][y]-T[nx-2][y])/dx**2+(T[nx-1][y-1]-T[nx-1][y])/dx**2)+(-(T[nx-1][y]-T[nx-1][y-1])/dy**2+(T[nx-1][y+1]-T[nx-1][y])/dy**2))
			dTdt[0][y] += (-k_cat_A_T*EM*T[0][y])#/(KM_A_T+T[0][y])
			dTdt[nx-1][y] += (-k_cat_A_T*EM*T[nx-1][y])#/(KM_A_T+T[nx-1][y])
			#A
			#dAdt[0][y] = alpha*((-(A[0][y]-leftBoundaryX)/dx**2+(A[1][y]-A[0][y])/dx**2)+(-(A[0][y]-A[0][y-1])/dy**2+(A[0][y+1]-A[0][y])/dy**2))
			#dAdt[nx-1][y] = alpha*((-(A[nx-1][y]-A[nx-2][y])/dx**2+(rightBoundaryX-A[nx-1][y])/dx**2)+(-(A[nx-1][y]-A[nx-1][y-1])/dy**2+(A[nx-1][y+1]-A[nx-1][y])/dy**2))
			#insulated boundaries
			dAdt[0][y] = alpha*((-(A[0][y]-A[1][y])/dx**2+(A[1][y]-A[0][y])/dx**2)+(-(A[0][y]-A[0][y-1])/dy**2+(A[0][y+1]-A[0][y])/dy**2))
			dAdt[nx-1][y] = alpha*((-(A[nx-1][y]-A[nx-2][y])/dx**2+(A[nx-1][y-1]-A[nx-1][y])/dx**2)+(-(A[nx-1][y]-A[nx-1][y-1])/dy**2+(A[nx-1][y+1]-A[nx-1][y])/dy**2))
			dAdt[0][y] += (k_cat_A_T*EM*T[0][y])#/(KM_A_T+T[0][y])
			dAdt[nx-1][y] += (k_cat_A_T*EM*T[nx-1][y])#/(KM_A_T+T[nx-1][y])
			
			#Corners
			#T
			#dTdt[0][0] = alpha*((-(T[0][0]-leftBoundaryX)/dx**2+(T[1][0]-T[0][0])/dx**2)+(-(T[0][0]-leftBoundaryY)/dy**2+(T[0][1]-T[0][0])/dy**2))
			#dTdt[nx-1][0] = alpha*((-(T[nx-1][0]-T[nx-2][0])/dx**2+(rightBoundaryX-T[nx-1][0])/dx**2)+(-(T[nx-1][0]-leftBoundaryY)/dy**2+(T[nx-1][1]-T[nx-1][0])/dy**2))
			#dTdt[0][ny-1] = alpha*((-(T[0][ny-1]-rightBoundaryX)/dx**2+(T[1][ny-1]-T[0][ny-1])/dx**2)+(-(T[0][ny-1]-T[0][ny-2])/dy**2+(rightBoundaryY-T[0][ny-1])/dy**2))
			#dTdt[nx-1][ny-1] = alpha*((-(T[nx-1][ny-1]-T[nx-2][ny-1])/dx**2+(leftBoundaryX-T[nx-1][ny-1])/dx**2)+(-(T[nx-1][ny-1]-T[nx-1][ny-2])/dy**2+(leftBoundaryY-T[nx-1][ny-1])/dy**2))
			#insulated boundaries
			dTdt[0][0] = alpha*((-(T[0][0]-T[1][0])/dx**2+(T[1][0]-T[0][0])/dx**2)+(-(T[0][0]-T[0][1])/dy**2+(T[0][1]-T[0][0])/dy**2))
			dTdt[nx-1][0] = alpha*((-(T[nx-1][0]-T[nx-2][0])/dx**2+(T[nx-2][0]-T[nx-1][0])/dx**2)+(-(T[nx-1][0]-T[nx-1][1])/dy**2+(T[nx-1][1]-T[nx-1][0])/dy**2))
			dTdt[0][ny-1] = alpha*((-(T[0][ny-1]-T[1][ny-1])/dx**2+(T[1][ny-1]-T[0][ny-1])/dx**2)+(-(T[0][ny-1]-T[0][ny-2])/dy**2+(T[0][ny-2]-T[0][ny-1])/dy**2))
			dTdt[nx-1][ny-1] = alpha*((-(T[nx-1][ny-1]-T[nx-2][ny-1])/dx**2+(T[nx-2][ny-1]-T[nx-1][ny-1])/dx**2)+(-(T[nx-1][ny-1]-T[nx-1][ny-2])/dy**2+(T[nx-1][ny-2]-T[nx-1][ny-1])/dy**2))
			dTdt[0][0] += (-k_cat_A_T*EM*T[0][0])#/(KM_A_T+T[0][0])
			dTdt[nx-1][0] += (-k_cat_A_T*EM*T[nx-1][0])#/(KM_A_T+T[nx-1][0])
			dTdt[0][ny-1] += (-k_cat_A_T*EM*T[0][ny-1])#/(KM_A_T+T[0][ny-1])
			dTdt[nx-1][ny-1] += (-k_cat_A_T*EM*T[nx-1][ny-1])#/(KM_A_T+T[nx-1][ny-1])
			#A
			#dAdt[0][0] = alpha*((-(A[0][0]-leftBoundaryX)/dx**2+(A[1][0]-A[0][0])/dx**2)+(-(A[0][0]-leftBoundaryY)/dy**2+(A[0][1]-A[0][0])/dy**2))
			#dAdt[nx-1][0] = alpha*((-(A[nx-1][0]-A[nx-2][0])/dx**2+(rightBoundaryX-A[nx-1][0])/dx**2)+(-(A[nx-1][0]-leftBoundaryY)/dy**2+(A[nx-1][1]-A[nx-1][0])/dy**2))
			#dAdt[0][ny-1] = alpha*((-(A[0][ny-1]-rightBoundaryX)/dx**2+(A[1][ny-1]-A[0][ny-1])/dx**2)+(-(A[0][ny-1]-A[0][ny-2])/dy**2+(rightBoundaryY-A[0][ny-1])/dy**2))
			#dAdt[nx-1][ny-1] = alpha*((-(A[nx-1][ny-1]-A[nx-2][ny-1])/dx**2+(leftBoundaryX-A[nx-1][ny-1])/dx**2)+(-(A[nx-1][ny-1]-A[nx-1][ny-2])/dy**2+(leftBoundaryY-A[nx-1][ny-1])/dy**2))
			dAdt[0][0] = alpha*((-(A[0][0]-A[1][0])/dx**2+(A[1][0]-A[0][0])/dx**2)+(-(A[0][0]-A[0][1])/dy**2+(A[0][1]-A[0][0])/dy**2))
			dAdt[nx-1][0] = alpha*((-(A[nx-1][0]-A[nx-2][0])/dx**2+(A[nx-2][0]-A[nx-1][0])/dx**2)+(-(A[nx-1][0]-A[nx-1][1])/dy**2+(A[nx-1][1]-A[nx-1][0])/dy**2))
			dAdt[0][ny-1] = alpha*((-(A[0][ny-1]-A[1][ny-1])/dx**2+(A[1][ny-1]-A[0][ny-1])/dx**2)+(-(A[0][ny-1]-A[0][ny-2])/dy**2+(A[0][ny-2]-A[0][ny-1])/dy**2))
			dAdt[nx-1][ny-1] = alpha*((-(A[nx-1][ny-1]-A[nx-2][ny-1])/dx**2+(A[nx-2][ny-1]-A[nx-1][ny-1])/dx**2)+(-(A[nx-1][ny-1]-A[nx-1][ny-2])/dy**2+(A[nx-1][ny-2]-A[nx-1][ny-1])/dy**2))
			dAdt[0][0] += (k_cat_A_T*EM*T[0][0])#/(KM_A_T+T[0][0])
			dAdt[nx-1][0] += (k_cat_A_T*EM*T[nx-1][0])#/(KM_A_T+T[nx-1][0])
			dAdt[0][ny-1] += (k_cat_A_T*EM*T[0][ny-1])#/(KM_A_T+T[0][ny-1])
			dAdt[nx-1][ny-1] += (k_cat_A_T*EM*T[nx-1][ny-1])#/(KM_A_T+T[nx-1][ny-1])

		T = np.add(T, dTdt*dt)
		A = np.add(A, dAdt*dt)
		integrals[0].append(integrateCartesian(T, dx, dy))
		integrals[1].append(integrateCartesian(A, dx, dy))
		#Plotting the simulation to heatmap
		if not j:# % 10:
			fig, ax = plt.subplots()
			c = ax.pcolormesh(X,Y, T, cmap='RdBu_r')#, vmin=tMin, vmax=tMax)
			ax.set_title('T')
			ax.axis([dx/2, Lx-dx/2, dy/2, Ly-dy/2])
			fig.colorbar(c, ax=ax)
			plt.show()
			
			fig, ax = plt.subplots()
			c = ax.pcolormesh(X,Y, A, cmap='RdBu_r')#, vmin=tMin, vmax=tMax)
			ax.set_title('A')
			ax.axis([dx/2, Lx-dx/2, dy/2, Ly-dy/2])
			fig.colorbar(c, ax=ax)
			plt.show()
			#outputName = '2D_outputs/output'+"{:05d}".format(j)+'.jpg'
			#plt.savefig(outputName, bbox_inches='tight')
			#plt.close()
	
	color = "%06x" % random.randint(0, 0xFFFFFF)
	plt.plot(t, integrals[0], "#"+str(color), label="T", linewidth=4)
	color = "%06x" % random.randint(0, 0xFFFFFF)
	plt.plot(t, integrals[1], "#"+str(color), label="A", linewidth=4)

	plt.xlabel('Time')
	plt.ylabel('Quantity')
	plt.xlim(-0.01)
	plt.ylim(-0.01)
	plt.grid()
	plt.legend(loc='best')
	plt.show()	
#########################################################################################################
################################ POLAR MODEL ############################################################
#########################################################################################################
#Model with 3 chemicals and one reaction
#T + A -> P + A
def FTCSmethod2DPolar(R, numberOfNodes, numberOfSegments, initMatrix, boundary, thermalConst, finalTime):
	#Getting the parameters
	n = numberOfNodes
	nS = numberOfSegments
	#Temperature / chemical #1 initialization
	T = initMatrix
	#Temperature flow / chemical intake
	TS = boundary
	#Chemical #2 initialization
	A = applyToMatrix(np.zeros((nS,n+1)),0,0)
	P = np.zeros((nS,n+1))
	dR = float(R)/float(n)
	#rad = np.linspace(dR/2, R-dR/2, n)
	rad = np.linspace(0, R, n+1)
	azm = np.linspace(0, 2*np.pi, nS)
	r, th = np.meshgrid(rad, azm)
	t_final = finalTime
	D_T = thermalConst
	dPhi = (2*np.pi)/numberOfSegments

	#Calculating suitable dt
	dt = calcTimeStepPolar(R, nS, n, float(D_T))

	#Reaction parameters
	kcat_A_T = 2.2
	KM_A_T = 0.052
	
	
	#Finding minimum and maximum values of temperature
	#They're either in initMatrix or boundary
	tMin = initMatrix[0][0]
	for i in initMatrix:
		tMin = tMin if tMin < min(i) else min(i)
	tMin = tMin if tMin < boundary else boundary

	tMax = initMatrix[0][0]
	for i in initMatrix:
		tMax = tMax if tMax > max(i) else max(i)
	tMax = tMax if tMax > boundary else boundary


	t = np.arange(0, t_final, dt)
	dTdt = np.empty((nS, n+1))
	dAdt = np.empty((nS, n+1))
	dPdt = np.empty((nS, n+1))
	
	integrals = [[],[]]
	integrals[0].append(integratePolar(T, dR, nS))
	integrals[1].append(integratePolar(P, dR, nS))

	#Running the simulation
	#For Neumann boundary	
	#boundT = np.ones(nS)*2
	for j in range(1,len(t)):
		time = t[j]
		for i in range(0,nS):
			phi = 2*np.pi*float(i)/float(nS-1)
			for k in range(1, n):
				#Inner part of the slice
###################################################################
				#TEMPERATURE / CHEMICAL T
###################################################################
				dTdt[i][k] = D_T*(((T[(i+1)%nS][k]-2*T[i][k]+T[(i-1)%nS][k])/((dPhi**2)*(k*dR)**2))+((T[i][k+1]-2*T[i][k]+T[i][k-1])/dR**2)+((T[i][k+1]-T[i][k-1])/(dR*k*dR)))
				#Reaction modifier A+T
				dTdt[i][k] += ((-kcat_A_T*A[i][k]*T[i][k])/(KM_A_T+T[i][k]))
###################################################################
				#CHEMICAL A
###################################################################
				#dAdt[i][k] = D_T*(((A[(i+1)%nS][k]-2*A[i][k]+A[(i-1)%nS][k])/((dPhi**2)*(k*dR)**2))+((A[i][k+1]-2*A[i][k]+A[i][k-1])/dR**2)+((A[i][k+1]-A[i][k-1])/(dR*k*dR)))
###################################################################
				#PRODUCT P 
###################################################################
				dPdt[i][k] = D_T*(((P[(i+1)%nS][k]-2*P[i][k]+P[(i-1)%nS][k])/((dPhi**2)*(k*dR)**2))+((P[i][k+1]-2*P[i][k]+P[i][k-1])/dR**2)+((P[i][k+1]-P[i][k-1])/(dR*k*dR)))
				#Reaction modifier A + T
				dPdt[i][k] += ((kcat_A_T*A[i][k]*T[i][k])/(KM_A_T+T[i][k]))
		#Zero point end of slice takes values from points in + directions 			
		#Angles:
		#90 degrees - nS/4
		#180 degrees - nS/2
		#270 degrees - 3*nS/2
######### Temperature/Chemical T ##########################################################
			#dTdt[i][0] = D_T*((T[0][1]*calculateS(1, dR, nS, n)-4*T[i][0]*calculateS(0, dR, nS, n)+T[int(nS/4)][1]*calculateS(1, dR, nS, n)+T[int(nS/2)][1]*calculateS(1, dR, nS, n)+T[int(3*nS/4)][1]*calculateS(1, dR, nS, n))/dR**2)/calculateS(0, dR, nS, n)
			dTdt[i][0] = D_T*((T[0][1]-4*T[i][0]+T[int(nS/4)][1]+T[int(nS/2)][1]+T[int(3*nS/4)][1])/dR**2)
		#Reaction modifier A + T
			dTdt[i][0] += ((-kcat_A_T*A[i][0]*T[i][0])/(KM_A_T+T[i][0]))
######### Chemical A ######################################################################
			#dAdt[i][0] = D_T*((A[0][1]-4*A[i][0]+A[nS/4][1]+A[nS/2][1]+A[3*nS/4][1])/dR**2)
######### Product P #######################################################################
			dPdt[i][0] = D_T*((P[0][1]-4*P[i][0]+P[nS/4][1]+P[nS/2][1]+P[3*nS/4][1])/dR**2)
		#Reaction modifier A + T
			dAdt[i][0] += ((kcat_A_T*A[i][0]*T[i][0])/(KM_A_T+T[i][0]))
		#Boundary end of slice takes values from boundary
		#Dirichlet boundaries
			#boundT = 2
			boundA = 0
			boundP = 0
		#Sine boundary
		#bound = 20*sin(2*np.pi*(float(i)/float(nS-1)))
		#Neumann boundary
			#if j == 0:
			#	boundT = 0.5
			#else:
			#	boundT = 0.5*(cos(phi))*2*dR+T[i][n-1]
			boundT = T[i][n-3]
			boundP = P[i][n-3]
			#print boundT
		#Boundary point end of slice takes values from boundary
############# Temperature / Chemical T ############################################################################
			dTdt[i][n] = D_T*(((T[(i+1)%nS][n]-2*T[i][n]+T[(i-1)%nS][n])/((dPhi**2)*(n*dR)**2))+((boundT-2*T[i][n]+T[i][n-1])/dR**2)+((boundT-T[i][n-1])/(dR*n*dR)))
		#Reaction modifier A + T
			dTdt[i][n] += ((-kcat_A_T*A[i][n]*T[i][n])/(KM_A_T+T[i][n]))
############# Chemical A ##########################################################################################
			#dAdt[i][n] = D_T*(((A[(i+1)%nS][n]-2*A[i][n]+A[(i-1)%nS][n])/((dPhi**2)*(n*dR)**2))+((boundA-2*A[i][n]+A[i][n-1])/dR**2)+((boundA-A[i][n-1])/(dR*n*dR)))
############# Product P ###########################################################################################
			dPdt[i][n] = D_T*(((P[(i+1)%nS][n]-2*P[i][n]+P[(i-1)%nS][n])/((dPhi**2)*(n*dR)**2))+((boundP-2*P[i][n]+P[i][n-1])/dR**2)+((boundP-P[i][n-1])/(dR*n*dR)))
		#Reaction modifier A + T
			dPdt[i][n] += ((kcat_A_T*A[i][n]*T[i][n])/(KM_A_T+T[i][n]))
			
		#Changes to T
		T = np.add(T, dTdt*dt)
		#print T
		#Changes to A
		#A = np.add(A, dAdt*dt)
		#Changes to P
		P = np.add(P, dPdt*dt)
		#Plotting the simulation to heatmap
		#if (j >= int(len(t))/2 -3 and j <= int(len(t))/2+3) or j >= len(t)-3:
		integrals[0].append(integratePolar(T, dR, nS))
		integrals[1].append(integratePolar(P, dR, nS))

		if j == 100:#>= len(t)-1:
			#if j > 10:
			#	return 0
			#Plot T
			#print T
			plt.subplot(projection="polar")
			c = plt.pcolormesh(th,r,T)#,vmin=tMin, vmax=tMax)
			plt.plot(azm,r,color='k',ls='none')
			plt.colorbar(c)
			plt.title("Chemical A", loc="left")
			#plt.grid()
			plt.show()
			
			#Plot A
			#print A
			#plt.subplot(projection="polar")
			#c = plt.pcolormesh(th,r,A, vmin=0, vmax=10)
			#plt.plot(azm,r,color='k',ls='none')
			#plt.colorbar(c)
			#plt.title("Chemical A", loc="left")
			#plt.grid()
			#plt.show()
		
			#Plot P
			#print P
			plt.subplot(projection="polar")
			c = plt.pcolormesh(th,r,P,vmin=0, vmax=3)
			plt.plot(azm,r,color='k',ls='none')
			plt.colorbar(c)
			plt.title("Product P", loc="left")
			plt.grid()
			plt.show()	
			
			#outputName = '2D_outputs/output'+"{:05d}".format(j)+'.jpg'
			#plt.savefig(outputName, bbox_inches='tight')
			#plt.close()

	color = "%06x" % random.randint(0, 0xFFFFFF)
	plt.plot(t, integrals[0], "#"+str(color), label="A", linewidth=4)
	color = "%06x" % random.randint(0, 0xFFFFFF)	
	plt.plot(t, integrals[1], "#"+str(color), label="P", linewidth=4)
	#color = "%06x" % random.randint(0, 0xFFFFFF)
	#plt.plot(t, integrals[1], "#"+str(color), label="A", linewidth=4)

	plt.xlabel('Time')
	plt.ylabel('Quantity')
	plt.xlim(-0.01)
	plt.ylim(-0.01)
	plt.grid()
	plt.legend(loc='best')
	plt.show()	

#########################################################################################################
################################ POLAR MODEL ############################################################
#########################################################################################################
#Model with 10 chemicals A B C D E F G H T P
#4 reactions:
#T + A -> P + A
#P + E + D -> F + E
#B + C -> D
#D + G -> H
def FTCSmethod2DPolarExtended(R, finalTime, names, inits, diffConsts, reactionMods, consts, bounds):
	#Getting the parameters
	nS = len(init[0])
	n = len(init[0][0])-1
	
	dR = float(R)/float(n)
	rad = np.linspace(0, R, n+1)
	azm = np.linspace(0, 2*np.pi, nS)
	r, th = np.meshgrid(rad, azm)
	t_final = finalTime
	dPhi = (2*np.pi)/nS

	dt = finalTime
	#Calculating suitable dt
	for dConst in diffConsts:
		newDt = calcTimeStepPolar(R, nS, n, float(dConst))
		dt = dt if dt < newDt else newDt 
	
	#Dictionaries
	U = dict(zip(names,init))
	diffC = dict(zip(names,diffConsts))
	rMods = dict(zip(names,reactionMods))
	
	dUmatrices = []
	for matrix in inits:
		dUmatrices.append(np.empty((nS, n+1)))

	dUdt = dict(zip(names,dUmatrices))
	
	boundType = []
	boundVal = []
	
	for b in bounds:
		boundType.append(b[0])
		boundVal.append(b[1])
	
	bt = dict(zip(names,boundType))
	bv = dict(zip(names,boundVal))
	#time range
	t = np.arange(0, t_final, dt)
	#Running the simulation
	for j in range(0,len(t)):
		time = t[j]
		for i in range(0,nS):
			phi = np.pi*float(i)/float(nS-1)
			for k in range(1, n):
				for name in names:
					#Diffusion modifier
					dUdt[name][i][k] = diffC[name]*(((U[name][(i+1)%nS][k]-2*U[name][i][k]+U[name][(i-1)%nS][k])/((dPhi**2)*(k*dR)**2))+((U[name][i][k+1]-2*U[name][i][k]+U[name][i][k-1])/dR**2)+((U[name][i][k+1]-U[name][i][k-1])/(dR*k*dR)))
					#Reaction modifier
					dUdt[name][i][k] += eval(rMods[name], locals(), consts)

			for name in names:
				#Center modifications
				k = 0
				#Center diffusion modifier
				dUdt[name][i][k] = diffC[name]*((U[name][0][k+1]-4*U[name][i][k]+U[name][nS/4][k+1]+U[name][nS/2][k+1]+U[name][3*nS/4][k+1])/dR**2)
				#Center reaction modifier
				dUdt[name][i][k] += eval(rMods[name], locals(), consts)

				#Boundary modifications
				k = n
				#Boundary diffusion modifier
				if bt[name] == "f" or bt[name] == "F":
					#Dirichlet boundary or an expression ( f(phi)
					dUdt[name][i][k] = diffC[name]*(((U[name][(i+1)%nS][k]-2*U[name][i][k]+U[name][(i-1)%nS][k])/((dPhi**2)*(k*dR)**2))+((eval(bv[name], globals(), locals())-2*U[name][i][k]+U[name][i][k-1])/dR**2)+((eval(bv[name], globals(), locals())-U[name][i][k-1])/(2*dPhi*k*dR)))
				elif bt[name] == "n" or bt[name] == "N":
					#Neumann boundary
					boundNeumann = (bv[name]/diffC[name]-(dUdt[name][(i+1)%n][n]-2*dUdt[name][i][n]+dUdt[name][(i-1)%n][n])/dPhi**2)*dR**2+2*dUdt[name][i][n]+dUdt[name][i][n-1]
					dUdt[name][i][k] = diffC[name]*(((U[name][(i+1)%nS][k]-2*U[name][i][k]+U[name][(i-1)%nS][k])/((dPhi**2)*(k*dR)**2))+((boundNeumann-2*U[name][i][k]+U[name][i][k-1])/dR**2)+((boundNeumann-U[name][i][k-1])/(2*dPhi*k*dR)))

				#Boundary reaction modifier
				dUdt[name][i][k] += eval(rMods[name], locals(), consts)
		
		#Apply changes
		for name in names:
			U[name] = np.add(U[name], dUdt[name]*dt)

		#Plotting the simulation to heatmap
		if (j):# >= int(len(t))/2 -3 and j <= int(len(t))/2+3) or j >= len(t)-3:
			#if j > 10:
			#	return 0
			for name in names:
				#Plot
				print U[name]
				plt.subplot(projection="polar")
				c = plt.pcolormesh(th,r,U[name])
				plt.plot(azm,r,color='k',ls='none')
				plt.colorbar(c)
				plt.title(name, loc="left")
				#plt.grid()
				plt.show()
			
			#outputName = '2D_outputs/output'+"{:05d}".format(j)+'.jpg'
			#plt.savefig(outputName, bbox_inches='tight')
			#plt.close()


def initFunction1(r):
	return np.e**(-(r**2)*10)


def applyToMatrix(matrix, funct, consts):
	ret = matrix
	for i in range(0,numberOfNodes):
		for j in range(0, numberOfNodes+1):
			ret[i][j] = 0.012#100 if j == 0 else 0
			#initFunction1(length*j/numberOfNodes)
	return ret

#finiteVolumeMethod(leftBoundary, rightBoundary, finalTime, length, numberOfNodes, D0, energyS)
#FTCSexample(leftBoundary, rightBoundary, finalTime, length, numberOfNodes, D0, energyS)
#FTCSmethod2D(length, length, numberOfNodes, numberOfNodes, np.zeros((numberOfNodes, numberOfNodes)), leftBoundary, leftBoundary, leftBoundary, leftBoundary, 0.001, finalTime)
test = np.zeros((numberOfNodes, numberOfNodes+1))
for i in test:
	i[len(i)-20] = 10
FTCSmethod2DPolar(length, numberOfNodes, numberOfNodes, test, leftBoundary, 0.01, finalTime)
names = ["T", "A", "P"]#, "B", "C", "D", "E", "F", "G" ,"H"]
init = []
init.append(np.zeros((numberOfNodes, numberOfNodes+1)))
init.append(applyToMatrix(np.zeros((numberOfNodes, numberOfNodes+1)),0,0))
init.append(np.zeros((numberOfNodes, numberOfNodes+1)))
#init.append(np.zeros((numberOfNodes, numberOfNodes+1)))
#init.append(np.zeros((numberOfNodes, numberOfNodes+1)))
#init.append(np.zeros((numberOfNodes, numberOfNodes+1)))
#init.append(np.zeros((numberOfNodes, numberOfNodes+1)))
#init.append(np.zeros((numberOfNodes, numberOfNodes+1)))
#init.append(np.zeros((numberOfNodes, numberOfNodes+1)))
#init.append(np.zeros((numberOfNodes, numberOfNodes+1)))
constants = {
	#Reaction A + T -> P parameters
	"kcat_A_T": 2.2,
	"KM_T": 0.052,
	#"kcat_B_C": 1.8,
	#"KM_P": 0.052,
	#"kcat_P_E_D": 1.1
	#"kcat_D_G": 6,
	#"KM_D": 0.052
}
#Reaction modifiers are different for reactions with enzymes:
#One substrate:
# V = (V_max*E*S)/(K_M+S)
#Two substrates:
# V = (V_max*E*S_1*S_2)/((K_S1+S1)*(K_S2+S2))
reactMods = []
#T
reactMods.append("((-kcat_A_T*U[\"A\"][i][k]*U[\"T\"][i][k])/(KM_T+U[\"T\"][i][k]))")
#A
reactMods.append("0")
#P
reactMods.append("((kcat_A_T*U[\"A\"][i][k]*U[\"T\"][i][k])/(KM_T+U[\"T\"][i][k]))")#+((-kcat_P_E_D*U[\"P\"][i][k]*U[\"E\"][i][k]*U[\"D\"][i][k])/((KM_P+U[\"P\"][i][k])+(KM_D+U[\"D\"][i][k])))")
#B
#reactMods.append("")
#C
#reactMods.append("")
#D
#reactMods.append("++((-kcat_P_E_D*U[\"P\"][i][k]*U[\"E\"][i][k]*U[\"D\"][i][k])/((KM_P+U[\"P\"][i][k])+(KM_D+U[\"D\"][i][k])))")
#E
#reactMods.append("0")
#F
#reactMods.append("((kcat_P_E_D*U[\"P\"][i][k]*U[\"E\"][i][k]*U[\"D\"][i][k])/((KM_P+U[\"P\"][i][k])+(KM_D+U[\"D\"][i][k])))")
#G
#reactMods.append("0")
#H
#reactMods.append("((kcat_D_G*U[\"D\"][i][k]*U[\"G\"][i][k])/(KM_D+U[\"D\"][i][k]))")

diffConsts = [0.01, 0.01, 0.01]#, 0.025, 0.05, 0.03, 0.033, 0.06, 0.04, 0.045, 0.015 ,0.005]
#2+sin(phi)
bounds = [("f", "2"), ("f","0"),("f","0")]
#FTCSmethod2DPolarExtended(length, finalTime, names, init, diffConsts, reactMods, constants, bounds)
