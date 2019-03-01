import numpy as np
import matplotlib.pyplot as plt
import math
from math import sin, cos, tan, pi
length = 1
numberOfNodes = 50
leftBoundary = 2
rightBoundary = 0.5
finalTime = 15

D0 = np.ones(numberOfNodes)*0.1
energyS = np.zeros(numberOfNodes)
def forX(x):
	#return 2-1.5*x+sin(math.pi*x)
	return 0
	#return 2*math.exp(-(x-0.5)**2/4)
	#return 4*x-4*x*x

def calcTimeStep(dx, D):
	return 0.5*dx**2/D

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
	#Running the simulation
	for j in range(1,len(t)):
		#T1S = 2*sin(t[j])
		#T2S = 2*cos(t[j])
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


def FTCSmethod2D(lengthX, lengthY, numberOfNodesX, numberOfNodesY, initMatrix, leftBoundaryX, leftBoundaryY, rightBoundaryX, rightBoundaryY, thermalConst, finalTime, timeStep):
	#Getting the parameters
	Lx = lengthX
	Ly = lengthY
	nx = numberOfNodesX
	ny = numberOfNodesY
	T = initMatrix
	T1Sx = leftBoundaryX
	T1Sy = leftBoundaryY
	T2Sx = rightBoundaryX
	T2Sy = rightBoundaryY
	dx = Lx/nx
	dy = Ly/ny
	x = np.linspace(dx/2, Lx-dx/2, nx)
	y = np.linspace(dy/2, Ly-dy/2, ny)
	X,Y = np.meshgrid(x,y)
	t_final = finalTime
	dt = timeStep
	alpha = thermalConst
	
	#Finding minimum and maximum values of temperature
	#They're either in initVector or one boundaries
	tMin = min(initVector)
	tMin = tMin if tMin < leftBoundaryX else leftBoundaryX
	tMin = tMin if tMin < rightBoundaryX else rightBoundaryX
	tMin = tMin if tMin < leftBoundaryY else leftBoundaryY
	tMin = tMin if tMin < rightBoundaryY else rightBoundaryY

	tMax = max(initVector)
	tMax = tMax if tMax > leftBoundaryX else leftBoundaryX
	tMax = tMax if tMax > rightBoundaryX else rightBoundaryX
	tMax = tMax if tMax > leftBoundaryY else leftBoundaryY
	tMax = tMax if tMax > rightBoundaryY else rightBoundaryY

	t = np.arange(0, t_final, dt)
	dTdt = np.empty(nx, ny)
	temperature = []

	#Running the simulation
	for j in range(1,len(t)):
		for i in range(1,nx-1):
			for k in range(1, ny-1):
				dTdt[i][j] = alpha*((-(T[i][k]-T[i-1][k])/dx**2+(T[i+1][k]-T[i][k])/dx**2)+(-(T[i][k]-T[i][k-1])/dy**2+(T[i][k+1]-T[i][k])/dy**2))
			dTdt[i][0] = alpha*((-(T[i][0]-T[i-1][0])/dx**2+(T[i+1][0]-T[i][0])/dx**2)+(-(T[i][0]-leftBoundaryY)/dy**2+(T[i][1]-T[i][0])/dy**2))
			dTdt[i][ny-1] = alpha*((-(T[i][ny-1]-T[i-1][ny-1])/dx**2+(T[i+1][ny-1]-T[i][ny-1])/dx**2)+(-(T[i][ny-1]-T[i][ny-2])/dy**2+(rightBoundaryY-T[i][ny-1])/dy**2))

		for y in range(1, ny-1):
			dTdt[0][y] = alpha*((-(T[0][y]-leftBoundaryX)/dx**2+(T[1][y]-T[0][y])/dx**2)+(-(T[0][y]-T[0][y-1])/dy**2+(T[0][y+1]-T[0][y])/dy**2))
			dTdt[nx-1][y] = alpha*((-(T[nx-1][y]-T[nx-2][y])/dx**2+(rightBoundaryX-T[n-1][y])/dx**2)+(-(T[n-1][y]-T[n-1][y-1])/dy**2+(T[n-1][y+1]-T[n-1][y])/dy**2))

			dTdt[0][0] = alpha*((-(T[0][0]-leftBoundaryX)/dx**2+(T[1][0]-T[0][0])/dx**2)+(-(T[0][0]-leftBoundaryY)/dy**2+(T[0][1]-T[0][0])/dy**2))
			dTdt[nx-1][0] = alpha*((-(T[nx-1][0]-T[nx-2][0])/dx**2+(rightBoundaryX-T[nx-1][0])/dx**2)+(-(T[nx-1][0]-leftBoundaryY)/dy**2+(T[nx-1][1]-T[nx-1][0])/dy**2))
			dTdt[0][ny-1] = alpha*((-(T[0][ny-1]-rightBoundaryX)/dx**2+(T[1][ny-1]-T[0][ny-1])/dx**2)+(-(T[0][ny-1]-T[0][ny-2])/dy**2+(rightBoundaryY-T[0][ny-1])/dy**2))
			dTdt[nx-1][ny-1] = alpha*((-(T[nx-1][ny-1]-T[nx-2][ny-1])/dx**2+(leftBoundaryX-T[nx-1][ny-1])/dx**2)+(-(T[nx-1][ny-1]-T[nx-1][ny-2])/dy**2+(leftBoundaryY-T[nx-1][ny-1])/dy**2))

		T += dTdt*dt
		temperature.append(T.tolist())
		#Laiko ciklas
		#Cia daryti atvaizdavima
	
	#Plotting the simulation to heatmap
	fig, ax = plt.subplots()
	c = ax.pcolormesh(x, t, temperature, cmap='RdBu_r', vmin=tMin, vmax=tMax)
	ax.set_title('1D heat diffusion')
	ax.axis([dx/2, L-dx/2, 0, t_final])
	fig.colorbar(c, ax=ax)
	
	plt.show()
			

#finiteVolumeMethod(leftBoundary, rightBoundary, finalTime, length, numberOfNodes, D0, energyS)

FTCSexample(leftBoundary, rightBoundary, finalTime, length, numberOfNodes, D0, energyS)