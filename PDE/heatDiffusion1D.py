import numpy as np
import matplotlib.pyplot as plt
import math
from math import sin, cos, tan, pi
length = 1
numberOfNodes = 50
leftBoundary = 2
rightBoundary = 0.5
finalTime = 60

D0 = np.ones(numberOfNodes)*0.1
energyS = np.zeros(numberOfNodes)
def forX(x):
	#return 2-1.5*x+sin(math.pi*x)
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
	print 0.001*diameter/D
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


def FTCSmethod2D(lengthX, lengthY, numberOfNodesX, numberOfNodesY, initMatrix, leftBoundaryX, leftBoundaryY, rightBoundaryX, rightBoundaryY, thermalConst, finalTime):
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

	#Running the simulation
	for j in range(1,len(t)):
		temperature = []
		for i in range(1,nx-1):
			for k in range(1, ny-1):
				#Inner square
				dTdt[i][k] = alpha*((-(T[i][k]-T[i-1][k])/dx**2+(T[i+1][k]-T[i][k])/dx**2)+(-(T[i][k]-T[i][k-1])/dy**2+(T[i][k+1]-T[i][k])/dy**2))
			#Left and right borders
			dTdt[i][0] = alpha*((-(T[i][0]-T[i-1][0])/dx**2+(T[i+1][0]-T[i][0])/dx**2)+(-(T[i][0]-leftBoundaryY)/dy**2+(T[i][1]-T[i][0])/dy**2))
			dTdt[i][ny-1] = alpha*((-(T[i][ny-1]-T[i-1][ny-1])/dx**2+(T[i+1][ny-1]-T[i][ny-1])/dx**2)+(-(T[i][ny-1]-T[i][ny-2])/dy**2+(rightBoundaryY-T[i][ny-1])/dy**2))

		for y in range(1, ny-1):
			#Top and bottom borders
			dTdt[0][y] = alpha*((-(T[0][y]-leftBoundaryX)/dx**2+(T[1][y]-T[0][y])/dx**2)+(-(T[0][y]-T[0][y-1])/dy**2+(T[0][y+1]-T[0][y])/dy**2))
			dTdt[nx-1][y] = alpha*((-(T[nx-1][y]-T[nx-2][y])/dx**2+(rightBoundaryX-T[nx-1][y])/dx**2)+(-(T[nx-1][y]-T[nx-1][y-1])/dy**2+(T[nx-1][y+1]-T[nx-1][y])/dy**2))
			
			#Corners
			dTdt[0][0] = alpha*((-(T[0][0]-leftBoundaryX)/dx**2+(T[1][0]-T[0][0])/dx**2)+(-(T[0][0]-leftBoundaryY)/dy**2+(T[0][1]-T[0][0])/dy**2))
			dTdt[nx-1][0] = alpha*((-(T[nx-1][0]-T[nx-2][0])/dx**2+(rightBoundaryX-T[nx-1][0])/dx**2)+(-(T[nx-1][0]-leftBoundaryY)/dy**2+(T[nx-1][1]-T[nx-1][0])/dy**2))
			dTdt[0][ny-1] = alpha*((-(T[0][ny-1]-rightBoundaryX)/dx**2+(T[1][ny-1]-T[0][ny-1])/dx**2)+(-(T[0][ny-1]-T[0][ny-2])/dy**2+(rightBoundaryY-T[0][ny-1])/dy**2))
			dTdt[nx-1][ny-1] = alpha*((-(T[nx-1][ny-1]-T[nx-2][ny-1])/dx**2+(leftBoundaryX-T[nx-1][ny-1])/dx**2)+(-(T[nx-1][ny-1]-T[nx-1][ny-2])/dy**2+(leftBoundaryY-T[nx-1][ny-1])/dy**2))

		T = np.add(T, dTdt*dt)
		#temperature = T.tolist()
		#Plotting the simulation to heatmap
		if j % 10:
			fig, ax = plt.subplots()
			c = ax.pcolormesh(X,Y, T, cmap='RdBu_r', vmin=tMin, vmax=tMax)
			ax.set_title('2D heat diffusion')
			ax.axis([dx/2, Lx-dx/2, dy/2, Ly-dy/2])
			fig.colorbar(c, ax=ax)
		
			outputName = '2D_outputs/output'+"{:05d}".format(j)+'.jpg'
			plt.savefig(outputName, bbox_inches='tight')
			plt.close()
			
def FTCSmethod2DPolar(R, numberOfNodes, numberOfSegments, initMatrix, boundary, thermalConst, finalTime):
	#Getting the parameters
	n = numberOfNodes
	nS = numberOfSegments
	#Temperature / chemical #1 initialization
	T = initMatrix
	#Temperature flow / chemical intake
	TS = boundary
	#Chemical #2 initialization
	Ch2 = np.zeros((nS,n+1))
	dR = float(R)/float(n)
	#rad = np.linspace(dR/2, R-dR/2, n)
	rad = np.linspace(0, R, n+1)
	azm = np.linspace(0, 2*np.pi, nS)
	r, th = np.meshgrid(rad, azm)
	t_final = finalTime
	alpha = thermalConst
	dPhi = (2*np.pi)/numberOfSegments

	#Calculating suitable dt
	dt = calcTimeStepPolar(R, nS, n, float(alpha))
	
	
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

	#Running the simulation
	for j in range(0,len(t)):
		for i in range(0,nS):
			symNode = T[(i+numberOfSegments/2)%numberOfSegments][0]
			for k in range(1, n):
###################################################################
			#TEMPERATURE / CHEMICAL #1 DIFFUSION
###################################################################
				#Inner part of the slice
				dTdt[i][k] = alpha*(((T[(i+1)%nS][k]-2*T[i][k]+T[(i-1)%nS][k])/((dPhi**2)*(k*dR)**2))+((T[i][k+1]-2*T[i][k]+T[i][k-1])/dR**2)+((T[i][k+1]-T[i][k-1])/(2*dPhi*k*dR))) #dTdt = alpha*((d^2 T/d^2 phi * r^2)+(d^2 T/d^2 r)+(dT/dr * r))
			#Wide and narrow ends of the slice
			#Zero point end of slice takes values from points in + directions 			
			#dTdt[i][0] = alpha*((-(T[i][0]-T[(i-1)%n][0])/dPhi**2+(T[(i+1)%n][0]-T[i][0])/dPhi**2)+(-(T[i][0]-symNode)/dR**2+(T[i][1]-T[i][0])/dR**2))
			#dTdt[i][0] =  alpha*(((T[(i+1)%n][1]-2*T[i][0]+T[(i-1)%n][0])/((dPhi**2)*(k*dR)**2))+((T[i][k+1]-2*T[i][k]+T[i][k-1])/dR**2))
			#Angles:
			#90 degrees - nS/4
			#180 degrees - nS/2
			#270 degrees - 3*nS/2
			dTdt[i][0] = alpha*((T[0][1]-4*T[i][0]+T[nS/4][1]+T[nS/2][1]+T[3*nS/4][1])/dR**2)
			#Boundary end of slice takes values from boundary
			#Dirichlet boundary
			bound = 2
			#Sine boundary
			#bound = 20*sin(2*np.pi*(float(i)/float(nS-1)))
			#Neumann boundary
			#bound = (boundary/alpha-(T[(i+1)%n][n]-2*T[i][n]+T[(i-1)%n][n])/dPhi**2)*dR**2+2*T[i][n]+T[i][n-1]
			#print bound
			dTdt[i][n] = alpha*(((T[(i+1)%nS][n]-2*T[i][n]+T[(i-1)%nS][n])/((dPhi**2)*(n*dR)**2))+((bound-2*T[i][n]+T[i][n-1])/dR**2)+((bound-T[i][n-1])/(2*dPhi*n*dR)))
#######################################################################################
			#CHEMICAL #2 DIFFUSSION 
#######################################################################################
			
		T = np.add(T, dTdt*dt)
		#Plotting the simulation to heatmap
		if j >= len(t)-3:
			print T
			plt.subplot(projection="polar")
			c = plt.pcolormesh(th,r,T,vmin=tMin, vmax=tMax)
			print tMin
			print tMax
			plt.plot(azm,r,color='k',ls='none')
			plt.colorbar(c)
			#plt.thetagrids([theta * 15 for theta in range(360//15)])
			#plt.grid()
			plt.show()
			#fig, ax = plt.subplots()
			#c = ax.pcolormesh(X,Y, T, cmap='RdBu_r', vmin=tMin, vmax=tMax)
			#ax.set_title('2D heat diffusion')
			#ax.axis([dx/2, Lx-dx/2, dy/2, Ly-dy/2])
			
		
			#outputName = '2D_outputs/output'+"{:05d}".format(j)+'.jpg'
			#plt.savefig(outputName, bbox_inches='tight')
			#plt.close()

def initFunction1(r):
	return np.e**(-(r**2)*10)


#finiteVolumeMethod(leftBoundary, rightBoundary, finalTime, length, numberOfNodes, D0, energyS)
#FTCSexample(leftBoundary, rightBoundary, finalTime, length, numberOfNodes, D0, energyS)
#FTCSmethod2D(length, length, numberOfNodes, numberOfNodes, np.zeros((numberOfNodes, numberOfNodes)), leftBoundary, leftBoundary, leftBoundary, leftBoundary, 0.001, finalTime)

test = np.zeros((numberOfNodes, numberOfNodes+1))
#for i in range(0,numberOfNodes):
#	for j in range(0, numberOfNodes+1):
#		test[i][j] = 10 if j == 0 else 0#initFunction1(length*j/numberOfNodes) 
#test = np.zeros((numberOfNodes, numberOfNodes+1))
print test
FTCSmethod2DPolar(length, numberOfNodes, numberOfNodes, test, leftBoundary, 0.01, finalTime)
