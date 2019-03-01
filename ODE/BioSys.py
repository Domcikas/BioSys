#!/usr/bin/env python2
from scipy.integrate import odeint
import matplotlib.pyplot as plt
import sys
import numpy as np
import random
from math import sin, cos, tan
assert len(sys.argv) > 1, "Error in "+sys.argv[0]+": too few arguments."
filename = ""
toFile = False
toScreen = True
for arg in sys.argv[1:]:
	parameter = arg.split("-")
	if parameter[0] == "":
		if parameter[1] == "H" or parameter[1] == "h":
			print sys.argv[0]+" help page."
			print "Usage: ./BioSys.py [options] input_file\n"
			print "Options:"
			print "-f,	to print output to file of the same name as input."
			print "-s,	to show output in screen."
			print "	This option is default and only turned off if -f is used.\n"
			print "Input file format:"
			print "x y 	  , where x equals to number of constants"
			print "	  , y equals to number of equations."
			print "A z , where A is constant name and z ir constant value."
			print "Dy A*Dy m , where Dy ir equation name"
			print "	  , A*Dy is the equation"
			print "	  , m is number for initial value of Dy."
			print "t1 t2 t3 , where t1 is time starting point, t2 is end point and t3 is number of steps between t1 and t2."
			sys.exit()
		elif parameter[1] == "F" or parameter[1] == "f":
			toFile = True
			toScreen = False
		elif parameter[1] == "S" or parameter[1] == "s":
			toScreen = True
			toFile = False
	else:
		filename = parameter[0];

assert filename, "Error in "+sys.argv[0]+": input file name is missing."
F = open(filename,'r')
currentLine = F.readline().rstrip()
systemSizes = currentLine.split(" ")
constantNames = []
constants = [] 
compositorNames = []
compositorValues = []
rates = []
for i in range(0, int(systemSizes[0])):
	currentLine = F.readline().rstrip()
	temp = currentLine.split(" ")
	constantNames.append(temp[0])
	constants.append(float(eval(temp[1], dict(zip(constantNames,constants)))))

for i in range(0, int(systemSizes[1])):
	currentLine = F.readline().rstrip()
	temp = currentLine.split(" ")
	compositorNames.append(temp[0])
	compositorValues.append(float(temp[2]))
	rates.append(temp[1])

def Biosystem(y, t, diff,  compositors, const, constVal):
	equations = []
	dictionary = dict(zip(compositors, y))
	temp = dict(zip(const, constVal))
	temp["t"] = t
	dictionary.update(temp)
	for i in range(0, len(diff)):
		equations.append(eval(diff[i], globals(),dictionary))
	return equations

currentLine = F.readline().rstrip()
time = currentLine.split(" ")
t = np.linspace(int(time[0]), int(time[1]), int(time[2]))

sol = odeint(Biosystem, compositorValues, t, args=(rates, compositorNames, constantNames, constants))

for i in range(0,len(compositorNames)):
	color = "%06x" % random.randint(0, 0xFFFFFF)
	plt.plot(t, sol[:, i], "#"+str(color), label=compositorNames[i], linewidth=4)

plt.xlabel('Time')
plt.ylabel('Quantity')
plt.xlim(-0.01)
plt.ylim(-0.01)
plt.grid()
plt.legend(loc='best')

if toScreen:
	plt.show()
elif toFile:
	plt.savefig(filename.split(".")[0]+".pdf", bbox_inches='tight')
	print sys.argv[0]+": Saved output of "+filename+" to "+filename.split(".")[0]+".png"
F.close()
