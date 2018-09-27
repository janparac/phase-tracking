############# SIMULATOR OF THE COHERENT RECEIVER CHANNELS ############
####         for deterministic functions for delta, theta, phi   #####  


from numpy import *

import datagenerator as dg

from randomwalker import RandomWalk


###########----entries----###########

mypoints=100000
Ts=4*10**(-7)
Fs=1/Ts

myrand1=RandomWalk()
myrand2=RandomWalk()
myrand3=RandomWalk()

def fundel(n):
	return 2+1*cos(2*pi*40*Ts*n)+sqrt((2*pi*pi*10**3)/(Fs))*myrand1.funrand(n)
	
def funthe(n):
	return 1+1*cos(2*pi*40*(4*10**(-7))*n)+sqrt((2*pi*pi*10**3)/(Fs))*myrand2.funrand(n)

def funphi(n):
	return 2+1*cos(2*pi*10*(4*10**(-7))*n)+sqrt((2*pi*pi*10**3)/(Fs))*myrand3.funrand(n)

dg.datagen(fundel, funthe, funphi, mypoints,1,0,0)


