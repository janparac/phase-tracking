############# SIMULATOR OF THE COHERENT RECEIVER CHANNELS ############
####         for deterministic functions for delta, theta, phi   #####  


from numpy import *

import datagenerator as dg


###########----entries----###########

mypoints=100000

def fundel(n):
	return 1+1*sin(2*pi*10*(20*10**(-7))*n)
	
def funthe(n):
	return 1+1*cos(2*pi*15*(40*10**(-7))*n)

def funphi(n):
	return 1.5+1*cos(2*pi*10*(40*10**(-7))*n)

dg.datagen(fundel, funthe, funphi, mypoints,1,0,0)

