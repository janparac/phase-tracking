
############# SIMULATOR OF THE COHERENT RECEIVER CHANNELS ############
####         with random-walk functions for delta, theta, phi   #####  


import sys,getopt
from numpy import *
import numpy.random as npr

import datagenerator as dg

#####################################
###########----entries----###########

mypoints=100000


class RandomWalk:

	last=0
	
	def funrand(self,t):

		step = npr.randint(0,2)
		if step == 1:
			self.last+= (10**(-2))*1 #+ (10**(-5))*npr.normal())
		else:
			self.last-= (10**(-2))*1 #+ (10**(-5))*npr.normal())	
		return self.last

######################################

def optionparser():
	sh=0
	allarg=sys.argv[1:]
	
	try:
		opts, args = getopt.getopt(allarg, "p")
	except getopt.GetoptError:
		print ("option not valid. Available options:")
		print ("FUNDEF_rand -p")
		sys.exit(2)
	for opt, arg in opts:
		if opt=='-p':
			print("Showing previous data...")
			sh=1
	return sh


show=optionparser()
dobj=RandomWalk()
tobj=RandomWalk()
pobj=RandomWalk()



dg.datagen(dobj.funrand, tobj.funrand, pobj.funrand, mypoints,1,1,show)


