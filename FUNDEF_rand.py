
############# SIMULATOR OF THE COHERENT RECEIVER CHANNELS ############
####         with random-walk functions for delta, theta, phi   #####  



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



dobj=RandomWalk()
tobj=RandomWalk()
pobj=RandomWalk()


dg.datagen(dobj.funrand, tobj.funrand, pobj.funrand, mypoints,1)


