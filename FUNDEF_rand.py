
############# SIMULATOR OF THE COHERENT RECEIVER CHANNELS ############
####         with random-walk functions for delta, theta, phi   #####  


import sys,getopt
from numpy import *

from randomwalker import RandomWalk

import datagenerator as dg

#####################################
###########----entries----###########

mypoints=100000

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


