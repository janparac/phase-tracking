import numpy.random as npr



class RandomWalk:

	last=0
	
	def funrand(self,t):

		step = npr.normal(0,1)
		self.last+= step*(10**(-3)) #+ (10**(-5))*npr.normal())	
		return self.last
