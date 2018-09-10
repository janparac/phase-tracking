

from pylab import *
import sympy as sm # many conflicts with pylab, numpy,math,...
from scipy import signal
import time


#### SET THESE FLAGS (Default: 0,0)######

direct_plot_mode=0
phase_plot_mode=0


#####----reading data------#################
print("started")
#example from CR-23-01-2018 Phase_tracking_sin.nb
t,rex,imx,rey,imy=genfromtxt("eout.csv",delimiter='\t',unpack='True')




####----phases plotting----###########

if phase_plot_mode :

	f1=figure()
	plot(rex,label='rex')
	plot(imx,label='imx')
	f2=figure()
	plot(rey,label='rey')
	plot(imy,label='imy')

#
########################################
######## DSP ALGORITHM         #########
########################################


#### W Matrix creation ########

delta, theta,phi=sm.symbols("delta theta phi", real=True)

R1=sm.Matrix([[sm.cos(theta),sm.sin(theta)],[-sm.sin(theta),sm.cos(theta)]])
M=sm.Matrix([[sm.exp(sm.I*(-delta/2)),0],[0,sm.exp(sm.I*(delta/2))]])
R2=sm.transpose(R1)
Fib=R2*M*R1
Ein=sm.Matrix([sm.cos(sm.pi/4)*sm.exp(sm.I*(sm.pi/2)),sm.sin(sm.pi/4)])*sm.exp(sm.I*(phi))
Eout=Fib*Ein

rexS=sm.simplify(sm.re(Eout[0]))
imxS=sm.simplify(sm.im(Eout[0]))
reyS=sm.simplify(sm.re(Eout[1]))
imyS=sm.simplify(sm.im(Eout[1]))

vec=sm.Matrix([[rexS,imxS,reyS,imyS]])
var=sm.Matrix([delta,theta,phi])

W=vec.jacobian(var)
Prod=sm.simplify(sm.transpose(W)*W)


Wt=sm.lambdify((delta,theta,phi),transpose(W),'numpy')
Prodv=sm.lambdify((delta,theta,phi),Prod,'numpy')
Fun=sm.lambdify((delta,theta,phi),sm.transpose(vec),'numpy')


############ main loop ##############
########### LSM algorithm ###########


##### declaration section ###########

B=array([[0.5],[0.5],[0.5]])
deB=array([[0],[0],[0]])
Yt=array([[0],[0],[0],[0]])
mod=1
deY=array([[0],[0],[0],[0]])
B1=[]
B2=[]
B3=[]
dell=[]
thel=[]
phil=[]

detl=[]

phibe=zeros(2)#phase buffer for even discon (array)
phibune=[0,0]#phase buffer for even discon unwrapped (list)
phiexe=[] #exact phi even discon(direct mode)
phibo=zeros(2)
phibuno=[0,0]
phiexo=[]

deBb=[0,0] # buffer for delta increments in discon

def myY(t):
	a=array([[rex[t]],[imx[t]],[rey[t]],[imy[t]]])
	return a

def modx(x):
	b=sqrt(rex[x]**2 + imx[x]**2 + rey[x]**2 + imy[x]**2)
	return b
	


print ("start loop")
t1=time.time()

######### loop section ###############

for i in range(len(rex)):

	
	mod= modx(i)
	deY=(myY(i)/mod)-Yt
	
	phasee=arctan(imy[i]/rey[i]) #phase at even discon
	phaseo=2*B[1,0]-arctan(rey[i]/imy[i])  #phase at odd discon
	phibe[1]=phasee
	phibo[1]=phaseo
	phibune=unwrap(phibe*2)/2
	phibuno=unwrap(phibo*2)/2	
	de=phibune[1]-B[2,0]
	do=phibuno[1]-B[2,0]
	if abs(de)>(pi/4):
		phibune[1]-=trunc(de/(pi/4))*(pi/4)
	if abs(do)>(pi/4):
		phibuno[1]-=trunc(do/(pi/4))*(pi/4)
	phiexe.append(phibune[1])
	phiexo.append(phibuno[1])
	phibe[0]=phibune[1]
	phibo[0]=phibuno[1]
	
	Wn=Wt(B[0,0],B[1,0],B[2,0])

	u=Prodv(B[0,0],B[1,0],B[2,0])
	
	d=linalg.det(u)
	detl.append(d)
	
	if d>10**(-7): #det far from zero: LSM. The smaller the value is the closer you can go to the discon point AND the more are the used resources
		uu=inv(u)
		deB=uu.dot(Wn).dot(deY)
	else: # det close to zero: altern method	
		print("det!!",i,' : ',deB[0,0])
		if abs(deBb[0]*deBb[1]<0):
			print ("----reverse----")
			deB[0,0]=deB[0,0]*0.1

		deB[1,0]=0

		if round(asscalar(B[0,0])/pi) & 1: #bitwise operat to check even/odd
			deB[2,0]=(phibuno[1])-B[2,0]
		else :
			deB[2,0]=(phibune[1])-B[2,0]

	deBb[0]=deBb[1]
	deBb[1]=deB[0,0]


	B=B+deB
	
	B1=B[0,0].tolist()
	B2=B[1,0].tolist()
	B3=B[2,0].tolist()
	
	
	dell.append(B1)
	thel.append(B2)
	phil.append(B3)
	
	Yt=Fun(B[0,0],B[1,0],B[2,0])
	
	
print("loop finished")
print("elapsed time: ",time.time()-t1)
###########################################



########### plot section ##################

phior=genfromtxt("phidiff.csv",delimiter='\t',unpack='True')
residual=array(phior)-roll(array(phil),0)

f1=figure()
plot(dell,color='r', label='delta') #linestyle='--', marker='o')
plot(thel,color='orange', label='theta')
plot(phil,color='green', label='phi')

if direct_plot_mode :
	phiexare=array(phiexe)
	phiexaro=array(phiexo)
	plot(detl,color='black',linestyle='dashed')
	plot(phiexare,color='blue',linestyle='dashed')
	plot(phiexaro,color='pink',linestyle='dashed')
#yticks([3.14*n for n in arange(-5,5)])
grid()
legend()


f2=figure()
title("Residual")
plot(residual,color='green')
grid()
show()

input()
	




########################################################################








#######  spare code  #################
#############################
#unwrap buffer |_|_|_|_|_|_|
#--------------|p|p|0|1|2|i|

#phibe=zeros(6)#array
#phibune=[0,0,0,0,0]#list
#c=0 #counter
#ofs=0
#############################

### buffered unwrap

#	if c<3:
#		phibe[c+2]=phase
#		c+=1
#	else:
#		phibe[5]=phase
#		phibune=unwrap(phibe*2)/2
#		phiexe.extend(phibune[2:6])
#		phibe[0:2]=phibune[4:6]
#		c=0










