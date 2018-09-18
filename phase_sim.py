################  PARAMETER RECOVERY ######################
############# from simulated data #########################

from pylab import *
import sympy as sm # many conflicts with pylab, numpy,math,...
from scipy import signal
import time


#### SET THESE FLAGS (Default: 0,0)######

direct_plot_mode=0
phase_plot_mode=0
plot_details=0
residual_plot=1


#####----reading data------#################
print("started")
#example from CR-23-01-2018 Phase_tracking_sin.nb
t,rex,imx,rey,imy=genfromtxt("eout.csv",delimiter='\t',unpack='True')




####----phases plotting (in need be)----###########

if phase_plot_mode :

	f1=figure()
	plot(rex,label='rex')
	plot(imx,label='imx')
	f2=figure()
	plot(rey,label='rey')
	plot(imy,label='imy')


########################################
######## DSP ALGORITHM         #########
########################################


#### W Matrix creation ########
#### Symbolic calculation #####

delta, theta,phi=sm.symbols("delta theta phi", real=True)

##### implementation of Eo= R*M*R*Ein

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

#### fundamental variables ##########

B=array([[1],[2],[2.5]]) # Beta-point
deB=array([[0],[0],[0]]) # Beta-increment
Yt=array([[0],[0],[0],[0]]) # 4 values of the coherent receiver 
mod=1 # module of the Y vector
deY=array([[0],[0],[0],[0]]) 
B1=[] # auxiliary list
B2=[] # auxiliary list
B3=[] # auxiliary list
dell=[] #retrieved delta
thel=[] # retrieved theta
phil=[] # tetreived phi
detl=[] # determinat of (tW*W)
sigdel=[]
sigthe=[]
sigphi=[]

#weighted average
theb=array([])
phib=array([])
them=B[1,0]
phim=B[2,0]
l=20
thew=them
phiw=phim
theml=[]


def myY(t):
	a=array([[rex[t]],[imx[t]],[rey[t]],[imy[t]]])
	return a

def modx(x):
	b=sqrt(rex[x]**2 + imx[x]**2 + rey[x]**2 + imy[x]**2)
	return b
###################################################	

######### loop section ###############

print ("start loop")
t1=time.time()


for i in range(100000):

	
	mod= modx(i)
	deY=(myY(i)/mod)-Yt
	

	Wn=Wt(B[0,0],B[1,0],B[2,0])
	u=Prodv(B[0,0],B[1,0],B[2,0])
	d=linalg.det(u)
	detl.append(1/d)
	

	uu=inv(u)
	sigd=uu[0,0]
	sigt=uu[1,1]
	sigp=uu[2,2]
	deB=uu.dot(Wn).dot(deY)
	#print(sigt)
	B=B+deB
	
	lt=theb.size
	if lt==l :
		
		theb[0:l-1]=theb[1:(l)]
		theb[l-1]=thew
		phib[0:l-1]=phib[1:(l)]
		phib[l-1]=phiw
	else :
		theb=append(theb,thew)
		phib=append(phib,phiw)
	lt=theb.size
	them=theb.sum()/lt
	phim=phib.sum()/lt
	theml.append(them)
	
	#if 10000<i<15000 :
	#	print(i,") B=",B[1,0],"sigt=",sigt,"them=",them)
	#	print("thew=",thew)
	sigstar=10**(-10)
	thew=(B[1,0]*(1/(sigt)) + them*sigstar)/((1/(sigt))+sigstar)
	phiw=(B[2,0]*(1/(sigp)) + phim*sigstar)/((1/(sigp))+sigstar)
	

	B1=B[0,0]
	B2=thew
	B3=phiw

	#print("the medio:",them,"thew:",thew)
	dell.append(B1)
	thel.append(B2)
	phil.append(B3)
	sigdel.append(sigd)
	sigthe.append(sigt)
	sigphi.append(sigp)
	
	Yt=Fun(B1,B2,B3)
	
	
print("loop finished")
print("elapsed time: ",time.time()-t1)
###########################################



########### plot section ##################


#

f1=figure()
ax1=f1.add_subplot(1,1,1)
ax1.plot(dell,color=(1,0,0), label='delta')#linestyle='--', marker='o')
ax1.plot(thel,color=(1,0.65,0), label='theta')
ax1.plot(phil,color=(0,1,0), label='phi')
#ax1.plot(theml,color=(0,0.8,0.8), label='them',linestyle='--', marker='o')
ax1.legend(loc=2)
ax1.set_ylabel("phase (rad)")
ax1.grid(linestyle='--')

if plot_details :
	ax12=ax1.twinx()
	ax12.plot(sigdel,color=(0.6,0,0), label='sigdel')
	ax12.plot(sigthe,color=(0.8,0.45,0), label='sigthe')
	ax12.plot(sigphi,color=(0,0.5,0), label='sigphi')
	ax12.plot(detl,color=(0,0,0), label='determ')
	ax12.legend(loc=1)
	ax12.set_ylabel("param variance (rad^2)")
	ax12.ticklabel_format(style='sci', axis='y', scilimits=(0,0),useMathText = True)
	#ax12.yaxis.major.formatter._useMathText = True


if direct_plot_mode :
	phiexare=array(phiexe)
	phiexaro=array(phiexo)
	ax1.plot(detl,color='black',linestyle='dashed')
	ax1.plot(phiexare,color='blue',linestyle='dashed')
	ax1.plot(phiexaro,color='pink',linestyle='dashed')
#yticks([3.14*n for n in arange(-5,5)])


if residual_plot :
	phior=genfromtxt("phidiff.csv",delimiter='\t',unpack='True')
	residual=array(phior)-roll(array(phil),0)
	f2=figure()
	ax2=f2.add_subplot(1,1,1)
	ax2.set_title("Residual")
	ax2.plot(residual,color='green')
	ax2.grid()

show()

	




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










