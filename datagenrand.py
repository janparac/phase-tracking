
############# SIMULATOR OF THE COHERENT RECEIVER CHANNELS ############
####         with random-walk functions for delta, theta, phi   #####  


from numpy import *
import numpy.random as npr
from matplotlib.pyplot import *
from mpl_toolkits.mplot3d import Axes3D
import csv

#####################################
###########----entries----###########

points=arange(100000)

###########----end entries---########

def randomwalk(p):
	x=[0]
	for j in arange(p-1):
		step_x = npr.randint(0,2)
		if step_x == 1:
			x.append(x[j] + (10**(-3))*1) #+ (10**(-5))*npr.normal())
		else:
			x.append(x[j] - (10**(-3))*1) #+ (10**(-5))*npr.normal())
	return x
	
fundel=randomwalk(len(points))
funthe=randomwalk(len(points))
funphi=randomwalk(len(points))
	


### jordan matrix model R(-t)M(d)R(t)
def fibmod(delta,theta,phi):
	R1=array([[cos(theta),sin(theta)],[-sin(theta),cos(theta)]])
	M=array([[e**(1j*(-delta/2)),0],[0,e**(1j*(delta/2))]])
	R2=transpose(R1)
	Ein=array([cos(pi/4)*e**(1j*(pi/2)),sin(pi/4)])*e**(1j*(phi))
	return R2.dot(M.dot(R1.dot(Ein)))

#####################################


t,rex,imx,rey,imy,Sq,Su,Sv=[],[],[],[],[],[],[],[]

for i in arange(len(points)):
	Eout=fibmod(fundel[i],funthe[i],funphi[i])
	Exr=Eout[0].real
	Exi=Eout[0].imag
	Eyr=Eout[1].real
	Eyi=Eout[1].imag
	rex.append(Exr)
	imx.append(Exi)
	rey.append(Eyr)
	imy.append(Eyi)
	
	t.append(i*4*10**(-7))

	Sq.append(Exr**2 + Exi**2 - Eyr**2 - Eyi**2)
	Su.append(2*(Exr*Eyr + Exi*Eyi))
	Sv.append(2*(Exi*Eyr + Exr*Eyi))


###########writing file################
outlist=list(zip(t,rex,imx,rey,imy))
f=open("eout.csv",'w')
w=csv.writer(f, delimiter='\t')
w.writerows(outlist)
f.close()

outphilist=list(zip(funphi))
f2=open("phidiff.csv",'w')
w=csv.writer(f2, delimiter='\t')
w.writerows(outphilist)
f2.close()

####plots of the simulated functions#####

fig1=figure() #new independent window
plot(fundel,'r',label="delta")
plot(funthe,'orange', label="theta")
plot(funphi,'g',label="phi")
#yticks(arange(0,7,pi/2),[r"$-\frac{\pi}{2}$", r"$-\frac{\pi}{4}$", r"$0$", r"$+\frac{\pi}{4}$",   r"$+\frac{\pi}{2}$"], fontsize=20)
#xlim(0,150)
#ylabel(("delta(n)", "theta"),fontsize=12,color='r')
legend()

u = np.linspace(0, 2 * np.pi, 100)
v = np.linspace(0, np.pi, 100)
x = 1 * np.outer(np.cos(u), np.sin(v))
y = 1 * np.outer(np.sin(u), np.sin(v))
z = 1 * np.outer(np.ones(np.size(u)), np.cos(v))



fig2=figure()
ax = fig2.gca(projection='3d') #get current axis
ax.scatter(Sq,Su,Sv,s=1)
ax.set_axis_off()
ax.plot_surface(x, y, z, color='grey',alpha=0.3)
show()





input()

