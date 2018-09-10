


from numpy import *
from matplotlib.pyplot import *
import csv



###########----entries----###########
def fundel(n):
	return 1+5.2832*sin(2*pi*20*(4*10**(-7))*n)
	
def funthe(n):
	return 1+1*cos(2*pi*20*(4*10**(-7))*n)

def funphi(n):
	return 1.5+1*cos(2*pi*10*(4*10**(-7))*n)

points=arange(100000)

def fibmod(delta,theta,phi):
	R1=array([[cos(theta),sin(theta)],[-sin(theta),cos(theta)]])
	M=array([[e**(1j*(-delta/2)),0],[0,e**(1j*(delta/2))]])
	R2=transpose(R1)
	Ein=array([cos(pi/4)*e**(1j*(pi/2)),sin(pi/4)])*e**(1j*(phi))
	return R2.dot(M.dot(R1.dot(Ein)))

#####################################


t,rex,imx,rey,imy=[],[],[],[],[]

for i in range(1, len(points)+1):
	Eout=fibmod(fundel(i),funthe(i),funphi(i))
	rex.append(Eout[0].real)
	imx.append(Eout[0].imag)
	rey.append(Eout[1].real)
	imy.append(Eout[1].imag)
	
	t.append(i*4*10**(-7))


###########writing file################
outlist=list(zip(t,rex,imx,rey,imy))
f=open("eout.csv",'w')
w=csv.writer(f, delimiter='\t')
w.writerows(outlist)
f.close()

outphilist=list(zip(funphi(points)))
f2=open("phidiff.csv",'w')
w=csv.writer(f2, delimiter='\t')
w.writerows(outphilist)
f2.close()


####plots of the simulated functions#####

fig1=figure() #new independent window
plot(fundel(points),'r',label="delta")
plot(funthe(points),'orange', label="theta")
plot(funphi(points),'g',label="phi")
#yticks(arange(0,7,pi/2),[r"$-\frac{\pi}{2}$", r"$-\frac{\pi}{4}$", r"$0$", r"$+\frac{\pi}{4}$",   r"$+\frac{\pi}{2}$"], fontsize=20)
#xlim(0,150)
#ylabel(("delta(n)", "theta"),fontsize=12,color='r')
grid()
legend()
fig1.show()
input()

