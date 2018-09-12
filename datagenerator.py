
############# SIMULATOR OF THE COHERENT RECEIVER CHANNELS ############
####         main function input: 3 functions, 2 integers                    
#### 			   output: 4-channels data file, parameters plot,
####				   poincare' sphere 


from numpy import *
import csv
from matplotlib.pyplot import *
from mpl_toolkits.mplot3d import Axes3D

def datagen(fundel,funthe,funphi, points,sphere):

	### jordan matrix model R(-t)M(d)R(t)
	def fibmod(delta,theta,phi):
		R1=array([[cos(theta),sin(theta)],[-sin(theta),cos(theta)]])
		M=array([[e**(1j*(-delta/2)),0],[0,e**(1j*(delta/2))]])
		R2=transpose(R1)
		Ein=array([cos(pi/4)*e**(1j*(pi/2)),sin(pi/4)])*e**(1j*(phi))
		return R2.dot(M.dot(R1.dot(Ein)))

	#####################################

	t,rex,imx,rey,imy,Sq,Su,Sv,dellist,thelist,philist=[],[],[],[],[],[],[],[],[],[],[]
	dd,pp,tt=0,0,0


	for i in arange(points):
		dd=fundel(i)
		tt=funthe(i)
		pp=funphi(i)
		dellist.append(dd)
		thelist.append(tt)
		philist.append(pp)

		Eout=fibmod(dd,tt,pp)
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
		Sv.append(2*(Exi*Eyr - Exr*Eyi))


	###########writing file################
	outlist=list(zip(t,rex,imx,rey,imy))
	f=open("eout.csv",'w')
	w=csv.writer(f, delimiter='\t')
	w.writerows(outlist)
	f.close()

	outphilist=list(zip(philist))
	f2=open("phidiff.csv",'w')
	w=csv.writer(f2, delimiter='\t')
	w.writerows(outphilist)
	f2.close()

	####plots of the simulated functions#####

	fig1=figure() #new independent window
	plot(dellist,'r',label="delta")
	plot(thelist,'orange', label="theta")
	plot(philist,'g',label="phi")
	#yticks(arange(0,7,pi/2),[r"$-\frac{\pi}{2}$", r"$-\frac{\pi}{4}$", r"$0$", r"$+\frac{\pi}{4}$",   r"$+\frac{\pi}{2}$"], fontsize=20)
	#xlim(0,150)
	#ylabel(("delta(n)", "theta"),fontsize=12,color='r')
	legend()

	if sphere :
		u = np.linspace(0, 2 * np.pi, 100)
		v = np.linspace(0, np.pi, 100)
		x = 1 * np.outer(np.cos(u), np.sin(v))
		y = 1 * np.outer(np.sin(u), np.sin(v))
		z = 1 * np.outer(np.ones(np.size(u)), np.cos(v))



		fig2=figure()
		ax = fig2.add_subplot(111, projection='3d')
		 #get current axis
		ax.scatter(Sq,Su,Sv,s=1)
		#ax.set_axis_off()
		ax.plot_surface(x, y, z, rstride=4, cstride=4, color='grey',alpha=0.3)
	show()


