####################################################
######  study of random-noise fourier spectrum #####
####################################################




from pylab import *
from scipy import signal

########-----entries-------######

### case 1)
### deterministic pure sine function
Fs = 1000        # Sampling frequency                    
T = 1/Fs            # Sampling period       
L = 1500        # Length of signal
t = (arange(L))*T       # Time vector
#S = 1*sin(2*pi*50*t) + sin(2*pi*120*t)


filtering=0

### case 2)
### random noise
### launch "FUNDEF_rand.py" in the working directory
S=genfromtxt("phidiff.csv")
Fs=1

#########filtering################

if filtering :
	fc=70 #cut frequency
	b,a=signal.butter(3, fc/(Fs/2), 'low')
	Sf=signal.lfilter(b, a, S)
	Sfl=signal.filtfilt(b, a, S) #phaselinearfilt

	fig1=figure()
	w, h = signal.freqz(b, a)
	plot(w*(Fs/(2*pi)), abs(h))
	grid()
	Ysf = fft(Sf)
	Ysfl = fft(Sfl)
	Ysfmod,Ysflmod =abs(Ysf),abs(Ysfl)


####----FFT---------########

fig2=figure()
compare=((arange(len(S))).astype(float))**(-1)
Ys = fft(S)
Ysmod=abs(Ys)
maxmod=max(Ysmod)
fvec=(arange(len(S)))*(Fs/len(S)) 
plot(fvec,Ysmod/maxmod)
plot(fvec,compare,c=(1,0,0))
loglog()
print(compare)
if filtering :

	plot(fvec,Ysfmod)
	plot(fvec,Ysflmod,color='r')
grid()

show()



