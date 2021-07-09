import scipy

from scipy.fftpack import fft,fftfreq
from scipy.signal import fftconvolve

#import matplotlib
#matplotlib.use('Agg')

import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import interp1d
import bootstrap


#midwave = 6709.14 # central wavelength of line
# velocity of light in angstrom/s
#epsilon = 0.6 # linear limbdarkening parameter
#dlam = 0.08 # for Hectochelle spectrum resolution = 0.2 (A) at nyquist 0.1 (A) at sampling 0.08 (A)
#width=0.9 #line width


def Fourier(l,f,midwave,width,dlam,epsilon):
	midwave=midwave+bootstrap.errfunction(0, dlam, 1)
	width=width+bootstrap.errfunction(0, dlam, 1)
	epsilon=epsilon+bootstrap.errfunction(0, 0.1, 1)
	cc = 2.99792458e+18
	x=np.arange(0, width + dlam, dlam)[1:]  
	y=interp1d(l-midwave+(width/2),f)
	z=y(x)

	x1=min(x)
	x2=max(x)
	y1=z[0]
	y2=z[-1]

	slope=(y2-y1)/(x2-x1)

	spec_obs=z-(slope*x)
	spec_obs=spec_obs/(max(spec_obs))

	spec_to_transform = (1-spec_obs)
	new_n = 100*len(spec_to_transform)
	spec_fft = np.abs(fft(spec_to_transform,n=new_n))**2

	x_fft = fftfreq(len(spec_fft),d=dlam)
	keep = x_fft>=0
	x_fft, spec_fft = x_fft[keep], spec_fft[keep]

	neg_to_pos = (np.diff(spec_fft[:-1])<=0) & (np.diff(spec_fft[1:])>=0)
	minima = x_fft[1:-1][neg_to_pos]

	q1 = 0.610 + 0.062*epsilon + 0.027*epsilon**2 + 0.012*epsilon**3 + 0.004*epsilon**4

	vsini_first_zero=cc/midwave*q1/minima[0]/1e13

	#plt.figure()
	#plt.plot(x,1-spec_obs)
	#plt.xlabel(r'$\Delta \lambda (A)$')
	#plt.ylabel('1-Flux (Counts)')
	#plt.title('Inverted Line')
	#plt.figure()
	#plt.plot(l,1-f)
	#plt.xlabel('wavelength $(\AA)$')
	#plt.ylabel('1-Flux (Counts)')
	#plt.title('Inverted Spectra')
	#plt.figure()
	#plt.plot(x_fft,spec_fft/spec_fft.max())
	#plt.xlabel('$\sigma\,(\AA^{-1})$')
	#plt.ylabel('|Power FFT|^2')
	#plt.gca().set_yscale('log')
	#plt.xlim(0.0,5)
	#plt.ylim(1e-7,1e1)
	#plt.axvline(x=minima[:1], color='r', alpha=0.5)
	#plt.text(1,1,'$v\sin(i)=$%.2f km/s' %(vsini_first_zero))
	#plt.title('Power FFT')
	
	#plt.savefig("Fourier_%0.3f.png" %midwave)

	return vsini_first_zero


