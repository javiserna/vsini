#!/home/anaconda3/bin/env python3
# -*- coding: utf-8 -*-
import os
import sys
import pandas as pd
import matplotlib.pyplot as plt
import scipy
from scipy.fftpack import fft,fftfreq
from scipy.signal import fftconvolve
from scipy.interpolate import interp1d
from scipy.optimize import curve_fit
from datetime import date, datetime
from scipy import stats
from lmfit import Model

from PyQt5 import QtWidgets

import fourier
import numpy as np

from rotation import Ui_rotation


class MyApp(QtWidgets.QMainWindow, Ui_rotation):

	def __init__(self):

		QtWidgets.QMainWindow.__init__(self)
		Ui_rotation.__init__(self)
		self.setupUi(self)
		self.progressBar.setMaximum(100)
		self.progressBar.setValue(0)
		self.buttom1.clicked.connect(self.getCSV)
		self.buttom2.clicked.connect(self.plotspec)
		self.buttom4.clicked.connect(self.plotline)
		self.buttom3.clicked.connect(self.run)

	def Increase_step(self):

		self.progressBar.setValue(self.progressBar.value()+1)

	def getCSV(self):

		global filePath
		filePath, _ = QtWidgets.QFileDialog.getOpenFileName(self, 'Open file', './')
		if filePath != "":
			self.df = pd.read_csv(str(filePath))

	def plotspec(self, event):

		global x, y
		x=self.df['col1']
		y=self.df['col2']
		spec=plt.figure(1)
		plt.plot(x,y)
		plt.xlabel('wavelength $(\AA)$')
		plt.ylabel('Flux')

		L=np.arange(min(x), max(x), x[1]-x[0])  
		F=interp1d(x,y)
		Z=F(L)

		clicks=[]

		def gaussian_line(x, amp, cen, wid, slope, intercept):
			"""1-d gaussian: gaussian+line(amp, cen, wid, slope, intercept)"""
			return ((amp / (np.sqrt(2*np.pi) * wid)) * np.exp(-(x-cen)**2 / (2*wid**2)))+ slope*x+intercept

		def onclick(event):

			ix=event.xdata

			if (event.dblclick == 1) and (event.button == 1):

				clicks.append(ix)
				fig=plt.axvline(ix, color='r')
				plt.draw()

			if len(clicks) == 2:

				coord=np.array(clicks)
				width=abs(coord[1]-coord[0])
				center=min(coord)+(width/2)
				amp=(((F(coord[1])-F(coord[0]))/2)-F(center))
				inter=F(coord[0])
				slope=0
				window=((center-(width/2))<=x) & (x<=(center+(width/2)))
				wav=x[window]
				flux=y[window]

				gmodel = Model(gaussian_line)
				result = gmodel.fit(flux, x=wav, amp=amp, cen=center, wid=width/10, slope=slope, intercept=inter)
				dely = result.eval_uncertainty(sigma=3)
				#param=result.fit_report()
				#print(param)

				popt_gauss, pcov_gauss = scipy.optimize.curve_fit(gaussian_line, wav, flux, p0=[amp, center, width/10, slope, inter])
				perr_gauss = np.sqrt(np.diag(pcov_gauss))

				plt.plot(wav, result.best_fit, 'k-', label='best fit')
				plt.fill_between(wav, result.best_fit-dely, result.best_fit+dely, color="grey", alpha=0.2, label='3-$\sigma$ uncertainty band')
				self.linecenter.setText(str(popt_gauss[1])[:7])
				self.linewidth.setText(str(6*popt_gauss[2])[:4]) #gaussian fit to 3 sigma level
				spec.canvas.mpl_disconnect(plot)
				
		plot=spec.canvas.mpl_connect('button_press_event', onclick)
		plt.show()

	def run(self):

		f=open("Fourier.out", "a")
		count = len(open("Fourier.out", "r").readlines())

		if count==0:
			f.write("#file date time midwave width resolution epsilon vsini vsini_err\n")

		cc = 2.99792458e+18
		midwave=float(self.linecenter.text())
		epsilon=float(self.limbdarkening.text())
		width=float(self.linewidth.text())
		dlam=midwave/(2*float(self.resolution.text())) #nyquist criteria
		vsini=[]

		for k in range(1001):
			while True:
				try:
					vsini.append(fourier.Fourier(x, y, midwave, width, dlam, epsilon))
					self.progressBar.setValue(int(k/10.))
					break
				except:
					pass

		vsini=np.array(vsini)
		#vrot="vsini "+str(np.nanmedian(vsini))[:5]+"±"+str(np.nanmedian(np.abs(vsini-np.nanmedian(vsini))))[:4] +"  ("+str(np.nanmedian(np.absolute(vsini - np.nanmedian(vsini)))*100/np.nanmedian(vsini))[:4]+"%)"
		vrot="vsini "+str(np.nanmedian(vsini))[:5]+"±"+str(np.nanstd(vsini))[:4] +"  ("+str(np.nanstd(vsini)*100/np.nanmedian(vsini))[:4]+"%)"
		self.resultado.setText(vrot)
		plt.hist(vsini, bins="scott", facecolor='blue', alpha=0.5) #plot histogram of vsini measures
		plt.xlabel("vsini (km/s)")
		plt.show()

		now = datetime.now()
		dt_string = now.strftime("%d/%m/%Y %H:%M:%S")
	
		f.write("%s %s %.2f %.2f %f %.1f %.2f %.2f\n" %(os.path.basename(filePath), dt_string, midwave, width, float(self.resolution.text()), epsilon, np.nanmedian(vsini), np.nanstd(vsini)))
		f.close()

	def plotline(self):

		midwave=float(self.linecenter.text())
		epsilon=float(self.limbdarkening.text())
		width=float(self.linewidth.text())
		dlam=midwave/(2*float(self.resolution.text())) #nyquist criteria

		z=np.arange(0, width + dlam, dlam)[1:]  
		w=interp1d(x-midwave+(width/2),y)
		q=w(z)

		x1=min(z)
		x2=max(z)
		y1=q[0]
		y2=q[-1]

		slope=(y2-y1)/(x2-x1)

		spec_obs=q-(slope*z)
		spec_obs=spec_obs/(max(spec_obs))

		plt.plot(z,spec_obs)
		plt.xlabel(r'$\Delta \lambda (A)$')
		plt.ylabel('Flux Normalized')
		plt.show()

if __name__ == "__main__":
	app =  QtWidgets.QApplication(sys.argv)
	window = MyApp()
	window.show()
	sys.exit(app.exec_())

