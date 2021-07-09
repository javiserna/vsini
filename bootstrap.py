import numpy as np
import random 
import scipy.stats
from scipy.special import erfinv
import matplotlib.pyplot as plt

def errfunction(mu,sigma,points):
	x=np.zeros(points)
	for i in range(0,points):
		x[i]=mu+sigma*np.sqrt(2)*erfinv(2*(random.uniform(0,1))-1)
	output=np.array(x)
	return output

#y=errfunction(4,2,1000)

#print y, np.mean(y), np.std(y)
#plt.figure()
#values, bins, _ = plt.hist(y, bins=100, normed=False, range=(0,10), color='blue')
#area = sum(np.diff(bins)*values)
#print area

#plt.plot(z, f(z), color='black')
#plt.xlabel(r'$x$')
#plt.ylabel(r'$PDF(x)$')
#plt.plot(bins[1:]-0.05, values/area, 'ok')
#plt.ylim(0,1)
#plt.figure()
#plt.plot(bins[1:]-0.05, values/area-f(bins[1:]-0.1),'-ok')
#plt.ylim(-0.3,0.3)
plt.show()
