import matplotlib.pyplot as plt #for plotting 
import numpy as np #advanced math and othe useful stuff (arrays, matrixes)
from scipy.integrate import quad #numerical integration 
from lmfit import Model #good for regression 
from scipy import special #special functions
import sys #to exit if input parameters are incorrect

#a couple of functions to simplify the code below
def mean(x):
	result = 0 
	for dot in x: result += dot
	result = result / len(x)
	return result
	
def mean_weights(x, weights):
	result = 0
	for i in range(len(x) -1):
		result += weights[i] * (x[i+1] - x[i]) * x[i] 
	return result
	
def gauss(x, s, mu):
	result = 1 - special.erfc((x - mu) / s) / 2
	return result
	
def noised_model(x, s, n, mu_s, mu_n):
	k = abs(s * n) / ((s**2 + n**2)/2) ** 0.5
	mu = mu_s + mu_n
	v = (mu * s**2 + x * n**2) / (s**2 + n**2)
	A = - v ** 2 + (mu**2 * s**2 + x**2 * n**2) / (s**2 + n**2)
	result = k / ((np.pi)**(0.5) * abs(s * n)) * np.exp( - A / k**2) * special.erfc( (mu - v) / k) / 2
	return result

def noised_model_integral(a, b, s, n, mu_s, mu_n):
	result = quad(noised_model, a, b, args = (s, n, mu_s, mu_n))[0]
	return result
	
def model_tofit(x, s, n, mu_s, mu_n):
	result = np.zeros((len(x),))
	result[0] = 0
	for i in range(1,len(x)):
		result[i] += result[i-1] +  noised_model_integral(x[i-1], x[i], s=s, n=n, mu_s = mu_s, mu_n = mu_n)
	return result

#we can choose which year do we use
str16 = 'C2015-16-flicker.dat'	
str17 = 'C2016-17-flicker.dat'
sigma_2016 = 0.045
sigma_2017 = 0.093

print('choose a year:')
year = int(input())
if(year == 16):
	filename = str16
	sigma_shift = sigma_2016
elif(year == 17):
	filename = str17
	sigma_shift = sigma_2017
else:
	print('you mistyped something, oops')
	sys.exit()
	
#opening file
dots_phase = []
dots_mag = []
f_in = open(filename, 'r')
for line in f_in:
	dot = [float(j) for  j in line.strip().split()]
	dots_phase.append(dot[0])
	dots_mag.append(dot[1])
f_in.close()

#shifting data by sigma_shift
m = mean(dots_mag)
for i in range(len(dots_mag)):
	dots_mag[i] += -m+sigma_shift
N = len(dots_phase)

#sorting magnitides
flickering_sorted = np.sort(dots_mag)

probability = np.linspace(0. + 1. / len(flickering_sorted), 1., len(flickering_sorted))

#plotting probability(magnitude)
plt.plot(flickering_sorted, probability, 'ro', markersize = 1)
plt.show()

'''
#testing distribution

x = np.linspace(-3, 15, 1000)
s = 1. # noise
n = 1. # signal 
mu_s = 0.0 # noise 
mu_n = 0.0 # signal 
test = noised_model(x = x, s = s, n = n, mu_s = mu_s, mu_n = mu_n)
print(mean_weights(x = x, weights = test))
plt.plot(x, test)
plt.show()

integrated_test = model_tofit(x = x, s = s, n = n, mu_s = mu_s, mu_n = mu_n)
plt.plot(x, integrated_test)
plt.show()
'''


#fitting distribution to my points:
main_model = Model(model_tofit)
result = main_model.fit(probability, x = flickering_sorted, s = 0.05, n = 0.05, mu_s = 0.05, mu_n = -0.05)
print(result.fit_report())
plt.plot(flickering_sorted, result.best_fit)
plt.plot(flickering_sorted, probability, 'ro', markersize = 1)
plt.show()

N_chi2 = 180
chi2 = np.zeros((N_chi2,))
chi2_reduced = np.zeros((N_chi2,))
chi2_tobe = np.zeros((N_chi2,))
L = 4
N_start = 200 #number of points in one bin to start
#print(len(flickering_sorted))
DOF = []
chi_95 = [3.8415, 5.9915, 7.8147, 9.4877, 11.0705, 12.5916, 14.0671, 15.5073, 16.9190, 18.3070, 19.6751, 21.0261, 22.3620, 23.6848, 24.9958, 26.2962, 27.5871, 28.8693, 30.1435, 31.4104, 32.6706, 33.9244, 35.1725, 36.4150, 37.6525, 38.8851, 40.1133, 41.3371, 42.5570, 43.7730, 44.9853, 46.1943, 47.3999, 48.6024, 49.8018, 50.9985, 52.1923, 53.3835, 54.5722, 55.7585]
chi_50 = [0.4549, 1.3863, 2.3660, 3.3567, 4.3515, 5.3481, 6.3458, 7.3441, 8.3428, 9.3418, 10.3410, 11.3403, 12.3398, 13.3393, 14.3389, 15.3385, 16.3382, 17.3379, 18.3377, 19.3374, 20.3372, 21.3370, 22.3369, 23.3367, 24.3366, 25.3365, 26.3363, 27.3362, 28.3361, 29.3360, 30.3359, 31.3359, 32.3358, 33.3357, 34.3356, 35.3356, 36.3355, 37.3355, 38.3354, 39.3353]
n = np.arange(1, 41)


for j in range(N_chi2):
	#print(j)
	N_points = -j + N_start #number of points in one bin
	N_bins = int(len(flickering_sorted) / N_points) #number of bins for chi2
	DOF.append(N_bins - L - 1)
	#print(N_bins) 

	#calculating theoretical freqs
	freqs = np.zeros((N_bins,))
	for i in range(N_bins - 1):
		freqs[i] = len(flickering_sorted) * (result.best_fit[(i+1) * N_points - 1] - result.best_fit[i * N_points])
		
	#missing part
	freqs[N_bins - 1] = len(flickering_sorted) * (1 -  result.best_fit[N_points * (N_bins - 1)])
	N_missing = len(flickering_sorted) - N_points * (N_bins - 1)
	
	if (N_points == 50):
		plt.plot(freqs, 'ro')
		plt.plot(50 * np.ones(len(freqs)), 'b')
		plt.show()
	
	
	#chi2 itself
	for i in range(N_bins - 1):
		chi2[j] += (N_points - freqs[i])**2  / freqs[i]
	if (N_missing > 5):
		chi2[j] += (N_missing - freqs[N_bins - 1]) **2 / freqs[N_bins - 1]

#plotting chi2 (DOF)	
plt.plot(DOF, chi2, 'ro', markersize = 2, label = 'real data')
plt.plot(n, chi_95, 'b', label = '5%')
plt.plot(n, chi_50, 'g', label = '50%')
plt.title('chi^2 ')
plt.ylabel('chi^2')
plt.xlabel('DOF = N_bins - L - 1')
plt.legend()
plt.show()


#optimal chi2
if (year == 16):
	N = 128
	N_bins = 8
else:
	N = 251
	N_bins = 9
	
freqs = np.zeros((N_bins,))
for i in range(N_bins):
	freqs[i] = len(flickering_sorted) * (result.best_fit[(i+1) * N - 1] - result.best_fit[i * N])

plt.plot(freqs, 'ro')
plt.plot(N * np.ones(len(freqs)), 'b')
plt.show()

chi2 = 0
for i in range(N_bins):
		chi2 += (N - freqs[i])**2  / freqs[i]
		
print('optimal chi2 = ' + str(chi2))

