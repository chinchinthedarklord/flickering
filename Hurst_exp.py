import matplotlib.pyplot as plt #for plotting
import numpy as np #advanced math and othe useful stuff (arrays, matrixes)
from lmfit import Model #good for regression 
import sys

#choosing model, it's a straight line log R/S - log tau
def linear(x, A, B):
	return A * x + B

linear_model = Model(linear)

#a couple of functions to simplify the code below
def mean(x):
	result = 0 
	for dot in x: result += dot
	result = result / len(x)
	return result
	
def sigma(x):
	result = 0
	meanx = mean(x) 
	for dot in x: result += dot * dot
	result = result / len(x) - meanx * meanx
	return result

#we can choose which year do we use
str16 = 'C2015-16-flicker.dat'	
str17 = 'C2016-17-flicker.dat'
sigma_2016 = 0.045
sigma_2017 = 0.093

print('choose a year:')
n = int(input())
if(n == 16):
	filename = str16
	sigma_shift = sigma_2016
elif(n == 17):
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
	dots_mag[i] += (-m+sigma_shift)
N = len(dots_phase)

#test, we plot our data to see that everything is fine
'''
plt.plot(dots_phase, dots_mag, 'ro', markersize = 1)
plt.plot(dots_phase, mean(dots_mag) *np.ones(N), 'b', markersize = 1)
plt.show()
'''

tau_start = 0.11
tau_end = 0.50
tau_amount = int(100 * (tau_end - tau_start) +1)	
tau = np.linspace(tau_start, tau_end, tau_amount)
RS_mean = np.zeros((tau_amount - 1,))

#calculating R/S (tau)
i = 0
for t in tau:
	shift = 0.01 # shifting phase
	N_sets = int((tau_end - t) / shift) # amount of shifts
	RS_datasets = np.zeros((N_sets,))
	if (N_sets > 0):
		for j in range(N_sets):
			# setting the begin end the end of current dataset
			index_start = 0
			index_end = 0
			while (dots_phase[index_start] < j * shift): 
				index_start += 1
			while (dots_phase[index_end] < t + j * shift): 
				index_end += 1
			dataset = dots_mag[index_start : index_end]
			mean_dataset = mean(dataset)
			X = np.zeros((len(dataset),))
			X[0] = dataset[0] - mean(dataset)
			for k in range(1, len(dataset)):
				X[k] = X[k-1] + (dataset[k] - mean(dataset))
			R = max(X) - min(X)
			S = (sigma(dataset) * len(dataset) / t)** 0.5
			RS_datasets[j] = R / S
		RS_mean[i] = mean(RS_datasets)
		#print(RS_mean[i], t, i) 
		i += 1
		
#test, we plot our data to see that everything is fine

plt.plot(np.log(tau[: len(tau) - 1]), np.log(RS_mean), label = 'raw data' )	
plt.xlabel('log tau')
plt.ylabel('log R/S')
plt.title('log R/S - log tau diagram ')
plt.legend()
plt.legend()	
plt.show()
	
#cutting of part, where data is not a couple of points 		
index_max = 0
max_elem = RS_mean[0]
for i in range(len(RS_mean)):
	if (RS_mean[i] > max_elem):
		max_elem = RS_mean[i]
		index_max = i
gap = 5
		
log_t = np.log(tau[: index_max])
log_RS = np.log(RS_mean[: index_max])

#does it have some fissure? (is there different H on different scales?)
print('do we need more than one H? (enter 1)')
answer = int(input())
if (answer == 1):
	print('log tau of the fissure:')
	fissure = float(input())
	index_fissure = int( (-tau_start +np.exp(fissure)) / shift)
else:
	index_fissure = index_max

log_t1 = np.log(tau[: index_fissure])
log_RS1 = np.log(RS_mean[: index_fissure])

#fit and its results
result1 = linear_model.fit(log_RS1, x = log_t1, A = 0., B = 0.)	
print(result1.fit_report())
#value of H +/- dH
H1 = float( '{:.3f}'.format(result1.params.get('A').value) )
dH1 = float( '{:.3f}'.format(result1.params.get('A').stderr) )

#plotting
plt.plot(log_t1, log_RS1, 'b', label = 'raw data' )	
plt.plot(log_t1, result1.best_fit, 'k--', label = 'regression, H = ' + str(H1) + ' +/- ' + str(dH1))

if (index_fissure != index_max):
	log_t2 = np.log(tau[index_fissure-1: index_max- gap])
	log_RS2 = np.log(RS_mean[index_fissure-1: index_max - gap])
	
	#fit and its results
	result2 = linear_model.fit(log_RS2, x = log_t2, A = 0., B = 0.)	
	print(result2.fit_report())
	#value of H +/- dH
	H2 = float( '{:.3f}'.format(result2.params.get('A').value) )
	dH2 = float( '{:.3f}'.format(result2.params.get('A').stderr) )
	
	plt.plot(log_t2, log_RS2, 'b')	
	plt.plot(log_t2, result2.best_fit, 'k--', label = 'regression, H = ' + str(H2) + ' +/- ' + str(dH2))

#to make graph pretty
plt.xlabel('log tau')
plt.ylabel('log R/S')
plt.title('log R/S - log tau diagram ')
plt.legend()
plt.show()

		
	
			