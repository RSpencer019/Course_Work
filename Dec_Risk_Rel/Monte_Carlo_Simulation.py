import numpy as np
import scipy.stats


''' Log normal Distribution '''
Ki_median = 5 * 10 ** -8
Ki_del = 1.0					# coefficient of variation (c.o.v.)

Ki_lambda = np.log(Ki_median)
Ki_zeta = np.sqrt( np.log(1 + Ki_del**2) )

n_values = np.array([1,3,6])	# number of layers
limit = 10 ** -7



print ('\nPart (A-C)')

def g(x):						# function for equivalent hydraulic conductivity
	return (1 / ((1/n)*np.sum(1/x)))

def pvalue(limit, lambdas, zetas):
	z = (limit - lambdas) / zetas
	return scipy.stats.norm.sf(z)

def first_order_approx():
	COV = rho * np.dot( Ki_zetas.T, Ki_zetas )
	Keq_zeta = np.sqrt( np.sum( np.dot(partial_derivatives.T, partial_derivatives) * COV) )
	Keq_mean = np.exp(Keq_lambda + ( Keq_zeta**2 / 2 ) )

	print ('    Keq: Mean = ', Keq_mean)
	print ('    Keq: Std. Dev. = ', Keq_zeta)
	print ('    Probability of Exceedance: ', pvalue(limit, Keq_lambda, Keq_zeta))

	
	''' Part (c): # of Realizations for C.O.V. of 10% '''
	accuracy = 0.1   # 10% from mean
	z_10 = 1.96      # 95% CI
	Keq_del_10 = accuracy / z_10
	Keq_zeta_10 = np.sqrt( np.log(1 + Keq_del_10**2) )
	realizations = ( Keq_zeta / Keq_zeta_10 )**2
	print ('    Monte Carlo Estimated Realizations (+/-10%): ', int(realizations))



for n in n_values:
	print ('n = ',n,':')

	Ki_lambdas = np.empty([1,n])
	Ki_lambdas.fill(Ki_lambda)
	Ki_zetas = np.empty([1,n])
	Ki_zetas.fill(Ki_zeta)

	Keq_lambda = g(Ki_lambdas)

	perturbation = np.ones([1,n])
	perturbation[0,0] = 1.00001
	partial_derivative_dgdx = (g(perturbation*Ki_lambdas) - g(Ki_lambdas)) / (perturbation[0,0]*Ki_lambda - Ki_lambda)

	partial_derivatives = np.empty([1,n])
	partial_derivatives.fill(partial_derivative_dgdx)


	print ('  Case I:')

	rho = np.identity(n)		# correlation coefficient (statistically independent)
	first_order_approx()

	print ('  Case II:')					
	
	''' fill correlation coefficient of 0.5 '''
	for i in range(len(rho)):
	    for j in range(len(rho[i])):
	        if np.abs(i-j) == 1:	# if adjacent
	            rho[i,j] = 0.5
	first_order_approx()





print ('\nPart (D)')

realizations_total = 100
sets = 10
CI95_z = 1.96


def Monte_Carlo():
	Keq_lambdas = []
	Keq_zetas = []
	for set in range(sets):

		Keq_vals = []
		for realization in range(realizations_total):

			''' Generate random Ki values '''
			z_rand = np.random.multivariate_normal(np.zeros([n]), rho)	# Generate random z values for each layer
			Ki_vals_rand =  Ki_lambdas + z_rand * Ki_zetas
			Keq_vals.append( g(Ki_vals_rand) )

		''' Create list of means and std. dev. from realizations '''
		Keq_lambdas.append( np.mean( Keq_vals ) )
		Keq_zetas.append( np.sqrt( np.var( Keq_vals ) ) )
	
	''' Calculate the mean of all sample means '''
	Keq_lambda = np.mean( Keq_lambdas )	
	Keq_zeta = np.mean( Keq_zetas )
	Keq_mean = np.exp(Keq_lambda + ( Keq_zeta**2 / 2 ) )

	''' Generate 95% CI '''
	Keq_lower = Keq_lambda - CI95_z * Keq_zeta
	Keq_upper = Keq_lambda + CI95_z * Keq_zeta
	Keq_CI_mean = np.exp([Keq_lower, Keq_upper] + ( Keq_zeta**2 / 2 ) )

	print ('    Keq: Mean = ', Keq_mean)
	print ('    Keq: Std. Dev. = ', Keq_zeta)
	print ('    Probability of Exceedance: ', pvalue(limit, Keq_lambda, Keq_zeta))
	print ('    95% CI: ', Keq_CI_mean)
	

for n in n_values:
	print ('n = ',n,':')

	''' Initiate Ki matrices '''
	Ki_lambdas = np.empty([1,n])
	Ki_lambdas.fill(Ki_lambda)
	Ki_zetas = np.empty([1,n])
	Ki_zetas.fill(Ki_zeta)

	print ('  Case I:')
	rho = np.identity(n)		# correlation coefficient (statistically independent)
	Monte_Carlo()

	print ('  Case II:')					
	''' fill correlation coefficient of 0.5 for adjacent layers '''
	for i in range(len(rho)):
	    for j in range(len(rho[i])):
	        if np.abs(i-j) == 1:	# if adjacent
	            rho[i,j] = 0.5
	Monte_Carlo()












''' Estimate # of realizations '''
'''
Keq_lambda_zeta = np.sqrt(np.var(Keq_lambdas))
Keq_lambda_del = np.sqrt( np.exp( Keq_lambda_zeta ** 2 ) - 1 )
print ('    % from mean w/ 95%CI: ', Keq_lambda_del*100*1.96)
'''





