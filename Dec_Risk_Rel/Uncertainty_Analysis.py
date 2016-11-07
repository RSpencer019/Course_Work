import numpy as np
import scipy.stats

# Given Data

a = np.array([[
	1, 1, 1, 1, 1, -1, -1
	]])

median_x = np.array([[
	0.029, 0.1, 1.4, 235, 9, 70, 25550
	]])

del_x = np.array([[					# c.o.v.
	1, 2, 0.2, 0.5, 0.5, 0.25, 0.3
	]])

rho = np.array([					# correlation coefficient  (statistically independent)
	[1,0,0,0,0,0,0],
	[0,1,0,0,0,0,0],
	[0,0,1,0,0,0,0],
	[0,0,0,1,0,0,0],
	[0,0,0,0,1,0,0],
	[0,0,0,0,0,1,0],
	[0,0,0,0,0,0,1],
	])


# Log normal calculations

lambda_x = np.log(median_x)
zeta_x = np.sqrt( np.log(1 + del_x**2) )




#### (B) ####

# Random variable transformations

COV = rho * np.dot( zeta_x.T, zeta_x )

lambda_R = np.dot( a, lambda_x.T )[0,0]
zeta_R = np.sqrt( np.sum( COV * np.dot( a.T, a ) ) )

median_R = np.exp(lambda_R)
mean_R = np.exp(lambda_R + ( zeta_R**2 / 2 ) )
variance_R = zeta_R ** 2
del_R = np.sqrt( np.exp( zeta_R ** 2 ) - 1 )


# Output

print ('(B)')

print ('COV\n',COV)
print ('lambda_R\n    ',lambda_R)
print ('zeta_R\n    ',zeta_R)

print ('median_R\n    ',median_R)
print ('mean_R\n    ',mean_R)
print ('variance_R\n    ',variance_R)
print ('c.o.v.\n    ',del_R)




#### (C) ####

percent_variance = np.sum( COV * np.dot( a.T, a ) , axis=1) / np.sum( COV * np.dot( a.T, a ) )


# Output

print ('(C)')
print ('percent_variance\n    ',percent_variance)


#### (D) ####
vals = np.log( np.array([10**-6, 10**-5, 10**-4]) )
z = (vals - lambda_R) / zeta_R
p_values = scipy.stats.norm.sf(z)
print ('(D)')
print ('p_values\n    ',p_values)


#### (E) ####

# Given Data

rho = np.array([					# correlation coefficient (W, T modified)
	[1,0,0,0,0,0,0],
	[0,1,0,0,0,0,0],
	[0,0,1,0,0,0,0],
	[0,0,0,1,0,0,0],
	[0,0,0,0,1,0,0],
	[0,0,0,0,0,1,-0.6],
	[0,0,0,0,0,-0.6,1],
	])

# Random variable transformations

zeta_R_old = zeta_R
variance_R_old = variance_R
COV = rho * np.dot( zeta_x.T, zeta_x )
zeta_R = np.sqrt( np.sum( COV * np.dot( a.T, a ) ) )
variance_R = zeta_R ** 2
variance_perc_change = (variance_R-variance_R_old)/variance_R_old*100

# Output

print ('(E)')
print ('variance_R\n    ',variance_R)
print ('variance_R_difference\n    ',variance_perc_change,'%')




#### (F) ####								
print ('(F)')

# Given Data
n = 100
lambda_Rs = np.empty([1,n])
lambda_Rs.fill(lambda_R)
zeta_Rs = np.empty([1,n])
zeta_Rs.fill(zeta_R_old)			

a = np.empty([1,n])
a.fill(1/n)

# Statistically Independent
rho = np.identity(n)					# correlation coefficient (statistically independent)

# Random variable transformations
COV = rho * np.dot( zeta_Rs.T, zeta_Rs )

lambda_R_ave = np.dot( a, lambda_Rs.T )[0,0]
zeta_R_ave = np.sqrt( np.sum( COV * np.dot( a.T, a ) ) )

mean_R_ave = np.exp(lambda_R_ave + ( zeta_R_ave**2 / 2 ) )
del_R_ave = np.sqrt( np.exp( zeta_R_ave ** 2 ) - 1 )


print ('Statistically Independent')
print ('mean_R_ave\n    ',mean_R_ave)
print ('c.o.v.\n    ',del_R_ave)

# Percent Exceedance
vals = np.log( np.array([10**-6, 10**-5, 10**-4]) )
z = (vals - lambda_R_ave) / zeta_R_ave
p_values = scipy.stats.norm.sf(z)
print ('p_values\n    ',p_values)


# Perfectly Correlated
rho = np.ones([n,n])					# correlation coefficient (perfectly correlated)

# Random variable transformations
COV = rho * np.dot( zeta_Rs.T, zeta_Rs )

lambda_R_ave = np.dot( a, lambda_Rs.T )[0,0]
zeta_R_ave = np.sqrt( np.sum( COV * np.dot( a.T, a ) ) )

mean_R_ave = np.exp(lambda_R_ave + ( zeta_R_ave**2 / 2 ) )
del_R_ave = np.sqrt( np.exp( zeta_R_ave ** 2 ) - 1 )

print ('Perfectly Correlated')
print ('mean_R_ave\n    ',mean_R_ave)
print ('c.o.v.\n    ',del_R_ave)


# Percent Exceedance
vals = np.log( np.array([10**-6, 10**-5, 10**-4]) )
z = (vals - lambda_R_ave) / zeta_R_ave
p_values = scipy.stats.norm.sf(z)
print ('p_values\n    ',p_values)




#### (G) ####

n = 100
U_n = scipy.stats.norm.ppf(1-1/n) * zeta_R_old + lambda_R
alpha_n = n * ( 1 / ( np.sqrt(2*np.pi) * zeta_R_old ) ) * np.exp( - (1/2) * (((U_n - lambda_R) / zeta_R_old)**2))
p_values = 1 - np.exp(-np.exp(-alpha_n*(vals-U_n)))
print ('(G)')
print ('Statistically Independent')
print ('p_values\n    ',p_values)

print ('Perfectly Correlated')
z = (vals - lambda_R) / zeta_R_old
p_values = scipy.stats.norm.sf(z)
print ('p_values\n    ',p_values)
