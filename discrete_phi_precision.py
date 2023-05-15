import mpmath
import time
import sys

start_time = time.time()

'''
Bi-seasonal discrete time risk model with income rate two.
This code is used to count discrete time survival probability
'''

# Set the desired precision to 50 digits
mpmath.mp.dps = 50

# Poisson pmf:
def poisson_pmf(k, lambd, xi):
    if k<xi: return mpmath.mpf('0')
    return mpmath.power(lambd, k-xi) * mpmath.exp(-lambd) / mpmath.factorial(k-xi)

# Poisson cdf:
def poisson_cdf(k, lambd, xi):
    sum=0
    for i in range(k+1):
        sum += poisson_pmf(i, lambd, xi)
    return sum

# Define the parameters
x_lambda = mpmath.mpf('2')
x_xi = mpmath.mpf('1')
y_lambda = mpmath.mpf('1')
y_xi = mpmath.mpf('1')

# Define big enough value for aprox
t_max=50
# What is the maximum u?
u_max=50

# Initialize the arrays
X_U = [mpmath.mpf('0') for i in range(0, u_max+8)]
Y_U = [mpmath.mpf('0') for i in range(0, u_max+8)]

x_n = [mpmath.mpf('0') for i in range(0, u_max+8)]
y_n = [mpmath.mpf('0') for i in range(0, u_max+8)]
s_n = [mpmath.mpf('0') for i in range(0, u_max+8)]

# Perform the Poisson cdf, pmf calculations using mpmath:
for i in range(u_max+8):
    X_U[i] = poisson_cdf(i, x_lambda, x_xi)
    Y_U[i] = poisson_cdf(i, y_lambda, y_xi)

for i in range(u_max+8):
    x_n[i] = poisson_pmf(i, x_lambda, x_xi)
    y_n[i] = poisson_pmf(i, y_lambda, y_xi)
    s_n[i] = poisson_pmf(i, x_lambda + y_lambda, x_xi + y_xi)

# Initializing matrix Phi_U_T
PHI_U_T = mpmath.ones(t_max, u_max + 9)

# Counting Phi(u, 1) values:
for i in range(u_max+4):
    PHI_U_T[0, i] = X_U[i+1]

# Counting Phi(u, 2) values:
for i in range(u_max+4):
    vid=mpmath.mpf('0')
    for k in range(i+2):
        vid=vid+x_n[k]*Y_U[i+3-k]
    PHI_U_T[1, i] = vid

# Counting Phi(u, T), T>=3 values:
for j in range(2, t_max):
    for i in range(u_max+1):
        vid=mpmath.mpf('0')
        for k in range(i+4):
            vid=vid+PHI_U_T[j-2, i+4-k]*s_n[k]
            
        PHI_U_T[j, i] = vid-(x_n[i+2]*y_n[0]*PHI_U_T[j-2, 2])-\
                (x_n[i+2]*y_n[1]+x_n[i+3]*y_n[0])*PHI_U_T[j-2, 1]


PHI_U_T_rounded = mpmath.matrix(PHI_U_T.rows, PHI_U_T.cols-8)
# convert each element of A to a string with 3 decimal places
for i in range(PHI_U_T.rows):
    for j in range(PHI_U_T.cols-8):
        PHI_U_T_rounded[i, j] = mpmath.nstr(PHI_U_T[i, j], 3, min_fixed=-4)

# Write the rounded matrix to a text file
with open("discrete_time.txt", "w") as f:
    f.write(str(PHI_U_T_rounded))

# Code execution timer
end_time = time.time()
print("Time taken: {:.6f} seconds".format(end_time - start_time))