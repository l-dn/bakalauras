import mpmath
import time
import sys

start_time = time.time()

'''
Bi-seasonal discrete time risk model with income rate two.
This code is used to count ultimate time survival probability
only when !!! s0 = s1 = s2 = 0, s3 > 0 !!!
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
x_lambda = mpmath.mpf('1/2')
x_xi = mpmath.mpf('2')
y_lambda = mpmath.mpf('1/3')
y_xi = mpmath.mpf('1')

# Define big enough value for aprox
n_max=150
# What is the maximum u?
u_max=25

# Initialize the arrays
X_U = [mpmath.mpf('0') for i in range(0, 3)]
Y_U = [mpmath.mpf('0') for i in range(0, 3)]
S_U = [mpmath.mpf('0') for i in range(0, 3)]

x_n = [mpmath.mpf('0') for i in range(0, n_max+3)]
y_n = [mpmath.mpf('0') for i in range(0, n_max+3)]
s_n = [mpmath.mpf('0') for i in range(0, n_max+3)]

# Perform the Poisson cdf, pmf calculations using mpmath:
for i in range(3):
    X_U[i] = poisson_cdf(i, x_lambda, x_xi)
    Y_U[i] = poisson_cdf(i, y_lambda, y_xi)
    S_U[i] = poisson_cdf(i, x_lambda + y_lambda, x_xi + y_xi)

for i in range(n_max+3):
    x_n[i] = poisson_pmf(i, x_lambda, x_xi)
    y_n[i] = poisson_pmf(i, y_lambda, y_xi)
    s_n[i] = poisson_pmf(i, x_lambda + y_lambda, x_xi + y_xi)

# Check first 4 s(i) values
print('s(0) = ', s_n[0])
print('s(1) = ', s_n[1])
print('s(2) = ', s_n[2])
print('s(3) = ', s_n[3])

print()

# 3 Different scenarios:
Scenario_1=False
Scenario_2=False
Scenario_3=False
Scenario_4=False

if (x_n[0]==0 and y_n[0]==0 and x_n[1]==0 and y_n[1]>0 and x_n[2]>0):
    Scenario_1=True
elif (x_n[0]==0 and y_n[0]==0 and y_n[1]==0 and x_n[1]>0 and y_n[2]>0):
    Scenario_2=True
elif (x_n[0]==0 and y_n[0]>0 and x_n[1]==0 and x_n[2]==0 and x_n[3]>0):
    Scenario_3=True
elif (x_n[0]>0 and y_n[0]==0 and y_n[1]==0 and y_n[2]>0 and y_n[3]>0):
    Scenario_4=True
else:  
    print("make sure: s0 = s1 = s2 = 0, s3 > 0, else use different code")
    
    sys.exit()

# Initialize Phi_U vector:
Phi_U = mpmath.matrix(u_max + 1, 1)

# Count starting Phi_U values:
if (Scenario_2 or Scenario_4):
    # Count Phi_U(0)
    Phi_U[0] = mpmath.mpf('4')-(x_lambda+y_lambda+x_xi+y_xi)
    # Count Phi_U(1)
    Phi_U[1] = Phi_[0] / (x_n[1]*y_n[2] + x_n[0]*y_n[3])

elif (Scenario_1):
    # Count Phi_U(0)
    Phi_U[0] = mpmath.mpf('0')
    # Count Phi_U(1)
    Phi_U[1] = (mpmath.mpf('4')-(x_lambda+y_lambda+x_xi+y_xi)) / y_n[1]

elif (Scenario_3):
    # Count Phi_U(0)
    Phi_U[0] = mpmath.mpf('0')
    # Count Phi_U(1)
    Phi_U[1] = mpmath.mpf('0')
    # Count Phi_U(2)
    Phi_U[2] = (mpmath.mpf('4')-(x_lambda+y_lambda+x_xi+y_xi)) / y_n[0]

# Remaining Phi_U
if (Scenario_1 or Scenario_2 or Scenario_4):
    for i in range(2, u_max+1):
        vid=mpmath.mpf('0')
        for k in range(1, i):
            vid=vid+s_n[i+3-k]*Phi_U[k]
        Phi_U[i]=(1/s_n[3])*(Phi_U[i-1]+(x_n[i+2]*y_n[0]+x_n[i+1]*y_n[1])*Phi_U[1] \
            +x_n[i+1]*y_n[0]*Phi_U[2]-vid)

elif (Scenario_3):
    for i in range(3, u_max+1):
        vid=mpmath.mpf('0')
        for k in range(1, i):
            vid=vid+s_n[i+3-k]*Phi_U[k]
        Phi_U[i]=(1/s_n[3])*(Phi_U[i-1]+(x_n[i+2]*y_n[0]+x_n[i+1]*y_n[1])*Phi_U[1] \
            +x_n[i+1]*y_n[0]*Phi_U[2]-vid)

print(Phi_U)

PHI_U_rounded = mpmath.matrix(Phi_U.rows, Phi_U.cols)
# convert each element of Phi to a string with 3 decimal places
for i in range(Phi_U.rows):
    for j in range(Phi_U.cols):
        PHI_U_rounded[i, j] = mpmath.nstr(Phi_U[i, j], 3)

# Write the rounded matrix to a text file
with open("ultimate_time_s3.txt", "w") as f:
    f.write(str(PHI_U_rounded))

# Code execution timer
end_time = time.time()
print("Time taken: {:.6f} seconds".format(end_time - start_time))
