import time
import sympy
import mult_roots
import generating

#import sys

start_time = time.time()

'''
Bi-seasonal discrete time risk model with income rate two.
'''

# Poisson pmf:
def poisson_pmf(k, lambd, xi):
    if k<xi: return sympy.N('0')
    return sympy.Pow(lambd, k-xi) * sympy.N(sympy.exp(-lambd), 40) / \
           sympy.factorial(k-xi)

# Define the parameters
x_lambda = sympy.N('1', 40)
x_xi = sympy.N('1')
y_lambda = sympy.N('19/10', 40)
y_xi = sympy.N('0')

kappa=2
kN=kappa*2

ES=x_lambda+y_lambda+x_xi+y_xi

# What is the maximum u?
u_max=40

# Initialize the arrays
x_n = [sympy.N('0') for i in range(0, u_max+1)]
y_n = [sympy.N('0') for i in range(0, u_max+1)]

# Perform the Poisson pmf calculations using function:
for i in range(u_max+1):
    x_n[i] = poisson_pmf(i, x_lambda, x_xi)
    y_n[i] = poisson_pmf(i, y_lambda, y_xi)

# define the symbolic variables
s = sympy.symbols('s')

# define and solve the equation
eq = sympy.Eq(s**(x_xi+y_xi) * sympy.N(sympy.exp((x_lambda+y_lambda)*(s-1)), 40) - s**4, 0)
sol = sympy.solve(eq, s)

# convert the solutions to complex floating-point numbers
mp_sol = [sympy.N(val, 40, chop=True) for val in sol]

# Only instances in complex circle with radius 1 and !=1,0 are needed
mp_sol = [val for val in mp_sol if sympy.Abs(val) < 1 and sympy.Abs(val) > 0]

# Find first Phi_U values
if (len(mp_sol) == 3):
    print("3 saknys.")
    [Phi_U, m_1, m_2] = mult_roots.trys_saknys (mp_sol, kN, ES, x_n, y_n, y_lambda, kappa, y_xi)

if (len(mp_sol) == 2):
    print("2 saknys.")
    [Phi_U, m_1, m_2] = mult_roots.dvi_saknys (mp_sol, kN, ES, x_n, y_n, y_lambda, kappa, y_xi)

if (len(mp_sol) == 1):
    print("1 saknis.")
    [Phi_U, m_1, m_2] = mult_roots.viena_saknis (mp_sol, kN, ES, x_n, y_n, y_lambda, kappa, y_xi)
#print([sympy.N(val, 4) for val in Phi_U])
    

# Find Phi_U values using generating function:
[func, Phi_U_G] = generating.gener_f(x_lambda, y_lambda, x_xi, y_xi, x_n, y_n, m_1, m_2, u_max)
Phi_U_G[0]=Phi_U[0]

#print(func)
print([sympy.N(val, 4) for val in Phi_U_G])

# Write the rounded Phi values to a text file
with open("ultimate_time_phi.txt", "w") as f:
    for val in Phi_U_G:
        f.write(str('u = ' + str(Phi_U_G.index(val)) + ', Phi(u) = ' + str(sympy.N(val, 4)) + '.\n'))
            
# Code execution timer
end_time = time.time()
print("Time taken: {:.6f} seconds".format(end_time - start_time))