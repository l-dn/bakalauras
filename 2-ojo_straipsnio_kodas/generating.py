#generating.py
import sympy
import factor

'''
This file is used to make survival probability-generating function
and count Phi_U values using it.
'''

def gener_f(x_lambda, y_lambda, x_xi, y_xi, x_n, y_n, m_1, m_2, u_max):
    
    # find coeficients:
    s_1 = m_2[0]*x_n[0]
    s_2 = m_2[0]*(x_n[0]+x_n[1])+m_2[1]*x_n[0]
    
    e_1 = m_1[0]*y_n[0]
    e_2 = m_1[0]*(y_n[0]+y_n[1])+m_1[1]*y_n[0]
       
    # define the symbolic variables
    s = sympy.symbols('s')
    
    # Unnessecery factorisation
    func = factor.smpl(x_lambda, y_lambda, x_xi, y_xi, s_1, s_2, e_1, e_2)
    
    func_N = sympy.N(func, 4)
    
    # For faster computation we will expand the gener. func. to tailor series:
    taylor = func.series(s, x0=0, n=u_max+10)
        
    # Initialize Phi_U and then add values:
    Phi_U_G = [sympy.N('0') for i in range(u_max + 1)]
    
    for i in range (u_max):
        if i == 0:
            dfdx = taylor
        elif i == 1:
            dfdx = sympy.Derivative(taylor, s)
        else:
            dfdx = sympy.Derivative(dfdx, s)
        Phi_U_G[i+1] = dfdx.subs(s, 0).evalf()*(1/sympy.factorial(i))
    
    return func_N, Phi_U_G
