#factor.py
import sympy

'''
Simplify doesnt do factorisation so this function does it is a pain but
i must do it this way.
'''

'''
original functions, but factorisation of s**1 doesnt work idunno
func = ((s_1 + s_2*s)*s**2 + (e_1 + e_2*s)*s**(x_xi)*sympy.exp(x_lambda*(s-1))) / \
        (s**(x_xi+y_xi)*sympy.exp((x_lambda+y_lambda)*(s-1))-s**4)
'''

def smpl(x_lambda, y_lambda, x_xi, y_xi, s_1, s_2, e_1, e_2):
        
    # define the symbolic variables
    s = sympy.symbols('s')     
        
        
    if x_xi+y_xi == 0:
        func = ((s_1 + s_2*s)*s**2 + (e_1 + e_2*s)*sympy.exp(x_lambda*(s-1))) / \
                (sympy.exp((x_lambda+y_lambda)*(s-1))-s**4)
        
    elif x_xi+y_xi == 1:
        if x_xi==1:
            func = ((s_1 + s_2*s)*s**1 + (e_1 + e_2*s)*sympy.exp(x_lambda*(s-1))) / \
                    (sympy.exp((x_lambda+y_lambda)*(s-1))-s**3)
        elif y_xi==1:
            func = ((s_1 + s_2*s)*s**1 + (e_1 + e_2)*sympy.exp(x_lambda*(s-1))) / \
                    (sympy.exp((x_lambda+y_lambda)*(s-1))-s**3)
            
    elif x_xi+y_xi == 2:
        if x_xi==2: 
            func = ((e_1 + e_2*s)*sympy.exp(x_lambda*(s-1))) / \
                    (sympy.exp((x_lambda+y_lambda)*(s-1))-s**2)
        elif y_xi==2:
            func = ((s_1 + s_2*s)) / \
                    (sympy.exp((x_lambda+y_lambda)*(s-1))-s**2)
        elif x_xi==1 and y_xi==1:
            func = (s_2*s + e_2*sympy.exp(x_lambda*(s-1))) / \
                    (sympy.exp((x_lambda+y_lambda)*(s-1))-s**2)
           
    elif x_xi+y_xi == 3:
        if x_xi==2 and y_xi==1: 
            func = (e_2*sympy.exp(x_lambda*(s-1))) / \
                (sympy.exp((x_lambda+y_lambda)*(s-1))-s)
        elif y_xi==2 and x_xi==1:
            func = (s_2) / \
                (sympy.exp((x_lambda+y_lambda)*(s-1))-s)
        elif x_xi==3:
            func = ((e_1 + e_2*s)*sympy.exp(x_lambda*(s-1))) / \
                (sympy.exp((x_lambda+y_lambda)*(s-1))-s)
        elif y_xi==3:
            func = (s_1 + s_2*s) / \
                (s*sympy.exp((x_lambda+y_lambda)*(s-1))-s**2)
            
    return func