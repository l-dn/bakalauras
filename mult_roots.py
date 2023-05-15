# mult_roots.py
import sympy

'''
This file is used to make matrix equation, the solve it and find first
Phi_U, m_1, m_2 values. There are 3 functions for 3 different situations -
1, 2 or 3 different roots.
'''

def viena_saknis(mp_sol, kN, ES, x_n, y_n, y_lambda, kappa, y_xi):
    
    # Matrices:
    M_1_11=y_n[1]*mp_sol[0]
    M_1_21=y_n[1]

    M_2_11=x_n[1]*mp_sol[0]
    M_2_21=x_n[1]

    G_11=(sympy.N(sympy.exp(y_lambda*(mp_sol[0]-1)), 40)/mp_sol[0]**2)*(mp_sol[0]**y_xi)
    G_21=1
     
    M = sympy.Matrix ([[M_1_11, M_2_11*G_11],
                       [M_1_21, M_2_21*G_21]])
    rhs = sympy.Matrix([[0],
                        [kN-ES]])
    # Finding first m_1 ir m_2 values:
    m = sympy.linsolve((M, rhs))

    # Initialize m_1, m_2 and then add first values:
    m_1 = [sympy.N('0') for i in range(2)]
    m_2 = [sympy.N('0') for i in range(2)]

    m_list = list(m)
    m_1[0], m_2[0] = m_list[0]
    
    # Assign only real parts:
    for i in range (1):
        m_1[i] = sympy.re(m_1[i])
        m_2[i] = sympy.re(m_2[i])

    # Finding m_1, m_2[1, 2]:
    for i in range (2, 3):
        sum_1 = sympy.N('0')
        sum_2 = sympy.N('0')
        if i == kappa:
            sum_1 +=  m_2[0]*x_n[1]
            sum_2 +=  m_1[0]*y_n[1]

        for j in range (i-2+1):
            sum_1 +=  m_2[j]*x_n[i-j]
            sum_2 +=  m_1[j]*y_n[i-j]

        m_2[i-1] = (m_1[i-kappa] - sum_1) / x_n[1]
        m_1[i-1] = (m_2[i-kappa] - sum_2) / y_n[1]

    # Initialize Phi_U and then add first values:
    Phi_U = [sympy.N('0') for i in range(3)]
    
    for i in range (1, 3):
        sum_1=sympy.N('0')
        for j in range (i):
            sum_1 += m_1[j]

        Phi_U[i] = sum_1
 
    # Finding Phi_U(0) value:
    Phi_U[0] = x_n[1]*y_n[1]*Phi_U[2]+x_n[1]*y_n[2]*Phi_U[1]


    return Phi_U, m_1, m_2

def dvi_saknys(mp_sol, kN, ES, x_n, y_n, y_lambda, kappa, y_xi):
    
    # Matrices:
    M_1_11=y_n[1]*mp_sol[0] + y_n[0]*(mp_sol[0]+1)
    M_1_21=y_n[1]*mp_sol[1] + y_n[0]*(mp_sol[1]+1)
    M_1_31=y_n[1] + 2*y_n[0]

    M_1_12=y_n[0]*mp_sol[0]
    M_1_22=y_n[0]*mp_sol[1]
    M_1_32=y_n[0]
    
    M_2_11=x_n[1]*mp_sol[0] + x_n[0]*(mp_sol[0]+1)
    M_2_21=x_n[1]*mp_sol[1] + x_n[0]*(mp_sol[1]+1)
    M_2_31=x_n[1] + 2*x_n[0]

    M_2_12=x_n[0]*mp_sol[0]
    M_2_22=x_n[0]*mp_sol[1]
    M_2_32=x_n[0]

    G_11=(sympy.N(sympy.exp(y_lambda*(mp_sol[0]-1)), 40)/mp_sol[0]**2)*(mp_sol[0]**y_xi)
    G_21=(sympy.N(sympy.exp(y_lambda*(mp_sol[1]-1)), 40)/mp_sol[1]**2)*(mp_sol[1]**y_xi)
    G_31=1

    G_12=(sympy.N(sympy.exp(y_lambda*(mp_sol[0]-1)), 40)/mp_sol[0]**2)*(mp_sol[0]**y_xi)
    G_22=(sympy.N(sympy.exp(y_lambda*(mp_sol[1]-1)), 40)/mp_sol[1]**2)*(mp_sol[1]**y_xi)
    G_32=1
     
    M = sympy.Matrix ([[M_1_11, M_1_12, M_2_11*G_11, M_2_12*G_12],
                       [M_1_21, M_1_22, M_2_21*G_21, M_2_22*G_22],
                       [M_1_31, M_1_32, M_2_31*G_31, M_2_32*G_32]])
    rhs = sympy.Matrix([[0],
                        [0],
                        [kN-ES]])
    # Finding first m_1[0 1] ir m_2[0 1] values:
    m = sympy.linsolve((M, rhs))
    
    # Initialize m_1, m_2 and then add first values:
    m_1 = [sympy.N('0') for i in range(4)]
    m_2 = [sympy.N('0') for i in range(4)]
    m_list = list(m)
    m_1[0], m_1[1], m_2[0], m_2[1] = m_list[0]
    
    # Assign only real parts:
    for i in range (2):
        m_1[i] = sympy.re(m_1[i])
        m_2[i] = sympy.re(m_2[i])
    
    # first we find m_2[1] which we need:
    m_2[1] = (m_1[0] - (m_2[0]*x_n[1] + m_2[0]*x_n[2])) / x_n[1]
   
    # Finding m_1, m_2[2 3]:
    for i in range (2, 4):
        sum_1 = sympy.N('0')
        sum_2 = sympy.N('0')
        if i == kappa:
            for j in range (2):
                sum_3 =  sympy.N('0')
                sum_4 =  sympy.N('0')
                for k in range (2-j):
                    sum_3 += x_n[k]
                    sum_4 += y_n[k]

                sum_1 += m_2[j]*sum_3
                sum_2 += m_1[j]*sum_4

        for j in range (i-1+1):
            sum_1 +=  m_2[j]*x_n[i-j]
            sum_2 +=  m_1[j]*y_n[i-j]

        m_2[i] = (m_1[i-kappa] - sum_1) / x_n[0]
        m_1[i] = (m_2[i-kappa] - sum_2) / y_n[0]

    # Initialize Phi_U and then add first values:
    Phi_U = [sympy.N('0') for i in range(5)]
    
    for i in range (1, 5):
        sum_1=sympy.N('0')
        for j in range (i):
            sum_1 += m_1[j]

        Phi_U[i] = sum_1
 
    # Finding Phi_U(0) value:
    Phi_U[0] = x_n[0]*y_n[0]*Phi_U[4] + \
               (x_n[0]*y_n[1]+x_n[1]*y_n[0])*Phi_U[3] + \
               (x_n[0]*y_n[2]+x_n[1]*y_n[1])*Phi_U[2] + \
               (x_n[0]*y_n[3]+x_n[1]*y_n[2])*Phi_U[1]
    
    return Phi_U, m_1, m_2

def trys_saknys(mp_sol, kN, ES, x_n, y_n, y_lambda, kappa, y_xi):
    
    # Matrices:
    M_1_11=y_n[1]*mp_sol[0] + y_n[0]*(mp_sol[0]+1)
    M_1_21=y_n[1]*mp_sol[1] + y_n[0]*(mp_sol[1]+1)
    M_1_31=y_n[1]*mp_sol[2] + y_n[0]*(mp_sol[2]+1)
    M_1_41=y_n[1] + 2*y_n[0]

    M_1_12=y_n[0]*mp_sol[0]
    M_1_22=y_n[0]*mp_sol[1]
    M_1_32=y_n[0]*mp_sol[2]
    M_1_42=y_n[0]
    
    M_2_11=x_n[1]*mp_sol[0] + x_n[0]*(mp_sol[0]+1)
    M_2_21=x_n[1]*mp_sol[1] + x_n[0]*(mp_sol[1]+1)
    M_2_31=x_n[1]*mp_sol[2] + x_n[0]*(mp_sol[2]+1)
    M_2_41=x_n[1] + 2*x_n[0]

    M_2_12=x_n[0]*mp_sol[0]
    M_2_22=x_n[0]*mp_sol[1]
    M_2_32=x_n[0]*mp_sol[2]
    M_2_42=x_n[0]

    G_11=(sympy.N(sympy.exp(y_lambda*(mp_sol[0]-1)), 40)/mp_sol[0]**2)*(mp_sol[0]**y_xi)
    G_21=(sympy.N(sympy.exp(y_lambda*(mp_sol[1]-1)), 40)/mp_sol[1]**2)*(mp_sol[1]**y_xi)
    G_31=(sympy.N(sympy.exp(y_lambda*(mp_sol[2]-1)), 40)/mp_sol[2]**2)*(mp_sol[1]**y_xi)
    G_41=1

    G_12=(sympy.N(sympy.exp(y_lambda*(mp_sol[0]-1)), 40)/mp_sol[0]**2)*(mp_sol[0]**y_xi)
    G_22=(sympy.N(sympy.exp(y_lambda*(mp_sol[1]-1)), 40)/mp_sol[1]**2)*(mp_sol[1]**y_xi)
    G_32=(sympy.N(sympy.exp(y_lambda*(mp_sol[2]-1)), 40)/mp_sol[2]**2)*(mp_sol[1]**y_xi)
    G_42=1
     
    M = sympy.Matrix ([[M_1_11, M_1_12, M_2_11*G_11, M_2_12*G_12],
                       [M_1_21, M_1_22, M_2_21*G_21, M_2_22*G_22],
                       [M_1_31, M_1_32, M_2_31*G_31, M_2_32*G_32],
                       [M_1_41, M_1_42, M_2_41*G_41, M_2_42*G_42]])
    rhs = sympy.Matrix([[0],
                        [0],
                        [0],
                        [kN-ES]])
    # Finding first m_1[0 1] ir m_2[0 1] values:
    m = sympy.linsolve((M, rhs))

    # Initialize m_1, m_2 and then add first values:
    m_1 = [sympy.N('0') for i in range(4)]
    m_2 = [sympy.N('0') for i in range(4)]

    m_list = list(m)
    m_1[0], m_1[1], m_2[0], m_2[1] = m_list[0]
    
    # Assign only real parts:
    for i in range (2):
        m_1[i] = sympy.re(m_1[i])
        m_2[i] = sympy.re(m_2[i])
            
    # Finding m_1, m_2[2 3]:
    for i in range (2, 4):
        sum_1 = sympy.N('0')
        sum_2 = sympy.N('0')
        if i == kappa:
            for j in range (2):
                sum_3 =  sympy.N('0')
                sum_4 =  sympy.N('0')
                for k in range (2-j):
                    sum_3 += x_n[k]
                    sum_4 += y_n[k]

                sum_1 += m_2[j]*sum_3
                sum_2 += m_1[j]*sum_4

        for j in range (i-1+1):
            sum_1 +=  m_2[j]*x_n[i-j]
            sum_2 +=  m_1[j]*y_n[i-j]

        m_2[i] = (m_1[i-kappa] - sum_1) / x_n[0]
        m_1[i] = (m_2[i-kappa] - sum_2) / y_n[0]

    # Initialize Phi_U and then add first values:
    Phi_U = [sympy.N('0') for i in range(5)]
    
    for i in range (1, 5):
        sum_1=sympy.N('0')
        for j in range (i):
            sum_1 += m_1[j]

        Phi_U[i] = sum_1
 
    # Finding Phi_U(0) value:
    Phi_U[0] = x_n[0]*y_n[0]*Phi_U[4] + \
               (x_n[0]*y_n[1]+x_n[1]*y_n[0])*Phi_U[3] + \
               (x_n[0]*y_n[2]+x_n[1]*y_n[1])*Phi_U[2] + \
               (x_n[0]*y_n[3]+x_n[1]*y_n[2])*Phi_U[1]
    
    return Phi_U, m_1, m_2