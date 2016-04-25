from math import exp, sqrt, log, pi, sinh, cosh
from __future__ import division

L =
n0 =
p0 =
Gamma0 =
tn0 = 10**(-6)
tp0 = 10**(-6)
epsin =
epsip =
k = 1.38064852*(10**(-23))
q = -1.602*(10**(-19))
h = 6.62607*(10**-34)

def JoP(L):
    return (-2)*q*n0*L*(Dn/(L(tp,Dp)**2))*(Beta5()/(Beta1()**2))*(I1(Beta5())/I0(Beta5()))

def JoN(L):
    return (-2)*p0*q*L*(Dp/(L(tp,Dp)**2))*(Beta2()/(Beta1()**2))*((f1*K1(Beta2()) - f2*I1(Beta2()))/(f1*K0(Beta2()) + f2*I0(Beta2())))

def J1N:
    factor1 = K1(Beta2())*(f1() - Beta4()*I0(Beta2())) - I1(Beta2())*(f2() + Beta4()*K0(Beta2()))
    factor2 = f1()*K0(Beta2()) + f2()*I0(Beta2())
    return (-2)*q*Gamma0*(Beta2()/(Beta1()**2))*(factor1/factor2)*(1 - exp(-1*Beta3))

def J1P:
    return (-2)*q*Gamma0*((L(tn,Dn)**2)/(l(tp,Dp)**2))*(Beta5()/(Beta1()**2))*(I1(Beta5())/I0(Beta5()))*(1-exp(-1*Beta6()))

def JgDepP(d2,x4):
    return (-1)*q*Gamma0*((d2**2 - x4**2)/(R**2))*(1 - exp(-1*Beta6()))

def JgDepN(d2):
    return (-1)*q0*Gamma0*(((d2+x2(V))**2-(d2**2))/R**2)*(1-exp(-1*Beta3()))

def JrDep:
    return (-1)*q*L*Umax*((r2**2-r1**2)/R**2)

def Beta1:
    return R/Lp

def Beta2:
    return (R-x1)/Lp

def Beta3:
    return alpha1*L

def Beta4:
    return Lp*Sp/Dp

def Beta5:
    return x4/Ln

def Beta6:
    return alpha2*L

def I1(Beta):
    # Function definition that pertains to I1
    pass

def I0(Beta):
    # Function definition that pertains to I0
    pass

def f1:
    # function of Beta_1 and Beta_4
    # return the value
    return I1(Beta1()) + Beta_4*I0(Beta1())

def f2:
    # function of Beta_1 and Beta_4
    # return the value
    return K1(Beta1()) - Beta_4*K0(Beta4())

def L(t,D):
    return sqrt(t*D)

def Umax(V,T):
    return (ni/(sqrt(tn0*tp0)))*sinh((q*V)/(2*k*T))

def r(V):
    return x4 + (log(Na/nip)/log((Na*Nd)/(nip*nin)))*(x2(V) + x3(V))

def kappa(T,V):
    return (pi*k*T)/(q*(VbiminusV()))

def VbiminusV(Na,Nd,x4):
    term_1 = ((q*Nd)/(2*epsin))*((d2+x2(V))**2)*log((d2+x2(V))/d2)
    term_2 = ((q*Nd)/(2*epsin))*((d2**2)-(d2+x2(V))**2) + ((q*Na)/(2*epsip))*x4*log(x4/d2)
    term_3 = ((q*Na)*(d2**2-x4**2))/(4*epsip)
    return term_3 + term_2 + term_1

def r1(V):
    return r(V) - ((x2(V) + x3(V))/2)*kappa(T,V)

def r2(V):
    return r(V) + ((x2(V) + x3(V))/2)*kappa(T,V)

def x3(V):
    pass

def x2(V):
    pass
