from __future__ import division
from numpy import array
from math import exp, sqrt, log, pi, sinh, cosh
from PNJunction import *

# TODO Complete the documentation
k = 1.38064852*(10**(-23))      # Boltzmann constant; m^2.kg/s^2.K
e = 1.602*(10**(-19))           # charge of electron; C
h = 6.62607*(10**-34)           # Planck constant; m^2.kg/s
u_p = 0.045                     # mobility of holes; m^2/V.s
u_n = 0.14                      # mobility of electrons; m^2/V.s
t_p = 10**(-6)                  # relaxation time of holes; s
t_n = 10**(-6)                  # relaxation time of electrons; s
c = 3*(10**8)                   # speed of light; m/s

def HoleDiffusionCoefficient(T):
    """
    """
    return (k*T/e)*u_p


def ElectronDiffusionCoefficient(T):
    """
    """
    return (k*T/e)*u_n


def LengthHoleDiffusion(T):
    """
    """
    return sqrt(HoleDiffusionCoefficient(T)*t_p)


def LengthElectronDiffusion(T):
    """
    """
    return sqrt(ElectronDiffusionCoefficient(T)*t_n)


def PDarkCurrentDensity(T, V, Pno):
    """
    """
    return (e*HoleDiffusionCoefficient(T,u_p)*Pno)/LengthHoleDiffusion(T, u_p)


def NDarkCurrentDensity(T, V, Npo):
    """
    """
    return (e*ElectronDiffusionCoefficient(T,u_n)*Npo)/LengthElectronDiffusion(T, u_n)


def HoleSurfaceRecombVelocity(wavelength):
    """
    """
    return 50


def ElectronSurfaceRecombVelocity(wavelength):
    """
    """
    return 100


def Reflectivity(wavelength):
    """
    """
    return 0.3280328678


def PhotonsGenerated(wavelength, width, intensity):
    """
    """
    return (intensity*Area(width)*wavelength)/(h*c)


def Area(width):
    """
    """
    return width*width*5


def AbsorptionCoefficient(wavelength):
    """
    """
    return 1.18*pow(10,-3)


def CarrierDensityHole(T, wavelength, x, width, intensity):
    """
    """
    numer_1 = e*PhotonsGenerated(wavelength, width, intensity)*LengthHoleDiffusion(T,u_p)*AbsorptionCoefficient(wavelength)*(1 - Reflectivity(wavelength))
    denom_1 = (AbsorptionCoefficient(wavelength)**2)*(LengthHoleDiffusion(T, u_p)**2) - 1

    numer_2 = (HoleSurfaceRecombVelocity(wavelength)*LengthHoleDiffusion(T, u_p))/HoleDiffusionCoefficient(T,u_p)
    numer_3 = AbsorptionCoefficient(wavelength)*LengthHoleDiffusion(T,u_p)
    numer_4 = exp(-1*AbsorptionCoefficient(wavelength)*x)
    numer_5 = numer_2*cosh(x/LengthHoleDiffusion(T,u_p)) + sinh(x/LengthHoleDiffusion(T,u_p))
    denom_2 = numer_2*sinh(x/LengthHoleDiffusion(T,u_p)) + cosh(x/LengthHoleDiffusion(T,u_p))
    numer_6 = AbsorptionCoefficient(wavelength)*LengthHoleDiffusion(T, u_p)*exp((-1)*AbsorptionCoefficient(wavelength)*x)

    return (numer_1/denom_1)*((numer_2 + numer_3 - numer_4*numer_5)/denom_2) - numer_6


def CarrierDensityElectron(T, wavelength, x, W, H_, width, intensity):
    """
    """
    numer_1 = e*PhotonsGenerated(wavelength, width, intensity)*LengthElectronDiffusion(T,u_n)*AbsorptionCoefficient(wavelength)*(1 - Reflectivity(wavelength))
    denom_1 = (AbsorptionCoefficient(wavelength)**2)*(LengthElectronDiffusion(T, u_n)**2) - 1
    numer_2 = exp((-1)*AbsorptionCoefficient(wavelength)*(x + W))
    numer_3 = AbsorptionCoefficient(wavelength)*LengthElectronDiffusion(T, u_n)

    numer_4 = (ElectronSurfaceRecombVelocity(wavelength)*LengthElectronDiffusion(T, u_n))/ElectronDiffusionCoefficient(T, u_n)
    numer_5 = (cosh(H_/LengthElectronDiffusion(T, u_n)) - exp(AbsorptionCoefficient(wavelength)*H_))
    numer_6 = sinh(H_/LengthElectronDiffusion(T, u_n))
    numer_7 = AbsorptionCoefficient(wavelength)*LengthElectronDiffusion(T,u_n)*exp(AbsorptionCoefficient(wavelength)*H_)
    denom_2 = numer_4*sinh(H_/LengthElectronDiffusion(T,u_n)) + cosh(H_/LengthElectronDiffusion(T,u_n))

    return (numer_1*numer_2/denom_1)*numer_3 - (numer_4*numer_5 + numer_6 + numer_7)/denom_2


def LightGeneratedCarriers(T, V, Pno, Npo, wavelength, width, x, W, H_, intensity):
    """
    """
    val_1 = (PDarkCurrentDensity(T, V, Pno) + NDarkCurrentDensity(T, V, Npo))*(exp((e*V)/(k*T)) - 1)
    val_2 = (-1)*CarrierDensityHole(T, wavelength, x, width, intensity) + (-1)*CarrierDensityElectron(T, wavelength, x, W, H_, width, intensity)

    return val_1 + val_2
