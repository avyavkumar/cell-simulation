from __future__ import division
from numpy import array
from math import exp, sqrt, log, pi, sinh, cosh
from PNJunction import PNJunction

class GeneratedCarriers:
    # TODO Complete the documentation
    k = 1.38064852*(10**(-23))      # Boltzmann constant; m^2.kg/s^2.K
    e = 1.602*(10**(-19))           # charge of electron; C
    h = 6.62607*(10**-34)           # Planck constant; m^2.kg/s
    u_p =                           # mobility of holes
    u_n =                           # mobility of electrons
    t_p = 10**(-6)

    @staticmethod
    def HoleDiffusionCoefficient(T, u_p):
        """
        """
        return (k*T/e)*u_p

    @staticmethod
    def ElectronDiffusionCoefficient(T, u_n):
        """
        """
        return (k*T/e)*u_n

    @staticmethod
    def LengthHoleDiffusion(T, u_p):
        """
        """
        return sqrt(HoleDiffusionCoefficient(T, u_p)*t_p)

    @staticmethod
    def LengthElectronDiffusion(T, u_n):
        """
        """
        return sqrt(ElectronDiffusionCoefficient(T, u_n)*t_n)

    @staticmethod
    def PDarkCurrentDensity(T, V, Pno):
        """
        """
        return (e*HoleDiffusion(T,u_p)*Pno*[exp((e*V)/(k*T)) - 1])/LengthHoleDiffusion(T, u_p)

    @staticmethod
    def NDarkCurrentDensity(T, V, Npo):
        """
        """
        return ((-1)*e*ElectronDiffusion(T,u_n)*Npo*[exp((e*V)/(k*T)) - 1])/LengthElectronDiffusion(T, u_n)

    @staticmethod
    def HoleSurfaceRecombVelocity(lambda):
        """
        """


    @staticmethod
    def ElectronSurfaceRecombVelocity(lambda):
        """
        """


    @staticmethod
    def Reflectivity(lambda):
        """
        """


    @staticmethod
    def PhotonsGenerated(lambda):
        """
        """


    @staticmethod
    def AbsorptionCoeffient(lambda):
        """
        """


    @staticmethod
    def CarrierDensityHole(T, lambda, u_p, x):
        """
        """
        numer_1 = e*PhotonsGenerated(lambda)*LengthHoleDiffusion(T,u_p)*AbsorptionCoeffient(lambda)*(1 - Reflectivity(lambda))
        denom_1 = (AbsorptionCoeffient(lambda)**2)*(LengthHoleDiffusion(T, u_p)**2) - 1

        numer_2 = (HoleSurfaceRecombVelocity(lambda)*LengthHoleDiffusion(T, u_p))/HoleDiffusionCoefficient(T,u_p)
        numer_3 = AbsorptionCoeffient(lambda)*LengthHoleDiffusion(T,u_p)
        numer_4 = exp(-1*AbsorptionCoeffient(lambda)*x)
        numer_5 = numer_2*cosh(x/LengthHoleDiffusion(T,u_p)) + sinh(x/LengthHoleDiffusion(T,u_p))
        denom_2 = numer_2*sinh(x/LengthHoleDiffusion(T,u_p)) + cosh(x/LengthHoleDiffusion(T,u_p))
        numer_6 = AbsorptionCoeffient(lambda)*LengthHoleDiffusion(T, u_p)*exp((-1)*AbsorptionCoeffient(lambda)*x)

        return (numer_1/denom_1)*((numer_2 + numer_3 - numer_4*numer_5)/denom_2) - numer_6

    @staticmethod
    def CarrierDensityElectron(T, lambda, u_n, x, W, H_):
        """
        """
        numer_1 = e*PhotonsGenerated(lambda)*LengthElectronDiffusion(T,u_n)*AbsorptionCoeffient(lambda)*(1 - Reflectivity(lambda))
        denom_1 = (AbsorptionCoeffient(lambda)**2)*(LengthElectronDiffusion(T, u_n)**2) - 1
        numer_2 = exp((-1)*AbsorptionCoeffient(lambda)*(x + W))
        numer_3 = AbsorptionCoeffient(lambda)*LengthElectronDiffusion(T, u_n)

        numer_4 = (ElectronSurfaceRecombVelocity(lambda)*LengthElectronDiffusion(T, u_n))/ElectronDiffusionCoefficient(T, u_n)
        numer_5 = (cosh(H_/LengthElectronDiffusion(T, u_n)) - exp(AbsorptionCoefficent(lambda)*H_))
        numer_6 = sinh(H_/LengthElectronDiffusion(T, u_n))
        numer_7 = AbsorptionCoefficent(lambda)*LengthElectronDiffusion(T,u_n)*exp(AbsorptionCoefficent(lambda)*H_)
        denom_2 = numer_4*sinh(H_/LengthElectronDiffusion(T,u_n)) + cosh(H_/LengthElectronDiffusion(T,u_n))

        return (numer_1*numer_2/denom_1)*numer_3 - (numer_4*numer_5 + numer_6 + numer_7)/denom_2

    @staticmethod
    def LightGeneratedCarriers(T, V, Pno, Npo, lambda, u_p, u_n, x, W, H_):
        """
        """
        val_1 = (PDarkCurrentDensity(T, V, Pno) + NDarkCurrentDensity(T, V, Npo))*(exp(e*V/k*T) - 1)
        val_2 = (-1)*CarrierDensityHole(T, lambda, u_p, x) + (-1)*CarrierDensityElectron(T, lambda, u_n, x, W, H_)
        return val_1 + val_2
