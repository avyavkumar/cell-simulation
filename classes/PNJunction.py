from __future__ import division
from numpy import array
from math import exp, sqrt, log, pi

class PNJunction:
    """This class encapsulates the methods required to measure parameters of a nanorod based solar cell"""
    k = 1.38064852*(10**(-23))      # Boltzmann constant; m^2.kg/s^2.K
    e = -1.602*(10**(-19))          # charge of electron; C
    me = 0.5*9.11*(10**(-31))
    mh = 0.58*9.11*(10**(-31))
    h = 6.62607*(10**-34)           # Planck constant; m^2.kg/s

    @staticmethod
    def WorkFunctionPType(T, Na, Eg):
        """Returns the work function of the P type diode
        @param: T The temperature of operation
        @param: Na Density of acceptor type dopant added to the material
        @param: Eg Activation energy
        """
        return (k*T/e)*log(Na/ICarrierConc(T,Eg))

    @staticmethod
    def workFunctionNType(T, Nd, Eg):
        """Returns the work function of the P type diode
        @param: T The temperature of operation
        @param: Nd Density of the donor type dopant added to the material
        @param: Eg Activation energy
        """
        return (k*T/e)*log(Nd/ICarrierConc(T, Eg))

    @staticmethod
    def BuiltinVoltage(T, Na, Nd, Eg):
        """ Calculates the internal built in voltage of the cell
        @param: T The temperature of operation
        @param: Na Density of acceptor type dopant added to the material
        @param: Nd Density of the donor type dopant added to the material
        @param: Eg Activation energy
        """
        return phiFp(T,Na,ICarrierConc(T, Eg)) + phiFn(T,Nd,ICarrierConc(T, Eg))

    @staticmethod
    def Nno(T, Na, Nd, ni, Npo):
        # TODO Complete the documentation of this method
        return Npo*exp(e*BuiltinVoltage(T,Na,Nd,ni)/(k*T))

    @staticmethod
    def Ppo(T, Na, Nd, ni, Pno):
        # TODO Complete the documentation of this method
        return Pno*exp(e*BuiltinVoltage(T,Na,Nd,ni)/(k*T))

    @staticmethod
    def DensityConductionBand(T):
        """ Returns the effective density of states in the conduction band
        @param: T Temperature of operation
        """
        return 2*(((2*pi*me*k*T)/(h**2))**(3/2))

    @staticmethod
    def DensityValenceBand(T):
        """ Return the effective density of states in the valence band
        @param: T Temperature of operation
        """
        return 2*(((2*pi*mh*k*T)/(h**2))**(3/2))

    @staticmethod
    def IFermiEnergy(T, Ec, Ev):
        """ Returns the intrinsic Fermi energy
        @param: T Temperature of operation
        @param: Ec Energy level of the conduction band
        @param: Ev Energy level of the valence band
        """
        return (Ec + Ev)/2 - (k*T)*log(DensityConductionBand(T)/DensityValenceBand(T))

    @staticmethod
    def ICarrierConc(T, Eg):
        """ Returns the intrinsic particle concentration
        @param: T Temperature of operation
        @param: Eg Activation energy
        """
        return sqrt(DensityConductionBand(T)*DensityValenceBand(T))*exp(-1*Eg/(2*k*T))

    @staticmethod
    def DNTypeConc(T, Eg, Ef, Ec, Ev):
        """ Returns the concentration of the doped N type particles
        @param: T Temperature of operation
        @param: Eg Activation energy
        @param: Ef Input value of the Fermi level
        @param: Ec Energy level of the conduction band
        @param: Ev Energy level of the valence band
        """
        return ICarrierConc(T, Eg)*exp((Ef - IFermiEnergy(T, Ec, Ev))/(k*T))

    @staticmethod
    def DPTypeConc(T, Eg, Ef, Ec, Ev):
        """ Returns the concentration of the doped P type particles
        @param: T Temperature of operation
        @param: Eg Activation energy
        @param: Ef Input value of the Fermi level
        @param: Ec Energy level of the conduction band
        @param: Ev Energy level of the valence band
        """
        return ICarrierConc(T, Eg)*exp((IFermiEnergy(T, Ec, Ev)- Ef)/(k*T))
