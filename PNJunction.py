from numpy import array
from __future__ import division, multiplication
from math import exp, sqrt, log, pi

class PNJunction:
    """This class encapsulates the methods required to measure parameters of a nanorod based solar cell"""
    k = 1.38064852*(10**(-23))
    e = -1.602*(10**(-19))
    me = 0.5*9.11*(10**(-31))
    mh = 0.58*9.11*(10**(-31))
    h = 6.62607*(10**-34)

    @staticmethod
    def phiFp(T, Na, ni):
    # TODO Complete the documentation of this method
        return (k*T/e)*log(Na/ni)

    @staticmethod
    def phiFn(T, Nd, ni):
    # TODO Complete the documentation of this method
        return (k*T/e)*log(Nd/ni)

    @staticmethod
    def BuiltinVoltage(T, Na, Nd, ni):
    """ Calculates the internal built in voltage of the cell
    @param: T The temperature of operation
    @param: Na Density of acceptor type dopant added to the material
    @param: Nd Density of the donor type dopant added to the material
    @param: ni Intrinsic carrier concentration
    """
    # Where is ni obtained from? Generated from the function or input?
        return phiFp(T,Na,ni) + phiFn(T,Nd,ni)

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
    def FermiEnergy(T, Ec, Ev):
    """ Returns the intrinsic Fermi energy
    @param: T Temperature of operation
            Ec Energy level of the conduction band
            Ev Energy level of the valence band
    """
        return (Ec + Ev)/2 - (k*T)log(DensityConductionBand(T)/DensityValenceBand(T))

    @staticmethod
    def ICarrierConc(T, Eg):
    # TODO Complete the documentation of this method
        return sqrt(DensityConductionBand(T)*DensityValenceBand(T))*exp(-1*Eg/(2*k*T))
