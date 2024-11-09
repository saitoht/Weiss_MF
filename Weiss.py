import numpy as np
import matplotlib.pyplot as plt
import sys, os
from scipy import integrate

R = 8.3144569919459 # J/mol K
kB = 8.6173303e-5 # eV/K
T2J = 5.584937543138665 # J/mol T
mu_B = 5.7883817982e-05 # eV/T
eV2Jmol = 96.4853e3

def main():
    """ m0: magnetic moment at T=0k, B=0T from DFT (mu_B)
        temp: temperature (K)
        Tc: Curie temperature (K)
        ndiv: number of devision for temperature
        colors: colors to plot physical quantities, should be the same dimension as magfield
        nmag: number of magnetic atom in the unit cell
        lamb: lambda, molecular field const 
        magfield: magnetic field B (T) 
        struc: crystal structure, now BCC or FCC"""
    ndiv = 1000
    ndivf = 100
    magfield = 0.5
    Tc = 1043.
    nmag = 1
    struc = 'carbide'
    
    usage = """ USAGE: python Weiss.py 2.2 [-p] [-h] """
    args = sys.argv
    if ( len(args) < 2 or '-h' in args ):
        print(usage)
        sys.exit()
    m0 = float(args[1])
    print("*moment @ T=0K: {0} (mu_B)".format(m0))
    lamb = Tc2lambda(0.5*m0, Tc, nmag)
    print("*lambda: {0} (T)".format(lamb*mu_B**2))
    tempc0A = np.linspace(1e-3,1.,2*ndiv,endpoint=False)
    tempc0B = np.linspace(1.,0.25,ndiv)
    tempc1A = np.sqrt(tempc0A)
    tempc1B = 2. - np.sqrt(tempc0B)
    tempc = np.append(tempc1A.tolist(), tempc1B.tolist())
    temp = Tc * tempc    
    iTc = np.where(temp>Tc)[0][0]
    magf = np.linspace(0.,magfield,ndivf)
    moment = np.zeros((len(magf),3*ndiv))
    moment[:,0] = m0
    for j, mf in enumerate(magf):
        for i in range(len(temp)-1):
            moment[j,i+1] = mmom(temp[i], mf, m0, lamb, moment[j,i])

    dGm = np.zeros(len(temp))
    for i in range(len(temp)):
        dGm[i] = -mu_B*eV2Jmol*integrate.simps(moment[:,i],magf[:])
    
    if ('-p' in args):
        plt.xlabel("Temperature (K)")
        plt.ylabel("Magnetic moment ($\mu_B$)")
        plt.xlim(0.,1.5*Tc)
        plt.ylim(0.,1.1*m0)
        plt.plot(temp, moment[len(magf)-1,:], c="black", label="B={0}T".format(magfield))
        plt.legend()
        plt.savefig("Moment.pdf")
        plt.show()

        plt.xlabel("Temperature (K)")
        plt.ylabel("$\Delta G_m$ (J/mol)")
        plt.plot(temp, dGm, c="black", label="B={0}T".format(magfield))
        plt.legend()
        plt.savefig("dGm.pdf")
        plt.show()

def coth(x):
    if (np.abs(x) < 1e-10):
        return 0.
    else:
        return 1./np.tanh(x)

def Brillouin(j, alpha):
    """ define Brillouin function   Eq. (S4) """
    cf = 0.5*((2*j+1)/j)
    cfp = 0.5/j
    return cf*coth(cf*alpha) - cfp*coth(cfp*alpha)

def mmom(temp, magfield, mmom0, lamb, mmomb):
    """ return magnetic moment at T, B    Eq. (S2) """
    alpha = mmom0*mu_B * (magfield + lamb*mmomb*mu_B)/(kB*temp)  ### Eq. (S3)
    return mmom0 * Brillouin(0.5*mmom0, alpha)

def Tc2lambda(J, Tc, n, g=2.):
    return 3.*kB*Tc/(J*(J+1.)*n*(g**2)*(mu_B**2))  # T/eV

if __name__=='__main__':
    main()
