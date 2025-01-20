import numpy as np
import matplotlib.pyplot as plt
import sys, os
from scipy import integrate

R = 8.3144569919459 # J/mol K
kB = 8.6173303e-5 # eV/K
T2J = 5.584937543138665 # J/mol T
mu_B = 5.7883817982e-05 # eV/T
eV2Jmol = 96.4853e3

plt.rcParams['font.family'] = 'Helvetica'
plt.rcParams["xtick.labelsize"]=14.0
plt.rcParams["ytick.labelsize"]=14.0
plt.rcParams["xtick.major.pad"] = 5
plt.rcParams["ytick.major.pad"] = 5
plt.rcParams["axes.labelsize"] = 14.0
plt.rcParams["axes.linewidth"] = 1.0
plt.rcParams["axes.labelpad"] = 6
plt.rcParams["xtick.direction"] = "out" 
plt.rcParams["ytick.direction"] = "out"
plt.rcParams["xtick.major.width"] = 1.0
plt.rcParams["ytick.major.width"] = 1.0
plt.rcParams["xtick.minor.width"] = 0.5
plt.rcParams["ytick.minor.width"] = 0.5
plt.rcParams["xtick.major.size"] = 4.5
plt.rcParams["ytick.major.size"] = 4.5
plt.rcParams["xtick.minor.size"] = 3.0
plt.rcParams["ytick.minor.size"] = 3.0
plt.rcParams["legend.edgecolor"] = 'black'
plt.rcParams["legend.fancybox"] = False

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
    m0 = [2.2, 2.2, 1.7, 1.7]
    magfield = [0.5,10.0,0.5,10.0]
    colors = ["magenta","red","cyan","blue"]
    Tc = [1043.,1043.0,800.,800.0]
    nmag = [1,1,1,1]
    
    usage = """ USAGE: python Weiss.py [-p] [-h] """
    args = sys.argv
    if ( '-h' in args ):
        print(usage)
        sys.exit()
    lengs = [len(m0),len(magfield),len(colors),len(Tc),len(nmag)]
    if ( not lengs[1:] == lengs[:-1] ):
        print("*ERROR: You should provide the same number of parameters in Weiss.py!")
        sys.exit()
    lamb = [Tc2lambda(0.5*m0[i], Tc[i], nmag[i]) for i in range(len(m0))]
    for i in range(len(m0)):
        print("*moment @ T=0K: {0} (mu_B)".format(m0[i]))
        print("*lambda: {0} (T)".format(lamb[i]*mu_B**2))
    m0max = np.max(m0)
    Tcmax = np.max(Tc)
    tempc0A = np.linspace(1e-3,1.,2*ndiv,endpoint=False)
    tempc0B = np.linspace(1.,0.25,ndiv)
    tempc1A = np.sqrt(tempc0A)
    tempc1B = 2. - np.sqrt(tempc0B)
    tempc = np.append(tempc1A.tolist(), tempc1B.tolist())
    temp = Tcmax * tempc
    iTc = [np.where(temp>Tci)[0][0] for Tci in Tc]
    magf = []
    for mag in magfield:
        magf.append(np.linspace(0.,mag,ndivf).tolist())
    magf = np.array(magf)
    moment = np.zeros((len(magf),len(magf[0]),3*ndiv))
    for i in range(len(magf)):
        moment[i,:,0] = m0[i]
    for k, magfi in enumerate(magf):
        for j, mf in enumerate(magfi):
            for i in range(len(temp)-1):
                moment[k,j,i+1] = mmom(temp[i], mf, m0[k], lamb[k], moment[k,j,i])

    dGm = np.zeros((len(magf),len(temp)))
    for j in range(len(magf)):
        for i in range(len(temp)):
            dGm[j,i] = -mu_B*eV2Jmol*integrate.simps(moment[j,:,i],magf[j,:])
    
    if ('-p' in args):
        plt.xlabel("Temperature (K)")
        plt.ylabel("Magnetic moment ($\mu_B$)")
        plt.xlim(0.,1.5*Tcmax)
        plt.ylim(0.,1.1*m0max)
        for i in range(len(magf)):
            plt.plot(temp, moment[i,len(magf[i])-1,:], c=colors[i], label="B={0}T".format(magfield[i]))
        plt.legend(fontsize=12)
        plt.minorticks_on()
        plt.savefig("Moment.pdf")
        plt.show()

        plt.xlabel("Temperature (K)")
        plt.ylabel("$\Delta G_m$ (J/mol)")
        plt.plot([0.,1.5*Tcmax],[0.,0.],c="black",linestyle="dashed",lw=0.5)
        for i in range(len(magf)):
            plt.plot(temp, dGm[i,:], c=colors[i], label="B={0}T".format(magfield[i]))
        plt.legend(fontsize=12)
        plt.xlim(0.,1.5*Tcmax)
        plt.minorticks_on()
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
