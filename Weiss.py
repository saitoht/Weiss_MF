import numpy as np
import matplotlib.pyplot as plt
import sys, os

R = 8.3144569919459 # J/mol K
kB = 8.6173303e-5 # eV/K
T2J = 5.584937543138665 # J/mol T
mu_B = 5.7883817982e-05 # eV/T

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
    magfield = [0.0, 10.0]
    colors = ['red', 'blue']
    Tc = 1000.
    nmag = 1
    struc = 'BCC'
    
    usage = """ USAGE: python Weiss.py 2.2 [-p] [-h] """
    args = sys.argv
    if ( len(args) < 2 or '-h' in args ):
        print(usage)
        sys.exit()
    m0 = float(args[1])
    print("*moment @ T=0K: {0} (mu_B)".format(m0))
    lamb = Tc2lambda(0.5*m0, Tc, nmag)
    print("*lambda: {0} (T)".format(lamb*mu_B**2))
    temp = np.linspace(1e-3,1.5*Tc,ndiv)
    moment = np.zeros((len(magfield),ndiv))
    moment[:,0] = m0
    for j, mf in enumerate(magfield):
        for i in range(len(temp)-1):
            moment[j,i+1] = mmom(temp[i], mf, m0, lamb, moment[j,i])

    print("*Calculate the specific heat Cm(T,B)")
    tempc = temp/Tc
    Cm = np.array([Cmcalc(tempc[i],mf,moment[j,i],struc,p=2.) for j, mf in enumerate(magfield) for i in range(len(tempc))]).reshape(len(magfield),ndiv)
    print(Cm)
    CmdT = np.array([Cm[j,i]/temp[i] for j in range(len(magfield)) for i in range(len(temp))]).reshape(len(magfield),ndiv)
    dHm = np.array([Simpson(temp,Cm[j,:],i) for j in range(len(magfield)) for i in range(len(temp)-2)]).reshape(len(magfield),ndiv-2)
    dSm = np.array([Simpson(temp,CmdT[j,:],i) for j in range(len(magfield)) for i in range(len(temp)-2)]).reshape(len(magfield),ndiv-2)

    if ('-p' in args):
        plt.xlabel("Temperature (K)")
        plt.ylabel("Magnetic moment ($\mu_B$)")
        for i in range(len(magfield)):
            plt.plot(temp, moment[i,:], c=colors[i], label="B={0}T".format(magfield[i]))
        plt.savefig("Moment.pdf")
        plt.legend()
        plt.show()

        plt.xlabel("Temperature (K)")
        plt.ylabel("$C_m$")
        for i in range(len(magfield)):
            plt.plot(temp, Cm[i,:], c=colors[i], label="B={0}T".format(magfield[i]))
        plt.savefig("Cm.pdf")
        plt.legend()
        plt.show()
        
        plt.xlabel("Temperature (K)")
        plt.ylabel("$\Delta H_m$")
        for i in range(len(magfield)):
            plt.plot(temp[0:ndiv-2], dHm[i,:], c=colors[i], label="B={0}T".format(magfield[i]))
        plt.savefig("dHm.pdf")
        plt.legend()
        plt.show()
        
        plt.xlabel("Temperature (K)")
        plt.ylabel("$\Delta S_m T$")
        mixS = np.array([temp[i]*dSm[j,i] for j in range(len(magfield)) for i in range(len(temp)-2)]).reshape(len(magfield),ndiv-2)
        for i in range(len(magfield)):
            plt.plot(temp[0:ndiv-2], mixS[i,:], c=colors[i], label="B={0}T".format(magfield[i]))
        plt.savefig("dSm.pdf")
        plt.legend()
        plt.show()

        plt.xlabel("Temperature (K)")
        plt.ylabel("$\Delta G_m$")
        for i in range(len(magfield)):
            plt.plot(temp[0:ndiv-2], dHm[i,:]-mixS[i,:], c=colors[i], label="B={0}T".format(magfield[i]))
        plt.savefig("dGm.pdf")
        plt.legend()
        plt.show()


def coth(x):
    if (np.abs(x) < 1e-20):
        return 0.
    else:
        return 1./np.tanh(x)

def Brillouin(j, alpha):
    """ define Brillouin function   Eq. (S4) """
    cf = 0.5*((2*j+1)/j)
    cfp = 0.5/j
    Bf = cf*coth(cf*alpha) - cfp*coth(cfp*alpha)
    return Bf

def mmom(temp, magfield, mmom0, lamb, mmonb):
    """ return magnetic moment at T, B    Eq. (S2) """
    alpha = mmom0*mu_B * (magfield + lamb*mmonb*mu_B)/(kB*temp)  ### Eq. (S3)
    return mmom0 * Brillouin(0.5*mmom0, alpha)

def fs(struc):
    if ( struc == 'BCC' ):
        return 0.285
    elif ( struc == 'FCC' ):
        return 0.105 

def kf(mm, struc):
    cf = 4.*(1.-fs(struc))/(1.-np.exp(-4.))
    return cf*kB*np.log(mm+1.)

def kp(mm, struc, p=2.):
    return 8.*p*fs(struc)*kB*np.log(mm+1.)

def Cmcalc(tempc, magfield, mm, struc, p=2.):
    if (tempc >= 1.):
        return kp(mm,struc,p=p)*tempc*np.exp(8.*p*(1.-tempc))
    else:
        return kf(mm,struc)*tempc*np.exp(-4.*(1.-tempc))

def Simpson(x, f, n):
    return sum([(x[i+1]-x[i])*(f[i]+4.*f[i+1]+f[i+2])/3. for i in range(0,n)])

def Tc2lambda(J, Tc, n, g=2.):
    return 3.*kB*Tc/(J*(J+1.)*n*(g**2)*(mu_B**2))  # T/eV

if __name__=='__main__':
    main()
