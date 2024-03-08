import numpy
from dadi import Numerics, PhiManip, Integration, Spectrum






def secondary_contact01(params, ns, pts): 
    nu1,nu2,T,Tm,m12,m21 = params
    
    xx = Numerics.default_grid(pts)
                        
    phi = PhiManip.phi_1D(xx)
    phi = PhiManip.phi_1D_to_2D(xx, phi)
    phi = Integration.two_pops(phi, xx, T, nu1, nu2, m12=0, m21=0)
    phi = Integration.two_pops(phi, xx, Tm, nu1, nu2, m12=m12, m21=m21)
    
    fs = Spectrum.from_phi(phi, ns, (xx,xx)) 
    return fs

def secondary_contact02_exp(params, ns, pts): 
    s,nu1,nu2,T,Tm,m12,m21 = params
    
    xx = Numerics.default_grid(pts)
                        
    phi = PhiManip.phi_1D(xx)
    phi = PhiManip.phi_1D_to_2D(xx, phi)
    
    # First, they split, no migration and they are allowed to change Ne(t)
    #
    nu1_func = lambda t: s * (nu1/s)**(t/T)
    nu2_func = lambda t: (1-s) * (nu2/(1-s))**(t/T)
    phi = Integration.two_pops(phi, xx, T, nu1_func, nu2_func, m12=0, m21=0)
    
    # Then, migration is on (secondary contact), but Ne is constant
    #
    phi = Integration.two_pops(phi, xx, Tm, nu1, nu2, m12=m12, m21=m21)
    
    fs = Spectrum.from_phi(phi, ns, (xx,xx)) 
    return fs


def growth(params, ns, pts):
    nu, T = params
    
    xx = Numerics.default_grid(pts)
    phi = PhiManip.phi_1D(xx)
    
    nu_func = lambda t: numpy.exp(numpy.log(nu) * t/T)
    phi = Integration.one_pop(phi, xx, T, nu_func)
    
    fs = Spectrum.from_phi(phi, ns, (xx,))
    return fs

def one_pop(params, ns, pts):
    nu1,T = params
    
    xx = Numerics.default_grid(pts)
    
    phi = PhiManip.phi_1D(xx)
    phi = Integration.one_pop(phi, xx, nu=nu1, T=T)
    
    fs = Spectrum.from_phi(phi, ns, (xx,))
    return fs

def pop_split(params, ns, pts):
    N1,N2,T = params
    """
    
    """
    xx = Numerics.default_grid(pts)

    phi = PhiManip.phi_1D(xx)
    phi = PhiManip.phi_1D_to_2D(xx, phi)
    phi = Integration.two_pops(phi, xx, T, N1, N2, m12=0, m21=0)

    fs = Spectrum.from_phi(phi, ns, (xx,xx))
    return fs

def fixTheta(params, ns, pts):
    nuA,nu2,T = params
    theta1=18000
    """
    
    """
    xx = Numerics.default_grid(pts)

    phi = PhiManip.phi_1D(xx, nu=nuA, theta0=theta1)
    phi = PhiManip.phi_1D_to_2D(xx, phi)
    phi = Integration.two_pops(phi, xx, T, nu1=1, nu2=nu2, theta0=theta1)

    fs = Spectrum.from_phi(phi, ns, (xx,xx))
    return fs

def pop_splitExp(params, ns, pts):
    s,nu1,nu2,T = params
    """
    
    """
    xx = Numerics.default_grid(pts)

    phi = PhiManip.phi_1D(xx)
    phi = PhiManip.phi_1D_to_2D(xx, phi)
    
    nu1_func = lambda t: s * (nu1/s)**(t/T)
    nu2_func = lambda t: (1-s) * (nu2/(1-s))**(t/T)
    
    phi = Integration.two_pops(phi, xx, T, nu1_func, nu2_func, m12=0, m21=0) 
    
    fs = Spectrum.from_phi(phi, ns, (xx,xx)) 
    return fs



def bottleneck(params, ns, pts):
    nu1,nuB,nuF,Tb,Tf = params
    
    xx = Numerics.default_grid(pts)

    phi = PhiManip.phi_1D(xx)
    phi = PhiManip.phi_1D_to_2D(xx, phi)
    
    phi = Integration.two_pops(phi, xx, Tb, nu1, nuB, m12=0, m21=0)
    phi = Integration.two_pops(phi, xx, Tf, nu1, nuF, m12=0, m21=0)
    
    fs = Spectrum.from_phi(phi, ns, (xx,xx))
    return fs


def bottleneck_twoSides(params, ns, pts):
    nu1,nuB1,nuB,nuF,Tb,Tf = params
    
    xx = Numerics.default_grid(pts)

    phi = PhiManip.phi_1D(xx)
    phi = PhiManip.phi_1D_to_2D(xx, phi)
    
    phi = Integration.two_pops(phi, xx, Tb, nuB1, nuB, m12=0, m21=0)
    phi = Integration.two_pops(phi, xx, Tf, nu1, nuF, m12=0, m21=0)
    
    fs = Spectrum.from_phi(phi, ns, (xx,xx))
    return fs


def bottleneck2(params, ns, pts):
    nu1,nuB,nuF,T,Tb,Tf = params
    
    xx = Numerics.default_grid(pts)
    
    phi = PhiManip.phi_1D(xx)
    phi = PhiManip.phi_1D_to_2D(xx, phi)
    
    phi = Integration.two_pops(phi, xx, T, nu1, nuF, m12=0, m21=0)
    
    phi = Integration.two_pops(phi, xx, Tb, nu1, nuB, m12=0, m21=0)
    phi = Integration.two_pops(phi, xx, Tf, nu1, nuF, m12=0, m21=0)
    
    fs = Spectrum.from_phi(phi, ns, (xx,xx))
    return fs

def IM(params, ns, pts): 
    s,nu1,nu2,T,m12,m21 = params
    
    xx = Numerics.default_grid(pts)
    
    phi = PhiManip.phi_1D(xx)
    phi = PhiManip.phi_1D_to_2D(xx, phi)

    nu1_func = lambda t: s * (nu1/s)**(t/T)
    nu2_func = lambda t: (1-s) * (nu2/(1-s))**(t/T)

    phi = Integration.two_pops(phi, xx, T, nu1_func, nu2_func,m12=m12, m21=m21)
    
    fs = Spectrum.from_phi(phi, ns, (xx,xx)) 
    return fs



def threePops_1(params, ns, pts):
    n1, n2, n3, T1, T2 = params
    """
    """
    xx = Numerics.default_grid(pts)

    phi = PhiManip.phi_1D(xx)

    phi = PhiManip.phi_1D_to_2D(xx, phi) 
    phi = Integration.two_pops(phi, xx, T1, nu1=n1, nu2=n2, m12=0, m21=0)
    
    phi = PhiManip.phi_2D_to_3D_split_2(xx, phi)
    phi = Integration.three_pops(phi, xx, T2, nu1=n1, nu2=n2, nu3=n3)

    fs = Spectrum.from_phi(phi, ns, (xx,xx,xx))
    return fs


def threePops_1_reverse(params, ns, pts):
    n1, n2, n3, T1, T2 = params
    """
    """
    xx = Numerics.default_grid(pts)

    phi = PhiManip.phi_1D(xx)

    phi = PhiManip.phi_1D_to_2D(xx, phi) 
    phi = Integration.two_pops(phi, xx, T1, nu1=n1, nu3=n3, m12=0, m21=0)
    
    phi = PhiManip.phi_2D_to_3D_split_3(xx, phi)
    phi = Integration.three_pops(phi, xx, T2, nu1=n1, nu2=n2, nu3=n3)

    fs = Spectrum.from_phi(phi, ns, (xx,xx,xx))
    return fs





def threePops_bot1(params, ns, pts):
    n1, n2, n3, n2b, n3b, T1, T2, Tb1, Tb2 = params
    """
    """
    xx = Numerics.default_grid(pts)

    phi = PhiManip.phi_1D(xx)

    phi = PhiManip.phi_1D_to_2D(xx, phi) 
    phi = Integration.two_pops(phi, xx, Tb1, nu1=n1, nu2=n2b, m12=0, m21=0)
    phi = Integration.two_pops(phi, xx, T1, nu1=n1, nu2=n2, m12=0, m21=0)
    
    phi = PhiManip.phi_2D_to_3D_split_2(xx, phi)
    phi = Integration.three_pops(phi, xx, Tb2, nu1=n1, nu2=n2, nu3=n3b)
    phi = Integration.three_pops(phi, xx, T2, nu1=n1, nu2=n2, nu3=n3)

    fs = Spectrum.from_phi(phi, ns, (xx,xx,xx))
    return fs



def threePops_bot1_reverse(params, ns, pts):
    n1, n2, n3, n2b, n3b, T1, T2, Tb1, Tb2 = params
    """
    """
    xx = Numerics.default_grid(pts)

    phi = PhiManip.phi_1D(xx)

    phi = PhiManip.phi_1D_to_2D(xx, phi) 
    phi = Integration.two_pops(phi, xx, Tb1, nu1=n1, nu3=n3b)
    phi = Integration.two_pops(phi, xx, T1, nu1=n1, nu3=n3)
    
    phi = PhiManip.phi_2D_to_3D_split_3(xx, phi)
    phi = Integration.three_pops(phi, xx, Tb2, nu1=n1, nu2=n2b, nu3=n3)
    phi = Integration.three_pops(phi, xx, T2, nu1=n1, nu2=n2, nu3=n3)

    fs = Spectrum.from_phi(phi, ns, (xx,xx,xx))
    return fs


def threePops_bot2(params, ns, pts):
    n1, n2, n3, n2b, n3b, n1b, T1, T2, Tb1, Tb2, Tb3 = params
    """
    """
    xx = Numerics.default_grid(pts)

    phi = PhiManip.phi_1D(xx)

    phi = PhiManip.phi_1D_to_2D(xx, phi) 
    phi = Integration.two_pops(phi, xx, Tb1, nu1=n1, nu2=n2b, m12=0, m21=0)
    phi = Integration.two_pops(phi, xx, T1, nu1=n1, nu2=n2, m12=0, m21=0)
    
    phi = PhiManip.phi_2D_to_3D_split_2(xx, phi)
    phi = Integration.three_pops(phi, xx, Tb2, nu1=n1, nu2=n2, nu3=n3b)
    phi = Integration.three_pops(phi, xx, T2, nu1=n1, nu2=n2, nu3=n3)
    phi = Integration.three_pops(phi, xx, Tb1, nu1=n1b, nu2=n2, nu3=n3)

    fs = Spectrum.from_phi(phi, ns, (xx,xx,xx))
    return fs






def threePops_botExp(params, ns, pts):
    n1, n2, n2b, n21, n3, n3b, T1, Tb1, T2, Tb2 = params
    """
    """
    xx = Numerics.default_grid(pts)

    phi = PhiManip.phi_1D(xx)
    phi = PhiManip.phi_1D_to_2D(xx, phi) 
    
    phi = Integration.two_pops(phi, xx, Tb1, nu1=n1, nu2=n2b)
    phi = Integration.two_pops(phi, xx, T1, nu1=n1, nu2=n21)
    
    n2_func = lambda t: n21 * (n2/n21)**(t/T2)
    n3_func = lambda t: n3b * (n3/n3b)**(t/T2)
    
    phi = PhiManip.phi_2D_to_3D_split_2(xx, phi)
    phi = Integration.three_pops(phi, xx, Tb2, nu1=n1, nu2=n21, nu3=n3b)
    phi = Integration.three_pops(phi, xx, T2, nu1=n1, nu2=n2_func, nu3=n3_func)

    fs = Spectrum.from_phi(phi, ns, (xx,xx,xx))
    return fs















def simple_canary((nuIb0, nuIb, nuMo0, nuMo, nuCa0, TIbMo, TMoCa, TF),
    ns, pts):
    """
    nuIb0: starting Iberian population size
    nuIb: ending Iberian population size
    nuMo0: starting Moroccan population size
    nuMo: ending Morrocan population size
    nuCa0: starting Canary island population size
    TIbMo: split time for Iberian and Moroccan populations
    TMoCa: split time for Moroccan and Canary island populations
    TF: simulation ending time
    """
    xx = Numerics.default_grid(pts)

    phi = PhiManip.phi_1D(xx)
    phi = Integration.one_pop(phi, xx, TIbMo, nu=nuIb0)

    phi = PhiManip.phi_1D_to_2D(xx, phi)
    phi = Integration.two_pops(phi, xx, TMoCa, nu1=nuIb, nu2=nuMo0)

    phi = PhiManip.phi_2D_to_3D_split_2(xx, phi)
    phi = Integration.three_pops(phi, xx, TF, nu1=nuIb, nu2=nuMo, nu3=nuCa0)


    fs = Spectrum.from_phi(phi, ns, (xx,xx,xx))

    return fs

def simple_canary_migration((nuIb0, nuIb, nuMo0, nuMo, nuCa0, TIbMo, TMoCa, mIbCa, TF),
    ns, pts):
    """
    nuIb0: starting Iberian population size
    nuIb: ending Iberian population size
    nuMo0: starting Moroccan population size
    nuMo: ending Morrocan population size
    nuCa0: starting Canary island population size
    TIbMo: split time for Iberian and Moroccan populations
    TMoCa: split time for Moroccan and Canary island populations
    mIbCa: migration rate from Ib to Ca
    TF: simulation ending time
    """
    xx = Numerics.default_grid(pts)

    phi = PhiManip.phi_1D(xx)
    phi = Integration.one_pop(phi, xx, TIbMo, nu=nuIb0)

    phi = PhiManip.phi_1D_to_2D(xx, phi)
    phi = Integration.two_pops(phi, xx, TMoCa, nu1=nuIb, nu2=nuMo0)

    phi = PhiManip.phi_2D_to_3D_split_2(xx, phi)
    phi = Integration.three_pops(phi, xx, TF, nu1=nuIb, nu2=nuMo, nu3=nuCa0, m13=mIbCa)

    fs = Spectrum.from_phi(phi, ns, (xx,xx,xx))

    return fs

def two_step((nMo, nCa, T1, T2, T3), ns, pts):
    """
    nMo: Moroccan constant population size
    nCa: Canary Island constant population size
    T1: time for first split (Iberia to Iberia and Morocco)
    T2: time for second split (Morocco to Morocco and Canary)
    """
    xx = Numerics.default_grid(pts)

    phi = PhiManip.phi_1D(xx)
    phi = Integration.one_pop(phi, xx, T1, nu=1)

    phi = PhiManip.phi_1D_to_2D(xx, phi)
    phi = Integration.two_pops(phi, xx, T2, nu1=1, nu2=nMo)

    phi = PhiManip.phi_2D_to_3D_split_2(xx, phi)
    phi = Integration.three_pops(phi, xx, T3, nu1=1, nu2=nMo, nu3=nCa)

    fs = Spectrum.from_phi(phi, ns, (xx, xx, xx))

    return fs

def two_step_migration((nMo, nCa, mIbCa, T1, T2), ns, pts):
    """
    nMo: Moroccan constant population size
    nCa: Canary Island constant population size
    T1: time for first split (Iberia to Iberia and Morocco)
    T2: time for second split (Morocco to Morocco and Canary)
    """
    xx = Numerics.default_grid(pts)

    phi = PhiManip.phi_1D(xx)

    phi = PhiManip.phi_1D_to_2D(xx, phi)
    phi = Integration.two_pops(phi, xx, T1, nu1=1, nu2=nMo)

    phi = PhiManip.phi_2D_to_3D_split_2(xx, phi)
    phi = Integration.three_pops(phi, xx, T2, nu1=1, nu2=nMo, nu3=nCa, m13=mIbCa)

    fs = Spectrum.from_phi(phi, ns, (xx, xx, xx))

    return fs

def Mo_Ib((nIb, nMo1, nMo2, T1), ns, pts):
    """
    nIb: Iberian constant population size
    nMo1: starting Moroccan size
    nMo2: ending Moroccan size
    T1: Time to first split
    This model doesn't consider the Canary Island population size or split time.
    """

    xx = Numerics.default_grid(pts)

    phi = PhiManip.phi_1D(xx)

    phi = PhiManip.phi_1D_to_2D(xx, phi)
    phi = Integration.two_pops(phi, xx, T1, nu1=nIb, nu2=nMo1)

    nu_funcMo = lambda t: (nMo2-1) * (t/T1) + 1

    phi = PhiManip.phi_2D_to_3D_split_2(xx, phi)
    phi = Integration.three_pops(phi, xx, 1, nu1=nIb, nu2=nu_funcMo, nu3=1, m13=0.2)

    fs = Spectrum.from_phi(phi, ns, (xx, xx, xx))

    return fs

def moroccoIberia_migration_Andrea((nuIb, nuMo, TIbMo, mIbMo, mMoIb),
    ns, pts):
    """
    nuIb: Iberian population size
    nuMo: Moroccan population size
    TIbMo: split time for Iberian and Moroccan populations
    mIbmo: migration rate from Ib to Mo
    mMoIb: migration rate from Mo to Ib
    """
    xx = Numerics.default_grid(pts)

    phi = PhiManip.phi_1D(xx)

    phi = PhiManip.phi_1D_to_2D(xx, phi)
    phi = Integration.two_pops(phi, xx, TIbMo, nu1=nuIb, nu2=nuMo, m12=mIbMo, m21=mMoIb)

    fs = Spectrum.from_phi(phi, ns, (xx,xx))

    return fs

def moroccoIberia_migration((nuIb, nuMo, nuCa, TIbMo, TMoCa, mIbMo, mMoIb, mMoCa),
    ns, pts):
    """
    nuIb: Iberian population size
    nuMo: Moroccan population size
    nuCa: Canary population size
    TIbMo: split time for Iberian and Moroccan populations
    TMoCa: split time for Moroccan and Canary populations
    mIbmo: migration rate from Ib to Mo
    mMoIb: migration rate from Mo to Ib
    mMoCa: migration rate from Mo to Ca
    """
    xx = Numerics.default_grid(pts)

    phi = PhiManip.phi_1D(xx)

    phi = PhiManip.phi_1D_to_2D(xx, phi)
    phi = Integration.two_pops(phi, xx, TIbMo, nu1=nuIb, nu2=nuMo, m12=mIbMo, m21=mMoIb)

    phi = PhiManip.phi_2D_to_3D_split_2(xx, phi)
    phi = Integration.three_pops(phi, xx, TMoCa, nu1=nuIb, nu2=nuMo, nu3=nuCa, m13=mMoCa)

    fs = Spectrum.from_phi(phi, ns, (xx,xx,xx))

    return fs
