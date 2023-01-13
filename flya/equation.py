#a simple equation to find the gas in diffuse Lyman alpha fraction

import numpy as np
import astropy.constants as const
import astropy.units as u

# constants
e =  4.8032e-10 # in esu
f_alpha = 0.416 # oscillator strength for Lyman alpha
m_e = const.m_e.to(u.g).value # electron mass 9.1093837015e-28 grams
m_p = const.m_p.to(u.g).value # proton mass 9.1093837015e-28 grams
lya = 1.21567e-5 # cm Lyman alpha wavelength
c = const.c.to(u.cm/u.s).value # c in cm/s
alpha_A = 4.2e-13 # recombination coefficient of H at T= 10^4 K
Mpc_to_cm = u.Mpc.to(u.cm) # 3.085677581467192e+24
G = const.G.to(u.cm**3/u.g/u.s**2).value # G value

# few values
rho_crit = 3* (100* u.km.to(u.cm)/Mpc_to_cm)**2/(8*np.pi*G)
print('Rho crit', rho_crit, '* h^2')
nu_lya =  c/lya # Hz Lyman alpha frequency

A =  np.pi * e**2 * f_alpha * alpha_A * rho_crit**2 / (m_e*nu_lya*m_p**2)
print('A = ', A)

# checks
#val = A*4**6*1.1667*(1-0.25)**2*0.0125**2/1e-12/(100* u.km.to(u.cm)/Mpc_to_cm)
#print(val) # expected to be around 0.9 from croft et al. 1998 paper at z = 3
#val = A*1.1667*(1-0.2454)**2*0.0224**2/1e-14/(100* u.km.to(u.cm)/Mpc_to_cm)
#print(val) # similar as above but at z = 0 with new cosmology

#---- function to calculte the H(z)
def Hz_flat(O_m, O_lambda, H0, z = 0.1):
    """
    H0 in std units km/s/Mpc
    """
    hz = H0*u.km.to(u.cm)/Mpc_to_cm *np.sqrt(O_lambda + O_m*(1+z)**3)
    hz_100 = 100*hz/H0 # in units of h

    return hz, hz_100


# assumed prams
z = 0.1
Gamma_HI = 10e-14
T0 = 4038 # K
gamma = 1.53
tau_avg = 0.0545

# cosmo params assumed
O_m = 0.3089
O_lambda = 0.6911
O_bh2 = 0.0223
H0 = 67.74
y_p = 0.24

# electron density correction
kHe = (2-y_p)/(2-2*y_p)
# beta
beta = 2 - 0.7*(gamma-1)

# formula
hz, hz_100 = Hz_flat(O_m = O_m, O_lambda = O_lambda, H0 = H0, z=z)

fLya = (Gamma_HI*hz/(A*kHe))**0.5 * tau_avg**(beta/2) *(T0/10000)**(0.35) / ( O_bh2 * (1-y_p) * (1+z)**3 )

print('f_Lya = ', fLya)







