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
nu_lya =  c/lya

A =  np.pi * e**2 * f_alpha * alpha_A * rho_crit**2 / (m_e*nu_lya*m_p**2)
print('A = ', A)

# check
val = A*4**6*1.1667*(1-0.25)**2*0.0125**2/1e-12/(100* u.km.to(u.cm)/Mpc_to_cm)
print(val)
val = A*1.1667*(1-0.2454)**2*0.0224**2/1e-14/(100* u.km.to(u.cm)/Mpc_to_cm)
print(val)

# formula


