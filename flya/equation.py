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
rho_crit = 3* (100* u.km.to(u.cm)/Mpc_to_cm)**2/(8*np.pi*G)
print('Rho crit', rho_crit, '* h^2')

# formula