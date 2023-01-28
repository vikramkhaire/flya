#a simple equation to find the gas in diffuse Lyman alpha fraction
import numpy as np
import astropy.constants as const
import astropy.units as u
import astropy.table as tab

#---- function to calculte the H(z)
def Hz_flat(O_m, O_lambda, H0, z = 0.03):
    """
    H0 in std units km/s/Mpc
    """
    Mpc_to_cm = u.Mpc.to(u.cm)  # 3.085677581467192e+24

    hz = H0*u.km.to(u.cm)/Mpc_to_cm *np.sqrt(O_lambda + O_m*(1+z)**3)
    hz_100 = 100*hz/H0 # in units of h

    return hz, hz_100

def diffuse_lya_fraction(taufile, Gamma_HI = None, simname = 'tng', tau_limits =[0.1, 4], z_data=0.03):

    # constants
    e = 4.8032e-10  # in esu
    f_alpha = 0.416  # oscillator strength for Lyman alpha
    m_e = const.m_e.to(u.g).value  # electron mass 9.1093837015e-28 grams
    m_p = const.m_p.to(u.g).value  # proton mass 9.1093837015e-28 grams
    lya = 1.21567e-5  # cm Lyman alpha wavelength
    c = const.c.to(u.cm / u.s).value  # c in cm/s
    Mpc_to_cm = u.Mpc.to(u.cm)  # 3.085677581467192e+24
    alpha_A = 4.2e-13  # recombination coefficient of H at T= 10^4 K
    G = const.G.to(u.cm ** 3 / u.g / u.s ** 2).value  # G value

    # few values
    rho_crit = 3 * (100 * u.km.to(u.cm) / Mpc_to_cm) ** 2 / (8 * np.pi * G)
    #print('Rho crit', rho_crit, '* h^2')
    nu_lya = c / lya  # Hz Lyman alpha frequency

    A = np.pi * e ** 2 * f_alpha * alpha_A * rho_crit ** 2 / (m_e * nu_lya * m_p ** 2)
    #print('A = ', A)

    # read cosmology parameters
    cosmo = tab.Table.read(taufile, hdu = 1)

    z = cosmo['z'][0]
    print('z', z)
    O_m = cosmo['Om0'][0]
    O_lambda = cosmo['Ode0'][0]
    h = cosmo['lit_h'] [0]
    O_bh2 = h**2 * cosmo['Ob0'] [0]
    H0 = h *100
    y_p = 1 - cosmo['X'] [0]
    # electron density correction
    kHe = (2 - y_p) / (2 - 2 * y_p)

    if Gamma_HI == None:
        Gamma_HI = 1e-12* cosmo['GAMMA'][0]
    #print('Gamma_HI', Gamma_HI, 'Obh2', O_bh2, 'yp', y_p, 'lambda', O_lambda)

    if z_data==0.03:
        if simname=='tng':
            T0 =  3852 # K
            gamma = 1.53
        if simname=='ill':
            T0 = 3888  # K
            gamma = 1.55

    if z_data == 0.1:
        if simname=='tng':
            T0 = 4010  # K
            gamma = 1.527
        if simname=='ill':
            T0 = 4231  # K
            gamma = 1.56


    # note the recombination coeffient has been taken to scale with T^{-0.7}
    beta = 2 - 0.7 * (gamma - 1)

    # read tau data
    data = tab.Table.read(taufile, hdu = 2)
    tau = data['TAU']

    # sort
    tau[tau<tau_limits[0]] = 0
    tau[tau>tau_limits[1]] = tau_limits[1]

    tau_avg = np.mean(tau**(1/beta))

    #formula
    hz, hz_100 = Hz_flat(O_m = O_m, O_lambda = O_lambda, H0 = H0, z=z)
    fLya = (Gamma_HI*hz/(A*kHe))**0.5 * tau_avg**(beta/2) *(T0/10000)**(0.35) / ( O_bh2 * (1-y_p) * (1+z)**3 )

    print(fLya, tau_avg, taufile)

    return fLya


def prep_input(taufile, Gamma12,  simname, z_data = 0.03):
    print('for', simname, 'Gamma_12', Gamma12)
    tau_array = np.arange(51)*0.001

    # tau max = 4
    taumax= 4
    file_name = 'taumin_{}_Gamma_{:0.3f}_taumax_4_z_{:0.2f}.fits'.format(simname, Gamma12, z_data)
    flya_array = []
    for tau_low in tau_array:
        fraction = diffuse_lya_fraction(taufile=taufile, simname= simname, tau_limits=[tau_low, taumax], z_data=z_data)
        flya_array.append(fraction)

    res = tab.Table([tau_array, flya_array], names=('taumin', 'frac'))
    res.write(file_name, overwrite=True)
    print(res)

    # tau max = 5
    taumax= 5
    file_name = 'taumin_{}_Gamma_{:0.3f}_taumax_5_z_{:0.2f}.fits'.format(simname, Gamma12, z_data)
    flya_array = []
    for tau_low in tau_array:
        fraction = diffuse_lya_fraction(taufile=taufile, simname= simname, tau_limits=[tau_low, taumax], z_data=z_data)
        flya_array.append(fraction)

    res = tab.Table([tau_array, flya_array], names=('taumin', 'frac'))
    res.write(file_name, overwrite=True)
    print(res)

    return


z_data = 0.03
zstr = 'z003'
seed = 1

simname = 'ill'
path = '/mnt/quasar/vikram/Illustris_z003/old/get_Gamma_HI'

taufile = path + '/' +'ran_skewers_'+zstr+'_random_OVT_tau_Gamma_0.01000_Nran_010000_seed_{}.fits'.format(seed)
prep_input(taufile= taufile, Gamma12= 0.004, simname= simname, z_data=z_data)

taufile = path + '/' +'ran_skewers_'+zstr+'_random_OVT_tau_Gamma_0.05000_Nran_010000_seed_{}.fits'.format(seed)
prep_input(taufile= taufile, Gamma12= 0.007, simname= simname, z_data=z_data)

taufile = path + '/' +'ran_skewers_'+zstr+'_random_OVT_tau_Gamma_0.07500_Nran_010000_seed_{}.fits'.format(seed)
prep_input(taufile= taufile, Gamma12= 0.04, simname= simname, z_data=z_data)

taufile = path + '/' +'ran_skewers_'+zstr+'_random_OVT_tau_Gamma_0.10000_Nran_010000_seed_{}.fits'.format(seed)
prep_input(taufile= taufile, Gamma12= 0.07, simname= simname, z_data=z_data)


simname = 'tng'
path = '/mnt/quasar/vikram/Illustris_z003/get_Gamma_HI'

taufile = path + '/' +'ran_skewers_'+zstr+'_random_OVT_tau_Gamma_0.01000_Nran_010000_seed_{}.fits'.format(seed)
prep_input(taufile= taufile, Gamma12= 0.004, simname= simname, z_data=z_data)

taufile = path + '/' +'ran_skewers_'+zstr+'_random_OVT_tau_Gamma_0.05000_Nran_010000_seed_{}.fits'.format(seed)
prep_input(taufile= taufile, Gamma12= 0.007, simname= simname, z_data=z_data)

taufile = path + '/' +'ran_skewers_'+zstr+'_random_OVT_tau_Gamma_0.07500_Nran_010000_seed_{}.fits'.format(seed)
prep_input(taufile= taufile, Gamma12= 0.04, simname= simname, z_data=z_data)

taufile = path + '/' +'ran_skewers_'+zstr+'_random_OVT_tau_Gamma_0.10000_Nran_010000_seed_{}.fits'.format(seed)
prep_input(taufile= taufile, Gamma12= 0.07, simname= simname, z_data=z_data)

