#a simple equation to find the gas in diffuse Lyman alpha fraction
import numpy as np
import astropy.constants as const
import astropy.units as u
import astropy.table as tab

#---- function to calculte the H(z)
def Hz_flat(O_m, O_lambda, H0, z = 0.1):
    """
    H0 in std units km/s/Mpc
    """
    Mpc_to_cm = u.Mpc.to(u.cm)  # 3.085677581467192e+24

    hz = H0*u.km.to(u.cm)/Mpc_to_cm *np.sqrt(O_lambda + O_m*(1+z)**3)
    hz_100 = 100*hz/H0 # in units of h

    return hz, hz_100


def get_tau_avg (beta, forward_file, tau_limits, dz_limit = None, find_error = True):
    """
    :param beta:
    :param forward_file:
    :param tau_limits:
    :param dz_limit:  a fraction of total dz (i.e 0.1 is 10%)
    :param find_error:
    :return:
    """

    # read tau data
    data = tab.Table.read(forward_file)

    if dz_limit is None:
        # fwd model
        flux = data['Flux']
        flux[flux < 0.0] = 0.0001  # a lower value
        tau = -np.log(flux)
        # sort
        tau[tau < tau_limits[0]] = 0
        tau[tau > tau_limits[1]] = tau_limits[1]
        tau_avg = np.mean(tau ** (1 / beta))

    else:
        max_ind = int(len(data)* dz_limit)
        flux = data['Flux'][:dz_limi+1]

        flux[flux < 0.0] = 0.0001  # a lower value
        tau = -np.log(flux)
        # sort
        tau[tau < tau_limits[0]] = 0
        tau[tau > tau_limits[1]] = tau_limits[1]
        tau_avg = np.mean(tau ** (1 / beta))

    if find_error:
        flux = data['Flux']

        if dz_limit is not None:
            max_ind = int(len(data) * dz_limit)
        else:
            max_ind = int(len(data) * 0.1)

        boot_means = []
        for _ in range(100):
            bootsample = np.random.choice(len(data), size=max_ind, replace=True)
            new_flux = []
            for i in bootsample:
                new_flux.append(flux[i])

            tau = -np.log(new_flux)
            # sort
            tau[tau < tau_limits[0]] = 0
            tau[tau > tau_limits[1]] = tau_limits[1]
            tau_avg = (np.mean(tau ** (1 / beta)))**(beta/2)
            boot_means.append(tau_avg)

        bootmean = np.mean(boot_means)
        bootmean_std = np.std(boot_means)


    # no noise no LSF
    flux  = data['Flux_nonoise_infres']
    flux [flux< 1e-4] = 0.0001 # a lower value
    tau = -np.log(flux)
    # sort
    tau[tau<tau_limits[0]] = 0
    tau[tau>tau_limits[1]] = tau_limits[1]
    tau_avg_perfect = np.mean(tau**(1/beta))

    return tau_avg, tau_avg_perfect, bootmean, bootmean_std



def diffuse_lya_fraction_forward(taufile, forward_file, Gamma_HI = None, simname = 'tng', tau_def = 'Dave'):

    if tau_def == 'Dave':
        tau_limits = [0.03, 4]
        #print('Assuming Dave et al 2010 definitoin')
    else: # Smith
        tau_limits = [0.015, 4]
        #print('Assuming Shull and Smith definition')

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
    #print('z', z)
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

    if simname == 'tng':
        T0 = 4038  # K
        gamma = 1.53
    if simname == 'ill':
        T0 = 4239  # K
        gamma = 1.57

    # note the recombination coeffient has been taken to scale with T^{-0.7}
    beta = 2 - 0.7 * (gamma - 1)

    tau_avg, tau_avg_perfect, bootmean, boot_std  = get_tau_avg(beta=beta, forward_file=forward_file, tau_limits=tau_limits)


    #formula
    hz, hz_100 = Hz_flat(O_m = O_m, O_lambda = O_lambda, H0 = H0, z=z)

    fLya = (Gamma_HI*hz/(A*kHe))**0.5 * tau_avg**(beta/2) *(T0/10000)**(0.35) / ( O_bh2 * (1-y_p) * (1+z)**3 )
    fLya_perfect = (Gamma_HI*hz/(A*kHe))**0.5 * tau_avg_perfect**(beta/2) *(T0/10000)**(0.35) / ( O_bh2 * (1-y_p) * (1+z)**3 )

    fLya_mean = (Gamma_HI*hz/(A*kHe))**0.5 * bootmean *(T0/10000)**(0.35) / ( O_bh2 * (1-y_p) * (1+z)**3 )
    fLya_std = (Gamma_HI*hz/(A*kHe))**0.5 * boot_std *(T0/10000)**(0.35) / ( O_bh2 * (1-y_p) * (1+z)**3 )
    #print(fLya, tau_avg, taufile)

    return fLya, fLya_perfect, fLya_mean, fLya_std



Gamma_12 = 0.05
SN_array = np.arange(21)*5+30


sim = 'tng'
path = '/mnt/quasar/vikram/Illustris_z01/get_Gamma_HI'
tau_file = path + '/' + 'ran_skewers_01_random_OVT_tau_Gamma_{:0.5f}_Nran_010000_seed_42.fits'.format(Gamma_12)

for SN in SN_array:
    fwd_file = path + '/flya' + '/forward_model_igmSN_{:0.0f}_res_cos_LP1.fits'.format(SN)
    flya, flya_perfect, mean, std = diffuse_lya_fraction_forward(taufile=tau_file, forward_file=fwd_file)
    print(flya, flya_perfect, mean, std, sim, 'SN', SN)


sim = 'ill'
path = '/mnt/quasar/vikram/Illustris_z01/old_Illustris/get_Gamma_HI'
tau_file = path + '/' + 'ran_skewers_01_random_OVT_tau_Gamma_{:0.5f}_Nran_010000_seed_42.fits'.format(Gamma_12)


for SN in SN_array:
    fwd_file = path + '/flya' + '/forward_model_igmSN_{:0.0f}_res_cos_LP1.fits'.format(SN)
    flya, flya_perfect, mean, std = diffuse_lya_fraction_forward(taufile=tau_file, forward_file=fwd_file)
    print(flya, flya_perfect, mean, std, sim, 'SN', SN)
    