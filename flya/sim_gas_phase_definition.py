import astropy.table as tab
import numpy as np


def get_gasphase(T, oden, density_cut = 120, Temp_cut = 1e5 ):
    # total baryons
    totb=np.sum(oden)

    #WHIM
    whim=oden[T>Temp_cut]
    whim=whim[whim<density_cut]
    gwhim=np.sum(whim)/totb
    print('rho-T', density_cut, Temp_cut, gwhim, 'WHIM')

    #Diffuse gas
    diff=oden[T<Temp_cut]
    diff=diff[diff<density_cut]
    gdiff=np.sum(diff)/totb
    print('rho-T', density_cut, Temp_cut, gdiff, 'diffuse-Lya')

    # Condensed gas
    cond=oden[T<Temp_cut]
    cond=cond[cond>density_cut]
    gcond=np.sum(cond)/totb
    print(gcond, 'condensed')

    # hot halo
    hh=oden[T>Temp_cut]
    hh=hh[hh>density_cut]
    ghh=np.sum(hh)/totb
    print(ghh,  'hot-halo')

    return gdiff, gwhim, ghh, gcond


def prep_input(ovt_file, simname):
    print(ovt_file)
    data = tab.Table.read(ovt_file, hdu=2)
    T = np.array(data['T']).flatten()
    oden = np.array(data['ODEN']).flatten()

    # Dave defination (Delta < 120, T< 1e5)
    print('Dave defination (Delta < 120, T< 1e5)')
    diff, whim, hh, cond = get_gasphase(T=T, oden=oden, density_cut=120)

    print('Smith & Shull definition (Delta < 1000, T< 1e5)')
    diff, whim, hh, cond = get_gasphase(T=T, oden=oden, density_cut=1000)

    return

tng_ovt = '/mnt/quasar/vikram/Illustris_z003/igm/ran_skewers_z003_random_OVT_tau.fits'
prep_input(ovt_file=tng_ovt, simname= 'tng')

ill_ovt = '/mnt/quasar/vikram/Illustris_z003/old/igm/ran_skewers_z003_random_OVT_tau.fits'
prep_input(ovt_file=ill_ovt, simname= 'ill')

