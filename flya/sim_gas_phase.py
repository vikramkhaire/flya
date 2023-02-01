import astropy.table as tab
import numpy as np


def get_gasphase(T, oden, density_cut = 100, Temp_cut = 1e5 ):


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
    #print(gcond, 'condensed')

    # hot halo
    hh=oden[T>Temp_cut]
    hh=hh[hh>density_cut]
    ghh=np.sum(hh)/totb
    #print(ghh,  'hot-halo')

    return gdiff, gwhim, ghh, gcond



def prep_input(ovt_file, simname):
    data = tab.Table.read(ovt_file, hdu=2)
    T = np.array(data['T']).flatten()
    oden = np.array(data['ODEN']).flatten()

    Teme_list = [5e4, 1e5, 5e5, 1e6, 5e6]
    #Teme_list = [1e5]


    for Temp in Teme_list:

        output_file = 'sim_gas_phase_{}_T_{:.0f}.fits'.format(simname, Temp)

        rho_array = np.arange(31) * 10 + 50

        diff_array = []
        whim_array = []
        hh_array = []
        cond_array = []
        for rho in rho_array:
            diff, whim, hh, cond = get_gasphase(T=T, oden=oden, Temp_cut=Temp, density_cut=rho)
            diff_array.append(diff)
            whim_array.append(whim)
            hh_array.append(hh)
            cond_array.append(cond)

        res = tab.Table([rho_array, [Temp] * len(rho_array), diff_array, whim_array, hh_array, cond_array],
                        names=('oden', 'T', 'diff', 'whim', 'hh', 'cond'))

        res.write(output_file, overwrite=True)
        print('for sim ----> ', simname)
        print(res)


    return

tng_ovt = '/mnt/quasar/vikram/Illustris_z01/max_skewers_cut/igm/ran_skewers_z01_random_OVT.fits'
prep_input(ovt_file=tng_ovt, simname= 'tng')

ill_ovt = '/mnt/quasar/vikram/Illustris_z01/old_Illustris/max_skewers_cut/igm/ran_skewers_z01_random_OVT.fits'
prep_input(ovt_file=ill_ovt, simname= 'ill')

