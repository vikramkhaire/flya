import astropy.table as tab
import numpy as np


def get_gasphase(T, oden, density_cut = 100, Temp_cut = 1e5 ):


    # total baryons
    totb=np.sum(oden)

    #WHIM
    whim=oden[T>Temp_cut]
    whim=whim[whim<density_cut]
    gwhim=np.sum(whim)/totb
    print(gwhim, 'WHIM')

    #Diffuse gas
    diff=oden[T<Temp_cut]
    diff=diff[diff<density_cut]
    gdiff=np.sum(diff)/totb
    print(gdiff, 'diffuse-Lya')

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


def calculate_and_write(density_cut, temperature_cut):
    # the code
    # ovt files
    ill_data = '/mnt/quasar/vikram/Illustris_z01/max_skewers_cut/igm/ran_skewers_z01_random_OVT.fits'
    old_ill_data = '/mnt/quasar/vikram/Illustris_z01/old_Illustris/max_skewers_cut/igm/ran_skewers_z01_random_OVT.fits'

    cut_statement = 'for the density {} and temperature {}'.format(density_cut, temperature_cut)
    # Ill old
    diff, whim, hot_halo, condensed = get_gasphase(old_ill_data, density_cut=density_cut, Temp_cut=temperature_cut)
    out_statement1 = r'for Old Illustris (diff, whim, hot_halo, condensed): {} {} {} {}'.format(diff, whim, hot_halo,
                                                                                                condensed)
    print(out_statement1)

    # TNG old
    diff, whim, hot_halo, condensed = get_gasphase(ill_data, density_cut=density_cut, Temp_cut=temperature_cut)
    out_statement2 = r'for TNG Illustris (diff, whim, hot_halo, condensed): {} {} {} {}'.format(diff, whim, hot_halo,
                                                                                                condensed)
    print(out_statement2)

    file.write(cut_statement)
    file.write("\n")
    file.write(out_statement1)
    file.write("\n")
    file.write(out_statement2)
    file.write("\n")

    return



def prep_input(ovt_file, simname):
    data = tab.Table.read(ovt_file, hdu=2)
    T = np.array(data['T']).flatten()
    oden = np.array(data['ODEN']).flatten()

    Temp = 1e5
    output_file =  'sim_gas_phase_{}_T1e5.fits'.format(simname)

    rho_array = np.arange(11)*10+50

    diff_array  = []
    whim_array = []
    hh_array = []
    cond_array = []
    for rho in rho_array:
        diff, whim, hh, cond = get_gasphase(T=T, oden=oden, Temp_cut=Temp, density_cut=rho)
        diff_array.append(diff)
        whim_array.append(whim)
        hh_array,append(hh)
        cond_array.append(cond)

    res = tab.Table([rho_array, [Temp]*len(rho_array), diff_array, whim_array, hh_array, cond_array],
                    names=('T', 'oden', 'diff', 'whim', 'hh', 'cond'))

    res.write(output_file, overwrite =  True)
    print(res)

    return

tng_ovt = '/mnt/quasar/vikram/Illustris_z01/max_skewers_cut/igm/ran_skewers_z01_random_OVT.fits'
prep_input(ovt_file=tng_ovt, simname= 'tng')

"""
ill_data = '/mnt/quasar/vikram/Illustris_z01/max_skewers_cut/igm/ran_skewers_z01_random_OVT.fits'
old_ill_data = '/mnt/quasar/vikram/Illustris_z01/old_Illustris/max_skewers_cut/igm/ran_skewers_z01_random_OVT.fits'

data = tab.Table.read(ovt_file, hdu=2)

T = np.array(data['T']).flatten()
oden = np.array(data['ODEN']).flatten()


#file = open("output_gas_phase_fraction_in_sims.txt", "w+")
file = open("output_gas_phase_fraction_in_sims_for_T_Delta.txt", "w+")

density_cut = 100
temperature_cut = 1e5
calculate_and_write(density_cut, temperature_cut)

density_cut = 120
temperature_cut = 5e4
calculate_and_write(density_cut, temperature_cut)


density_cut = 1
temperature_cut = 3e4
calculate_and_write(density_cut, temperature_cut)

file.close()

"""

