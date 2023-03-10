# find minimum tau for getting the Lyman alpha fraction

import numpy as np
import astropy.table as tab
from scipy.interpolate import interp1d

def find_dlyaf(simname, mintau_file, z=0.03):

    taumin_dave = 0.02
    #taumin_smith = 0.015

    # refer to notes
    if z ==0.03:
        if simname=='tng':
            dlyaf_dave = 0.3759719791811993
            dlyaf_smith = 0.3882154260965749

        if simname=='ill':
            dlyaf_dave = 0.22038954375541137
            dlyaf_smith = 0.2421547492733484

    if z == 0.1:
        if simname=='tng':
            dlyaf_dave = 0.3846925721047914
            dlyaf_smith = 0.397219265153281

        if simname=='ill':
            dlyaf_dave = 0.23242597144345933
            dlyaf_smith = 0.2558800721245848



    data = tab.Table.read(mintau_file)

    f = interp1d( data['taumin'], data['frac'])

    frac_dave_est =  f(taumin_dave)
    print(dlyaf_dave, frac_dave_est, 'diff = ', (frac_dave_est-dlyaf_dave)*100)

    #frac_smith_est = f(taumin_smith)
    #print(dlyaf_smith, frac_smith_est, 'diff = ', (frac_smith_est-dlyaf_smith)*100)


    return #frac_smith_est, frac_smith_est



def find_mintau(simname, mintau_file, z=0.03):

    # refer to notes
    if z==0.03:
        if simname=='tng':
            dlyaf_dave = 0.3759719791811993
            dlyaf_smith = 0.3882154260965749

        if simname=='ill':
            dlyaf_dave = 0.22038954375541137
            dlyaf_smith = 0.2421547492733484

    if z==0.1:
        if simname=='tng':
            dlyaf_dave = 0.3846925721047914
            dlyaf_smith = 0.397219265153281

        if simname=='ill':
            dlyaf_dave = 0.23242597144345933
            dlyaf_smith = 0.2558800721245848


    data = tab.Table.read(mintau_file)

    f = interp1d(data['frac'], data['taumin'], fill_value =  'extrapolate')

    tau_min_dave =  f(dlyaf_dave)
    tau_min_smith = f(dlyaf_smith)

    return tau_min_dave, tau_min_smith


def prep_input_z01(simname):
    print('for', simname)
    path = '/home/vikram/flya/flya/data'

    Gamma_12_list = [0.01, 0.05, 0.075, 0.10]

    """
    for taumax in [4, 5]:
        for Gamma12 in Gamma_12_list:
            # tau max = 4
            file_name = path + '/' + 'taumin_{}_Gamma_{:0.3f}_taumax_{:0.0f}.fits'.format(simname, Gamma12, taumax)
            dave, smith = find_mintau(simname=simname, mintau_file=file_name)

            print('taumin', dave, smith, 'Gamma12=', Gamma12, 'taumax = ', taumax)
    """

    # -------- for validataion
    for taumax in [4, 5]:
        for Gamma12 in Gamma_12_list:
            # tau max = 4
            file_name = path + '/' + 'taumin_{}_Gamma_{:0.3f}_taumax_{:0.0f}.fits'.format(simname, Gamma12, taumax)
            find_dlyaf(simname=simname, mintau_file=file_name)

    return



def prep_input_z003(simname):

    z = 0.03
    print('for', simname)
    path = '/home/vikram/flya/flya/data'

    Gamma_12_list = [0.007, 0.04, 0.07]

    """
        for taumax in [4, 5]:
        for Gamma12 in Gamma_12_list:
            # tau max = 4
            file_name = path + '/' + 'taumin_{}_Gamma_{:0.3f}_taumax_{:0.0f}_z_{:0.02f}.fits'.format(simname, Gamma12, taumax, z)
            dave, smith = find_mintau(simname=simname, mintau_file=file_name, z=z)

            print('taumin', dave, smith, 'Gamma12=', Gamma12, 'taumax = ', taumax)
    """

    # -------- for validataion
    for taumax in [4]: #, 5]:
        for Gamma12 in Gamma_12_list:
            # tau max = 4
            file_name = path + '/' + 'taumin_{}_Gamma_{:0.3f}_taumax_{:0.0f}_z_{:0.02f}.fits'.format(simname, Gamma12, taumax, z)
            find_dlyaf(simname=simname, mintau_file=file_name, z=0.03)
            


    return

prep_input_z003('tng')
prep_input_z003('ill')

"""
at z= 0.1
taumin_dave = 0.03
taumin_smith = 0.015
These are final values so that
Dave-definition tau_min = 0.03  (Values within 2.6 % for both sims)
For Smith defination tau_min = 0.015  (less than 2% for illustris, 2.4 % for tng)

TODO : Repeat above analysis for z= 0.03
Check with few Nyx simulations 
 ------------------- done!-------------

at z =0.03
taumin_dave = 0.02 (less than 2 % difference in tng and 2.6% in ill)
taumin_smith = 0.015 (less than 2.5 % difference in tng and 1.3% in ill)

"""