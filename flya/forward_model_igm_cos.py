#------------------
# for reducing  numpy threads
import os
os.environ["NUMEXPR_NUM_THREADS"] = "20"
os.environ["OMP_NUM_THREADS"] = "20"
os.environ["VECLIB_MAXIMUM_THREADS"] = "20"
os.environ["OPENBLAS_NUM_THREADS"] = "20"
os.environ["MKL_NUM_THREADS"] = "20"
#-----------------
import matplotlib as mpl
mpl.use("Agg")
import astropy.table as tab
from enigma.whim.forward_model import forward_model
from astropy.cosmology import FlatLambdaCDM
import logging
import numpy as np
import astropy.constants as const
import astropy.units as u
os.nice(10)

def run_forward_z01(simname, outfileFirstName = 'igm', zsim =0.1, dz_total = 1):

    # constants
    c = const.c.to(u.km / u.s).value
    lya = 1215.67
    lya *= u.AA

    # the photoinoization rate in units 10^-12 s^-1 # see notes for choosing this value
    Gamma_12 = 0.05

    if simname == 'tng':
        # simulation file
        hdf5file = '/mnt/quasar/sims/IllustrisTNG/TG100-1/z01/gadget_format/grid/Ill_TG100-1_z0.10.h5'
        # path to catalogfile
        path = '/mnt/quasar/vikram/Illustris_z01'

    if simname == 'ill':
        # simulation file
        hdf5file = '/mnt/quasar/sims/Illustris/Illustris1/z0.1/gadget_format/grid/Ill1_z0.10.h5'
        # path to catalogfile
        path = '/mnt/quasar/vikram/Illustris_z01/old_Illustris'


    outpath = path + '/get_Gamma_HI'

    file = outpath + '/' + 'ran_skewers_01_random_OVT_tau_Gamma_{:0.5f}_Nran_010000_seed_42.fits'.format(Gamma_12)
    print('from:', file)

    params = tab.Table.read(file, hdu = 1)
    #------------ find number of spectra
    # assigning box and cosmology parameters
    ncell = params['Ng'][0]
    OmegaM = params['Om0'][0]
    OmegaL = params['Ode0'][0]
    OmegaB = params['Ob0'][0]
    h = params['lit_h'][0]
    Lbox = params['Lbox'][0]
    Lbox *= u.Mpc
    # not using the astropy units
    if zsim is None:
        z = params['z'][0]
    else:
        z = zsim

    # Computing Velocity Data Points
    cosmo = FlatLambdaCDM(H0=100.0 * h, Om0=OmegaM, Ob0=OmegaB)
    hubblez = cosmo.H(z)
    distance = Lbox / h
    dz = (distance * hubblez / c).value # total dz across one side of the sim

    nmodel = int(dz_total//dz  +1)
    print ('need {} models to get total dz {}'.format(nmodel, dz_total))
    #--------------------------------------------

    model, res = forward_model.model_readin(taufile=file)

    SN_array = np.arange(21)*5+ 30
    for SN in SN_array:
        print('running for SN', SN)
        model_all_effects, model_in_lengthened_only, res_array, ind = forward_model.forward_model(
            model=model, data_resolution=res, nmodel_spectra=nmodel, SN=SN)

        # storing  the forward models
        # setting names for output dir
        outpathname = outpath + '/flya'  # + '/' + 'linefitting_forward_model_SN_{:03d}'.format(int(SN))

        if not os.path.exists(outpathname):
            print('creating: ', outpathname)
            os.mkdir(outpathname)

        outfilename = outpathname + '/forward_model_' + outfileFirstName + 'SN_{:0.0f}_res_cos_LP1.fits'.format(SN)
        model_all_effects.write(outfilename, overwrite=True)

        """
        print('stroring in dir: ', outpathname)
        if not os.path.isdir(outpathname):
            os.mkdir(outpathname)
            print('creating dir: ', outpathname)
            forward_model.store_spectra_for_vpfitting(model_all_effects=model_all_effects, res_array=res_array,
                outpath=outpathname, outfile=outfileFirstName)
        else:
            if  os.path.exists(outfilename):
                print('file exists {} : doing no calculations'.format(outfilename))
            else:
                forward_model.store_spectra_for_vpfitting(model_all_effects=model_all_effects, res_array=res_array,
                    outpath=outpathname, outfile=outfileFirstName)
        """

    return

def run_forward_z003(simname, outfileFirstName = 'igm', zsim =0.03, dz_total = 1):

    # constants
    c = const.c.to(u.km / u.s).value
    lya = 1215.67
    lya *= u.AA

    # the photoinoization rate in units 10^-12 s^-1 # see notes for choosing this value
    Gamma_12 = 0.04

    if simname == 'tng':
        # simulation file
        hdf5file = '/mnt/quasar/sims/IllustrisTNG/TG100-1/z01/gadget_format/grid/Ill_TG100-1_z0.10.h5'
        # path to catalogfile
        path = '/mnt/quasar/vikram/Illustris_z003'

    if simname == 'ill':
        # simulation file
        hdf5file = '/mnt/quasar/sims/Illustris/Illustris1/z0.1/gadget_format/grid/Ill1_z0.10.h5'
        # path to catalogfile
        path = '/mnt/quasar/vikram/Illustris_z003/old'


    outpath = path + '/get_Gamma_HI'

    file = outpath + '/' + 'ran_skewers_z003_random_OVT_tau_Gamma_{:0.5f}_Nran_010000_seed_42.fits'.format(Gamma_12)
    print('from:', file)

    params = tab.Table.read(file, hdu = 1)
    #------------ find number of spectra
    # assigning box and cosmology parameters
    ncell = params['Ng'][0]
    OmegaM = params['Om0'][0]
    OmegaL = params['Ode0'][0]
    OmegaB = params['Ob0'][0]
    h = params['lit_h'][0]
    Lbox = params['Lbox'][0]
    Lbox *= u.Mpc
    # not using the astropy units
    if zsim is None:
        z = params['z'][0]
    else:
        z = zsim

    # Computing Velocity Data Points
    cosmo = FlatLambdaCDM(H0=100.0 * h, Om0=OmegaM, Ob0=OmegaB)
    hubblez = cosmo.H(z)
    distance = Lbox / h
    dz = (distance * hubblez / c).value # total dz across one side of the sim

    nmodel = int(dz_total//dz  +1)
    print ('need {} models to get total dz {}'.format(nmodel, dz_total))
    #--------------------------------------------

    model, res = forward_model.model_readin(taufile=file)

    SN_array = np.arange(21)*5+ 40
    for SN in SN_array:
        print('running for SN', SN)
        model_all_effects, model_in_lengthened_only, res_array, ind = forward_model.forward_model(
            model=model, data_resolution=res, nmodel_spectra=nmodel, SN=SN)

        # storing  the forward models
        # setting names for output dir
        outpathname = outpath + '/flya'  # + '/' + 'linefitting_forward_model_SN_{:03d}'.format(int(SN))

        if not os.path.exists(outpathname):
            print('creating: ', outpathname)
            os.mkdir(outpathname)

        outfilename = outpathname + '/forward_model_' + outfileFirstName + 'SN_{:0.0f}_res_cos_LP1.fits'.format(SN)
        model_all_effects.write(outfilename, overwrite=True)


    return


SimNameArray= ['ill', 'tng']
#--------------multiprocessing routine
import multiprocessing as mp
pool = mp.Pool(processes=2)
results = [pool.apply_async(run_forward_z003, args=(simulation, )) for  simulation in SimNameArray]
output = [p.get() for p in results]
