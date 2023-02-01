import os
os.environ["NUMEXPR_NUM_THREADS"] = "1"
os.environ["OMP_NUM_THREADS"] = "1"
os.environ["VECLIB_MAXIMUM_THREADS"] = "1"
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"
#-----------------
import matplotlib as mpl
mpl.use("Agg")
import astropy.table as tab
import numpy as np
from pypowerspec import forward_model
from pypowerspec import model_readin
from pypowerspec import data_readin
import multiprocessing as mp
os.nice(10)

def do_danforth_fwd_z003(path, Gamma_12, SN =60, nskew = 1000):
    datadir = '/mnt/quasar/vikram/Danforth_data/mw_format/'


    modelfile = path + '/' + 'ran_skewers_z003_random_OVT_tau_Gamma_{:0.5f}_Nran_010000_seed_1.fits'.format(Gamma_12)

    print("performing forward models in", modelfile)
    param = tab.Table.read(modelfile, hdu =1)
    print('expected Gamma_12', Gamma_12, 'in file', param['GAMMA'])

    forward_model_filename = path + '/flya/danforth_sn/forward_model_igm_danforth_sn_{:.0f}.fits'.format(SN)

    dirname = os.path.dirname(forward_model_filename)
    if not os.path.exists(dirname):
        print('creating: ', dirname)
        os.mkdir(dirname)


    model = model_readin.model_readin(filename=modelfile, nskew = nskew)
    # use nskew = for selecting number of skewers to fit

    data, resolution_array, z_qso_array = data_readin.data_readin(
        use_metalmasking=False, min_z=0.005, max_z=0.06, dataset='COS_data', minsn =60, path=datadir, fill_with_noise=True, use_emissionmasking=True)

    model_all_effects, model_in_lengthened_only, res_array, ind = forward_model.forward_model(model,data=data, SN = SN, data_resolution=resolution_array)

    model_all_effects.write(forward_model_filename, overwrite = True )

    return

# for more qso's
SN = 34

# for tng
path ='/mnt/quasar/vikram/Illustris_z003/get_Gamma_HI'
Gamma_12 = 0.04
do_danforth_fwd_z003(path= path, Gamma_12 = Gamma_12, SN=SN)

# for ill
path ='/mnt/quasar/vikram/Illustris_z003/old/get_Gamma_HI'
Gamma_12 = 0.04
do_danforth_fwd_z003(path= path, Gamma_12 = Gamma_12, SN=SN)
