import astropy.table as tab
import numpy as np
from scipy import interpolate
np.random.seed(42)


def rebin_flux (flux, wave, velocity, dv_to_rebin = 3):
    """
    rebinning the model spectrum if the forward models are not based on the the observational data
    :param flux: flux array
    :param wave: wavelength array
    :param velocity: velocity array
    :param dv_to_rebin: rebin velocity width in km/s
    :return: rebinned flux wavelength and velocity arrays
    """

    # rebinned velocity arraya
    new_velocity_array = np.arange(velocity[0], velocity[-1], dv_to_rebin)

    # rebin flux
    flux_from_velocity = interpolate.interp1d(velocity, flux)
    flux_new = flux_from_velocity(new_velocity_array)

    # rebin_wavelength
    wavelength_from_velocity = interpolate.interp1d (velocity, wave)
    new_wave_array = wavelength_from_velocity (new_velocity_array)

    # replace zero or negative by 1e-310 number
    flux_new[flux_new <= 0] =  1e-310

    return flux_new, new_velocity_array, new_wave_array



def uniform_grid(data, wave_to_grid):
    """
    # first read the data, remove nan values, write interpolation funciton for various versions of flux
    # then use fill value = nan for the points where wavelength is not covered
    # figure out the mask issue
    'Flux',
    'Noise',
    'Mask',
    'Flux_nonoise',
    'Flux_nonoise_infres',
    'corresponding_data_flux',

    Args:
        data:
        wave_to_grid:

    Returns
    """
    flux = interpolate.interp1d(data['Wave'][i], data['Flux'][i], fill_value = NaN)
    flux_nonoise = interpolate.interp1d(data['Wave'][i], data['Flux_nonoise'][i])
    flux_nonoise_infres = interpolate.interp1d(data['Wave'][i], data['Flux_nonoise_infres'][i])
    data_flux = interpolate.interp1d(data['Wave'][i], data['corresponding_data_flux'][i])
    noise = interpolate.interp1d(data['Wave'][i], data['Noise'][i])
    mask = interpolate.interp1d(data['Wave'][i], data['Mask'][i])


    return uniform_data


def find_indices(list_to_check, item_to_find):
    array = np.array(list_to_check)
    indices = np.where(array == item_to_find)[0]
    return list(indices)



filename = 'forward_model_igm_danforth_sn_34.fits'
d = tab.Table.read (filename)

# find the unique QSOs in fwd models
object_list = list(d['corresponding_data_obj'])
obj = list(set(object_list))
print('set of qso:', obj)

# find minimum number of instances of for all QSOs
# and to find the largest chunk of common wavelength for stacking
number_of_objects = []
wave_min = []
wave_max = []
wave_diff = []
for i in obj:
    number_of_objects.append(len (d[d['corresponding_data_obj'] == i]))
    print(number_of_objects)
    index = object_list.index(i)
    wave = d['Wave'][index]
    wave = wave[~np.isnan(wave)]
    wave_min.append(np.min(wave))
    wave_max.append(np.max(wave))
    diff = wave[1]-wave[0]
    wave_diff.append(diff)

min_number_of_objects = np.min(number_of_objects)
print('minimum number of fwd models per qso is', min_number_of_objects)
start_wavelength= np.min(wave_min)
end_wavelength = np.max(wave_max)
wavelength_scale =  np.min(wave_diff)
print('wavelength range for stacking is ', start_wavelength, end_wavelength, ' with Delta wavelength', wavelength_scale)

# create wavelenght array
final_wave = np.arange(start_wavelength, end_wavelength+0.1*wavelength_scale, wavelength_scale)
# 0.1* scale is for making sure that the last element is taken into account

# find sets of all QSO dataset
set_of_indices = np.zeros((len(obj), min_number_of_objects), dtype=np.int8) # integer array
for qso, i in zip(obj, np.arange(len(obj))):
    obj_ind_list = find_indices(object_list, qso)
    # choose random min_number of obj
    selected_indices = np.random.choice(obj_ind_list, min_number_of_objects, replace = False)
    #print(selected_indices)
    set_of_indices[i] =  np.array(selected_indices)

print(set_of_indices)


# stacking

for j in range(min_number_of_objects):
    qso_ind = selected_indices[:, j]     # indices of qso's to stack

    data = tab.Table()
    # create a data table
    for index in qso_ind:
        data = tab.vstack(data, d[index])

    data_to_stack = uniform_grid(data=data, wave_to_grid=final_wave)

    for l in final_wave:






