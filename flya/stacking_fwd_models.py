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
    Args:
        data:
        wave_to_grid:

    Returns
    """

    # rewrite the data table
    new_data = tab.Table()
    for i in range(len(data)):
        flux = interpolate.interp1d(data['Wave'][i], data['Flux'][i], fill_value=np.nan, bounds_error=False)
        new_flux = flux(wave_to_grid)
        flux_nonoise = interpolate.interp1d(data['Wave'][i], data['Flux_nonoise'][i], fill_value=np.nan, bounds_error=False)
        new_flux_nonoise = flux_nonoise(wave_to_grid)
        flux_nonoise_infres = interpolate.interp1d(data['Wave'][i], data['Flux_nonoise_infres'][i], fill_value=np.nan, bounds_error=False)
        new_flux_nonoise_infres = flux_nonoise_infres(wave_to_grid)
        data_flux = interpolate.interp1d(data['Wave'][i], data['corresponding_data_flux'][i], fill_value=np.nan, bounds_error=False)
        new_data_flux = data_flux(wave_to_grid)
        noise = interpolate.interp1d(data['Wave'][i], data['Noise'][i], fill_value=np.nan, bounds_error=False)
        new_noise = noise(wave_to_grid)
        # masks
        mask = interpolate.interp1d(data['Wave'][i], data['Mask'][i], fill_value=-1, bounds_error=False) # -1 for out of bound
        new_masks = np.array(mask(wave_to_grid), dtype = np.int16)
        tab_line = tab.Table([[wave_to_grid], [new_flux], [new_noise], [new_masks], [new_flux_nonoise], [new_flux_nonoise_infres], [new_data_flux]],
                             names = ('Wave', 'Flux', 'Noise', 'Mask', 'Flux_nonoise', 'Flux_nonoise_infres', 'corresponding_data_flux'))
        new_data = tab.vstack([new_data, tab_line])


    return new_data


def find_indices(list_to_check, item_to_find):
    array = np.array(list_to_check)
    indices = np.where(array == item_to_find)[0]

    return list(indices)


#------------- main code --------

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
set_of_indices = np.zeros((len(obj), min_number_of_objects), dtype=np.int16) # integer array
for qso, i in zip(obj, np.arange(len(obj))):
    print(qso, i )
    obj_ind_list = find_indices(object_list, qso)
    # choose random min_number of obj
    selected_indices = np.random.choice(obj_ind_list, min_number_of_objects, replace = False)
    #print(selected_indices)
    set_of_indices[i, :] =  selected_indices

print(set_of_indices)


# stacking

stacked_data = tab.Table()

#min_number_of_objects = 2
for j in range(min_number_of_objects):
    j= 0
    #qso_ind = selected_indices[:, j]     # indices of qso's to stack
    qso_ind = set_of_indices[:, j]  # indices of qso's to stack

    data = tab.Table()
    # create a data table
    for index in qso_ind:
        print(index)
        #print(d[index])
        data = tab.vstack([data, d[index]])

    data_to_stack = uniform_grid(data=data, wave_to_grid=final_wave)

    stack_wave = []
    stack_flux = []
    stack_noise = []
    stack_flux_nonoise = []
    stack_flux_nonoise_infres = []
    stack_data_flux =[]
    num_count = []

    # this is for the mean method
    # TODO code for other methods (median and SN weighted)
    new_mask_array = []
    new_data = tab.Table()
    for k in range(len(final_wave)):
        #np.count_nonzero(np.isnan(data))
        stack_wave.append(np.nanmean(data_to_stack['Wave'][:, k]))
        stack_flux.append(np.nanmean(data_to_stack['Flux'][:, k]))
        stack_noise.append(np.nanmean(data_to_stack['Noise'][:, k]))
        stack_flux_nonoise.append(np.nanmean(data_to_stack['Flux_nonoise'][:, k]))
        stack_flux_nonoise_infres.append(np.nanmean(data_to_stack['Flux_nonoise_infres'][:, k]))
        stack_data_flux.append(np.nanmean(data_to_stack['corresponding_data_flux'][:, k]))
        # number of no nan elements
        count = np.count_nonzero(np.isnan(data_to_stack['Wave'][:, k]))
        num_count.append(count)
        if -1 in data_to_stack['Mask'][:, k]:
            new_mask_array.append(-1)
        else:
            new_mask_array.append(int(np.ceil(np.nanmean(data_to_stack['Mask'][:, k]))))

    # find ways to deal with masks () --> fix new_masks
    tab_line = tab.Table(
        [[stack_wave], [stack_flux], [stack_noise], [new_mask_array], [stack_flux_nonoise], [stack_flux_nonoise_infres], [stack_data_flux]],
        names=('Wave', 'Flux', 'Noise', 'Mask', 'Flux_nonoise', 'Flux_nonoise_infres', 'corresponding_data_flux'))
    new_data = tab.vstack([new_data, tab_line])


new_data.write('test_stack.fits', overwrite = True)
