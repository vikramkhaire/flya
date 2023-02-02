import astropy.table as tab
import numpy as np

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
    number_of_objects.append(len (d[object_list== i]))
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
print('wavelength range for stacking is ', start_wavelength, end_wavelength)

print(wave_diff)
