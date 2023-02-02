import astropy.table as tab
import numpy as np

filename = 'forward_model_igm_danforth_sn_34.fits'
d = tab.Table.read (filename)

# find the unique QSOs in fwd models
obj = list(set(d['corresponding_data_obj']))
print('set of qso:', obj)

# find minimum number of instances of for all QSOs
number_of_objects = []
for i in obj:
    number_of_objects.append(len (d[d['corresponding_data_obj']== i]))

min_number_of_objects = np.min(number_of_objects)
print('minimum number of fwd models per qso is', min_number_of_objects)

