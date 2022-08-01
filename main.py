from qpcr_plate_layout import *
from qpcr_data_analysis import *

print('Creating plate \n')
plate = main()

print('Loading reaction data from sheet \n')
load_reaction_timeseries(plate, results_df)

print('Fitting logistic curve to reaction data for each well \n')
compute_ct_from_rn(plate)

print('Generating standard curve plots')
plot_standard_curve(plate)