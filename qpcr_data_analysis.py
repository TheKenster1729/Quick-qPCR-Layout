import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from params import *
from qpcr_plate_layout import *

results_df = pd.read_excel(RESULTS_SPREADSHEETS_FILENAME, sheet_name = "Amplification Data", usecols = "A:E", skiprows = 42)

def f(x, a, b, c, d):
    try:
        result = a/(1. + np.exp(-c*(x - d))) + b
    except RuntimeWarning:
        pass
    return result

def logit_curve_fit(rns, num_cycles = 40):
    import scipy.optimize as opt

    cycles = np.arange(1, num_cycles + 1)

    try:
        (a_, b_, c_, d_), _ = opt.curve_fit(f, cycles, rns, method = 'trf')
    except RuntimeError:
        (a_, b_, c_, d_), _ = opt.curve_fit(f, cycles, rns, method = 'lm')

    rn_crit = f(d_, a_, b_, c_, d_)
    # d is the inflection point

    return a_, b_, c_, d_, rn_crit

def compute_ct_from_rn(plate):
    from scipy.optimize import leastsq

    calculate_average_crit_rn(plate)

    mean_rn_crit = plate.mean_rn_crit

    def f_to_solve(x, a, b, c, d, mean_rn_crit):
        return abs(a/(1. + np.exp(-c*(x - d))) + b - mean_rn_crit)
    
    for col in plate.plate.iloc[1:, 1:]:
        for element in plate.plate[col]:
            if type(element) is Well and element.sample_name != '' and element.sample_name != 'Blank':
                a, b, c, d = element.logit_results[:4]
                ct = leastsq(f_to_solve, x0 = d, args = (a, b, c, d, mean_rn_crit))
                element.ct = ct[0][0]

def plot_standard_curve(plate, extra_info = False):
    from sklearn.linear_model import LinearRegression

    # TODO: error bars for highest and lowest values
    for family in SAMPLES:
        family_on_plate = plate.view_plate_sample_names().apply(lambda x: x.str.contains(family)).reset_index(drop = True)
        family_wells = plate.plate.iloc[1:, 1:].reset_index(drop = True)[family_on_plate].dropna(how = 'all').dropna(axis = 1, how = 'all')
        copies = family_wells.apply(lambda x: x.apply(lambda y: y.quantity))
        log_copies = np.log10(copies)
        ct_vals = family_wells.apply(lambda x: x.apply(lambda y: y.ct))
        plt.scatter(log_copies, ct_vals)

        log_copies_for_reg = log_copies.melt(value_vars = log_copies.columns)['value']
        ct_vals_for_reg = ct_vals.melt(value_vars = ct_vals.columns)['value']
        model = LinearRegression()
        model_fitted = model.fit(log_copies_for_reg.to_numpy().reshape(-1, 1), ct_vals_for_reg.to_numpy().reshape(-1, 1))
        plt.plot(log_copies_for_reg, model_fitted.predict(log_copies_for_reg.to_numpy().reshape(-1, 1)))

        score = model_fitted.score(log_copies_for_reg.to_numpy().reshape(-1, 1), ct_vals_for_reg.to_numpy().reshape(-1, 1))*100
        efficiency = abs(model_fitted.coef_[0][0])/3.33*100

        if extra_info:
            # TODO: add significance information (confidence interval and so on)
            pass

        plt.title('Standard Curve for Sample %s \n Effiency = %3f,  R^2 = %3f' % (family, efficiency, score))
        plt.xlabel('Log Copies')
        plt.ylabel('Ct')
        plt.show()

def load_reaction_timeseries(plate, data):
    for well in plate.well_positions:
        if plate.check_well_occupied(well) and plate.access_well(well).sample_name != 'Blank':
            reaction_timeseries = data[data['Well Position'] == well][['Cycle', 'Rn']]
            plate.add_well_timeseries(well, reaction_timeseries)
            plate.access_well(well).logit_results = logit_curve_fit(reaction_timeseries['Rn'], num_cycles = len(reaction_timeseries['Cycle']))

def calculate_average_crit_rn(plate):
    rn_crits = plate.plate.iloc[1:, 1:].apply(lambda x: x.apply(lambda y: y.logit_results[4]))
    mean_rn_crit = rn_crits.mean(skipna = True).mean(skipna = True)
    plate.mean_rn_crit = mean_rn_crit