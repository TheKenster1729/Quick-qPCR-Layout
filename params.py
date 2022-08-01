# PAREMETERS
# Insert all N sample family names here
# i.e. [sample_family_1_name, sample_family_2_name, ..., sample_family_N_name]
SAMPLES = ["D5: T1_QS2_004_BC26", "F9: T1_QS2_007_BC17"]

# Specify identification (only 'dil', for dilution series, supported right now)
IDENTIFICATION = ['dil', 'dil']

# Insert starting well for each sample family here (please account for replicates)
# i.e. [sample_family_1_starting_well, sample_family_2_starting_well, ..., sample_family_N_starting_well]
STARTING_WELLS = ["A1", "E1"]

# Insert number of replicates for each sample family here
# i.e. [sample_family_1_replicates, sample_family_2_replicates, ..., sample_family_N_replicates]
NUM_REPLICATES = [3, 3]

# Specify initial concentrations (in molecules/μL)
# i.e. [sample_family_1_initial_conc, sample_family_2_initial_conc, ..., sample_family_N_initial_conc]
INITIAL_CONCS = [10**10, 10**10]

# Specify sizes of dilution series (not including NTC or blank)
# i.e. [sample_family_1_dilution_series_size, sample_family_2_dilution_series_size, ..., sample_family_N_dilution_series_size]
DILUTION_SERIES_SIZES = [7, 7]

# Specify dilution factors (amount each subsample is diluted relative to previous)
# i.e. [sample_family_1_dilution_factor, sample_family_2_dilution_factor, ..., sample_family_N_dilution_factor]
DILUTION_FACTORS = [10, 10]

# Specify blank/ntc conditions
# i.e. [(blank_sample_family_1, ntc_sample_family_1), (blank_sample_family_2, ntc_sample_family_2), ...,(blank_sample_family_N, ntc_sample_family_N)]
BLANK_NTC_CONDITIONS = [(True, True), (True, True)]

# Specify whether to add a spacer well (to reduce risk of NTC contamination)
SPACER_WELL = [True, True]

RESULTS_SPREADSHEETS_FILENAME = "2022-07-25_JS_KC_d5_f9_qpcr_barcodes.xls" # add filename (including .xls) of results spreadsheet

all_variables = dir()