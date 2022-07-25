from random import sample
import pandas as pd
import numpy as np
from http.server import SimpleHTTPRequestHandler, HTTPServer
from IPython.display import display
from itertools import product
import csv

# contains metadata for each of the wells of the plate
# very useful when exporting to QuantStudio plate template
class Well:
    def __init__(self, sample_name: str, quantity: float, biogroup_name = '', target_name = '', task = '', 
        comments = '', reporter = 'FAM', quencher = 'NFQ-MGB', biogroup_color = '', units = 'Copies/mL'):
        self.sample_name = sample_name
        self.biogroup_name = biogroup_name
        self.target_name = sample_name if target_name == '' else target_name
        self.task = task
        self.reporter = reporter
        self.quencher = quencher
        self.quantity = quantity
        self.comments = comments
        self.biogroup_color = biogroup_color
        self.units = units
        if self.sample_name == 'NTC':
            self.task = 'NTC'

class qPCR_Plate:
    def __init__(self, num_wells = 96):
        if num_wells == 96:
            self.num_wells = 96
            base_layout = pd.DataFrame(data = np.zeros((9, 13)))
            base_layout.iloc[:, 0] = ['', 'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H']
            base_layout.iloc[0, :] = ['', '1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12']
            base_layout.iloc[1:, 1:] = Well('', 0)
            self.plate = base_layout
        # TODO: add preset for 384-well plate

    def initialize_sample_family(self, family_name, starting_quantity, 
        num_dilutions, dilution_factor = 10, blanks = True, ntcs = True):
        if blanks == True and ntcs == True:
            offset = 2
            wells = np.array([Well('Blank', '')] + [Well('NTC', '')] + [Well(family_name + '_' + '1', starting_quantity, task = 'STANDARD')])
        elif blanks == True and ntcs == False:
            offset = 1
            wells = np.array([Well('Blank', '')] + [Well(family_name + '_' + '1', starting_quantity, task = 'STANDARD')])
        elif blanks == False and ntcs == True:
            offset = 1
            wells = np.array([Well('NTC', '')] + [Well(family_name + '_' + '1', starting_quantity, task = 'STANDARD')])
        else:
            offset = 0
            wells = np.array([Well(family_name + '_' + '1', starting_quantity, task = 'STANDARD')])

        for dil in range(1, num_dilutions):
            previous_well_quantity = wells[dil + offset - 1].quantity
            wells = np.append(wells, Well(family_name + '_{}'.format(dil + 1), 
            previous_well_quantity/dilution_factor, task = 'STANDARD'))
        return wells
        
    def view_plate_sample_names(self):
        plate_with_sample_names = self.plate.iloc[1:, 1:].apply(lambda x: x.apply(lambda y: y.sample_name), axis = 1)
        plate_with_sample_names.index = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H']
        return plate_with_sample_names
        # print(self.plate.to_string(index = False, header = False))

    def view_plate_quantities(self):
        plate_with_quantities = self.plate.iloc[1:, 1:].apply(lambda x: x.apply(lambda y: y.quantity), axis = 1)
        plate_with_quantities.index = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H']
        return plate_with_quantities

    def change_well(self, well_coords, new_value):
        row_location = self.plate[0][self.plate[0] == well_coords[0]].index.values[0]
        self.plate.iloc[row_location, int(well_coords[1:])] = new_value

    def check_row_addition_validity(self, position, wells_to_add):
        try:
            row_location = self.plate[0][self.plate[0] == position[0]].index.values[0]
            column_location = int(position[1:])
            self.plate.iloc[row_location, column_location:column_location + len(wells_to_add)] = wells_to_add
            return True
        except ValueError:
            print('Error, check that proposed wells fit on this row. No modification performed')
    
    def check_replicate_addition_validity(self, position, num_replicates, wells_to_add):
        try:
            row_location = self.plate[0][self.plate[0] == position[0]].index.values[0]
            test_plate = self.plate
            well_block = [wells_to_add] * num_replicates
            column_location = int(position[1:])
            test_plate.iloc[row_location:row_location + num_replicates - 1, column_location:column_location + len(wells_to_add)] = well_block
            return True
        except ValueError:
            print('Error, check that proposed replicates can be accomodated. No modification performed')

    def add_samples(self, position: str, family_name: str, starting_quantity: float, num_subsamples: int, 
        dilution_factor = 10, num_replicates = 1, identifier = 'dil'):
        """Efficiently name an entire row. \n
        It is often the case that each row contains subsamples from the same 
        parent sample, e.g. 'qPCR_barcode_1', and that each of the subsamples
        are just dilutions of the parent sample (so the subsamples would be 
        named qPCR_barcode_1_1, qPCR_barcode_1_2, and so on). With just the 
        name of the parent sample, and the number of subsamples, this 
        function will rename the entire row accordingly.

        Args:
            row(str): the row to be renamed ('A'-'H' on 96-well plate, 
            'A'-'P' on 384-well plate)
            sample_name(str): name of parent sample
            num_subsaples(int): number of subsamples to include
            identifier(str): how the subsamples differ, e.g. the default, 'dil', 
            specifies that the subsamples are dilutions
        """
        row_location = self.plate[0][self.plate[0] == position[0]].index.values[0]
        column_location = int(position[1:])

        if identifier == 'dil':
            wells_to_add = self.initialize_sample_family(family_name, starting_quantity, num_subsamples, dilution_factor)
        elif identifier == None:
            # TODO: add additional presets
            pass
        else:
            # TODO: add additional presets
            pass

        if num_replicates == 1:
            self.plate.iloc[row_location, column_location:column_location + len(wells_to_add)] = wells_to_add
        elif num_replicates > 1:
            self.plate.iloc[row_location:row_location + num_replicates, column_location:column_location + len(wells_to_add)] = wells_to_add
        else:
            raise ValueError("Please input a valid number of replicates to add")
    
    def display_readable_format(self):
        # helper function
        def combine_columns(column_1, column_2):
            joint_column = []
            for x, y in zip(column_1, column_2):
                joint_string ='\n'.join([str(x), str(y)])
                joint_column.append(joint_string)
            if self.num_wells == 96:
                joint_series = pd.Series(data = joint_column, index = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H'])
            if self.num_wells == 384:
                # TODO: implement this when 384 well plate is done
                pass
            return joint_series

        names = self.view_plate_sample_names()
        quantities = self.view_plate_quantities()
        readable_df = names.combine(quantities, combine_columns).to_html()

        return readable_df

    def export_template(self, filename):
        to_export = pd.DataFrame(columns = ['Well', 'Well Position', 'Sample Name', 
        'Biogroup Name', 'Biogroup Color', 'Target Name', 'Task', 'Reporter',
        'Quencher', 'Quantity', 'Comments'])
        pivot_frame = pd.melt(self.plate.iloc[1:, 1:].transpose())
        if self.num_wells == 96:
            to_export['Well'] = [i for i in range(1, 97)]
            letters = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H']
            nums = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12']
            well_positions = [a[0] + a[1] for a in [i for i in product(letters, nums)]]
            to_export['Well Position'] = well_positions
            to_export['Sample Name'] = pivot_frame['value'].apply(lambda x: x.sample_name)
            to_export['Biogroup Name'] = pivot_frame['value'].apply(lambda x: x.biogroup_name)
            to_export['Biogroup Color'] = pivot_frame['value'].apply(lambda x: x.biogroup_color)
            to_export['Target Name'] = pivot_frame['value'].apply(lambda x: x.target_name)
            to_export['Task'] = pivot_frame['value'].apply(lambda x: x.task)
            to_export['Reporter'] = pivot_frame['value'].apply(lambda x: x.reporter)
            to_export['Quencher'] = pivot_frame['value'].apply(lambda x: x.quencher)
            to_export['Quantity'] = pivot_frame['value'].apply(lambda x: x.quantity)
            to_export['Comments'] = pivot_frame['value'].apply(lambda x: x.comments)

        if self.num_wells == 384:
            to_export['Well'] = [i for i in range(1, 385)]

        with open('{}.csv'.format(filename), 'w', newline = '') as f:
            writer = csv.writer(f, delimiter = ',')
            to_export_list = [to_export.columns] + to_export.to_numpy().tolist()
            writer.writerow(['[Sample Setup]'])
            writer.writerows(to_export_list)
            f.close()
        
def pre_processing(mol_mass, sample_conc, desired_conc, sample_vol = 1):
    molecular_conc = sample_conc * 1/(mol_mass * 1.66054e-15)
    final_volume = (sample_vol * molecular_conc)/desired_conc
    water_to_add = final_volume - sample_vol
    print("To achieve an initial concentration of", desired_conc,"molecules/μL in a final volume of",
        f"{final_volume:.2f}", "μL, add", f"{water_to_add:.2f}", "μL of water to", 
        sample_vol, "μL of sample.", "\n")

# PAREMETERS
# Insert all N sample family names here
# i.e. [sample_family_1_name, sample_family_2_name, ..., sample_family_N_name]
SAMPLES = ["qPCR_001_replicate", "qPCR_002_replicate"]

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
DILUTION_SERIES_SIZES = [10, 10]

# Specify dilution factors (amount each subsample is diluted relative to previous)
# i.e. [sample_family_1_dilution_factor, sample_family_2_dilution_factor, ..., sample_family_N_dilution_factor]
DILUTION_FACTORS = [10, 10]

# Specify blank/ntc conditions
# i.e. [(blank_sample_family_1, ntc_sample_family_1), (blank_sample_family_2, ntc_sample_family_2), ...,(blank_sample_family_N, ntc_sample_family_N)]
BLANK_NTC_CONDITIONS = [(True, True), (True, True)]

if __name__ == "__main__":
    plate = qPCR_Plate()

    plate.add_samples('A1', 'qPCR_001_replicate', 10**10, 10, num_replicates = 3)
    plate.add_samples('F1', 'qPCR_002_replicate', 10**10, 10, num_replicates = 3)

print(plate.view_plate_sample_names())
print(plate.view_plate_quantities())

plate.export_template('qPCR_001_replicate_template')
# st.write(plate.display_readable_format())

# pre_processing(78346, 15.2, 1e10, sample_vol = 1)