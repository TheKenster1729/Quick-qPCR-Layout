from prettytable import PrettyTable
import pandas as pd
import numpy as np

class qPCR_Plate:
    def __init__(self, num_wells = 96, include_blanks = True, include_NTCs = True):
        if num_wells == 96:
            base_layout = pd.DataFrame(data = np.zeros((9, 13)))
            base_layout.iloc[:, 0] = ['', 'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H']
            base_layout.iloc[0, :] = ['', '1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12']
            self.plate = base_layout
            if include_blanks == True:
                # TODO: set up blank wells and include some way to set how many
                pass
            if include_NTCs == True:
                # TODO: set up NTC wells and include some way to set how many
                pass

    def view_plate(self):
        print(self.plate.to_string(index = False, header = False))

    def change_well(self, well_coords, new_value):
        row_location = self.plate[0][self.plate[0] == well_coords[0]].index.values[0]
        self.plate.iloc[row_location, int(well_coords[1:])] = new_value

    def change_row(self, row: str, sample_name: str, num_subsamples: int, identifier = 'dil'):
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
            # TODO: need to check that there are enough rows to include all subsamples
            identifier(str): how the subsamples differ, e.g. the default, 'dil', 
            specifies that the subsamples are dilutions
        """
        
plate = qPCR_Plate()
plate.change_well('H12', 'S7')
plate.view_plate()

# generate the underlying base layout
base_layout = pd.DataFrame(data = np.zeros((9, 13)))
base_layout.iloc[:, 0] = ['', 'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H']
base_layout.iloc[0, :] = ['', '1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12']


# base_layout[:, 0] = np.array(['', 'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H'])

myTable = PrettyTable(['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12'])
myTable.add_column('', ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H'])


with open('table_page.html', 'w') as f:
    f.write(myTable.get_html_string())
    f.close()