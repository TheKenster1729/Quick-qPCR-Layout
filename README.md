# Quick qPCR Layout
Quickly generate plate layouts for qPCR

The program reads in data from the params.py file, the only file the user will (hopefully) have to edit, and uses it to generate a representation of a qPCR plate. This representation is then exported as a .csv to be used as the plate template for a ThermoFischer QuantStudio 5 compatible qPCR machine. The output data from the machine can then be loaded in to the program for further analysis.

## qPCR Analysis Pipeline
1. Clone the repo to your computer.
2. To create a layout for a plate, enter the relevant information for your samples in the 'params.py' file, and save. If you are only creating a plate layout template, you may set the parameter 'RESULTS_SPREADSHEETS_FILENAME' to an empty string ("").
3. To generate the .csv template to load into the qPCR machine: open a command line terminal, cd to the local copy of the repo, and enter the following command: 'python qpcr_plate_layout.py'.
4. To analyze output data from the machine, save the output file in the local repo folder. Enter the relevant parameters for the sample in 'params.py', save, and run the following command: 'python main.py'. Standard curve plots should be displayed automatically.

## Example
