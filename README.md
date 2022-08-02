# Quick qPCR Layout
Quickly generate plate layouts for qPCR.

The program reads in data from the 'params.py' file, the only file the user will (hopefully) have to edit, and uses it to generate a representation of a qPCR plate. This representation is then exported as a .csv to be used as the plate template for a ThermoFischer QuantStudio 5 compatible qPCR machine. The output data from the machine can then be loaded in to the program for further analysis.

## qPCR Analysis Pipeline
1. Clone the repo to your computer.
2. To create a layout for a plate, enter the relevant information for your samples in the 'params.py' file, and save. If you are only creating a plate layout template, you may set the parameter 'RESULTS_SPREADSHEETS_FILENAME' to an empty string ("").
3. To generate the .csv template to load into the qPCR machine: open a command line terminal, cd to the local copy of the repo, and enter the following command: 'python qpcr_plate_layout.py'.
4. To analyze output data from the machine, save the output file in the local repo folder. Enter the relevant parameters for the sample in 'params.py', save, and run the following command: 'python main.py'. Standard curve plots should be displayed automatically.

## Example
To create the following plate layout:
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>1</th>
      <th>2</th>
      <th>3</th>
      <th>4</th>
      <th>5</th>
      <th>6</th>
      <th>7</th>
      <th>8</th>
      <th>9</th>
      <th>10</th>
      <th>11</th>
      <th>12</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>A</th>
      <td>Blank</td>
      <td>NTC</td>
      <td></td>
      <td>D5: T1_QS2_004_BC26_1</td>
      <td>D5: T1_QS2_004_BC26_2</td>
      <td>D5: T1_QS2_004_BC26_3</td>
      <td>D5: T1_QS2_004_BC26_4</td>
      <td>D5: T1_QS2_004_BC26_5</td>
      <td>D5: T1_QS2_004_BC26_6</td>
      <td>D5: T1_QS2_004_BC26_7</td>
      <td></td>
      <td></td>
    </tr>
    <tr>
      <th>B</th>
      <td>Blank</td>
      <td>NTC</td>
      <td></td>
      <td>D5: T1_QS2_004_BC26_1</td>
      <td>D5: T1_QS2_004_BC26_2</td>
      <td>D5: T1_QS2_004_BC26_3</td>
      <td>D5: T1_QS2_004_BC26_4</td>
      <td>D5: T1_QS2_004_BC26_5</td>
      <td>D5: T1_QS2_004_BC26_6</td>
      <td>D5: T1_QS2_004_BC26_7</td>
      <td></td>
      <td></td>
    </tr>
    <tr>
      <th>C</th>
      <td>Blank</td>
      <td>NTC</td>
      <td></td>
      <td>D5: T1_QS2_004_BC26_1</td>
      <td>D5: T1_QS2_004_BC26_2</td>
      <td>D5: T1_QS2_004_BC26_3</td>
      <td>D5: T1_QS2_004_BC26_4</td>
      <td>D5: T1_QS2_004_BC26_5</td>
      <td>D5: T1_QS2_004_BC26_6</td>
      <td>D5: T1_QS2_004_BC26_7</td>
      <td></td>
      <td></td>
    </tr>
    <tr>
      <th>D</th>
      <td></td>
      <td></td>
      <td></td>
      <td></td>
      <td></td>
      <td></td>
      <td></td>
      <td></td>
      <td></td>
      <td></td>
      <td></td>
      <td></td>
    </tr>
    <tr>
      <th>E</th>
      <td>Blank</td>
      <td>NTC</td>
      <td></td>
      <td>F9: T1_QS2_007_BC17_1</td>
      <td>F9: T1_QS2_007_BC17_2</td>
      <td>F9: T1_QS2_007_BC17_3</td>
      <td>F9: T1_QS2_007_BC17_4</td>
      <td>F9: T1_QS2_007_BC17_5</td>
      <td>F9: T1_QS2_007_BC17_6</td>
      <td>F9: T1_QS2_007_BC17_7</td>
      <td></td>
      <td></td>
    </tr>
    <tr>
      <th>F</th>
      <td>Blank</td>
      <td>NTC</td>
      <td></td>
      <td>F9: T1_QS2_007_BC17_1</td>
      <td>F9: T1_QS2_007_BC17_2</td>
      <td>F9: T1_QS2_007_BC17_3</td>
      <td>F9: T1_QS2_007_BC17_4</td>
      <td>F9: T1_QS2_007_BC17_5</td>
      <td>F9: T1_QS2_007_BC17_6</td>
      <td>F9: T1_QS2_007_BC17_7</td>
      <td></td>
      <td></td>
    </tr>
    <tr>
      <th>G</th>
      <td>Blank</td>
      <td>NTC</td>
      <td></td>
      <td>F9: T1_QS2_007_BC17_1</td>
      <td>F9: T1_QS2_007_BC17_2</td>
      <td>F9: T1_QS2_007_BC17_3</td>
      <td>F9: T1_QS2_007_BC17_4</td>
      <td>F9: T1_QS2_007_BC17_5</td>
      <td>F9: T1_QS2_007_BC17_6</td>
      <td>F9: T1_QS2_007_BC17_7</td>
      <td></td>
      <td></td>
    </tr>
    <tr>
      <th>H</th>
      <td></td>
      <td></td>
      <td></td>
      <td></td>
      <td></td>
      <td></td>
      <td></td>
      <td></td>
      <td></td>
      <td></td>
      <td></td>
      <td></td>
    </tr>
  </tbody>
</table>
