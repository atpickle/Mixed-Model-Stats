# Mixed-Model-Stats

Microglia Area_RMixedFormated.csv is example formating of the input data sheet.

Use Convert_Excel_to_Mixed_Models_CSV.ipynb to convert an unorganized excel sheet consisting of the naming convention described below in the first column, and the raw data in the second column.

Naming convention as separated by underscores:
Tri_1ugLPS_01_1_C_0.tiff
0 - Tri (Culture Type)
1 - 1ugLPS (Treatment)
2 - 01 (Well number)
3 - 1 (Image number)
4 - C (Channel)
5 - 0 (Channel # where 0 is microglia)

Only items 1, 2, 3 are grabbed by the code to organize into columns.
