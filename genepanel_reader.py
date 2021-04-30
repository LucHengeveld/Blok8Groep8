"""
Margo Raijmakers
28-04-2021
"""
import pandas as pd


def excel_reader(file_name):
    """This function converts the table in the excel file to a table in
    Python.

    :param file_name: the name of the file
    :return: the table
    """
    return pd.read_excel(file_name)


def main():
    file_name = "GenPanelOverzicht_DG-3.1.0_HAN.xlsx"
    gp_table = excel_reader(file_name)
    print(gp_table['GenePanels_Symbol'].tolist())


main()
