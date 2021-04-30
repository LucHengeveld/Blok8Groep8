"""
Margo Raijmakers
28-04-2021
"""
import pandas as pd


def excel_reader(file_name):
    """This function converts the table in the excel file to a table in
    Python.

    :param file_name: the name of the file
    :return: the gene panels table
    """
    return pd.read_excel(file_name)


def get_column(gp_table, column_name):
    """This function gets the values from a column from the gene panels
    table and returns it in a list.

    :param gp_table: the gene panels table
    :param column_name: the name of the column
    :return: the values from a column in a list
    """
    return gp_table[column_name].tolist()


def main():
    file_name = "GenPanelOverzicht_DG-3.1.0_HAN.xlsx"
    gp_table = excel_reader(file_name)
    genes = get_column(gp_table, 'GenePanels_Symbol')
    gene_panels = get_column(gp_table, 'GenePanel')
    print(len(genes))
    print(len(gene_panels))

main()
