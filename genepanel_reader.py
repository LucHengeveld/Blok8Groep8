"""
Margo Raijmakers
28-04-2021
"""
import pandas as pd
import re


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


def make_genepanel_list_set(gene_panels_list):
    """This function makes a list and set of all the gene panels with
    the (AD, etc.) removed.

    :param gene_panels_list: the values from a column in a list
    :return: the list and set of all the gene panels
    """
    gps_list = []
    gps_set = []
    for gene_panels in gene_panels_list:
        gp_list = []
        gp_list.clear()
        for gene_panel in re.sub('[a-zA-Z]+;', ",", gene_panels).split(";"):
            gp = re.sub('( (.*))', "", gene_panel)
            gp_list.append(gp)
            if gp in gps_set:
                pass
            else:
                gps_set.append(gp)
        gps_list.append(gp_list[:])
    return gps_list, gps_set


def main():
    file_name = "GenPanelOverzicht_DG-3.1.0_HAN.xlsx"
    gp_table = excel_reader(file_name)
    genes = get_column(gp_table, "GenePanels_Symbol")
    gene_panels_list = get_column(gp_table, "GenePanel")
    gps_list, gps_set = make_genepanel_list_set(gene_panels_list)

    # gene_panel_dict = {}
    # for gene_panel in gps_set:
    #     print(gene_panel)
    #     for gene_panels in gps_list:
    #         if gene_panel in gene_panels:
    #             print(gene_panel)

    pos_count = 0
    for gene_panels in gps_list:
        pos_count += 1
        print(pos_count)
        for gene_panel in gps_set:
            print(gene_panel)
            if gene_panel in gene_panels:
                print(gene_panel)
    # todo key:value = genpanel:[gen, gen, gen]
    # todo posities in de gene_panels_list die genpanel bevatten -> zelfde positie van genes, die genen in een lijst
    # todo dict(zip(keys, values))

main()



