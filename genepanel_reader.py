import re
import pandas as pd


def excel_reader(file_name):
    """
    This function converts the table in the excel file to a table in
    Python.

    :param file_name: the name of the file

    :return: the gene panels table
    """
    return pd.read_excel(file_name)


def get_column(gp_table, column_name):
    """
    This function gets the values from a column from the gene panels
    table and returns it in a list.

    :param gp_table: the gene panels table
    :param column_name: the name of the column

    :return: the values from a column in a list
    """
    return gp_table[column_name].tolist()


def make_genepanel_list_set(gene_panels_list):
    """
    This function makes a list and set of all the gene panels with
    the (AD, etc.) removed.

    :param gene_panels_list: the values from a column in a list

    :return: the list and set of all the gene panels
    """
    gps_list = []
    gps_set = []
    for gene_panels in gene_panels_list:
        gp_list = []
        gp_list.clear()
        # Remove the (AD, etc.)
        for gene_panel in re.sub('[a-zA-Z]+;', ",", gene_panels).split(";"):
            gp = re.sub('( (.*))', "", gene_panel)
            gp_list.append(gp)
            if gp in gps_set:
                pass
            else:
                gps_set.append(gp)
        gps_list.append(gp_list[:])
    return gps_list, gps_set


def make_gene_dict(genes, gps_list):
    """
    This function makes a dictionary with the genes as keys and the
    gene panels in a list as values.

    :param genes: the list of all the genes
    :param gps_list: the list of all the gene panels

    :return: the dictionary with the genes as keys and the gene panels
    as values
    """
    return dict(zip(genes, gps_list))


def make_gene_panel_dict(gps_set, genes_dict):
    """
    This function makes a dictionary with the gene panels as keys
    and the genes as values.

    :param gps_set: the set of all the gene panels
    :param genes_dict: the dictionary with the genes as keys and the
    gene panels as values

    :return gene_panel_dict: the dictionary with the gene panels as keys
    and the genes as values
    """
    gene_panel_dict = {}
    for gene_panel in gps_set:
        gene_panel_dict[gene_panel] = [k for k, v in genes_dict.items() if
                                       gene_panel in v]
    return gene_panel_dict


def gene_synonyms(synonyms, gps_list, genes_dict, gene_panel_dict):
    """
    This function adds all synonyms from the excel file to the
    genes_dict and the gene_panel_dict.

    :param synonyms: the list of all the genes
    :param gps_list: the list of all the gene panels
    :param genes_dict: the dictionary with the genes as keys and the
    gene panels as values
    :param gene_panel_dict: the dictionary with the gene panels as keys
    and the genes as values

    :return genes_dict: the dictionary with the genes/gene synonyms as
    keys and the gene panels as values
    :return gene_panel_dict: the dictionary with the gene panels as keys
    and the genes and gene synonyms as values
    """
    # Loops through the synonyms and adds them and the corresponding
    # genepanels to the genes_dict and gene_panel_dict
    for row in range(len(synonyms)):
        synonyms[row] = str(synonyms[row]).upper().split("|")
        for gene in range(len(synonyms[row])):
            genes_dict[synonyms[row][gene]] = gps_list[row]
            for genpanel in gps_list[row]:
                gene_panel_dict[genpanel].append(synonyms[row][gene])
    return genes_dict, gene_panel_dict

