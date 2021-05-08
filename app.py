import re
from flask import Flask, render_template, request
import pandas as pd
from Bio import Entrez, Medline
import requests
from datetime import datetime

app = Flask(__name__)


@app.route('/')
@app.route('/home.html', methods=["POST", "GET"])
# 1. Filters ophalen
# 2. genpanels
# 3. retrieve data / making query
# 4. test.py
# 5. resultaten weergeven
def get_input():
    if request.method == 'POST':

        email = request.form.get("email", "")

        or_filter = request.form.get("or_filter", "")
        or_list = request.form.getlist('or_list', "")
        or_list.insert(0, or_filter)

        and_filter = request.form.get("and_filter", "")
        not_filter = request.form.get("not_filter", "")
        genepanel_file = request.form.get("genepanel_file", "")
        gene_filter = request.form.get("gene_filter", "")
        date_filter = request.form.get("date_filter", "")
        genepanel_filter = request.form.get("genepanel_filter", "")

        print(or_list)
        print(and_filter)
        print(not_filter)
        print(genepanel_file)
        print(gene_filter)
        print(date_filter)
        print(genepanel_filter)
        print(email)

        Entrez.email = email

        gp_table = excel_reader(genepanel_file)
        genes = get_column(gp_table, "GenePanels_Symbol")
        gene_panels_list = get_column(gp_table, "GenePanel")
        gps_list, gps_set = make_genepanel_list_set(gene_panels_list)
        genes_dict = make_gene_dict(genes, gps_list)
        gene_panel_dict = make_gene_panel_dict(gps_set, genes_dict)
        or_list, and_filter, not_filter, gene_filter = \
            retrieve_data(or_list, and_filter, not_filter, gene_filter)
        query = making_query(or_list, and_filter, not_filter, gene_filter)
        id_list = get_pubmed_ids(query, date_filter)
        pubtator_link = get_pubtator_link(id_list)
        results = read_pubtator_file(pubtator_link)
        results = pubmed_hyperlink(results)
        results = publication_date(results)

        return render_template("home.html",
                               email=email,
                               or_filter=or_filter,
                               or_list=or_list,
                               and_filter=and_filter,
                               not_filter=not_filter,
                               gene_filter=gene_filter,
                               date_filter=date_filter,
                               genepanel_file=genepanel_file,
                               genepanel_filter=genepanel_filter)
    else:
        return render_template("home.html",
                               email="",
                               or_filter="",
                               or_list="",
                               and_filter="",
                               not_filter="",
                               gene_filter="",
                               date_filter="",
                               genepanel_file="",
                               genepanel_filter="")


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
    """This function makes a dictionary with the genes as keys and the
    gene panels in a list as values.

    :param genes: the list of all the genes
    :param gps_list: the list of all the gene panels
    :return: the dictionary with the genes as keys and the gene panels
    as values
    """
    return dict(zip(genes, gps_list))


def make_gene_panel_dict(gps_set, genes_dict):
    """This function makes a dictionary with the gene panels as keys
    and the genes as values.

    :param gps_set: the set of all the gene panels
    :param genes_dict: the dictionary with the genes as keys and the
    gene panels as values
    :return: the dictionary with the gene panels as keys and the genes
    as values
    """
    gene_panel_dict = {}
    for gene_panel in gps_set:
        gene_panel_dict[gene_panel] = [k for k, v in genes_dict.items() if
                                       gene_panel in v]
    return gene_panel_dict


def retrieve_data(or_list, and_filter, not_filter, gene_filter):
    """ This function will prepare the different parts of the query.

    :param or_list: The input from get_input.
    :param and_filter: The input from get_input.
    :param not_filter: The input from get_input.
    :param gene_filter: The input from get_input.
    :return or_list: List with OR search terms.
    :return and_fitler: List with AND search terms.
    :return not_filter: List with NOT search terms.
    :return gene_filter: List with gene filter search terms.
    """

    try:
        # The or_list will be edited here
        if or_list is not None:
            or_list = str(or_list)
            if "['" in or_list:  # [] gets replaced by ()
                list_or = or_list.replace("['", "(").replace("']",
                                                             " [tiab])")
                # print("List: ", list_or)
            if "" not in list_or:  # Checks if list is empty
                return
            or_list = list_or.replace(",", " [tiab] OR")
            if ' OR ' in or_list:  # some OR's gets replaced by AND
                or_list = or_list.replace("' [tiab] OR '", " [tiab]) "
                                                           "AND (")
            # print("OR search: ", or_list)
        else:
            or_list = str(or_list)

        # The and_filter list will be edited here
        if and_filter is not None:
            and_filter = str(and_filter)
            # Commas get replaced by OR
            and_filter = and_filter.replace(",", " [tiab] OR")
            # print("AND search: ", and_filter)
        else:
            and_filter = str(and_filter)

        # The not_filter list will be edited here
        if not_filter is not None:
            not_filter = str(not_filter)
            # Commas get replaced by NOT
            not_filter = not_filter.replace(",", " [tiab] NOT")
            # print("NOT search: ", not_filter)
        else:
            not_filter = str(not_filter)

        # # The gene_filter list will be edited here
        if gene_filter is not None:
            gene_filter = str(gene_filter)
            # Commas get replaced by OR
            gene_filter = gene_filter.replace(",", " [tiab] OR")
            # print("Gene filter search: ", gene_filter)
        else:
            gene_filter = str(gene_filter)

        return or_list, and_filter, not_filter, gene_filter

    except ValueError:
        print("Error: something went wrong. Please check the info "
              "page.")


def making_query(or_list, and_filter, not_filter, gene_filter):
    """ This function combine's the 4 filters into a query. This
    query can be used for searching the PubMed.

    :param or_list: List with OR search terms.
    :param and_filter: List with AND search terms.
    :param not_filter: List with NOT search terms.
    :param gene_filter: List with gene filter search terms.
    :return query: Combining the or_list, and_filter, not_filter and
    gene_filter.
    """

    try:
        query = []  # Creating empty list

        if or_list != "":  # Query of the or_list added to empty
            # query list
            query_or = or_list
            query.append(query_or)
            # print("Query or: ", query_or)
        else:
            pass

        if and_filter != "":  # Query of the and_filter added to empty
            # query list
            query_and = " AND (", and_filter
            query_and = str(query_and).replace(" ', '", " ").replace(
                "'", "").replace(", ", "").replace(")", " [tiab])")
            query.append(query_and)
            # print("Query and: ", query_and)
        else:
            pass

        if not_filter != "":  # Query of the not_filter added to empty
            # query list
            query_not = " AND (NOT", not_filter
            query_not = str(query_not).replace(" ', '", "").replace(
                "'", "").replace(",", "").replace(")", " [tiab])")
            query.append(query_not)
            # print("Query not: ", query_not)
        else:
            pass

        if gene_filter != "":  # Query of the gene_filter added to
            # empty query list
            query_gene = " AND ", gene_filter
            query_gene = str(query_gene).replace("'", "").replace(", "
                        "", "(").replace(")", " [tiab])")
            query.append(query_gene)
            # print("Query gene: ", query_gene)
        else:
            pass

        query = str(query).replace("', '( ", " ").replace("['", "(") \
            .replace("']", ")")
        print("Query: ", query)
        return query

    except ValueError:
        print("Error: something went wrong. Please check the info "
              "page.")


def get_pubmed_ids(query, date_filter):
    """This function enters the query into an esearch and retrieves the
    pubmed id's from the found articles.

    :param query: Pubmed search query.
    :param date_filter: Date input from the webapplication.
    :return id_list: List with all found pubmed id's.
    """
    # Retrieves the input date from the webapplication.
    date = date_filter.replace("-", "/")

    # Retrieves today's date.
    datetoday = str(datetime.date(datetime.now())).replace("-", "/")

    # Enters the query, input date and today's date in an esearch.
    handle = Entrez.esearch(db="pubmed", term=query, mindate=date,
                            maxdate=datetoday)

    # Reads the esearch and saves the ID's in id_list
    record = Entrez.read(handle)
    handle.close()
    id_list = (record["IdList"])

    return id_list


def get_pubtator_link(id_list):
    """This function uses the id's from the id_list to create a
    pubtator link which filters the genes, mutations and diseases out
    of the title and abstract.

    :param id_list: List with all found pubmed id's.
    :return pubtator_link: Pubtator link with the title, abstact, genes,
    diseases and mutations of each article.
    """
    # Standard pubtator link format:
    link = "https://www.ncbi.nlm.nih.gov/research/pubtator-api/publications" \
           "/export/pubtator?pmids=idvalues&concepts=gene,mutation,disease"

    # Creates a string with all the ID's. Separated by a comma.
    id_string = ""
    for i in range(len(id_list)):
        id_string += id_list[i] + ","
    id_string = id_string[:-1]

    # Adds the ID values to the standard link.
    pubtator_link = link.replace("idvalues", id_string)

    return pubtator_link


def read_pubtator_file(pubtator_link):
    """This function reads the pubtator link as a text file and
    retrieves the genes, diseases and mutations out of each article.

    :param pubtator_link: Pubtator link with the title, abstact, genes,
    diseases and mutations of each article.
    :return results: Dictionary with as key the article ID and as value
    a list with the structure [title, abstract, genelist, diseaselist,
    mutationlist]
    """
    # Retrieves the pubtator link with the article ID's in text format
    pubtator_text = requests.get(pubtator_link).text

    # Splits the text in lines.
    lines = pubtator_text.split("\n")

    # Checks if the line is the title, the abstract, a gene, a disease,
    # a mutation and adds the information to lists.
    results = {}
    genelist = []
    diseaselist = []
    mutationlist = []
    for i in range(len(lines)):
        if "|t|" in lines[i]:
            article_id = lines[i].split("|t|")[0]
            title = lines[i].split("|t|")[1]

        elif "|a|" in lines[i]:
            abstract = lines[i].split("|a|")[1]

        elif lines[i] != "":
            if lines[i] != "":
                if "Gene" in lines[i]:
                    gene = lines[i].split("\t")[3] + " " + \
                           lines[i].split("\t")[-1]
                    if gene not in genelist:
                        genelist.append(gene)
                elif "Disease" in lines[i]:
                    disease = lines[i].split("\t")[3] + " " + \
                              lines[i].split("\t")[-1]
                    if disease not in diseaselist:
                        diseaselist.append(disease)
                elif "Mutation" in lines[i]:
                    mutation = lines[i].split("\t")[3] + " " + \
                               lines[i].split("\t")[-1]
                    if mutation not in mutationlist:
                        mutationlist.append(mutation)

        # if the line is empty, which means its the end of an article,
        # it will add the title, abstract, genelist, diseaselist and
        # mutationlist to the results dictionary.
        else:
            if genelist:
                results[article_id] = [title, abstract, genelist, diseaselist,
                                       mutationlist]
                genelist = []
                diseaselist = []
                mutationlist = []

    return results


def pubmed_hyperlink(results):
    """This function creates the hyperlink to the pubmed article by
    using the key's (article ID's) from the results dictionary. The
    hyperlink will be added to the values of each key.

    :param results: Dictionary with as key the article ID and as value
    a list with the structure [title, abstract, genelist, diseaselist,
    mutationlist]
    :return results: Dictionary with as key the article ID and as value
    a list with the structure [title, abstract, genelist, diseaselist,
    mutationlist, hyperlink]
    """
    # Standard pubmed hyperlink
    standard_hyperlink = "https://pubmed.ncbi.nlm.nih.gov/id/"

    # Adds the article ID's to the standard pubmed hyperlink and appends
    # this to the results dictionary.
    for key in results:
        hyperlink = standard_hyperlink.replace("id", key)
        results[key].append(hyperlink)

    return results


def publication_date(results):
    """This function retrieves the publication date of the article by
    using the key's of results (article ID's) in efetch. The publication
    date will be added to the values of each key.

    :param results: Dictionary with as key the article ID and as value
    a list with the structure [title, abstract, genelist, diseaselist,
    mutationlist, hyperlink]
    :return results: Dictionary with as key the article ID and as value
    a list with the structure [title, abstract, genelist, diseaselist,
    mutationlist, hyperlink, publication date]
    """
    # Adds all ID's to a string, separated by a comma.
    id_string = ""
    for key in results:
        id_string += key + ","
    id_string = id_string[:-1]

    # Enters the id_string in efetch to retrieve the publication date
    # of each article.
    handle = Entrez.efetch(db="pubmed", id=id_string, retmode="text",
                           rettype="medline")
    records = Medline.parse(handle)
    records = list(records)

    # Adds the publication date of every article to the results
    # dictionary
    for record in records:
        results[record["PMID"]].append(record["DP"])

    return results


@app.route('/info.html', methods=["POST", "GET"])
def info():
    return render_template('info.html')


if __name__ == '__main__':
    app.run()
