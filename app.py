import re
from flask import Flask, render_template, request, send_file
import pandas as pd
from Bio import Entrez, Medline
import requests
from datetime import datetime
import csv
from xlsxwriter.workbook import Workbook
import ast

app = Flask(__name__)


@app.route('/')
@app.route('/home.html', methods=["POST", "GET"])
def get_input():
    if request.method == 'POST':
        try:
            or_filter = request.form["or_filter"]
            or_list = request.form.getlist('or_list')
            or_list.insert(0, or_filter)

            and_filter = request.form["and_filter"]
            not_filter = request.form["not_filter"]
            gene_filter = request.form["gene_filter"]
            date_filter = request.form["date_filter"]

            try:
                genepanel_file = request.files["genepanel_file"]
                genepanel_file_name = genepanel_file.filename
                genepanel_file.save(genepanel_file_name)
            except:
                genepanel_file = ""
                genepanel_file_name = "No file selected."

            genepanel_filter = request.form["genepanel_filter"]

            email = request.form["email"]

            try:
                use_co_occurence = request.form["occurence"]
            except:
                use_co_occurence = "Not selected"

            print(email)
            print(or_filter)
            print(or_list)
            print(and_filter)
            print(not_filter)
            print(genepanel_file_name)
            print(gene_filter)
            print(date_filter)
            print(genepanel_filter)
            print(use_co_occurence)

            Entrez.email = email

            # todo module genpanel reader
            if genepanel_file:
                gp_table = excel_reader(genepanel_file)
                genes = get_column(gp_table, "GenePanels_Symbol")
                synonyms = get_column(gp_table, "Aliases")
                gene_panels_list = get_column(gp_table, "GenePanel")
                gps_list, gps_set = make_genepanel_list_set(
                    gene_panels_list)
                genes_dict = make_gene_dict(genes, gps_list)
                gene_panel_dict = make_gene_panel_dict(gps_set, genes_dict)
                genes_dict, gene_panel_dict = gene_synonyms(synonyms, gps_list,
                                                            genes_dict,
                                                            gene_panel_dict)
            else:
                gene_panel_dict = {}
                genes_dict = {}

            # todo module pubmed/pubtator
            or_list2, and_filter2, not_filter2, gene_filter2 = retrieve_data(
                or_list, and_filter, not_filter, gene_filter)

            query = making_query(or_list2, and_filter2, not_filter2,
                                 gene_filter2)
            # query = "((ABC transporter [tiab] OR transporter [tiab] OR transport [" \
            #         "tiab]) AND (disease [tiab] OR mutation [tiab] OR mutations [" \
            #         "tiab] OR liver disease [tiab]) AND (lipids [tiab] OR " \
            #         "cholesterol [tiab] OR bile salts [tiab] OR canalicular membrane " \
            #         "[tiab] OR phosphatidylcholine [tiab] OR PC [tiab]) AND (ABCB4 [" \
            #         "tiab] OR ABCB4 deficiency [tiab])) "

            id_list = get_pubmed_ids(query, date_filter)
            pubtator_link = get_pubtator_link(id_list)

            results = read_pubtator_file(pubtator_link, gene_panel_dict,
                                         genepanel_filter, gene_filter)
            results = pubmed_hyperlink(results)
            results = publication_date(results)
            results = genepanel_results(results, genes_dict)

            titlepoints = 10
            sentencepoints = 5
            abstractpoints = 3
            articlepoints = 1
            diseasepoints = co_occurrence(results, articlepoints,
                                          abstractpoints, sentencepoints,
                                          titlepoints, 3)
            results = add_co_occurrence_to_results(results, diseasepoints)

            # todo module relevance score
            filters = get_values_for_relevance(or_list, and_filter,
                                               gene_filter)
            relevance_score = get_relevance_score(results, filters)
            relevance_score = sort_relevance_score(relevance_score)
            # print(relevance_score)

            return render_template("results.html",
                                   or_list=or_list,
                                   and_filter=and_filter,
                                   not_filter=not_filter,
                                   gene_filter=gene_filter,
                                   date_filter=date_filter,
                                   genepanel_file_name=genepanel_file_name,
                                   genepanel_filter=genepanel_filter,
                                   email=email,
                                   use_co_occurence=use_co_occurence,
                                   results=results,
                                   relevance_score=relevance_score)
        except:
            or_filter = request.form["or_filter"]
            or_list = request.form.getlist('or_list')
            or_list.insert(0, or_filter)

            and_filter = request.form["and_filter"]
            not_filter = request.form["not_filter"]
            gene_filter = request.form["gene_filter"]
            date_filter = request.form["date_filter"]

            try:
                genepanel_file = request.files["genepanel_file"]
                genepanel_file_name = genepanel_file.filename
                genepanel_file.save(genepanel_file_name)
            except:
                genepanel_file = ""
                genepanel_file_name = "No file selected."

            genepanel_filter = request.form["genepanel_filter"]

            email = request.form["email"]

            try:
                use_co_occurence = request.form["occurence"]
            except:
                use_co_occurence = "Not selected"

            return render_template("home_error.html",
                                   or_list=or_list,
                                   and_filter=and_filter,
                                   not_filter=not_filter,
                                   gene_filter=gene_filter,
                                   date_filter=date_filter,
                                   genepanel_file_name=genepanel_file_name,
                                   genepanel_filter=genepanel_filter,
                                   email=email,
                                   use_co_occurence=use_co_occurence)

    else:
        return render_template("home.html",
                               or_list="",
                               and_filter="",
                               not_filter="",
                               gene_filter="",
                               date_filter="",
                               genepanel_file_name="",
                               genepanel_filter="",
                               email="",
                               use_co_occurence="",
                               results="")


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


def gene_synonyms(synonyms, gps_list, genes_dict, gene_panel_dict):
    for row in range(len(synonyms)):
        synonyms[row] = str(synonyms[row]).split("|")
        for gene in range(len(synonyms[row])):
            genes_dict[synonyms[row][gene]] = gps_list[row]
            for genpanel in gps_list[row]:
                gene_panel_dict[genpanel].append(synonyms[row][gene])

    return genes_dict, gene_panel_dict


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
            or_list2 = str(or_list)
            if "['" in or_list2:  # [] gets replaced by ()
                list_or2 = or_list2.replace("['", "(").replace("']",
                                                               " [tiab])")
                # print("List: ", list_or)
            if "" not in list_or2:  # Checks if list is empty
                return
            or_list2 = list_or2.replace(",", " [tiab] OR")
            if ' OR ' in or_list2:  # some OR's gets replaced by AND
                or_list2 = or_list2.replace("' [tiab] OR '", " [tiab]) "
                                                             "AND (")
            # print("OR search: ", or_list)
        else:
            or_list2 = str(or_list)

        # The and_filter list will be edited here
        if and_filter is not None:
            and_filter2 = str(and_filter)
            # Commas get replaced by OR
            and_filter2 = and_filter2.replace(",", " [tiab] OR")
            # print("AND search: ", and_filter)
        else:
            and_filter2 = str(and_filter)

        # The not_filter list will be edited here
        if not_filter is not None:
            not_filter2 = str(not_filter)
            # Commas get replaced by NOT
            not_filter2 = not_filter.replace(",", " [tiab] NOT")
            # print("NOT search: ", not_filter)
        else:
            not_filter2 = str(not_filter)

        # # The gene_filter list will be edited here
        if gene_filter is not None:
            gene_filter2 = str(gene_filter)
            # Commas get replaced by OR
            gene_filter2 = gene_filter.replace(",", " [tiab] OR")
            # print("Gene filter search: ", gene_filter)
        else:
            gene_filter2 = str(gene_filter)

        return or_list2, and_filter2, not_filter2, gene_filter2

    except ValueError:
        print("Error: something went wrong. Please check the info "
              "page.")


def making_query(or_list2, and_filter2, not_filter2, gene_filter2):
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
        if or_list2 != "":  # Query of the or_list added to empty
            # query list
            query_or = or_list2
            query.append(query_or)
            # print("Query or: ", query_or)
        else:
            pass
        if and_filter2 != "":  # Query of the and_filter added to empty
            # query list
            query_and = " AND (", and_filter2
            query_and = str(query_and).replace(" ', '", " ").replace(
                "'", "").replace(", ", "").replace(")", " [tiab])")
            query.append(query_and)
            # print("Query and: ", query_and)
        else:
            pass
        if not_filter2 != "":  # Query of the not_filter added to empty
            # query list
            query_not = " AND (NOT", not_filter2
            query_not = str(query_not).replace(" ', '", "").replace(
                "'", "").replace(",", "").replace(")", " [tiab])")
            query.append(query_not)
            # print("Query not: ", query_not)
        else:
            pass
        if gene_filter2 != "":  # Query of the gene_filter added to
            # empty query list
            query_gene = " AND ", gene_filter2
            query_gene = str(query_gene).replace("'", "").replace(
                ", """, "(").replace(")", " [tiab])")
            query.append(query_gene)
            # print("Query gene: ", query_gene)
        else:
            pass
        query = str(query).replace("', '( ", " ").replace("['", "(") \
            .replace("']", ")")
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
    pubtator link which filters the genes and diseases out
    of the title and abstract.

    :param id_list: List with all found pubmed id's.
    :return pubtator_link: Pubtator link with the title, abstact, genes
    and diseases of each article.
    """
    # Standard pubtator link format:
    link = "https://www.ncbi.nlm.nih.gov/research/pubtator-api/publications" \
           "/export/pubtator?pmids=idvalues&concepts=gene,disease"

    # Creates a string with all the ID's. Separated by a comma.
    id_string = ""
    for i in range(len(id_list)):
        id_string += id_list[i] + ","
    id_string = id_string[:-1]

    # Adds the ID values to the standard link.
    pubtator_link = link.replace("idvalues", id_string)

    return pubtator_link


def read_pubtator_file(pubtator_link, gene_panel_dict, genepanel_filter,
                       gene_filter):
    """This function reads the pubtator link as a text file and
    retrieves the title, abstract, genes and diseases out of each article.

    :param pubtator_link: Pubtator link with the title, abstact, genes
    and diseases of each article.
    :return results: Dictionary with as key the article ID and as value
    a list with the structure [title, abstract, genelist, diseaselist]
    """
    # Retrieves the pubtator link with the article ID's in text format
    pubtator_text = requests.get(pubtator_link).text

    # Splits the text in lines.
    lines = pubtator_text.split("\n")

    # Checks if the line is the title, the abstract, a gene, a disease
    # and adds the information to lists.
    results = {}
    genelist = []
    diseaselist = []
    genepanel_filter_lijst = []
    if gene_filter != "":
        genefilter_lijst = gene_filter.split(", ")
    else:
        genefilter_lijst = []
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
                    genefilter_boolean = False
                    if not " " in lines[i].split("\t")[3]:
                        for j in range(len(genefilter_lijst)):
                            if genefilter_lijst and lines[i].split("\t")[3] == \
                                    genefilter_lijst[j]:
                                genefilter_boolean = True
                        if genefilter_boolean == True or not genefilter_lijst:
                            if gene.upper() not in genelist:
                                if gene.lower() not in genelist:
                                    if gene not in genelist:
                                        if genepanel_filter != "":
                                            genepanelboolean = False
                                            if "," in genepanel_filter.replace(" ", ""):
                                                genepanel_filter_lijst = genepanel_filter.upper().replace(" ", "").split(",")
                                            else:
                                                genepanel_filter_lijst.append(
                                                    genepanel_filter.upper())
                                            for j in genepanel_filter_lijst:
                                                if j in gene_panel_dict.keys():
                                                    if lines[i].split("\t")[
                                                        3].upper() \
                                                            in gene_panel_dict[j]:
                                                        genepanelboolean = True
                                            if genepanelboolean == False:
                                                if gene not in genelist:
                                                    genelist.append(gene)
                                        else:
                                            genelist.append(gene)

                elif "Disease" in lines[i]:
                    disease = lines[i].split("\t")[3] + " " + \
                              lines[i].split("\t")[-1]
                    if disease not in diseaselist:
                        diseaselist.append(disease)

        # if the line is empty, which means its the end of an article,
        # it will add the title, abstract, genelist and diseaselist to the results dictionary.
        else:
            if genelist:
                results[article_id] = [title, abstract, genelist, diseaselist]
                genelist = []
                diseaselist = []

    return results


def pubmed_hyperlink(results):
    """This function creates the hyperlink to the pubmed article by
    using the key's (article ID's) from the results dictionary. The
    hyperlink will be added to the values of each key.

    :param results: Dictionary with as key the article ID and as value
    a list with the structure [title, abstract, genelist, diseaselist]
    :return results: Dictionary with as key the article ID and as value
    a list with the structure [title, abstract, genelist, diseaselist,
    hyperlink]
    """
    # Standard pubmed hyperlink
    standard_hyperlink = "https://www.ncbi.nlm.nih.gov/research/pubtator/?view=docsum&query=id"

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
    hyperlink]
    :return results: Dictionary with as key the article ID and as value
    a list with the structure [title, abstract, genelist, diseaselist,
    hyperlink, publication date]
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


def co_occurrence(results, articlepoints, abstractpoints, sentencepoints,
                  titlepoints, pos):
    """This function gives points to every combination of genes and
    diseases in the PubMed articles (title and/or abstract).
    :param results: Dictionary with as key the article ID and as value
    a list with the structure [title, abstract, genelist, diseaselist,
     hyperlink]
    :param articlepoints: The amount of points for a gene and
    disease in the same article
    :param abstractpoints: The amount of points for a gene and
    disease in the same abstract
    :param sentencepoints: The amount of points for a gene and
    disease in the same sentence
    :param titlepoints: The amount of points for a gene and
    disease in the same title
    :param pos: The position in the results dictionary. 3 = disease
    :return: The list with all the points per gene and disease
    combination per article
    """
    points = {}
    for key in results:
        pointsperid = {}
        for gene in results[key][2]:
            gene = gene.rsplit(" ", 1)[0]
            pointspergene = []
            valuespergene = []
            for value in results[key][pos]:
                value = value.rsplit(" ", 1)[0]
                valuespergene.append(value)
                count = 0

                # If gene and disease are in the same article,
                # but not in the title or abstract
                if gene in results[key][0] and value in results[key][1] \
                        or gene in results[key][1] and value in \
                        results[key][0]:
                    count += articlepoints

                # If gene and disease are both in the abstract
                if gene in results[key][1] and value in results[key][1]:
                    count += abstractpoints

                # If gene and disease are both in the same line
                # of the abstract
                for line in results[key][1].split(". "):
                    if gene in line and value in line:
                        count += sentencepoints

                # If gene and disease are both in the title
                if gene in results[key][0] and value in results[key][0]:
                    count += titlepoints

                pointspergene.append(count)
            pointsperid[gene] = zip(pointspergene, valuespergene)
        points[key] = pointsperid
    return points


def add_co_occurrence_to_results(results, diseasepoints):
    """This function adds the 3 diseases with the highest co-occurrence
    score to the results.

    :param results: Dictionary with as key the article ID and as value
    a list with the structure [title, abstract, genelist, diseaselist,
    hyperlink]
    :param diseasepoints: The list with all the points per gene and
    disease combination per article
    :return: the results dictionary
    """
    for key, value in results.items():
        diseasesperarticle = []
        for key2, value2 in diseasepoints.get(key).items():
            diseasespergene = []
            for disease in sorted(value2, reverse=True)[:3]:
                diseasespergene.append(disease[1])
            diseasesperarticle.append(diseasespergene)
        results[key].append(diseasesperarticle)
    return results


def get_values_for_relevance(or_list, and_filter, gene_filter):
    """This function gets all values to calculate the relevance score
    with.

    :param or_list: a list of OR search terms
    :param and_filter: a list of AND search terms
    :param gene_filter: a list of gene filter search terms
    :return: a list of filter terms
    """
    filters = []
    for or_filter in or_list:
        for or_fil in or_filter.split(", "):
            if or_fil == "":
                pass
            else:
                filters.append(or_fil)
    for and_fil in and_filter.split(", "):
        if and_fil == "":
            pass
        else:
            filters.append(and_fil)
    for gene_fil in gene_filter.split(", "):
        if gene_fil == "":
            pass
        else:
            filters.append(gene_fil)
    return filters


def get_relevance_score(results, filters):
    """This function calculates the relevance score of every article
    with the searched query.

    :param results: Dictionary with as key the article ID and as value
    a list with the structure [title, abstract, genelist, diseaselist,
    hyperlink]
    :param filters: a list of filter terms
    :return: a 2D list of every article with the relevance score in a
    list
    """
    relevance_score = []
    for key, value in results.items():
        rel_score = 0
        for filter in filters:
            # Add +2 for every filter in the title
            rel_score += 2 * len(re.findall(filter.lower(), value[0].lower()))
            # Add +1 for every filter in the abstract
            rel_score += len(re.findall(filter.lower(), value[1].lower()))
        relevance_score.append([key, rel_score])
    return relevance_score


def sort_relevance_score(relevance_score):
    """This function sorts the articles by relevance score.

    :param relevance_score: a 2D list of every article with the
    relevance score in a list
    :return: a sorted 2D list of every article with the relevance score
    in a list
    """
    return sorted(relevance_score, key=lambda l: l[1], reverse=True)


def genepanel_results(results, genes_dict):
    for key in results:
        genepanel = []
        for gene in results[key][2]:
            if gene.rsplit(" ", 1)[0] in genes_dict.keys():
                genepanel.append(genes_dict[gene.rsplit(" ", 1)[0]])
            else:
                genepanel.append("-")
        results[key].append(genepanel)

    return results


@app.route('/results.html', methods=["POST"])
def save_results():
    results = request.form['results']
    results = ast.literal_eval(results)
    try:
        selected_extension = request.form["file_extension"]
    except:
        selected_extension = "tsv"

    if selected_extension == "txt":
        output_file = "results.txt"
    else:
        output_file = "results.tsv"

    with open(output_file, 'w', newline='') as out_file:
        tsv_writer = csv.writer(out_file, delimiter='\t')
        tsv_writer.writerow(
            ["Gene name", "Gene ID", "Gene Panels", "Pubmed ID",
             "Pubmed Hyperlink", "Publication Date", "Possible diseases"])
        for pubmed_id in results:
            for gene_index in range(len(results[pubmed_id][2])):
                gene_name = results[pubmed_id][2][gene_index].rsplit(" ", 1)[0]
                gene_id = results[pubmed_id][2][gene_index].rsplit(" ", 1)[1]
                genepanelstring = ""
                genepanels = results[pubmed_id][6][gene_index]
                for i in genepanels:
                    genepanelstring += i + ";"
                genepanelstring = genepanelstring[:-1]
                hyperlink = results[pubmed_id][4]
                date = results[pubmed_id][5]
                try:
                    diseases = results[pubmed_id][-1][gene_index]
                except IndexError:
                    diseases = ["-"]
                diseasestring = ""
                for j in diseases:
                    if j != diseases[-1]:
                        diseasestring += j + ","
                    else:
                        diseasestring += j

                tsv_writer.writerow(
                    [gene_name, gene_id, genepanelstring, pubmed_id, hyperlink,
                     date, diseasestring])

    out_file.close()

    if selected_extension == "xlsx":
        xlsx_file = 'results.xlsx'
        workbook = Workbook(xlsx_file, {'strings_to_numbers': True})
        worksheet = workbook.add_worksheet()
        tsv_reader = csv.reader(open(output_file, 'rt'), delimiter='\t')
        for row, data in enumerate(tsv_reader):
            worksheet.write_row(row, 0, data)
        workbook.close()

    if selected_extension == "xlsx":
        output = "results.xlsx"
    elif selected_extension == "txt":
        output = "results.txt"
    else:
        output = "results.tsv"

    if selected_extension == "xlsx":
        return send_file(attachment_filename=output, filename_or_fp=output,
                         mimetype="xlsx", as_attachment=True)
    else:
        return send_file(attachment_filename=output, filename_or_fp=output,
                         mimetype="text/tsv", as_attachment=True)


@app.route('/info.html', methods=["POST", "GET"])
def info():
    return render_template('info.html')


if __name__ == '__main__':
    app.run()
