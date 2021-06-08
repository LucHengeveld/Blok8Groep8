import requests
from Bio import Entrez, Medline


def get_pubtator_link(id_list):
    """
    This function uses the id's from the id_list to create a
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
    """
    This function reads the pubtator link as a text file and
    retrieves the title, abstract, genes and diseases out of each
    article.

    :param gene_filter: given gene filter
    :param genepanel_filter: given gene panel filter
    :param gene_panel_dict: the dictionary with the gene panels as keys
    and the genes as values
    :param pubtator_link: Pubtator link with the title, abstact, genes
    and diseases of each article

    :return results: Dictionary with as key the article ID and as value
    a list with the structure [title, abstract, genelist, diseaselist]
    """
    # Retrieves the pubtator link with the article ID's in text format
    pubtator_text = requests.get(pubtator_link).text

    # Splits the text in lines.
    lines = pubtator_text.split("\n")

    # Creates an empty dictionary and 3 empty lists
    results = {}
    genelist = []
    diseaselist = []
    genepanel_filter_lijst = []

    # Check if the user entered a gene filter
    if gene_filter != "":
        genefilter_lijst = gene_filter.split(", ")
    else:
        genefilter_lijst = []

    # Loops through the PubTator text lines
    for i in range(len(lines)):

        # If "|t| is in the line, it saves it as the ID and title
        if "|t|" in lines[i]:
            article_id = lines[i].split("|t|")[0]
            title = lines[i].split("|t|")[1]

        # If "|a| is in the line, it saves it as the abstract
        elif "|a|" in lines[i]:
            abstract = lines[i].split("|a|")[1]

        # If the line isnt emtpy it checks if the line contains "Gene"
        elif lines[i] != "":
            if "Gene" in lines[i]:
                # Saves the gene ID and gene name to a variable
                gene = lines[i].split("\t")[3] + " " + \
                       lines[i].split("\t")[-1]
                genefilter_boolean = False

                # Continues if the gene name doesnt have a space or -
                # This filters out most false positives
                if " " not in lines[i].split("\t")[3] and "-" not in \
                        lines[i].split("\t")[3]:

                    # If a gene filter is entered by the user, it
                    # compares the PubTator genes to the gene filter
                    for j in range(len(genefilter_lijst)):
                        if genefilter_lijst and lines[i].split("\t")[3] \
                                == genefilter_lijst[j]:
                            genefilter_boolean = True

                    # Checks if the gene is already added to the list
                    if genefilter_boolean or not genefilter_lijst:
                        if gene.upper() not in genelist and gene.lower() not \
                                in genelist and gene not in genelist:

                            # Continues if genepanel filter isn't empty
                            if genepanel_filter != "":
                                genepanelboolean = False

                                # Checks if multiple genepanels are
                                # entered by the user
                                if "," in genepanel_filter.replace(" ", ""):
                                    genepanel_filter_lijst = genepanel_filter.\
                                        upper().replace(" ", "").split(",")
                                else:
                                    genepanel_filter_lijst.append(
                                        genepanel_filter.upper())

                                # Loops through the genepanel filters
                                # and leaves the corresponding genes out
                                # of the results
                                for j in genepanel_filter_lijst:
                                    if j in gene_panel_dict.keys():
                                        if lines[i].split("\t")[3].upper() in gene_panel_dict[j]:
                                            genepanelboolean = True
                                if not genepanelboolean:
                                    if gene not in genelist:
                                        genelist.append(gene.upper())
                            else:
                                # If no genepanel filter is entered, it
                                # adds all found genes to the gene list
                                genelist.append(gene.upper())

            # If the line isnt emtpy it checks if the line contains "Gene"
            elif "Disease" in lines[i]:
                disease = lines[i].split("\t")[3] + " " + \
                          lines[i].split("\t")[-1]
                if disease not in diseaselist:
                    diseaselist.append(disease)

        # if the line is empty, which means its the end of an article,
        # it will add the title, abstract, genelist and diseaselist to
        # the results dictionary.
        else:
            if genelist:
                results[article_id] = [title, abstract, genelist, diseaselist]
                genelist = []
                diseaselist = []

    return results


def pubtator_hyperlink(results):
    """
    This function creates the hyperlink to the PubTator article by
    using the key's (article ID's) from the results dictionary. The
    hyperlink will be added to the values of each key.

    :param results: dictionary with as key the article ID and as value
    a list with the structure [title, abstract, genelist, diseaselist]

    :return results: dictionary with as key the article ID and as value
    a list with the structure [title, abstract, genelist, diseaselist,
    hyperlink]
    """
    # Standard PubTator hyperlink
    standard_hyperlink = "https://www.ncbi.nlm.nih.gov/research/pubtator/?" \
                         "view=docsum&query=id"

    # Adds the article ID's to the standard pubmed hyperlink and appends
    # this to the results dictionary.
    for key in results:
        hyperlink = standard_hyperlink.replace("id", key)
        results[key].append(hyperlink)

    return results


def publication_date(results):
    """
    This function retrieves the publication date of the article by
    using the key's of results (article ID's) in efetch. The publication
    date will be added to the values of each key.

    :param results: dictionary with as key the article ID and as value
    a list with the structure [title, abstract, genelist, diseaselist,
    hyperlink]

    :return results: dictionary with as key the article ID and as value
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


def genepanel_results(results, genes_dict):
    """
    This function adds the genepanels to a gene in the results
    dictionary.

    :param results: dictionary with as key the article ID and as value
    a list with the structure [title, abstract, genelist, diseaselist,
    hyperlink, publication date]
    :param genes_dict: the dictionary with the genes/gene synonyms as
    keys and the gene panels as values

    :return results: dictionary with as key the article ID and as value
    a list with the structure [title, abstract, genelist, diseaselist,
    hyperlink, publication date, genepanels]
    """
    # Loops through the genes and adds the genepanels to them
    for key in results:
        genepanel = []
        for gene in results[key][2]:
            if gene.rsplit(" ", 1)[0] in genes_dict.keys():
                genepanel.append(genes_dict[gene.rsplit(" ", 1)[0]])
            else:
                genepanel.append("-")
        results[key].append(genepanel)

    return results
