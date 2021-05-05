from Bio import Entrez, Medline
import requests
from datetime import datetime


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
    """This function uses the id's from the id_list to create a pubtator link
    which filters the genes, mutations and diseases out of the title and
    abstract.

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
    """This function reads the pubtator link as a text file and retrieves the
    genes, diseases and mutations out of each article.

    :param pubtator_link: Pubtator link with the title, abstact, genes,
    diseases and mutations of each article.
    :return results: Dictionary with as key the article ID and as value a list
    with the structure [title, abstract, genelist, diseaselist, mutationlist]
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
    """This function creates the hyperlink to the pubmed article by using the
    key's (article ID's) from the results dictionary. The hyperlink will be
    added to the values of each key.

    :param results: Dictionary with as key the article ID and as value a list
    with the structure [title, abstract, genelist, diseaselist, mutationlist]
    :return results: Dictionary with as key the article ID and as value a list
    with the structure [title, abstract, genelist, diseaselist, mutationlist,
    hyperlink]
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
    """This function retrieves the publication date of the article by using
    the key's of results (article ID's) in efetch. The publication date will be
    added to the values of each key.

    :param results: Dictionary with as key the article ID and as value a list
    with the structure [title, abstract, genelist, diseaselist, mutationlist,
    hyperlink]
    :return results: Dictionary with as key the article ID and as value a list
    with the structure [title, abstract, genelist, diseaselist, mutationlist,
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


def main():
    date_filter = "2020-01-01"
    query = "((ABC transporter [tiab] OR transporter [tiab] OR transport [" \
            "tiab]) AND (disease [tiab] OR mutation [tiab] OR mutations [" \
            "tiab] OR liver disease [tiab]) AND (lipids [tiab] OR " \
            "cholesterol [tiab] OR bile salts [tiab] OR canalicular membrane "\
            "[tiab] OR phosphatidylcholine [tiab] OR PC [tiab]) AND (ABCB4 [" \
            "tiab] OR ABCB4 deficiency [tiab])) "

    id_list = get_pubmed_ids(query, date_filter)
    pubtator_link = get_pubtator_link(id_list)
    results = read_pubtator_file(pubtator_link)
    results = pubmed_hyperlink(results)
    results = publication_date(results)


main()
