import re


def get_values_for_relevance(or_list, and_filter, gene_filter):
    """
    This function retrieves all values out of the or_list, and_filter
    and gene_filter and puts them in a list.

    :param or_list: list of OR search terms
    :param and_filter: list of AND search terms
    :param gene_filter: list of gene filter search terms

    :return filters: list with all values from the or_list, and_filter
    and gene_filter
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
    """
    This function calculates the relevance score of every article
    with the searched query.

    :param results: dictionary with as key the article ID and as value
    a list with the structure [title, abstract, genelist, diseaselist,
    hyperlink, publication date, genepanels, diseases-co-occurence]
    :param filters: list with all values from the or_list, and_filter
    and gene_filter

    :return relevance_score: 2D list of every article with the relevance
    score
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
    """
    This function sorts the articles by relevance score.

    :param relevance_score: 2D list of every article with the relevance
    score

    :return: A sorted 2D list (high to low score) of every article with
    the relevance score
    """
    return sorted(relevance_score, key=lambda l: l[1], reverse=True)
