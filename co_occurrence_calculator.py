"""
This module calculates the co-occurrence score.
"""


def co_occurrence(results, pos):
    """
    This function gives points to every combination of genes and
    diseases in the PubMed articles (title and/or abstract).

    :param results: dictionary with as key the article ID and as value
    a list with the structure [title, abstract, genelist, diseaselist,
    hyperlink, publication date, genepanels]
    :param pos: the position in the results dictionary. 3 = disease

    :return points: the list with all the points per gene and disease
    combination per article
    """
    # The amount of points for a gene and disease in the same title
    titlepoints = 10

    # The amount of points for a gene and disease in the same sentence
    sentencepoints = 5

    # The amount of points for a gene and disease in the same abstract
    abstractpoints = 3

    # The amount of points for a gene and disease in the same article
    articlepoints = 1

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
    """
    This function adds the 3 diseases with the highest co-occurrence
    score to the results.

    :param results: dictionary with as key the article ID and as value
    a list with the structure [title, abstract, genelist, diseaselist,
    hyperlink]
    :param diseasepoints: the list with all the points per gene and
    disease combination per article

    :return results: dictionary with as key the article ID and as value
    a list with the structure [title, abstract, genelist, diseaselist,
    hyperlink, publication date, genepanels]
    """
    for key, value in results.items():
        diseasesperarticle = []
        for key2, value2 in diseasepoints.get(key).items():
            diseasespergene = []
            # Add the 3 diseases with the highest score to the results
            for disease in sorted(value2, reverse=True)[:3]:
                diseasespergene.append(disease[1])
            diseasesperarticle.append(diseasespergene)
        results[key].append(diseasesperarticle)
    return results
