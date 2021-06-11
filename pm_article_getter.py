"""
This module retrieves the PubMed articles.
"""


from datetime import datetime
from Bio import Entrez


def get_pubmed_ids(query, date_filter):
    """
    This function enters the query into an esearch and retrieves the
    pubmed id's from the found articles.

    :param query: Pubmed search query
    :param date_filter: date input from the web application

    :return id_list: List with all found pubmed id's
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
