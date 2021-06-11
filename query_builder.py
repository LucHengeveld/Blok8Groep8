"""
This module makes the query.
"""


def make_query(or_list, and_filter, not_filter):
    """
    This function creates and returns a query out of the entered
    keywords.

    :param or_list: list with all entered or keywords
    :param and_filter: string with the entered and keywords
    :param not_filter: string with the entered not keywords

    :return query: string with a pubmed search query
    """
    # Creates empty string
    query = ""

    # If the or_list is not empty it adds the values to the query string
    # Also adds brackets and [tiab]
    if or_list:
        for i in range(len(or_list)):
            if ", " in or_list[i]:
                keywords = or_list[i].split(", ")
                for j in range(len(keywords)):
                    keywords[j] = "(" + str(keywords[j]) + " [tiab])"
                    if keywords[j] != keywords[-1]:
                        query += keywords[j] + " OR "
                    else:
                        if or_list[i] != or_list[-1]:
                            query += keywords[j] + " AND "
                        else:
                            query += keywords[j]
            else:
                if or_list[i] != or_list[-1]:
                    query += "(" + or_list[i] + " [tiab]) AND "
                else:
                    query += "(" + or_list[i] + " [tiab])"

    # If the and_filter is not empty it adds the values to the query
    # string. Also adds brackets and [tiab]
    if and_filter:
        if query:
            if ", " in and_filter:
                keywords = and_filter.split(", ")
                for i in range(len(keywords)):
                    query += " AND (" + keywords[i] + " [tiab])"
            else:
                query += " AND (" + and_filter + " [tiab])"
        else:
            if ", " in and_filter:
                keywords = and_filter.split(", ")
                for i in range(len(keywords)):
                    if i == 0:
                        query += "(" + keywords[i] + " [tiab])"
                    else:
                        query += " AND (" + keywords[i] + " [tiab])"
            else:
                query += "(" + and_filter + " [tiab])"

    # If the not_filter is not empty it adds the values to the query
    # string. Also adds brackets and [tiab]
    if not_filter:
        if query:
            if ", " in not_filter:
                keywords = not_filter.split(", ")
                for i in range(len(keywords)):
                    query += " NOT (" + keywords[i] + " [tiab])"
            else:
                query += " NOT (" + not_filter + " [tiab])"
        else:
            if ", " in not_filter:
                keywords = not_filter.split(", ")
                for i in range(len(keywords)):
                    if i == 0:
                        query += "NOT (" + keywords[i] + " [tiab])"
                    else:
                        query += " NOT (" + keywords[i] + " [tiab])"
            else:
                query += "NOT (" + not_filter + " [tiab])"

    # Returns the created query
    return query
