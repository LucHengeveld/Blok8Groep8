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

        data = retrieve_data(or_list, and_filter, not_filter, gene_filter)

        return render_template("home.html",
                               email=email,
                               or_filter=or_filter,
                               or_list=or_list,
                               and_filter=and_filter,
                               not_filter=not_filter,
                               gene_filter=gene_filter,
                               date_filter=date_filter,
                               genepanel_file=genepanel_file,
                               genepanel_filter=genepanel_filter,
                               retrieve_data=data)
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
                               genepanel_filter="",
                               retrieve_data="")


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
            # Comma's get replaced by OR
            and_filter = and_filter.replace(",", " [tiab] OR")
            # print("AND search: ", and_filter)
        else:
            and_filter = str(and_filter)


        # The not_filter list will be edited here
        if not_filter is not None:
            not_filter = str(not_filter)
            # Comma's get replaced by NOT
            not_filter = not_filter.replace(",", " [tiab] NOT")
            # print("NOT search: ", not_filter)
        else:
            not_filter = str(not_filter)

        # # The gene_filter list will be edited here
        if gene_filter is not None:
            gene_filter = str(gene_filter)
            # Comma's get replaced by OR
            gene_filter = gene_filter.replace(",", " [tiab] OR")
            # print("Gene filter search: ", gene_filter)
        else:
            gene_filter = str(gene_filter)

        making_query(or_list, and_filter, not_filter, gene_filter)
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


@app.route('/info.html', methods=["POST", "GET"])
def info():
    return render_template('info.html')


if __name__ == '__main__':
    app.run()
