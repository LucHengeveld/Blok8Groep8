def var():
    """
    Handmatig de variabelen voor de query beheren. Deze functie kan
    later weg, want deze functie wordt gebruikt om te kijken hoe de
    query in elkaar komt te zitten.
    De or_list is verplicht!
    """
    or_list = ['ABC transporters, transporter, transport', 'disease, mutations, mutation, liver disease',
               'Lipids, cholesterol, bile salts, canalicular '
               'membrane, phosphatidylcholine, PC ']
    and_filter = "ABCB4, ABCB4 deficiency"
    not_filter = "human, monkey"
    gene_filter = "celline, gene"

    return or_list, and_filter, not_filter, gene_filter


def retrieve_data(or_list, and_filter, not_filter, gene_filter):
    """
    In deze functie worden de onderdelen voor de query klaargemaakt.
    Als eerste wordt er gecontroleerd of de lijsten leeg zijn of niet.
    Als ze niet leeg zijn, dan worden komma's en de haakjes
    gereplaced met de juiste onderdelen (bv. OR en AND)
    """
    try:
        # De or_list wordt hier bewerkt
        if or_list is not None:
            or_list = str(or_list)
            if "['" in or_list:  # [] wordt vervangen door ()
                list = or_list.replace("['", "(").replace("']", ")")
                # print("List: ", list)
            if "" in list:  # komma wordt vervangen door OR
                or_list = list.replace(",", " OR")
                if ' OR ' in or_list:  # sommige OR's worden
                    # vervangen door AND
                    or_data = or_list.replace("' OR '", ") AND (")
                print("OR: ", or_data)
        else:
            or_list = str(or_list)
            print(or_list)

        # De and_filter lijst wordt hier bewerkt
        if and_filter is not None:
            and_filter = str(and_filter)
            # De komma's worden vervangen door OR
            and_filter = and_filter.replace(",", " OR")
            print("AND: ", and_filter)
        else:
            and_filter = str(and_filter)
            print(and_filter)

        # De not_filter lijst wordt hier bewerkt
        if not_filter is not None:
            not_filter = str(not_filter)
            # De komma's worden vervangen door NOT
            not_filter = not_filter.replace(",", " NOT")
            print("NOT: ", not_filter)
        else:
            not_filter = str(not_filter)
            print(not_filter)

        # De gene_filter lijst wordt hier bewerkt
        if gene_filter is not None:
            gene_filter = str(gene_filter)
            # De komma's worden vervangen door OR
            gene_filter = gene_filter.replace(",", " OR")
            print("Gene filter: ", gene_filter)
        else:
            gene_filter = str(gene_filter)
            print(gene_filter)

        return or_data, and_filter, not_filter, gene_filter
    except ValueError:
        print("Error: something went wrong. Please check the info "
              "page.")


def making_query(or_list, and_filter, not_filter, gene_filter):
    """
    In deze functie wordt de query in elkaar gezet. Hierbij wordt
    gekeken naar welke van de 4 velden zijn ingevuld en vervolgens
    wordt het bijbehorende stuk aan de lege lijst query toegevoegd.
    de or_list is verplicht.
    """
    try:
        query = []  # lege lijst wordt aangemaakt

        if or_list != "":  # query voor de or_list wordt aan de
            # lege lijst toegevoegd
            query_or = or_list
            query.append(query_or)
            # print("Query or: ", query_or)
        else:
            pass

        if and_filter != "":  # query voor de and_filter wordt aan
            # de lijst toegevoegd
            query_and = " AND (", and_filter
            query_and = str(query_and).replace(" ', '",
                " ").replace("'", "").replace(", ","")
            query.append(query_and)
            # print("Query and: ", query_and)
        else:
            pass

        if not_filter != "":  # query voor de not_filter wordt
            # aan de lijst toegevoegd
            query_not = " AND (NOT", not_filter
            query_not = str(query_not).replace(" ', '",
                "").replace("'", "").replace(",", "")
            query.append(query_not)
            # print("Query not: ", query_not)
        else:
            pass

        if gene_filter != "":  # query voor de
            # gene_filter wordt aan de lijst toegevoegd
            query_gene = " AND ", gene_filter
            query_gene = str(query_gene).replace("'",
                "").replace(",", "(")
            query.append(query_gene)
            # print("Query gene: ", query_gene)
        else:
            pass

        query = str(query).replace("', '( ", " ").replace("['", "(").replace("']", ")")
        print("Query lijst: ", query)
        return query

    except ValueError:
        print("Error: something went wrong. Please check the info "
              "page.")


def main():
    or_list, and_filter, not_filter, gene_filter = var()
    or_list, and_filter, not_filter, gene_filter = retrieve_data(
        or_list, and_filter, not_filter, gene_filter)
    making_query(or_list, and_filter, not_filter, gene_filter)


main()
