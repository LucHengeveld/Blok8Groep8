from Bio import Entrez
import requests
from datetime import datetime

Entrez.email = "luchengeveld@hotmail.com"


def search_query():
    query = "((ABC transporter [tiab] OR transporter [tiab] OR transport [" \
            "tiab]) AND (disease [tiab] OR mutation [tiab] OR mutations [" \
            "tiab] OR liver disease [tiab]) AND (lipids [tiab] OR " \
            "cholesterol [tiab] OR bile salts [tiab] OR canalicular membrane "\
            "[tiab] OR phosphatidylcholine [tiab] OR PC [tiab]) AND (ABCB4 [" \
            "tiab] OR ABCB4 deficiency [tiab])) "
    return query


def get_pubmed_ids(query, date_filter):
    date = date_filter.replace("-", "/")
    datetoday = str(datetime.date(datetime.now())).replace("-", "/")
    handle = Entrez.esearch(db="pubmed", term=query, mindate=date,
                            maxdate=datetoday)
    record = Entrez.read(handle)
    handle.close()
    id_list = (record["IdList"])

    get_pubtator_link(id_list)


def get_pubtator_link(id_list):
    link = "https://www.ncbi.nlm.nih.gov/research/pubtator-api/publications" \
           "/export/pubtator?pmids=idvalues&concepts=gene,mutation,disease"
    id_string = ""
    for i in range(len(id_list)):
        id_string += id_list[i] + ","
    id_string = id_string[:-1]
    pubtator_link = link.replace("idvalues", id_string)
    print(pubtator_link)

    read_pubtator_file(pubtator_link)


def read_pubtator_file(pubtator_link):
    pubtator_text = requests.get(pubtator_link).text
    lines = pubtator_text.split("\n")
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

        else:
            if genelist:
                # print(article_id)
                # print(title)
                # print(abstract)
                # print(genelist)
                # print(diseaselist)
                # print(mutationlist)
                # print("")
                results[article_id] = [title, abstract, genelist, diseaselist,
                                       mutationlist]
                genelist = []
                diseaselist = []
                mutationlist = []

    # for key in results:
    #     print(key + "\t" + str(results[key]))
    return results


def main():
    # output van website: jaar-maand-dag bv 2021-05-02 is 2 mei 2021
    date_filter = "2020-01-01"
    query = search_query()
    get_pubmed_ids(query, date_filter)


main()
