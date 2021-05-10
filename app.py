import re
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
        or_list = request.form.getlist('or_list')
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
        genepanel_file = "C:/Users/luche/Documents/HAN/Leerjaar_2/Informatica Jaar 2/Blok 8/GenPanelOverzicht_DG-3.1.0_HAN.xlsx"
        gp_table = excel_reader(genepanel_file)
        genes = get_column(gp_table, "GenePanels_Symbol")
        gene_panels_list = get_column(gp_table, "GenePanel")
        gps_list, gps_set = make_genepanel_list_set(gene_panels_list)
        genes_dict = make_gene_dict(genes, gps_list)
        gene_panel_dict = make_gene_panel_dict(gps_set, genes_dict)
        or_list2, and_filter2, not_filter2, gene_filter2 = \
            retrieve_data(or_list, and_filter, not_filter, gene_filter)
        query = making_query(or_list, and_filter, not_filter, gene_filter, or_list2, and_filter2, not_filter2, gene_filter2)
        query = "((ABC transporter [tiab] OR transporter [tiab] OR transport [" \
                "tiab]) AND (disease [tiab] OR mutation [tiab] OR mutations [" \
                "tiab] OR liver disease [tiab]) AND (lipids [tiab] OR " \
                "cholesterol [tiab] OR bile salts [tiab] OR canalicular membrane " \
                "[tiab] OR phosphatidylcholine [tiab] OR PC [tiab]) AND (ABCB4 [" \
                "tiab] OR ABCB4 deficiency [tiab])) "
        id_list = get_pubmed_ids(query, date_filter)
        pubtator_link = get_pubtator_link(id_list)
        results = read_pubtator_file(pubtator_link)
        results = pubmed_hyperlink(results)
        results = publication_date(results)
        results = genepanel_results(results, genes_dict)

        # Results: ID[title, abstract, genelist, diseaselist, mutationlist, hyperlink, publication date, genepanel_list]
        #               0       1       2           3           4               5           6               7
        # results = {'32942997': ['Whole-exome sequencing reveals ANO8 as a genetic risk factor for intrahepatic cholestasis of pregnancy.', "BACKGROUND: Intrahepatic cholestasis of pregnancy (ICP) is characterized by pruritus and cholestasis in late pregnancy and results in adverse pregnancy outcomes, including preterm delivery and birth weight, which are affected by the genetic and environmental background. However, until now, the genetic architecture of ICP has remained largely unclear. METHODS: Twenty-six clinical data points were recorded for 151 Chinese ICP patients. The data generated from whole-exome sequencing (WES) using the BGISEQ-500 platform were further analyzed by Burrows-Wheeler Aligner (BWA) software, Genome Analysis Toolkit (GATK), ANNOVAR tool, etc. R packages were used to conduct t-test, Fisher's test and receiver operating characteristic (ROC) curve analyses. RESULTS: We identified eighteen possible pathogenic loci associated with ICP disease in known genes, covering ABCB4, ABCB11, ATP8B1 and TJP2. The loci Lys386Gln, Gly527Gln and Trp708Ter in ABCB4, Leu589Met, Gln605Pro and Gln1194Ter in ABCB11, and Arg189Ser in TJP2 were novel discoveries. In addition, WES analysis indicated that the gene ANO8 involved in the transport of bile salts is newly identified as associated with ICP. The functional network of the ANO8 gene confirmed this finding. ANO8 contained 8 rare missense mutations that were found in eight patients among the 151 cases and were absent from 1029 controls. Out of the eight SNPs, 3 were known, and the remaining five are newly identified. These variants have a low frequency, ranging from 0.000008 to 0.00001 in the ExAC, gnomAD - Genomes and TOPMED databases. Bioinformatics analysis showed that the sites and their corresponding amino acids were both highly conserved among vertebrates. Moreover, the influences of all the mutations on protein function were predicted to be damaging by the SIFT tool. Combining clinical data, it was found that the mutation group (93.36 micromol/L) had significantly (P = 0.038) higher total bile acid (TBA) levels than the wild-type group (40.81 micromol/L). CONCLUSIONS: To the best of our knowledge, this is the first study to employ WES technology to detect genetic loci for ICP. Our results provide new insights into the genetic basis of ICP and will benefit the final identification of the underlying mutations.", ['ANO8 57719', 'ABCB4 5244', 'ABCB11 8647', 'ATP8B1 5205', 'TJP2 9414'], ['gallstones MESH:D042882', 'obesity MESH:D009765', 'ABCB4 deficiency of the hepatic lecithin MESH:C535935', 'gallstone MESH:D042882', 'acute cholecystitis MESH:D041881', 'cholangitis and biliary pancreatitis MESH:D008105', 'intrahepatic cholestasis MESH:D002780', 'Intrahepatic cholestasis MESH:D002780', 'ICP MESH:D001932', 'pruritus and cholestasis MESH:D011537', 'ICP disease MESH:D001932'], [], 'https://pubmed.ncbi.nlm.nih.gov/32942997/', '2020 Sep 17'], '33231505': ['Hypothyroidism increases cholesterol gallstone prevalence in mice by elevated hydrophobicity of primary bile acids.', 'BACKGROUND: Thyroid hormone (TH) deficiency has been associated with increased cholesterol gallstone prevalence. Hypothyroidism impacts hepatic lipid homeostasis, biliary secretion, gallbladder motility and gallstone (LITH) gene expression, i.e. potential factors contributing to cholesterol gallstone disease (CGD). However, how TH deficiency may propel gallstone formation is still poorly understood. Therefore, we performed molecular studies in a CGD mouse model under lithogenic conditions and modulation of TH status. METHODS: Male three months old C57BL/6 mice were randomly divided into a control (eu), a hypothyroid (hypo), a gallstone (litho) and a gallstone+hypothyroid (litho+hypo) group and treated for two, four and six weeks (n=8/treatment period). Gallstone prevalence, biliary composition and cholesterol crystals, hepatic expression of genes participating in cholesterol, bile acid (BA) and phosphatidylcholine synthesis (Hmgcr, Cyp7a1, Pcyt1a) and canalicular transport (Abcg5, Bsep, Abcb4) were investigated. RESULTS: Increased cholesterol gallstone prevalence was observed in hypothyroid mice under lithogenic diet after four and six weeks of treatment (4 weeks: 25% vs. 0%; 6 weeks: 75% vs. 37.5%). Interestingly, neither the composition of the three main biliary components, cholesterol, BAs and phosphatidylcholine nor the hepatic expression of genes involved in synthesis and transport could explain the differences in cholesterol gallstone formation in the mice. However, TH deficiency resulted in significantly increased hydrophobicity of primary BAs in bile. Furthermore downregulation of hepatic sulfonation enzymes Papss2 and Sult2a8 as well as diminished biliary BA sulfate concentrations in mice, were observed under hypothyroid conditions all contributing to a lithogenic biliary milieu as evidenced by microscopic cholesterol crystals and macroscopic gallstone formation. CONCLUSIONS: We describe a novel pathogenic link between TH deficiency and CGD and suggest that the increased hydrophobic character of biliary BAs due to the diminished expression of hepatic detoxification enzymes promotes cholesterol crystal precipitation and finally enhances cholesterol gallstone formation in the bile of hypothyroid mice.', ['Thyroid hormone 21823', 'TH 21823', 'Hmgcr 15357', 'Cyp7a1 13122', 'Pcyt1a 13026', 'Abcg5 27409', 'Bsep 27413', 'Abcb4 18670', 'Papss2 23972', 'Sult2a8 76971'], ['Hypothyroidism increases cholesterol gallstone MESH:D042882', 'deficiency MESH:D007153', 'cholesterol gallstone MESH:D042882', 'Hypothyroidism impacts hepatic lipid homeostasis MESH:D007037', 'gallbladder motility and gallstone MESH:D042882', 'cholesterol gallstone disease MESH:D042882', 'CGD MESH:D006105', 'TH deficiency MESH:D007153', 'gallstone MESH:D042882', 'hypothyroid MESH:D007037', 'gallstone+hypothyroid MESH:D042882'], [], 'https://pubmed.ncbi.nlm.nih.gov/33231505/', '2021 Jan 5'], '33256620': ['Case report: progressive familial intrahepatic cholestasis type 3 with compound heterozygous ABCB4 variants diagnosed 15 years after liver transplantation.', "BACKGROUND: Progressive familial intrahepatic cholestasis (PFIC) type 3 is an autosomal recessive disorder arising from mutations in the ATP-binding cassette subfamily B member 4 (ABCB4) gene. This gene encodes multidrug resistance protein-3 (MDR3) that acts as a hepatocanalicular floppase that transports phosphatidylcholine from the inner to the outer canalicular membrane. In the absence of phosphatidylcholine, the detergent activity of bile salts is amplified and this leads to cholangiopathy, bile duct loss and biliary cirrhosis. Patients usually present in infancy or childhood and often progress to end-stage liver disease before adulthood. CASE PRESENTATION: We report a 32-year-old female who required cadaveric liver transplantation at the age of 17 for cryptogenic cirrhosis. When the patient developed chronic ductopenia in the allograft 15 years later, we hypothesized that the patient's original disease was due to a deficiency of a biliary transport protein and the ductopenia could be explained by an autoimmune response to neoantigen that was not previously encountered by the immune system. We therefore performed genetic analyses and immunohistochemistry of the native liver, which led to a diagnosis of PFIC3. However, there was no evidence of humoral immune response to the MDR3 and therefore, we assumed that the ductopenia observed in the allograft was likely due to chronic rejection rather than autoimmune disease in the allograft. CONCLUSIONS: Teenage patients referred for liver transplantation with cryptogenic liver disease should undergo work up for PFIC3. An accurate diagnosis of PFIC 3 is key for optimal management, therapeutic intervention, and avoidance of complications before the onset of end-stage liver disease.", ['ABCB4 5244', 'ATP-binding cassette subfamily B member 4 5244', 'multidrug resistance protein-3 5244', 'MDR3 5244', 'PFIC3 5244', 'PFIC 3 5244'], ['familial intrahepatic cholestasis MESH:D002780', 'PFIC MESH:C535933', 'autosomal recessive disorder MESH:D030342', 'cholangiopathy, bile duct loss and biliary cirrhosis MESH:D008105', 'liver disease MESH:D008107', 'cirrhosis MESH:D005355', 'chronic ductopenia MESH:D002908', 'deficiency MESH:D007153', 'ductopenia ', 'autoimmune disease MESH:D001327'], [], 'https://pubmed.ncbi.nlm.nih.gov/33256620/', '2020 Nov 30'], '33340584': ['Synthetic human ABCB4 mRNA therapy rescues severe liver disease phenotype in a BALB/c.Abcb4-/- mouse model of PFIC3.', 'BACKGROUND&AIMS: Progressive familial intrahepatic cholestasis type 3 (PFIC3) is a rare lethal autosomal recessive liver disorder caused by loss-of-function variations of the ABCB4 gene, encoding a phosphatidylcholine transporter (ABCB4/MDR3). Currently, no effective treatment exists for PFIC3 outside of liver transplantation. METHODS: We have produced and screened chemically- and genetically-modified mRNA variants encoding human ABCB4 (hABCB4 mRNA) encapsulated in lipid nanoparticles (LNPs). We examined their pharmacological effects in a cell-based and in a new in vivo mouse model resembling human PFIC3 due to homozygous disruption of the Abcb4 gene in fibrosis-susceptible BALB/c.Abcb4-/- mice. RESULTS: We show that treatment with liver-targeted hABCB4 mRNA resulted in de novo expression of functional hABCB4 protein and restored phospholipid transport in cultured cells and in PFIC3 mouse livers. Importantly, repeated injections of the hABCB4 mRNA effectively rescued the severe disease phenotype in young Abcb4-/- mice, with rapid and dramatic normalization of all clinically relevant parameters such as inflammation, ductular reaction and liver fibrosis. Synthetic mRNA therapy also promoted favorable hepatocyte-driven liver regeneration to restore normal homeostasis, including liver weight, body weight, liver enzymes, and portal vein blood pressure. CONCLUSION: Our data provide strong preclinical proof-of-concept for hABCB4 mRNA therapy as a potential treatment option for patients with PFIC3.', ['ABCB4 5244', 'Abcb4 5244', 'ABCB4 18670', 'MDR3 18671', 'hABCB4 5244', 'Abcb4 18670'], ['liver disease MESH:D008107', 'familial intrahepatic cholestasis MESH:D002780', 'autosomal recessive liver disorder MESH:D017093', 'fibrosis MESH:D005355', 'inflammation MESH:D007249', 'liver fibrosis MESH:D008103'], [], 'https://pubmed.ncbi.nlm.nih.gov/33340584/', '2020 Dec 17'], '33546617': ['Whole-exome sequencing identifies novel mutations in ABC transporter genes associated with intrahepatic cholestasis of pregnancy disease: a case-control study.', 'BACKGROUND: Intrahepatic cholestasis of pregnancy (ICP) can cause premature delivery and stillbirth. Previous studies have reported that mutations in ABC transporter genes strongly influence the transport of bile salts. However, to date, their effects are still largely elusive. METHODS: A whole-exome sequencing (WES) approach was used to detect novel variants. Rare novel exonic variants (minor allele frequencies: MAF < 1%) were analyzed. Three web-available tools, namely, SIFT, Mutation Taster and FATHMM, were used to predict protein damage. Protein structure modeling and comparisons between reference and modified protein structures were performed by SWISS-MODEL and Chimera 1.14rc, respectively. RESULTS: We detected a total of 2953 mutations in 44 ABC family transporter genes. When the MAF of loci was controlled in all databases at less than 0.01, 320 mutations were reserved for further analysis. Among these mutations, 42 were novel. We classified these loci into four groups (the damaging, probably damaging, possibly damaging, and neutral groups) according to the prediction results, of which 7 novel possible pathogenic mutations were identified that were located in known functional genes, including ABCB4 (Trp708Ter, Gly527Glu and Lys386Glu), ABCB11 (Gln1194Ter, Gln605Pro and Leu589Met) and ABCC2 (Ser1342Tyr), in the damaging group. New mutations in the first two genes were reported in our recent article. In addition, compared to the wild-type protein structure, the ABCC2 Ser1342Tyr-modified protein structure showed a slight change in the chemical bond lengths of ATP ligand-binding amino acid side chains. In placental tissue, the expression level of the ABCC2 gene in patients with ICP was significantly higher (P < 0.05) than that in healthy pregnant women. In particular, the patients with two mutations in ABC family genes had higher average values of total bile acids (TBA), aspartate transaminase (AST), direct bilirubin (DBIL), total cholesterol (CHOL), triglycerides (TG) and high-density lipoprotein (HDL) than the patients who had one mutation, no mutation in ABC genes and local controls. CONCLUSIONS: Our present study provide new insight into the genetic architecture of ICP and will benefit the final identification of the underlying mutations.', ['ABC transporter 9429', 'ABC 10058', 'ABCB4 5244', 'ABCB11 8647', 'ABCC2 1244', 'aspartate transaminase 26503', 'AST 26503'], ['intrahepatic cholestasis of pregnancy disease MESH:C535932', 'Intrahepatic cholestasis MESH:D002780', 'stillbirth MESH:D050497', 'CHOL '], [], 'https://pubmed.ncbi.nlm.nih.gov/33546617/', '2021 Feb 5']}
        for key in results: #28587926
            print(results[key][-1])
        return render_template("homeresults.html",
                               or_list=or_list,
                               and_filter=and_filter,
                               not_filter=not_filter,
                               gene_filter=gene_filter,
                               date_filter=date_filter,
                               genepanel_file=genepanel_file,
                               genepanel_filter=genepanel_filter,
                               email=email,
                               results=results)
    else:
        return render_template("home.html",
                               or_list="",
                               and_filter="",
                               not_filter="",
                               gene_filter="",
                               date_filter="",
                               genepanel_file="",
                               genepanel_filter="",
                               email="",
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


def making_query(or_list, and_filter, not_filter, gene_filter, or_list2, and_filter2, not_filter2, gene_filter2):
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
            query_gene = str(query_gene).replace("'", "").replace(", """, "(")\
                .replace(")", " [tiab])")
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
    pubtator link which filters the genes, mutations and diseases out
    of the title and abstract.

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
    """This function reads the pubtator link as a text file and
    retrieves the genes, diseases and mutations out of each article.

    :param pubtator_link: Pubtator link with the title, abstact, genes,
    diseases and mutations of each article.
    :return results: Dictionary with as key the article ID and as value
    a list with the structure [title, abstract, genelist, diseaselist,
    mutationlist]
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
    """This function creates the hyperlink to the pubmed article by
    using the key's (article ID's) from the results dictionary. The
    hyperlink will be added to the values of each key.

    :param results: Dictionary with as key the article ID and as value
    a list with the structure [title, abstract, genelist, diseaselist,
    mutationlist]
    :return results: Dictionary with as key the article ID and as value
    a list with the structure [title, abstract, genelist, diseaselist,
    mutationlist, hyperlink]
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
    """This function retrieves the publication date of the article by
    using the key's of results (article ID's) in efetch. The publication
    date will be added to the values of each key.

    :param results: Dictionary with as key the article ID and as value
    a list with the structure [title, abstract, genelist, diseaselist,
    mutationlist, hyperlink]
    :return results: Dictionary with as key the article ID and as value
    a list with the structure [title, abstract, genelist, diseaselist,
    mutationlist, hyperlink, publication date]
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
    for key in results:
        genepanel = []
        for gene in results[key][2]:
            if gene.rsplit(" ", 1)[0] in genes_dict.keys():
                genepanel.append(genes_dict[gene.rsplit(" ", 1)[0]])
            else:
                genepanel.append("-")
        results[key].append(genepanel)

    return results


@app.route('/homeresults.html', methods=["POST"])
def save_results():
    print("test")
    return render_template("home.html")


@app.route('/info.html', methods=["POST", "GET"])
def info():
    return render_template('info.html')


if __name__ == '__main__':
    app.run()
