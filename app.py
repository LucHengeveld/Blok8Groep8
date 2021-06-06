from flask import Flask, render_template, request, send_file
from Bio import Entrez
import csv
from xlsxwriter.workbook import Workbook
import ast

import genepanel_reader as gpr
import query_builder as qb
import pm_article_getter as pmag
import results_getter as rg
import co_occurrence_calculator as coc
import relevance_score_calculator as rsc

app = Flask(__name__)


@app.route('/')
@app.route('/home.html', methods=["POST", "GET"])
def get_input():
    """
    This function retrieves all filters from the webapplication, calls
    all other functions and renders the results page.

    :return render template: results.html with all entered filters and
    results
    """
    if request.method == 'POST':
        try:
            or_filter = request.form["or_filter"]
            or_list = request.form.getlist('or_list')
            or_list.insert(0, or_filter)

            and_filter = request.form["and_filter"]
            not_filter = request.form["not_filter"]
            gene_filter = request.form["gene_filter"]
            date_filter = request.form["date_filter"]

            try:
                genepanel_file = request.files["genepanel_file"]
                genepanel_file_name = genepanel_file.filename
                genepanel_file.save(genepanel_file_name)
            except:
                genepanel_file = ""
                genepanel_file_name = "No file selected."

            genepanel_filter = request.form["genepanel_filter"]

            email = request.form["email"]

            try:
                use_co_occurence = request.form["occurence"]
            except:
                use_co_occurence = "Not selected"

            Entrez.email = email

            if genepanel_file:
                gp_table = gpr.excel_reader(genepanel_file)
                genes = gpr.get_column(gp_table, "GenePanels_Symbol")
                synonyms = gpr.get_column(gp_table, "Aliases")
                gene_panels_list = gpr.get_column(gp_table, "GenePanel")
                gps_list, gps_set = gpr.make_genepanel_list_set(
                    gene_panels_list)
                genes_dict = gpr.make_gene_dict(genes, gps_list)
                gene_panel_dict = gpr.make_gene_panel_dict(gps_set, genes_dict)
                genes_dict, gene_panel_dict = gpr.gene_synonyms(synonyms,
                                                                gps_list,
                                                                genes_dict,
                                                                gene_panel_dict
                                                                )
            else:
                gene_panel_dict = {}
                genes_dict = {}

            query = qb.make_query(or_list, and_filter, not_filter)

            if not query:
                return render_template("home_error.html",
                                       or_list=or_list,
                                       and_filter=and_filter,
                                       not_filter=not_filter,
                                       gene_filter=gene_filter,
                                       date_filter=date_filter,
                                       genepanel_file_name=genepanel_file_name,
                                       genepanel_filter=genepanel_filter,
                                       email=email,
                                       use_co_occurence=use_co_occurence)

            id_list = pmag.get_pubmed_ids(query, date_filter)
            pubtator_link = rg.get_pubtator_link(id_list)

            results = rg.read_pubtator_file(pubtator_link, gene_panel_dict,
                                            genepanel_filter, gene_filter)
            results = rg.pubtator_hyperlink(results)
            results = rg.publication_date(results)
            results = rg.genepanel_results(results, genes_dict)

            diseasepoints = coc.co_occurrence(results, 3)
            results = coc.add_co_occurrence_to_results(results, diseasepoints)

            filters = rsc.get_values_for_relevance(or_list, and_filter,
                                                   gene_filter)
            relevance_score = rsc.get_relevance_score(results, filters)
            relevance_score = rsc.sort_relevance_score(relevance_score)

            return render_template("results.html",
                                   or_list=or_list,
                                   and_filter=and_filter,
                                   not_filter=not_filter,
                                   gene_filter=gene_filter,
                                   date_filter=date_filter,
                                   genepanel_file_name=genepanel_file_name,
                                   genepanel_filter=genepanel_filter,
                                   email=email,
                                   use_co_occurence=use_co_occurence,
                                   results=results,
                                   relevance_score=relevance_score)
        except:
            or_filter = request.form["or_filter"]
            or_list = request.form.getlist('or_list')
            or_list.insert(0, or_filter)

            and_filter = request.form["and_filter"]
            not_filter = request.form["not_filter"]
            gene_filter = request.form["gene_filter"]
            date_filter = request.form["date_filter"]

            try:
                genepanel_file = request.files["genepanel_file"]
                genepanel_file_name = genepanel_file.filename
                genepanel_file.save(genepanel_file_name)
            except:
                genepanel_file = ""
                genepanel_file_name = "No file selected."

            genepanel_filter = request.form["genepanel_filter"]

            email = request.form["email"]

            try:
                use_co_occurence = request.form["occurence"]
            except:
                use_co_occurence = "Not selected"

            return render_template("home_error.html",
                                   or_list=or_list,
                                   and_filter=and_filter,
                                   not_filter=not_filter,
                                   gene_filter=gene_filter,
                                   date_filter=date_filter,
                                   genepanel_file_name=genepanel_file_name,
                                   genepanel_filter=genepanel_filter,
                                   email=email,
                                   use_co_occurence=use_co_occurence)

    else:
        return render_template("home.html",
                               or_list="",
                               and_filter="",
                               not_filter="",
                               gene_filter="",
                               date_filter="",
                               genepanel_file_name="",
                               genepanel_filter="",
                               email="",
                               use_co_occurence="",
                               results="")


@app.route('/results.html', methods=["POST"])
def save_results():
    """
    This function will write and return the results in a file when the
    "Save results" button is pressed on the webapplication. There are 3
    options for file extensions, .txt, .tsv and .xlsx

    :return send_file: downloads the file with the results to the user's
    pc
    """
    results = request.form['results']
    results = ast.literal_eval(results)
    try:
        selected_extension = request.form["file_extension"]
    except:
        selected_extension = "tsv"

    if selected_extension == "txt":
        output_file = "results.txt"
    else:
        output_file = "results.tsv"

    # Write results to file
    with open(output_file, 'w', newline='') as out_file:
        tsv_writer = csv.writer(out_file, delimiter='\t')
        tsv_writer.writerow(
            ["Gene name", "Gene ID", "Gene Panels", "Pubmed ID",
             "Pubmed Hyperlink", "Publication Date", "Possible diseases"])
        for pubmed_id in results:
            for gene_index in range(len(results[pubmed_id][2])):
                gene_name = results[pubmed_id][2][gene_index].rsplit(" ", 1)[0]
                gene_id = results[pubmed_id][2][gene_index].rsplit(" ", 1)[1]
                genepanelstring = ""
                genepanels = results[pubmed_id][6][gene_index]
                for i in genepanels:
                    genepanelstring += i + ";"
                genepanelstring = genepanelstring[:-1]
                hyperlink = results[pubmed_id][4]
                date = results[pubmed_id][5]
                try:
                    diseases = results[pubmed_id][-1][gene_index]
                except IndexError:
                    diseases = ["-"]
                diseasestring = ""
                for j in diseases:
                    if j != diseases[-1]:
                        diseasestring += j + ","
                    else:
                        diseasestring += j

                tsv_writer.writerow(
                    [gene_name, gene_id, genepanelstring, pubmed_id, hyperlink,
                     date, diseasestring])

    out_file.close()

    if selected_extension == "xlsx":
        xlsx_file = 'results.xlsx'
        workbook = Workbook(xlsx_file, {'strings_to_numbers': True})
        worksheet = workbook.add_worksheet()
        tsv_reader = csv.reader(open(output_file, 'rt'), delimiter='\t')
        for row, data in enumerate(tsv_reader):
            worksheet.write_row(row, 0, data)
        workbook.close()

    if selected_extension == "xlsx":
        output = "results.xlsx"
    elif selected_extension == "txt":
        output = "results.txt"
    else:
        output = "results.tsv"

    if selected_extension == "xlsx":
        return send_file(attachment_filename=output, filename_or_fp=output,
                         mimetype="xlsx", as_attachment=True)
    else:
        return send_file(attachment_filename=output, filename_or_fp=output,
                         mimetype="text/tsv", as_attachment=True)


@app.route('/info.html', methods=["POST", "GET"])
def info():
    """
    This function shows the info page when the user selects it in the
    menu bar on the webapplication. The info page contains information
    about the application

    :return render template: shows the info.html page to the user
    """
    return render_template('info.html')


if __name__ == '__main__':
    app.run()
