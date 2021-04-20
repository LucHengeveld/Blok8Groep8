"""
Margo Raijmakers
18-04-2021
"""

from flask import Flask, render_template, request

app = Flask(__name__)


@app.route('/', methods=["POST", "GET"])
def get_filters():
    """
    todo docstrings
    """
    if request.method == "POST":
        or_search = request.form.get("or_search", "")
        and_search = request.form.get("and_search", "")
        gene_search = request.form.get("gene_search", "")
        date_search = request.form.get("date_search", "")
        genepanel_search = request.form.get("or_search", "")
        print(or_search)
        print(and_search)
        print(gene_search)
        print(date_search)
        print(genepanel_search)
        return render_template("WebApp.html", or_search=or_search,
                               and_search=and_search,
                               gene_search=gene_search,
                               date_search=date_search,
                               genepanel_search=genepanel_search)
    else:
        return render_template("WebApp.html", or_search="",
                               and_search="",
                               gene_search="",
                               date_search="",
                               genepanel_search="")


if __name__ == '__main__':
    app.run()
