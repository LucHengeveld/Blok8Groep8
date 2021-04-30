from flask import Flask, render_template, request

app = Flask(__name__)

@app.route('/')
@app.route('/home.html', methods=["POST", "GET"])
def get_input():
    if request.method == 'POST':

        or_filter = request.form.get("or_filter", "")
        or_list = request.form.getlist('or_list')
        or_list.insert(0, or_filter)

        and_filter = request.form.get("and_filter", "")
        not_filter = request.form.get("not_filter", "")
        gene_filter = request.form.get("gene_filter", "")
        date_filter = request.form.get("date_filter", "")
        genepanel_filter = request.form.get("genepanel_filter", "")

        print(or_list)
        print(and_filter)
        print(not_filter)
        print(gene_filter)
        print(date_filter)
        print(genepanel_filter)

        return render_template("home.html",
                               or_filter=or_filter,
                               or_list=or_list,
                               and_filter=and_filter,
                               not_filter=not_filter,
                               gene_filter=gene_filter,
                               date_filter=date_filter,
                               genepanel_filter=genepanel_filter)
    else:
        return render_template("home.html",
                               or_filter="",
                               or_list="",
                               and_filter="",
                               not_filter="",
                               gene_filter="",
                               date_filter="",
                               genepanel_filter="")


@app.route('/info.html', methods=["POST", "GET"])
def info():
    return render_template('info.html')


if __name__ == '__main__':
    app.run()

