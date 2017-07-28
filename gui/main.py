from flask import Flask, render_template, request, session, redirect, url_for, jsonify, send_file
import parsers

app = Flask(__name__)

app.secret_key = "something-from-os.urandom(24)"


import logging
log = logging.getLogger('werkzeug')
log.setLevel(logging.ERROR)


def run(prefix):
    app.config['prefix'] = prefix
    app.config['data'] = parsers.read_data(prefix)
    app.run()


@app.route('/', methods=['GET', 'POST'])
def index():
    return redirect('regions')


@app.route('/regions', methods=['GET', 'POST'])
def regions():
    if request.method == 'POST':
        session['region'] = request.form['region']
        return request.form['region']
    else:
        return render_template('regions.html', results=app.config.get('data')['regions'])


@app.route('/profiles', methods=['GET', 'POST'])
def profiles():
    return render_template('profiles.html', region=session['region'])


@app.route('/genes', methods=['GET', 'POST'])
def genes():
    return render_template('genes.html', region=session['region'])


@app.route('/chromosomes', methods=['GET', 'POST'])
def chromosomes():
    return render_template('chroms.html', region=session['region'])


@app.route('/analysis', methods=['GET', 'POST'])
def analysis():
    return render_template('analysis.html')

