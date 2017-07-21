from flask import Flask, render_template, request, jsonify, send_file
import parsers

app = Flask(__name__)

import logging
log = logging.getLogger('werkzeug')
log.setLevel(logging.ERROR)


def run(prefix):
    app.config['prefix'] = prefix
    app.config['data'] = parsers.read_data(prefix)
    app.run()


@app.route('/', methods=['GET', 'POST'])
def index():
    return 'data:{}'.format(app.config.get('data'))


@app.route('/regions', methods=['GET', 'POST'])
def regions():
    record = {
        'region': 'ERCC5_15',
        'rc': '494',
        'medcov': '72',
        'mincov': '22',
        'medqcov': '69',
        'minqcov': '22',
        'maxflbq': '0.141',
        'maxflmq': '0.019',
        'pass': 'PASS'
    }
    results = []
    for _ in range(100):
        results.append(record)

    record = {
        'region': 'BRCA2_8',
        'rc': '494',
        'medcov': '72',
        'mincov': '22',
        'medqcov': '69',
        'minqcov': '22',
        'maxflbq': '0.141',
        'maxflmq': '0.019',
        'pass': 'PASS'
    }
    for _ in range(100):
        results.append(record)

    return render_template('regions.html', results=results)