from flask import Flask, render_template, request, jsonify, send_file

app = Flask(__name__)

import logging
log = logging.getLogger('werkzeug')
log.setLevel(logging.ERROR)


def run(prefix):
    app.config['prefix'] = prefix
    app.config['data'] = read_data(prefix)
    app.run()


def read_data(prefix):
    return 'content'


@app.route('/', methods=['GET', 'POST'])
def index():
    return 'Data from {}: {}'.format(app.config.get('prefix'), app.config.get('data'))
