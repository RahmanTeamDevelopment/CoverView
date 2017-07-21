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
