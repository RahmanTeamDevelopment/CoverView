from flask import Flask, render_template, request, session, redirect, url_for, jsonify, send_file
import parsers
from coverview import reference
import helper

app = Flask(__name__)

app.secret_key = "something-from-os.urandom(24)"

import logging
log = logging.getLogger('werkzeug')
log.setLevel(logging.ERROR)

def run(prefix, reffn):
    app.config['prefix'] = prefix
    app.config['data'] = parsers.read_data(prefix)
    app.config['regionlist'] = [x['region'] for x in app.config['data']['regions']]
    app.config['passedregions'] = [x['region'] for x in app.config['data']['regions'] if x['pass_or_fail'] == 'PASS']
    app.config['region'] = ''
    app.config['metadata'] = parsers.read_metadata(prefix)

    sequences = {}
    ref = reference.Reference(reffn)
    for region, coords in app.config['data']['region_coords'].iteritems():
        sequences[region] = ref.getSequence(coords[0], int(coords[1]), int(coords[2]))
    app.config['sequences'] = sequences

    app.config['genes_by_chrom'], app.config['failed_regions_stat'], app.config['failed_genes_stat'] \
        = helper.create_fail_statistics(app.config.get('data')['regions'], app.config['data']['region_coords'])

    app.run()


@app.route('/', methods=['GET', 'POST'])
def index():
    return redirect('regions')


@app.route('/regions', methods=['GET', 'POST'])
def regions():
    if request.method == 'POST':
        app.config['region'] = request.form['region']
        return request.form['region']
    else:
        return render_template('regions.html', region=app.config['region'], results=app.config.get('data')['regions'])


@app.route('/profiles', methods=['GET', 'POST'])
def profiles():
    if request.method == 'POST':
        region = request.json['region']
        return jsonify(app.config.get('data')['profiles'][region])
    else:
        return render_template(
            'profiles.html',
            region=app.config['region'],
            regionlist=app.config['regionlist'],
            passedregions=app.config['passedregions'],
            regioncoords=app.config['data']['region_coords'],
            sequences=app.config['sequences']
        )


@app.route('/genes', methods=['GET', 'POST'])
def genes():
    return render_template('genes.html', region=app.config['region'])


@app.route('/analysis', methods=['GET', 'POST'])
def analysis():
    return render_template(
        'analysis.html',
        metadata=app.config['metadata'],
        regions_stat=app.config['failed_regions_stat'],
        genes_stat=app.config['failed_genes_stat'],
        pass_def=create_pass_critera_string()
    )

def create_pass_critera_string():
    criteria = app.config['metadata']['config_opts']['pass']
    ret = []
    for k,v in criteria.iteritems():
        metrics = k[:-4]
        relation = '>=' if k[-3:] == 'MIN' else '<='
        ret.append('{} {} {}'.format(metrics, relation, v))
    return '; '.join(ret)


