from flask import Flask, render_template, request, redirect, jsonify
import parsers
from coverview_ import reference
import helper
from gevent.wsgi import WSGIServer
import signal
import sys
import webbrowser

app = Flask(__name__)

app.secret_key = "something-from-os.urandom(24)"


def run(prefix, reffn):
    app.config['prefix'] = prefix
    app.config['data'] = parsers.read_data(prefix)
    app.config['regionlist'] = [x['region'] for x in app.config['data']['regions'] if x['region'] in app.config['data']['region_coords'] and region_size(x['region']) > 1]
    app.config['passedregions'] = [x['region'] for x in app.config['data']['regions'] if x['pass_or_flag'] == 'PASS']
    app.config['region'] = ''
    app.config['metadata'] = parsers.read_metadata(prefix)
    app.config['gene'] = ''
    app.config['flagged'] = False
    app.config['togene'] = ''

    sequences = {}
    ref = reference.Reference(reffn)
    for region, coords in app.config['data']['region_coords'].iteritems():
        sequences[region] = ref.getSequence(coords[0], int(coords[1]), int(coords[2]))
    app.config['sequences'] = sequences

    app.config['all_genes'], app.config['genes_by_chrom'], app.config['flagged_regions_stat'], app.config['flagged_genes_stat'] \
        = helper.create_flag_statistics(app.config.get('data')['regions'], app.config['data']['region_coords'])

    signal.signal(signal.SIGINT, signal_handler)

    webbrowser.open('http://127.0.0.1:5000/', new=2)

    http_server = WSGIServer(('', 5000), app, log=None)
    http_server.serve_forever()




def signal_handler(signal, frame):
    print('\n\nCoverView GUI is closed. Bye!\n')
    sys.exit(0)


def region_size(region):
    x = app.config['data']['region_coords'][region]
    return int(x[2])-int(x[1])


@app.route('/', methods=['GET', 'POST'])
def index():
    return redirect('regions')


@app.route('/regions', methods=['GET', 'POST'])
def regions():
    if request.method == 'POST':
        app.config['region'] = request.json['region']
        app.config['gene'] = request.json['gene']
        app.config['flagged'] = request.json['flagged']
        return request.json['region']
    else:
        return render_template(
            'regions.html',
            region=app.config['region'],
            results=app.config.get('data')['regions'],
            pass_def=app.config['metadata']['config_opts']['pass'],
            gene=app.config['gene'],
            flagged=app.config['flagged'],
            regionlist=app.config['regionlist']
        )


@app.route('/profiles', methods=['GET', 'POST'])
def profiles():
    if request.method == 'POST':
        app.config['region'] = request.json['region']
        if request.json['region'] == '':
            return ''
        else:
            return jsonify(app.config.get('data')['profiles'][app.config['region']])
    else:
        return render_template(
            'profiles.html',
            region=app.config['region'],
            regionlist=app.config['regionlist'],
            passedregions=app.config['passedregions'],
            regioncoords=app.config['data']['region_coords'],
            sequences=app.config['sequences'],
            passdef=app.config['metadata']['config_opts']['pass']
        )


@app.route('/genes', methods=['GET', 'POST'])
def genes():
    if request.method == 'POST':
        app.config['togene'] = request.json['gene']
        return request.json['gene']
    else:
        return render_template(
            'genes.html',
            region=app.config['region'],
            all_genes=app.config['all_genes'],
            genes_by_chrom=app.config['genes_by_chrom'],
            summary=app.config['data']['summary'],
            togene=app.config['togene']
        )


@app.route('/analysis', methods=['GET', 'POST'])
def analysis():
    return render_template(
        'analysis.html',
        metadata=app.config['metadata'],
        regions_stat=app.config['flagged_regions_stat'],
        genes_stat=app.config['flagged_genes_stat'],
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


@app.route('/usage', methods=['GET', 'POST'])
def usage():
    return render_template('usage.html')