HEADERS=coverview/*.pxd
PYX=coverview/*.pyx
PY=coverview/*.py bamgen/*.py bin/CoverView.py
PEP8=pep8 --max-line-length=120
SCRIPTS=bin/coverview

.ONESHELL:
env:
	virtualenv -p python2.7 --no-site-packages --always-copy env
	source env/bin/activate
	pip install -U pip
	pip install -r requirements.txt --no-cache-dir --ignore-installed

.ONESHELL:
pep8:
	${PEP8} ${PY}

.ONESHELL:
clean:
	source env/bin/activate
	pip uninstall -y CoverView

.ONESHELL:
wheels: env
	source env/bin/activate
	pip wheel .

.ONESHELL:
env/bin/coverview: ${HEADERS} ${PYX} ${PY} ${SCRIPTS} env
	source env/bin/activate
	pip install . --no-cache-dir --ignore-installed
	touch env/bin/coverview

install: env/bin/coverview

.ONESHELL:
profile: install
	source env/bin/activate
	time python -m cProfile -s cumulative env/bin/CoverView.py --input ../Data/NA21144.mapped.ILLUMINA.bwa.GIH.exome.20121211.bam -b chrom20_exons.bed > profile.out

.ONESHELL:
regression_test: install pep8
	source env/bin/activate
	coverview --input test/16768_sorted_picard.bam -b test/TSCP_coverviewInput.bed -c test/CoverView_default.json

.ONESHELL:
unittest: pep8 install
	pytest test
