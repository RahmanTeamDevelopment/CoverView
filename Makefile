HEADERS=coverview/*.pxd
PYX=coverview/*.pyx
PY=coverview/*.py bamgen/*.py bin/CoverView.py
PEP8=pep8 --max-line-length=120
SCRIPTS=bin/coverview

env:
	virtualenv -p python2.7 --no-site-packages --always-copy env

pep8:
	${PEP8} ${PY}

clean:
	pip uninstall -y CoverView

wheels:
	pip wheel .

env/bin/coverview: ${HEADERS} ${PYX} ${PY} ${SCRIPTS}
	pip install .
	touch env/bin/coverview

install: env/bin/coverview

profile: libs
	source env/bin/activate; time python -m cProfile -s cumulative env/bin/CoverView.py --input ../Data/NA21144.mapped.ILLUMINA.bwa.GIH.exome.20121211.bam -b chrom20_exons.bed > profile.out

test: check libs
	source env/bin/activate; CoverView.py --input test/16768_sorted_picard.bam -b test/TSCP_coverviewInput.bed -c test/CoverView_default.json

unittest: pep8 install
	pytest test
