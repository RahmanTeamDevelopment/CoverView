HEADERS=coverview/*.pxd
PYX=coverview/*.pyx
PY=coverview/*.py bamgen/*.py
PEP8=pep8 --max-line-length=120

pep8:
	${PEP8} ${PY}

clean:
	rm -rf build
	rm -rf dist
	rm -rf CoverView.egg-info
	pip uninstall -y CoverView

wheels:
	pip wheel .

env/bin/coverview: ${HEADERS} ${PYX} ${PY}
	pip install .

install: env/bin/coverview

profile: libs
	source env/bin/activate; time python -m cProfile -s cumulative env/bin/CoverView.py --input ../Data/NA21144.mapped.ILLUMINA.bwa.GIH.exome.20121211.bam -b chrom20_exons.bed > profile.out

test: check libs
	source env/bin/activate; CoverView.py --input test/16768_sorted_picard.bam -b test/TSCP_coverviewInput.bed -c test/CoverView_default.json

unittest: pep8 install
	pytest test
