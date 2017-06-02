HEADERS=coverview/*.pxd
PYX=coverview/*.pyx
PY=coverview/*.py bamgen/*.py bin/CoverView.py testutils/*.py
PEP8=pep8 --max-line-length=120
SCRIPTS=bin/coverview

env:
	virtualenv -p python2.7 --no-site-packages --always-copy env
	pip install -U pip
	pip install -r requirements.txt --no-cache-dir --ignore-installed

pep8:
	${PEP8} ${PY}

clean:
	pip uninstall -y CoverView
	find . -name __pycache__ | xargs rm -rf

cleanAll:
	rm -rf env

wheels: env
	pip wheel .

env/bin/coverview: ${HEADERS} ${PYX} ${PY} ${SCRIPTS} env
	pip install -U .
	touch env/bin/coverview
	touch env/bin/CoverView.py

install: env/bin/coverview

profile: install
	time python -m cProfile -s cumulative env/bin/CoverView.py --input ../Data/NA21144.mapped.ILLUMINA.bwa.GIH.exome.20121211.bam -b chrom20_exons.bed > profile.out

regression_test: install pep8
	coverview --input test/16768_sorted_picard.bam -b test/TSCP_coverviewInput.bed -c ../regression_test_data_for_coverview/CoverView_default.json

unittest: install
	@echo ''
	@echo 'Running unit tests'
	@echo ''
	@pytest test/unit

acceptancetest: install
	@echo ''
	@echo 'Running acceptance tests'
	@echo ''
	@pytest test/acceptance

smoketest: install
	@echo ''
	@echo 'Running smoke tests'
	@echo ''
	./test/smoke/check_installation_succeeded.bash

test: pep8 smoketest unittest acceptancetest
	@echo ''
	@echo 'Finished running all tests'
	@echo ''

test_coverage:
	pytest --cov=coverview test
