HEADERS=coverview/*.pxd
PYX=coverview/*.pyx
PY=coverview/*.py bamgen/*.py bin/CoverView.py testutils/*.py
FLAKE8=flake8 --max-line-length=120
SCRIPTS=bin/coverview

flake8:
	${FLAKE8} ${PY} test

clean:
	pip uninstall -y CoverView
	find . -name __pycache__ | xargs rm -rf
	cd docs/sphinx; make clean

cleanAll: clean cleanDocs
	rm -rf env

cleanDocs:
	cd docs/sphinx; make clean

.PHONY:
pdfdocs: docs/sphinx/*.rst
	cd docs/sphinx; make latexpdf

.PHONY:
docs: docs/sphinx/*.rst
	cd docs/sphinx; make html;
	cp -rf docs/sphinx/_build/html/* docs/

wheels:
	pip wheel .

env/bin/coverview: ${HEADERS} ${PYX} ${PY} ${SCRIPTS}
	./install.sh
	touch env/bin/coverview

.PHONY:
install: env/bin/coverview ;

profile: install
	time python -m cProfile -s cumulative env/bin/CoverView.py --input ../Data/NA21144.mapped.ILLUMINA.bwa.GIH.exome.20121211.bam -b chrom20_exons.bed > profile.out

regression_test: install
	coverview --input ../regression_test_data_for_coverview/16768_sorted_picard.bam -b ../regression_test_data_for_coverview/TSCP_coverviewInput.bed -c ../regression_test_data_for_coverview/CoverView_default.json

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

test: flake8 smoketest unittest acceptancetest
	@echo ''
	@echo 'Finished running all tests'
	@echo ''

test_coverage:
	pytest --cov=coverview test
