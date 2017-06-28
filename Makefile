HEADERS=coverview/*.pxd
PYX=coverview/*.pyx
PY=coverview/*.py bamgen/*.py bin/CoverView.py testutils/*.py
PEP8=pep8 --max-line-length=120
SCRIPTS=bin/coverview

pep8:
	${PEP8} ${PY}

clean:
	pip uninstall -y CoverView
	find . -name __pycache__ | xargs rm -rf

cleanAll: clean
	rm -rf env

wheels:
	pip wheel .

install: ${HEADERS} ${PYX} ${PY} ${SCRIPTS}
	./install.sh

profile: install
	time python -m cProfile -s cumulative env/bin/CoverView.py --input ../Data/NA21144.mapped.ILLUMINA.bwa.GIH.exome.20121211.bam -b chrom20_exons.bed > profile.out

regression_test: install pep8
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

test: pep8 smoketest unittest acceptancetest
	@echo ''
	@echo 'Finished running all tests'
	@echo ''

test_coverage:
	pytest --cov=coverview test
