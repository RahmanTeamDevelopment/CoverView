check:
	pep8 bamgen/bamgen.py

clean:
	rm -f *.bam
	rm -f *.bai
	rm -f coverview/*.c
	rm -f coverview/*.so
	rm -f coverview/*.pyc
	rm -rf build
	rm -rf dist
	rm -rf CoverView.egg-info

wheels:
	pip wheel .

libs:
	pip install .

profile: libs
	time python -m cProfile -s cumulative bin/CoverView.py --input ../Data/NA21144.mapped.ILLUMINA.bwa.GIH.exome.20121211.bam -b chrom20_exons.bed > profile.out

test: check libs
	python bin/CoverView.py --input test/16768_sorted_picard.bam -b test/TSCP_coverviewInput.bed -c test/CoverView_default.json
