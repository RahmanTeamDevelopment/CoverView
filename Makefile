check:
	pep8 bamgen/bamgen.py

clean:
	rm -f *.bam
	rm -f *.bai
	rm -f coverage/*.c
	rm -f coverage/*.so
	rm -rf build
	rm -rf dist

libs:
	pip install -e .

test: check
	./run_unit_tests.bash

profile: libs
	time python -m cProfile -s cumulative bin/CoverView.py --input ../Data/NA21144.mapped.ILLUMINA.bwa.GIH.exome.20121211.bam -b chrom20_exons.bed > profile.out
