check:
	pep8 bamgen/bamgen.py

clean:
	rm -f *.bam
	rm -f *.bai
	rm -f coverage/*.c
	rm -f coverage/*.so
	rm -rf coverage/build

libs:
	cd coverage; python setup.py build_ext --inplace

test: check
	./run_unit_tests.bash

profile:
	time python -m cProfile -s cumulative CoverView.py --input ../Data/NA21144.mapped.ILLUMINA.bwa.GIH.exome.20121211.bam -b chrom20_exons.bed > profile.out
