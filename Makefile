
check:
	pep8 bamgen/bamgen.py

clean:
	rm -f *.bam
	rm -f *.bai

test: check
	./run_unit_tests.bash
