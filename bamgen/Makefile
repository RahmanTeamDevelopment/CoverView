
check:
	pep8 bamgen.py

clean:
	rm -f *.bam
	rm -f *.bai

test: check
	python bamgen.py config.json
