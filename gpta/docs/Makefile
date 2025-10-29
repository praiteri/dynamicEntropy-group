.PHONY: all
all: clean
	-pdflatex manual.tex
	-pdflatex manual.tex
	make clean

.PHONY: clean
clean:
	-rm -f *.aux
	-rm -f *.log
	-rm -f *.idx
	-rm -f *.ind
	-rm -f *.ilg
	-rm -f *.toc
	-rm -f *.out

.PHONY: clean-all
clean-all: clean
	-rm -f manual.pdf

