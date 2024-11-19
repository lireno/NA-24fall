TEXFILE = *.tex
TEX = xelatex
BIB = bibtex

all: pdf

test:
	gcc -g 1.c -o 1
	gdb ./1

pdf: $(TEXFILE).pdf

$(TEXFILE).pdf: $(TEXFILE)
	$(TEX) $(TEXFILE)
	$(TEX) $(TEXFILE)
	rm -f *.aux *.bbl *.blg *.log *.out *.toc *.run.xml *-blx.bib *.synctex.gz

clean:
	rm -f *.aux *.bbl *.blg *.log *.out *.toc *.pdf *.run.xml *-blx.bib *.synctex.gz

.PHONY: all clean
