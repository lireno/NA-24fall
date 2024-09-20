# Variables
TEXFILE = answer
TEX = pdflatex
BIB = bibtex

all: pdf

pdf: $(TEXFILE).pdf

$(TEXFILE).pdf: $(TEXFILE).tex
	# First pdflatex pass
	$(TEX) $(TEXFILE).tex
	# Run bibtex
	$(BIB) $(TEXFILE)
	# Second and third pdflatex pass
	$(TEX) $(TEXFILE).tex
	$(TEX) $(TEXFILE).tex

clean:
	rm -f *.aux *.bbl *.blg *.log *.out *.toc *.pdf *.run.xml *-blx.bib *.synctex.gz
.PHONY: all clean
