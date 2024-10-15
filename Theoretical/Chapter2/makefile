# Makefile for LaTeX compilation using xelatex -> bibtex -> xelatex * 2

# Variables
TEXFILE = *.tex
TEX = xelatex
BIB = bibtex

all: pdf

pdf: $(TEXFILE).pdf

$(TEXFILE).pdf: $(TEXFILE)
	# First xelatex pass
	$(TEX) $(TEXFILE)
	# Run bibtex
	$(BIB) $(TEXFILE)
	# Second and third xelatex pass
	$(TEX) $(TEXFILE)
	$(TEX) $(TEXFILE)

clean:
	rm -f *.aux *.bbl *.blg *.log *.out *.toc *.pdf *.run.xml *-blx.bib
.PHONY: all clean
