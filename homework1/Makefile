# Makefile for LaTeX compilation using xelatex -> bibtex -> xelatex * 2

# Variables
TEXFILE = main
TEX = xelatex
BIB = bibtex

all: pdf

pdf: $(TEXFILE).pdf

$(TEXFILE).pdf: $(TEXFILE).tex
	# First xelatex pass
	$(TEX) $(TEXFILE).tex
	# Run bibtex
	$(BIB) $(TEXFILE)
	# Second and third xelatex pass
	$(TEX) $(TEXFILE).tex
	$(TEX) $(TEXFILE).tex

clean:
	rm -f *.aux *.bbl *.blg *.log *.out *.toc *.pdf *.run.xml *-blx.bib
.PHONY: all clean
