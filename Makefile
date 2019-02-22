ALL_TEX := $(shell find . -name \*.tex)

ALL_SVG_FIG := $(shell find . -name \*.svg)
#ALL_PDF_FIG := $(ALL_SVG_FIG:.svg=.pdf)
ALL_PDF_FIG := $(join $(dir $(ALL_SVG_FIG)),  $(addprefix ., $(notdir $(ALL_SVG_FIG:.svg=.pdf))))

LAST_SUBMISSION := 6a1b120db
	

all: thesis.pdf

thesis.pdf: $(ALL_TEX) $(ALL_PDF_FIG) Makefile thesis.bib
	@echo "Generating pdf"
	latexmk -pdf thesis.tex


thesis-changes.pdf: thesis.tex $(ALL_TEX)
	latexdiff-vc -r $(LAST_SUBMISSION) $<  --pdf  --force   
	mv thesis-diff$(LAST_SUBMISSION).pdf $@ 
	rm thesis-diff$(LAST_SUBMISSION).tex 


.%.pdf: %.svg
	@echo converting $<  to pdf
	inkscape -f $< -A $@
clean:
	latexmk -C
	rm *.{acn,tdo,ist,gl,aux,log,glo,glsdefs} -f
	rm -f $(ALL_PDF_FIG)
