PROJECT=context-dependence
BUILDTEX=pdflatex $(PROJECT).tex

all: peeling-schematic.pdf Tmer-dependency.pdf
	$(BUILDTEX)
	bibtex $(PROJECT)
	$(BUILDTEX)
	$(BUILDTEX)

clean:
	rm -f *.log *.bak *.aux *.bbl *.blg *.idx *.toc *.out *~

# save inkscape svg files as .ink.svg and this'll do the right thing
%.svg : %.ink.svg
	inkscape $< --export-plain-svg=$@

%.pdf : %.ink.svg
	inkscape $< --export-pdf=$@

%.svg : %.pdf
	inkscape $< --export-plain-svg=$@

%.png : %.pdf
	convert -density 300 $< -flatten $@

