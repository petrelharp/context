PROJECT=context-dependence
BUILDTEX=pdflatex $(PROJECT).tex

all: pruning-schematic.pdf Tmer-dependency.pdf overview.pdf
	$(BUILDTEX)
	biber $(PROJECT)
	$(BUILDTEX)
	$(BUILDTEX)

clean:
	rm -f *.log *.bak *.aux *.bbl *.blg *.idx *.toc *.out *.bcf *.run.xml *~

# save inkscape svg files as .ink.svg and this'll do the right thing
%.svg : %.ink.svg
	inkscape $< --export-plain-svg=$@

%.pdf : %.ink.svg
	inkscape $< --export-pdf=$@

%.svg : %.pdf
	inkscape $< --export-plain-svg=$@

%.png : %.pdf
	convert -density 300 $< -flatten $@

continuous:
	make
	while :; do inotifywait -e modify $(PROJECT).tex $(PROJECT).bib; make; done

diff-to-first-submission.tex : context-dependence.tex
	latexdiff-git --force -r b9f2c8a27be749f32546775019d3c17b85d61b56 $<
	mv context-dependence-diffb9f2c8a27be749f32546775019d3c17b85d61b56.tex $@
