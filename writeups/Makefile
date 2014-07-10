PROJECT=context-dependence
BUILDTEX=pdflatex $(PROJECT).tex

all:
	$(BUILDTEX)
	bibtex $(PROJECT)
	$(BUILDTEX)
	$(BUILDTEX)

clean:
	rm -f *.log *.bak *.aux *.bbl *.blg *.idx *.toc *.out *~
