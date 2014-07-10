PROJECT=context-dependence
BUILDTEX=pdflatex $(PROJECT).tex

all:
	$(BUILDTEX)
	bibtex $(PROJECT)
	$(BUILDTEX)
	$(BUILDTEX)
