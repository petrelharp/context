.PHONY : clean
FIGS = sequence-context-1.pdf sequence-context-2.pdf sequence-context-3.pdf sequence-context-4.pdf sequence-context-5.pdf sequence-context-6.pdf

hutch-july-2014.pdf : hutch-july-2014.tex $(FIGS)
	pdflatex $<

clean :
	rm -f *.out *.log *.aux *.bbl *.blg *.nav *.snm *.toc

sequence-context-%.pdf : sequence-context.ink.svg
	./export-layers-svg.sh $< layer0 layer$* >$@

%.pdf : %.ink.svg
	inkscape --without-gui --export-pdf=$@ $<
