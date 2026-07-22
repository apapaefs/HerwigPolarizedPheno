.PHONY: all clean

all: main.pdf

main.pdf: main.tex references.bib jheppub.sty JHEP.bst
	latexmk -pdf -interaction=nonstopmode -halt-on-error main.tex

clean:
	latexmk -C main.tex

