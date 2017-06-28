all: User_Guide.pdf

User_Guide.pdf: User_Guide.tex User_Guide.bib
	pdflatex User_Guide.tex
	bibtex User_Guide
	pdflatex User_Guide.tex
	pdflatex User_Guide.tex

clean:
	rm -rf *.out *.log *.aux *.bbl
	rm -rf User_Guide.pdf
