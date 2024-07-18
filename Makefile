all: html pdf epub

html: *.Rmd
	Rscript -e "bookdown::render_book('./', 'bookdown::gitbook', clean = FALSE)"
