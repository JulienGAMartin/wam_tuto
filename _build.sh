#!/bin/sh

Rscript -e "bookdown::render_book('index.Rmd', 'bookdown::gitbook', params = list(echo_sol = TRUE))"
Rscript -e "bookdown::render_book('index.Rmd', 'bookdown::pdf_book', params = list(echo_sol = FALSE))"
Rscript -e "bookdown::render_book('index.Rmd', 'bookdown:::epub_book')"
