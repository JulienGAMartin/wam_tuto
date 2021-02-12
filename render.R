#!/usr/bin/env Rscript
rm(list=ls());bookdown::render_book("index.Rmd", "bookdown::gitbook", params = list(echo_sol = TRUE))
rm(list=ls());bookdown::render_book("index.Rmd", "bookdown::pdf_book", params = list(echo_sol = FALSE))
rm(list=ls());bookdown::render_book("index.Rmd", 'bookdown::epub_book',params = list(echo_sol = FALSE))
bookdown::calibre('docs/wam_tuto.epub', 'mobi')
