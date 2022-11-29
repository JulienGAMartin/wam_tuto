#!/usr/bin/env Rscript
rm(list = ls())
bookdown::render_book("index.Rmd", "bookdown::gitbook", params = list(echo_sol = TRUE, html_pdf = TRUE))
rm(list = ls())
bookdown::render_book("index.Rmd", "bookdown::pdf_book", params = list(echo_sol = FALSE, html_pdf = FALSE))
rm(list = ls())
bookdown::render_book("index.Rmd", "bookdown::epub_book", params = list(echo_sol = FALSE, html_pdf = TRUE))
