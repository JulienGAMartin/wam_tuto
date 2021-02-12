#!/usr/bin/env Rscript
rm(list=ls());bookdown::render_book("index.Rmd", "bookdown::gitbook", params = list(echo_sol = TRUE))
rm(list=ls());bookdown::render_book("index.Rmd", "bookdown::pdf_book", params = list(echo_sol = FALSE))
rm(list=ls());bookdown::render_book("index.Rmd", 'bookdown::epub_book')
bookdown::calibre('docs/Labo_BIO4558.epub', 'mobi')


#bookdown::preview_chapter("06-anova_mult.Rmd", output_format = 'bookdown::html_chapters', output_dir = "docs/preview")
bookdown::preview_chapter("index.Rmd", output_dir = "docs/preview")
bookdown::preview_chapter("01-intro.Rmd", output_dir = "docs/preview")
bookdown::preview_chapter("02-Gpower.Rmd", output_dir = "docs/preview")
bookdown::preview_chapter("03-reg_lin.Rmd", output_dir = "docs/preview")
bookdown::preview_chapter("04-t_test.Rmd", output_dir = "docs/preview")
bookdown::preview_chapter("05-anova.Rmd", output_dir = "docs/preview")
bookdown::preview_chapter("06-anova_mult.Rmd", output_dir = "docs/preview")
bookdown::preview_chapter("07-reg_mult.Rmd", output_dir = "docs/preview")
bookdown::preview_chapter("08-ancova_glm.Rmd", output_dir = "docs/preview")
bookdown::preview_chapter("09-model_freq.Rmd", output_dir = "docs/preview")
