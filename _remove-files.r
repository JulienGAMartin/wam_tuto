lang <- Sys.getenv("QUARTO_PROFILE")
files <- list.files(lang, pattern = "*.[Rq]md")
file.remove(files)
