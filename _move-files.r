lang <- Sys.getenv("QUARTO_PROFILE")
files <- list.files(
  lang,
  pattern = "*.[Rq]md",
  full.names = TRUE
  )
file.copy(files, ".", overwrite = TRUE)
