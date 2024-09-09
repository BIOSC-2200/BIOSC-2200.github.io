exs <- c("GE-HWE/index.qmd",
         "GE-ME/index.qmd",
         "QT/index.qmd",
         "QT-Case-Study/index.qmd",
         "TGI/index.qmd",
         "TGI-Case-Study/index.qmd")

for (ex in exs) {
  syscall <- paste0("quarto render ", ex, " --to html")
  system(syscall)
}
