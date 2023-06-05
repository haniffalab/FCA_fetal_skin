library(rmd2jupyter)

cat("Enter in-file path: ")
r <- readLines(con = "stdin", n = 1)
if (!(grepl("/", r, fixed = TRUE))) {
    r <- file.path(getwd(), r)
}
stopifnot(endsWith(r, ".md"))
w <- gsub(".md", ".Rmd", r)
r <- readLines(r)
cat('---
title: "Untitled"
output: html_document
---
', file=w)  ## basic YAML header, may be customized
for (i in seq(r)) {cat(r[i], '\n', file=w, append=TRUE)}

rmd2jupyter(w)