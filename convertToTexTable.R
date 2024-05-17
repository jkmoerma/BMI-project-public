convertToTexTable <- function(tab, filename, caption=NULL, reflabel=NULL, rows.named=FALSE){
  sink(file=filename)
  cat("\\begin{table}[h]")
  cat("\n")
  cat("\\centering")
  cat("\n")
  cat("\\begin{tabular}")
  if (rows.named) {marg <- "c|"}
  else {marg <- ""}
  cat(paste(c("{||", marg, rep("c", times=ncol(tab)), "||}")))
  cat("\n")
  cat("\\hline\\hline")
  cat("\n")
  column_names <- gsub(pattern="_", replacement="\\_", x=colnames(tab), fixed=TRUE)
  if (rows.named) {cat(" & ")}
  cat(paste(column_names, collapse=" & "))
  cat(" \\\\ \n \\hline")
  
  for (i in 1:nrow(tab)) {
    cat("\n")
    row <- gsub(pattern="_", replacement="\\_", x=tab[i, ], fixed=TRUE)
    if (rows.named) {cat(paste0(rownames(tab)[i], " & "))}
    cat(paste(row, collapse=" & "))
    cat(" \\\\")
  }
  cat(" \n")
  cat("\\hline\\hline \n\\end{tabular} \n")
  if (!is.null(caption)) {cat(paste0("\\caption{", caption, "} \n"))}
  if (!is.null(reflabel)) {cat(paste0("\\label{", reflabel, "} \n"))}
  
  cat("\\end{table}")
  cat("\n")
  sink()
}

if (FALSE) {
  # a check of functionality for this function
  tab <- matrix(1:9, ncol=3)
  convertToTexTable(tab, "tab1.tex")
  convertToTexTable(tab, "tab2.tex", caption="This is a test", reflabel="test_matrix")
}