
#' Write a matrix or data frame to a .tex table, ready to include in your LaTeX report
#' 
#' @param tab matrix or data frame
#' @param filename Name of the file the matrix or data frame must be written to. Using a .tex extension is advised, suitability for other file extensions was not checked
#' @param caption Information and explanations about the table as a character string.
#' @param reflabel Label of the .tex table for referencing as a character string.
#' @param rows.named Logical, set to TRUE for matrices only. If TRUE, The rownames of the matrix are included in the .tex file
#' @param minipage Logical, set to TRUE if the .tex table is included in a "minipage" environment in LaTeX
#' @return 
#' @examples 
#' 
#' getwd()
#'
#' tab <- matrix(1:9, ncol=3)
#' colnames(tab) <- c("col1", "col2", "col3")
#' rownames(tab) <- c("row1", "row2", "row3")
#' 
#' convertToTexTable(tab, "tab.tex")
#' convertToTexTable(tab, "tabWithCaption.tex", caption="This is a test", reflabel="test_matrix")
#' convertToTexTable(tab, "tabWithRowsNamed.tex", caption="This is a test", reflabel="test_matrix", rows.named=TRUE)
#' convertToTexTable(tab, "tabMinipage.tex", minipage=TRUE)
#' 
convertToTexTable <- function(tab, filename, caption=NULL, reflabel=NULL, 
                              rows.named=FALSE, minipage=FALSE){
  sink(file=filename)
  if (!minipage) {
    cat("\\begin{table}[h]")
    cat("\n")
    cat("\\centering")
    cat("\n")
  }
  cat("\\begin{tabular}")
  if (rows.named) {marg <- "c|"}
  else {marg <- ""}
  cat(paste(c("{||", marg, rep("c", times=ncol(tab)), "||}")))
  cat("\n")
  cat("\\hline\\hline")
  cat("\n")
  column_names <- gsub(pattern="_", replacement="\\_", x=colnames(tab), fixed=TRUE)
  if (rows.named) {
    cat(" & ")
    row_names <- gsub(pattern="_", replacement="\\_", x=rownames(tab), fixed=TRUE)
  }
  cat(paste(column_names, collapse=" & "))
  cat(" \\\\ \n \\hline")
  
  for (i in 1:nrow(tab)) {
    cat("\n")
    row <- gsub(pattern="_", replacement="\\_", x=tab[i, ], fixed=TRUE)
    row <- gsub(pattern="%", replacement="\\%", x=tab[i, ], fixed=TRUE)
    if (rows.named) {cat(paste0(row_names[i], " & "))}
    cat(paste(row, collapse=" & "))
    cat(" \\\\")
  }
  cat(" \n")
  cat("\\hline\\hline \n\\end{tabular} \n")
  if (!is.null(caption)) {cat(paste0("\\caption{", caption, "} \n"))}
  if (!is.null(reflabel)) {cat(paste0("\\label{", reflabel, "} \n"))}
  
  if (!minipage) {
    cat("\\end{table}")
  }
  cat("\n")
  sink()
}
