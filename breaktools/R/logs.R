log = function(..., collapse=NULL) {
  msg = paste0(">>> ", paste0(..., collapse=collapse))
  message(msg)
  # cat(file=stderr(), paste0(msg, "\n"))
}