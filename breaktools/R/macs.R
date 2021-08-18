#' @export
macs2 = function(name, sample, effective_size, control=NULL, maxgap=NULL, qvalue=0.01, extsize=2000, slocal=50000, llocal=10000000, output_dir="data/macs2") {
  bed_sample = paste("-t", sample)
  bed_control = ifelse(is.null(control), "", paste("-c", control))
  maxgap = ifelse(is.null(maxgap), "", paste("--max-gap", sprintf("%0.0f", maxgap)))

  cmd = paste("macs2 callpeak ", bed_sample, bed_control, "--seed 123 -f BED --keep-dup all --nomodel --bdg --trackline", maxgap, "-g", effective_size, "-n", name, "--outdir", output_dir, "--slocal", sprintf("%0.0f", slocal), "--extsize", sprintf("%0.0f", extsize), "-q", qvalue, "--llocal", sprintf("%0.0f", llocal))
  # cmd = paste0("macs2 callpeak {bed_sample} {bed_control} --seed 123 {maxgap} -f BED -g {effsize} --keep-dup all -n {name} --outdir {output_dir} --nomodel --slocal {slocal} --extsize {extsize} -q {qvalue} --llocal {llocal} --bdg --trackline", bed_sample=bed_sample, bed_control=bed_control, name=name, output_dir=output_dir, extsize=extsize, qvalue=qvalue, maxgap=maxgap, llocal=sprintf("%0.0f", llocal), slocal=sprintf("%0.0f", slocal), effsize=effective_size)
  log(cmd)
  output = system(paste(cmd, " 2>&1"), intern = T)
  output = paste0(output, collapse="\n")
  log(output)

  readr::read_tsv(paste0(output_dir, "/", name, "_peaks.xls"), comment="#", col_names=names(macs_cols()$cols), col_types=macs_cols()) %>%
    dplyr::slice(-1) %>%
    dplyr::select(-macs_comment)
}
