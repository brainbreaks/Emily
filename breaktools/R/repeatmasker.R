#' @export
repeatmasker_read = function(path, columns=c("repeatmasker_chrom", "repeatmasker_start", "repeatmasker_end", "repeatmasker_strand", "repeatmasker_name", "repeatmasker_class", "repeatmasker_family")) {
  cache_path = paste0("tmp/", digest::digest(paste(path, paste(columns, collapse="")), algo="md5"), "_repeatmasker_df.rda")
  if(file.exists(cache_path)) {
    load(cache_path)
    return(repeatmasker_df %>% dplyr::mutate(repeatmasker_id=1:n()))
  }
  else {
    if(!dir.exists("tmp")) {
      dir.create("tmp", recursive=T)
    }

    repeatmasker_df = readr::read_tsv(path, col_names=names(repeatmasker_cols()$cols), col_types=repeatmasker_cols(), skip=1) %>%
      dplyr::filter(repeatmasker_chrom %in% paste0("chr", c(1:99, "X", "Y"))) %>%
      data.frame()
    repeatmasker_df = repeatmasker_df[columns]


    save(repeatmasker_df, file=cache_path)
    return(repeatmasker_df %>% dplyr::mutate(repeatmasker_id=1:n()))
  }
}

#' @export
repeatmasker_cols = function() {
  readr::cols(
    repeatmasker_bin=readr::col_double(),
    repeatmasker_score=readr::col_double(),
    repeatmasker_mismatches_per_kb=readr::col_double(),
    repeatmasker_deletions_per_kb=readr::col_double(),
    repeatmasker_insertions_per_kb=readr::col_double(),
    repeatmasker_chrom=readr::col_factor(),
    repeatmasker_start=readr::col_double(),
    repeatmasker_end=readr::col_double(),
    repeatmasker_genoLeft=readr::col_double(),
    repeatmasker_strand=readr::col_factor(),
    repeatmasker_name=readr::col_factor(),
    repeatmasker_class=readr::col_factor(),
    repeatmasker_family=readr::col_factor(),
    repeatmasker_repStart=readr::col_character(),
    repeatmasker_repEnd=readr::col_character(),
    repeatmasker_repLeft=readr::col_character(),
    repeatmasker_id=readr::col_factor()
  )
}

