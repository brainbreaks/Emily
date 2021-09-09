setwd("/mnt/e/Workspace/Emily")
# library(randomForest)
library(dplyr)
library(reshape2)
library(readr)
library(GenomicRanges)
library(GenomicFeatures)
library(ggplot2)
library(tidyr)
library(baseline)
library(smoother)
library(dbscan)
library(ggpmisc)
library(rtracklayer)
library(ggbeeswarm)
library(Rtsne)
devtools::load_all('breaktools')

liftOverRanges = function(ranges, chain) {
  ranges$ranges_id = 1:length(ranges)
  as.data.frame(unlist(rtracklayer::liftOver(ranges, chain))) %>%
    dplyr::group_by(ranges_id) %>%
    dplyr::mutate(start=min(start), end=max(end)) %>%
    dplyr::ungroup() %>%
    dplyr::distinct(ranges_id, .keep_all=T) %>%
    dplyr::select(-ranges_id)
}

repliseq_read = function(path) {
  repliseq_df = as.data.frame(t(readr::read_tsv(path, col_names=F))) %>%
    dplyr::rename(repliseq_chrom="V1", repliseq_start="V2", repliseq_end="V3") %>%
    dplyr::mutate(repliseq_start=as.numeric(repliseq_start), repliseq_end=as.numeric(repliseq_end)) %>%
    reshape2::melt(id.vars=c("repliseq_chrom", "repliseq_start", "repliseq_end"), variable.name="repliseq_fraction", value.name="repliseq_value") %>%
    dplyr::mutate(repliseq_fraction=as.numeric(gsub("V", "", repliseq_fraction))-3, repliseq_value=as.numeric(repliseq_value))
  repliseq_df.keep = repliseq_df %>%
    dplyr::arrange(repliseq_start) %>%
    dplyr::group_by(repliseq_chrom, repliseq_start, repliseq_end) %>%
    dplyr::summarize(repliseq_isna=sum(is.na(repliseq_value))>10) %>%
    dplyr::group_by(repliseq_chrom) %>%
    dplyr::mutate(first_value_position=min(which(!repliseq_isna)), last_value_position=max(which(!repliseq_isna)), keep_value=dplyr::between(1:n(), first_value_position, last_value_position)) %>%
    dplyr::filter(keep_value) %>%
    dplyr::ungroup() %>%
    dplyr::select(repliseq_chrom, repliseq_start, repliseq_end)
  repliseq_df.f = repliseq_df %>%
    dplyr::inner_join(repliseq_df.keep, by=c("repliseq_chrom", "repliseq_start", "repliseq_end"))

  repliseq_df.f
}

repliseq_summarize = function(repliseq_df, window=5) {
  th.repliseq_value_norm = 0.1
  repliseq_df = repliseq_df %>%
    dplyr::group_by(repliseq_chrom, repliseq_start, repliseq_end) %>%
    dplyr::mutate(repliseq_value_norm=((repliseq_value)/max(repliseq_value))) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(repliseq_chrom, repliseq_start, repliseq_end) %>%
    do((function(z) {
      zz<<-z

      i.max = z$repliseq_fraction[which.max(z$repliseq_value_norm)]
      lb = which(z$repliseq_value_norm < th.repliseq_value_norm & z$repliseq_fraction < i.max)
      lb = z$repliseq_fraction[ifelse(length(lb), lb[length(lb)]+1, 1)]
      ub = which(z$repliseq_value_norm < th.repliseq_value_norm & z$repliseq_fraction > i.max)
      ub = z$repliseq_fraction[ifelse(length(ub), ub[1]-1, nrow(z))]

      z$repliseq_value_in_scope = dplyr::between(z$repliseq_fraction, lb, ub)
      z
    })(.)) %>%
    dplyr::ungroup()

  repliseq_time_df = repliseq_df %>%
    dplyr::mutate(repliseqTime_chrom=repliseq_chrom, repliseqTime_start=repliseq_start, repliseqTime_end=repliseq_end) %>%
    dplyr::group_by(repliseqTime_chrom, repliseqTime_start, repliseqTime_end) %>%
    dplyr::mutate(repliseqTime_avg=weighted.mean(repliseq_fraction[repliseq_value_in_scope], repliseq_value_norm[repliseq_value_in_scope], na.rm=T)) %>%
    dplyr::mutate(lb=dplyr::between(repliseq_fraction, floor(repliseqTime_avg[1]-4), ceiling(repliseqTime_avg[1]-2)), repliseqTime_min=ifelse(any(lb), weighted.mean(repliseq_fraction[lb], repliseq_value_norm[lb]), NA_real_)) %>%
    dplyr::mutate(ub=dplyr::between(repliseq_fraction, floor(repliseqTime_avg[1]+2), ceiling(repliseqTime_avg[1]+4)), repliseqTime_max=ifelse(any(ub), weighted.mean(repliseq_fraction[ub], repliseq_value_norm[ub]), NA_real_)) %>%
    dplyr::summarize(repliseqTime_avg=repliseqTime_avg[1], repliseqTime_min=repliseqTime_min[1], repliseqTime_max=repliseqTime_max[1], repliseqTime_lb=min(which(repliseq_value_in_scope)), repliseqTime_ub=max(which(repliseq_value_in_scope))) %>%
    dplyr::group_by(repliseqTime_chrom) %>%
    dplyr::mutate(repliseqTime_avg=smoother::smth.gaussian(repliseqTime_avg, window=window), repliseqTime_min=smoother::smth.gaussian(repliseqTime_min, window=window), repliseqTime_max=smoother::smth.gaussian(repliseqTime_max, window=window)) %>%
    dplyr::mutate(repliseqTime_avg=zoo::na.fill(repliseqTime_avg, "extend"), repliseqTime_min=zoo::na.fill(repliseqTime_min, "extend"), repliseqTime_max=zoo::na.fill(repliseqTime_max, "extend")) %>%
    dplyr::ungroup()

  repliseq_time_df
}

repliseq_preprocess = function() {
  # Load repliseq data
  repliseq_ESC_df = repliseq_read("data/zhao_bmc_repliseq_2020/GSE137764_mESC_Gaussiansmooth_scaled_autosome.mat") %>% dplyr::mutate(repliseq_celltype="esc")
  repliseq_NPC_df = repliseq_read("data/zhao_bmc_repliseq_2020/GSE137764_mNPC_Gaussiansmooth_scaled_autosome.mat") %>% dplyr::mutate(repliseq_celltype="npc")
  repliseq_ESC2NPC_df = dplyr::bind_rows(repliseq_ESC_df, repliseq_NPC_df)
  readr::write_tsv(repliseq_ESC2NPC_df,"data/zhao_bmc_repliseq_2020/preprocessed/repliseq.tsv")

  repliseqTime_ESC_df = repliseq_summarize(repliseq_ESC_df) %>% dplyr::mutate(repliseqTime_celltype="esc")
  repliseqTime_NPC_df = repliseq_summarize(repliseq_NPC_df) %>% dplyr::mutate(repliseqTime_celltype="npc")
  repliseqTime_ESC2NPC_df = dplyr::bind_rows(repliseqTime_ESC_df, repliseqTime_NPC_df)
  readr::write_tsv(repliseqTime_ESC2NPC_df,"data/zhao_bmc_repliseq_2020/preprocessed/repliseqTime.tsv")


  repliseq_merged_df = repliseq_ESC2NPC_df %>%
    dplyr::group_by(repliseq_chrom, repliseq_start, repliseq_end, repliseq_fraction) %>%
    dplyr::summarize(repliseq_value=mean(repliseq_value, na.rm=T)) %>%
    dplyr::arrange(repliseq_chrom, repliseq_start, repliseq_end, repliseq_fraction)
  repliseqTime_merged_df = repliseq_summarize(repliseq_merged_df) %>% dplyr::mutate(repliseqTime_celltype="esc/npc")
  readr::write_tsv(repliseqTime_merged_df,"data/zhao_bmc_repliseq_2020/preprocessed/repliseqTime_merged.tsv")
}



main = function() {
  m9_m10 = rtracklayer::import.chain("data/mm9/mm9ToMm10.over.chain")

  # chromosomes_map_df = readr::read_tsv("data/chromosome_synonyms.tsv")

  # Load genome info
  chromsizes_cols = readr::cols(seqnames=col_character(), seqlengths=col_double())
  genome_info = with(readr::read_tsv("data/mm10/annotation/mm10.chrom.sizes", col_names=names(chromsizes_cols$cols), col_types=chromsizes_cols),
           GenomeInfoDb::Seqinfo(seqnames, seqlengths, isCircular=rep(F, length(seqnames)), genome=rep("mm10", length(seqnames))))
  genome_info = genome_info[paste0("chr", c(1:19, "X", "Y"))]

  genome_txdb = GenomicFeatures::makeTxDbFromGFF('data/mm10/annotation/mm10.refGene.gtf.gz', format="gtf")
  genes_df = as.data.frame(GenomicFeatures::genes(genome_txdb)) %>%
    dplyr::mutate(gene_chrom=as.character(seqnames), gene_start=start, gene_end=end, gene_strand=strand, gene_length=gene_end-gene_start)
  genes_ranges = GenomicRanges::makeGRangesFromDataFrame(genes_df, keep.extra.columns=T)

  genes_reduced_ranges = genes_ranges
  strand(genes_reduced_ranges) = "*"
  genes_reduced_ranges = GenomicRanges::reduce(genes_reduced_ranges)
  genes_reduced_ranges$gene_cluster = 1:length(genes_reduced_ranges)
  genes_df = as.data.frame(IRanges::mergeByOverlaps(genes_ranges, genes_reduced_ranges)) %>%
    dplyr::select(-dplyr::matches("_ranges\\.")) %>%
    dplyr::group_by(gene_chrom, gene_cluster) %>%
    dplyr::summarize(gene_cluster_size=n(), gene_dominated_length=max(gene_end-gene_start), gene_id=paste0(gene_id, collapse=","), gene_start=min(gene_start), gene_end=max(gene_end), gene_strand="*", gene_length=gene_end-gene_start) %>%
    dplyr::select(-gene_dominated_length)

  chromsizes_cols = readr::cols(seqnames=col_character(), seqlengths=col_double())
  genome_info_mm9 = with(readr::read_tsv("data/mm9/annotation/mm9.chrom.sizes", col_names=names(chromsizes_cols$cols), col_types=chromsizes_cols),
             GenomeInfoDb::Seqinfo(seqnames, seqlengths, isCircular=rep(F, length(seqnames)), genome=rep("mm9", length(seqnames))))
  genome_info_mm9 = genome_info_mm9[paste0("chr", c(1:19, "X", "Y"))]

  #
  # Read TLX
  #
  repeatmasker_df = repeatmasker_read("data/mm10/annotation/ucsc_repeatmasker.tsv")
  samples_df = readr::read_tsv("data/tlx_samples.tsv")
  tlx_df = tlx_read_many(samples_df)
  tlx_df = tlx_remove_rand_chromosomes(tlx_df)
  tlx_df = tlx_mark_bait_chromosome(tlx_df)
  tlx_df = tlx_mark_bait_junctions(tlx_df, 1.5e6)
  tlx_df = tlx_mark_repeats(tlx_df, repeatmasker_df)
  tlx_df = tlx_df %>%
    dplyr::filter(tlx_is_bait_chromosome & !tlx_is_bait_junction) %>%
    dplyr::select(-Seq) %>%
    dplyr::mutate(tlx_id=1:n()) %>%
    dplyr::ungroup()
  tlx_ranges = GenomicRanges::makeGRangesFromDataFrame(tlx_df %>% dplyr::mutate(seqnames=Rname, start=Junction, end=Junction), keep.extra.columns=T, ignore.strand=T)
  libsize_df = tlx_df %>%
    dplyr::group_by(tlx_group, tlx_group_i) %>%
    dplyr::summarize(sample_size=sum(!tlx_control), control_size=sum(tlx_control))

  #
  # Load repliseq
  #
  repliseq_df = readr::read_tsv("data/zhao_bmc_repliseq_2020/preprocessed/repliseq.tsv") %>% dplyr::filter(repliseq_celltype=="npc")
  repliseqTime_df = readr::read_tsv("data/zhao_bmc_repliseq_2020/preprocessed/repliseqTime.tsv") %>% dplyr::filter(repliseqTime_celltype=="npc") %>%
    dplyr::mutate(repliseqTime_id=1:n()) %>%
    dplyr::group_by(repliseqTime_celltype, repliseqTime_chrom) %>%
    dplyr::mutate(repliseqTime_smooth=smoother::smth.gaussian(repliseqTime_avg, window=50))

  # Calculate IZ
  starts = unique(repliseqTime_df$repliseqTime_start)
  repliseqIZ_df = repliseqTime_df %>%
    dplyr::rename(iz_chrom="repliseqTime_chrom", iz_celltype="repliseqTime_celltype") %>%
    dplyr::group_by(iz_chrom, iz_celltype) %>%
    dplyr::do((function(z){
      zz<<-z
      z.valleys = which(ggpmisc:::find_peaks(17-z$repliseqTime_smooth, ignore_threshold=4/16, span=9, strict=F))
      z.peaks = which(ggpmisc:::find_peaks(z$repliseqTime_smooth, ignore_threshold=4/16, span=9, strict=F))
      data.frame(iz_start=c(z.peaks, z.valleys), iz_type=rep(c("peak", "valley"), c(length(z.peaks), length(z.valleys)))) %>%
        dplyr::mutate(iz_time=z$repliseqTime_smooth[iz_start], iz_start=z$repliseqTime_start[iz_start])
    })(.)) %>%
    dplyr::arrange(iz_celltype, iz_chrom, iz_start) %>%
    dplyr::group_by(iz_chrom, iz_celltype) %>%
    dplyr::mutate(iz_type_group=rep(seq_along(rle(iz_type)$values), rle(iz_type)$lengths)) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(iz_celltype, iz_chrom, iz_type_group) %>%
    dplyr::do((function(z){
      if(nrow(z)==1) return(z)
      zz<<-z

      repliseqTime_df %>%
        dplyr::filter(repliseqTime_celltype==z$iz_celltype[1] & repliseqTime_chrom==z$iz_chrom[1] & dplyr::between(repliseqTime_start, min(z$iz_start), max(z$iz_start))) %>%
        dplyr::arrange(repliseqTime_smooth) %>%
        dplyr::slice(1) %>%
        dplyr::mutate(iz_type=z$iz_type[1], iz_type_group=z$iz_type_group[1]) %>%
        dplyr::select(iz_chrom=repliseqTime_chrom, iz_celltype=repliseqTime_celltype, iz_start=repliseqTime_start, iz_type, iz_time=repliseqTime_smooth, iz_type_group)
    })(.)) %>%
    dplyr::ungroup()

    # dplyr::arrange(iz_start) %>%
    # dplyr::group_by(iz_chrom, iz_type) %>%
    # dplyr::mutate(iz_cluster=dbscan::dbscan(as.matrix(iz_start), eps=50000, minPts=1)$cluster) %>%
    # dplyr::group_by(iz_chrom, iz_type, iz_cluster) %>%
    # dplyr::mutate(iz_dist=ifelse(iz_type=="peak", 17-iz_time, iz_time), iz_start=iz_start[which.min(iz_dist)]) %>%
    # # dplyr::mutate(iz_start=starts[which.min(abs(starts-mean(iz_start)))]) %>%
    # dplyr::ungroup() %>%
    # dplyr::distinct(iz_chrom, iz_start, iz_type, iz_celltype, .keep_all=T)

  x = repliseqTime_df %>% dplyr::filter(repliseqTime_chrom=="chr6" & dplyr::between(repliseqTime_start, 20e6, 120e6))
  y = repliseqIZ_df %>% dplyr::filter(iz_chrom=="chr6" & dplyr::between(iz_start, 20e6, 120e6)) %>% data.frame()
  ggplot(x) +
    geom_line(aes(x=repliseqTime_start, y=repliseqTime_avg)) +
    geom_line(aes(x=repliseqTime_start, y=repliseqTime_smooth, color="smooth")) +
    geom_vline(aes(xintercept=iz_start, color=iz_type), data=y)


  repliseqTimeRegions_df = repliseqTime_df %>%
    dplyr::mutate(repliseqTime_type=factor(repliseqTime_type), repliseqTime_chrom=factor(repliseqTime_chrom)) %>%
    dplyr::filter(c(T, diff(as.numeric(repliseqTime_type))!=0 | diff(as.numeric(repliseqTime_chrom))!=0)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(repliseqTime_type=as.character(repliseqTime_type), repliseqTime_chrom=as.character(repliseqTime_chrom)) %>%
    dplyr::select(repliseqTime_celltype, repliseqTime_chrom, repliseqTime_start, repliseqTime_end, repliseqTime_type)
  repliseqTimeRegions_df$repliseqTime_type = c(repliseqTimeRegions_df$repliseqTime_type[1], zoo::rollapply(repliseqTimeRegions_df$repliseqTime_type, width=3, align="center", FUN=function(z) ifelse(all(z==c("Late", "Middle", "Late")), "Late", z[1])), rev(repliseqTimeRegions_df$repliseqTime_type)[1])
  repliseqTimeRegions_df$repliseqTime_type = c(repliseqTimeRegions_df$repliseqTime_type[1], zoo::rollapply(repliseqTimeRegions_df$repliseqTime_type, width=3, align="center", FUN=function(z) ifelse(all(z==c("Early", "Middle", "Early")), "Early", z[1])), rev(repliseqTimeRegions_df$repliseqTime_type)[1])
  repliseqTimeRegions_df = repliseqTimeRegions_df %>%
    dplyr::mutate(repliseqTime_type=factor(repliseqTime_type), repliseqTime_chrom=factor(repliseqTime_chrom)) %>%
    dplyr::filter(c(T, diff(as.numeric(repliseqTime_type))!=0 | diff(as.numeric(repliseqTime_chrom))!=0)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(repliseqTime_type=as.character(repliseqTime_type), repliseqTime_chrom=as.character(repliseqTime_chrom))


    dplyr::inner_join(as.data.frame(genome_info) %>% tibble::rownames_to_column("repliseqTime_chrom"), by=c("repliseqTime_chrom")) %>%
    dplyr::group_by(repliseqTime_celltype, repliseqTime_chrom) %>%
    dplyr::mutate(repliseqTime_end=c(repliseqTime_start[2:n()], seqlengths[1])) %>%
    dplyr::mutate(repliseqTime_nextstart=c(repliseqTime_start[2:n()], NA_integer_), repliseqTime_nextend=c(repliseqTime_end[2:n()], NA_integer_)) %>%
    dplyr::mutate(repliseqTime_prevstart=c(NA_integer_, repliseqTime_start[1:(n()-1)]), repliseqTime_prevend=c(NA_integer_, repliseqTime_end[1:(n()-1)])) %>%
    dplyr::mutate(repliseqTime_prevtype=c(NA_character_, repliseqTime_type[1:(n()-1)]), repliseqTime_nexttype=c(repliseqTime_type[2:n()], NA_character_)) %>%
    dplyr::filter(repliseqTime_type %in% c("Early", "Late"))

  repliseqTimeRegions_df %>%
    dplyr::filter(repliseqTime_type=="Early") %>%
    dplyr::inner_join(repliseqTimeRegions_df %>% dplyr::filter(repliseqTime_type=="Middle"), by=c("repliseqTime_celltype", "repliseqTime_chrom", "repliseqTime_start"="repliseqTime_end"))

  table(prev=repliseqTimeRegions_df$repliseqTime_prevtype, type=repliseqTimeRegions_df$repliseqTime_type)

  repliseqTime_ranges = GenomicRanges::makeGRangesFromDataFrame(repliseqTime_df %>% dplyr::mutate(seqnames=repliseqTime_chrom, start=repliseqTime_start, end=repliseqTime_end), keep.extra.columns=T)


  table(repliseqTime_df$repliseqTime_type)


  #
  # Join TLX and repliseq
  #

  repliseqTime2breaks = as.data.frame(IRanges::mergeByOverlaps(repliseqTime_ranges, tlx_ranges))
repliseqTime2breaks$repliseqTime_ty

  #
  # calculate treatment/control fraction for each pair
  #
  repliseqTime2breaks_sum = repliseqTime2breaks %>%
    dplyr::group_by(repliseqTime_chrom, repliseqTime_start, repliseqTime_end, repliseqTime_avg, tlx_group, tlx_group_i) %>%
    dplyr::summarize(
      n_control=sum(!is.na(tlx_control) & tlx_control),
      n_sample=sum(!is.na(tlx_control) & !tlx_control),
      frac=n_sample/n_control) %>%
    dplyr::inner_join(libsize_df, by=c("tlx_group", "tlx_group_i")) %>%
    dplyr::mutate(frac_norm=frac*(control_size/sample_size)) %>%
    dplyr::mutate(frac_lognorm=log2(ifelse(frac_norm==0, min(frac_norm[is.finite(frac_norm)]), frac_norm))) %>%
    dplyr::group_by(repliseqTime_chrom, repliseqTime_start, repliseqTime_end, repliseqTime_avg, tlx_group) %>%
    dplyr::summarize(tlx_group_i=1, frac_lognorm=mean(frac_lognorm)) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(tlx_group, tlx_group_i, repliseqTime_chrom, repliseqTime_start, repliseqTime_end) %>%
    dplyr::group_by(repliseqTime_chrom, tlx_group)

  table(repliseqTime2breaks_sum$repliseqTime_chrom)
table(repliseqTime2breaks_sum$tlx_c)

  repliseqTime2breaks_sum_ranges = GenomicRanges::makeGRangesFromDataFrame(repliseqTime2breaks_sum %>% dplyr::mutate(seqnames=repliseqTime_chrom, start=repliseqTime_start, end=repliseqTime_end), keep.extra.columns=T)

  repliseqTime2breaks_sum %>%
    dplyr::group_by(repliseqTime_chrom, repliseqTime_start)


  #
  # Calculate MACS2
  #
  macs_df = tlx_macs2(tlx_df %>% dplyr::filter(tlx_group=="APH"), effective_size=1.87e9, extsize=1e4, maxgap=4e4, exttype="along", qvalue=0.05, pileup=5, slocal=1e4, llocal=1e7, exclude_bait_region=T, exclude_repeats=T) %>%
    dplyr::mutate(macs_length=as.numeric(macs_length))
  macs_ranges = GenomicRanges::makeGRangesFromDataFrame(macs_df %>% dplyr::mutate(seqnames=macs_chrom, start=macs_start-1e5, end=macs_end+1e5), keep.extra.columns=T)


  # Create random regions
  random_lengths = macs_df %>% dplyr::select(macs_chrom, macs_length)
  random_df = random_lengths %>%
    dplyr::inner_join(as.data.frame(genome_info) %>% tibble::rownames_to_column("macs_chrom"), by=c("macs_chrom")) %>%
    dplyr::group_by(macs_chrom) %>%
    dplyr::do((function(z){
      zz<<-z
      n = 5
      z.start = sample(z$seqlengths[1], nrow(z)*n)

      data.frame(macs_chrom=z$macs_chrom[1], seqlengths=z$seqlengths[1], start=z.start, macs_length=rep(z$macs_length, each=n)) %>%
        dplyr::mutate(macs_start=ifelse(start+macs_length>=seqlengths, start-macs_length, start), macs_end=start+macs_length, macs_group="Random positions") %>%
        dplyr::mutate(macs_cluster=paste0(macs_chrom, "_", 1:n()))
    })(.)) %>%
    dplyr::ungroup()

  hits_df = dplyr::bind_rows(random_df, macs_df)
  hits_ranges = GenomicRanges::makeGRangesFromDataFrame(hits_df %>% dplyr::select(-seqlengths) %>% dplyr::mutate(seqnames=macs_chrom, start=macs_start-1e5, end=macs_end+1e5), keep.extra.columns=T)

  x = as.data.frame(IRanges::mergeByOverlaps(hits_ranges, repliseqTime2breaks_sum_ranges)) %>%
    dplyr::filter(is.finite(frac_norm)) %>%
    dplyr::mutate(macs_name=paste0(macs_group, " ", macs_chrom, ":", macs_start, "-", macs_end))
  ggplot(x) +
    geom_line(aes(x=repliseqTime_start, y=frac_lognorm, color=tlx_group, group=paste(tlx_group, tlx_group_i))) +
    geom_line(aes(x=repliseqTime_start, y=repliseqTime_avg)) +
    coord_cartesian(ylim=c(-20,16)) +
    facet_wrap(~macs_name, scales="free_x")
  ggplot(x) +
    geom_boxplot(aes(x=tlx_group, y=repliseqTime_avg, fill=macs_group))
  x_sum = x %>%
    dplyr::group_by(macs_group, tlx_group, macs_name) %>%
    dplyr::summarize(frac_lognorm)


  x = repliseqTime2breaks_sum %>%
    dplyr::filter(repliseqTime_chrom=="chr6" & dplyr::between(repliseqTime_start, 68e6, 75e6)) %>%
    dplyr::filter(is.finite(frac_norm))
  ggplot(x) +
    geom_line(aes(x=repliseqTime_start, y=repliseqTime_avg)) +
    geom_line(aes(x=repliseqTime_start, y=log2(frac_norm), color=tlx_group, group=paste(tlx_group, tlx_group_i))) +
    facet_wrap(~tlx_group_i)





  dim(repliseqTime2breaks)
  repliseqTime_df$breaks_count = IRanges::countOverlaps(repliseqTime_ranges, tlx_ranges)
  repliseqTime_df$repliseqTime_cut = cut(repliseqTime_df$repliseqTime_avg, seq(0, 16, 2))

  repliseqTime_df %>%
    dplyr::filter(breaks_count>0) %>%
    dplyr::group_by(repliseqTime_cut) %>%
    dplyr::summarize(x=mean(breaks_count), n=n())

  ggplot(repliseqTime_df %>% dplyr::filter(breaks_count>0)) +
    geom_violin(aes(x=repliseqTime_cut, y=log2(breaks_count)))


  data_ranges = GenomicRanges::makeGRangesFromDataFrame(macs_df %>% dplyr::mutate(seqnames=macs_chrom, start=macs_start, end=macs_end), keep.extra.columns=T)
  data$hits = IRanges::countOverlaps(data_ranges, repliseqTime_ranges)





  data = dplyr::bind_rows(macs_df, random_df) %>% dplyr::select(macs_group, macs_chrom, macs_start, macs_end)
  data_ranges = GenomicRanges::makeGRangesFromDataFrame(data %>% dplyr::mutate(seqnames=macs_chrom, start=macs_start, end=macs_end), keep.extra.columns=T)
  data$hits = IRanges::countOverlaps(data_ranges, repliseqTime_ranges)




  ggplot(repliseqTime_df$repliseqTime_celltype %>% dplyr::filter(repliseqTime_chrom=="chr6" & celltype=="npc")) +
    geom_line(aes())

  x = repliseq_df %>% dplyr::filter(repliseq_chrom=="chr6" & dplyr::between(repliseq_start, 68e6, 75e6))
  ggplot(x) +
    geom_tile(aes(x=repliseq_start, y=repliseq_fraction, fill=repliseq_value))
  x = repliseqTime_df %>% dplyr::filter(repliseqTime_chrom=="chr6" & dplyr::between(repliseqTime_start, 68e6, 75e6))
  ggplot(x) +
    geom_line(aes(x=repliseqTime_start, y=repliseqTime_avg))




  # Calculate IZ
  starts = unique(repliseqTime_df$repliseqTime_start)
  repliseqIZ_df = repliseqTime_df %>%
    reshape2::melt(id.vars=c("repliseqTime_celltype", "repliseqTime_chrom", "repliseqTime_start", "repliseqTime_end"), value.name="repliseqTime_time", variable.name="repliseqTime_source") %>%
    tidyr::extract(repliseqTime_source, "repliseqTime_source", "repliseqTime_([^.]+)") %>%
    dplyr::filter(repliseqTime_source %in% c("avg")) %>%
    dplyr::rename(iz_chrom="repliseqTime_chrom", iz_source="repliseqTime_source", iz_celltype="repliseqTime_celltype") %>%
    dplyr::group_by(iz_chrom, iz_celltype, iz_source) %>%
    dplyr::do((function(z){
      zz<<-z
      data.frame(iz_start=which(ggpmisc:::find_peaks(17-z$repliseqTime_time, ignore_threshold = 1/16, span=9, strict=F))) %>%
      dplyr::mutate(iz_start=z$repliseqTime_start[iz_start])
    })(.)) %>%
    dplyr::arrange(iz_start) %>%
    dplyr::group_by(iz_chrom) %>%
    dplyr::mutate(iz_cluster=dbscan::dbscan(as.matrix(iz_start), eps=50000*5, minPts=1)$cluster) %>%
    dplyr::group_by(iz_chrom, iz_cluster, iz_source) %>%
    dplyr::mutate(
      iz_start=starts[which.min(abs(starts-mean(iz_start)))],
      iz_source=paste(unique(iz_source), collapse=",")) %>%
    dplyr::inner_join(repliseqTime_df %>% dplyr::select(repliseqTime_chrom, repliseqTime_start, repliseqTime_celltype, repliseqTime_avg), by=c("iz_chrom"="repliseqTime_chrom", "iz_start"="repliseqTime_start", "iz_celltype"="repliseqTime_celltype")) %>%
    dplyr::distinct(iz_chrom, iz_start, iz_source, iz_celltype, .keep_all=T)
  repliseqIZmerged_df = repliseqIZ_df %>%
    reshape2::dcast(iz_chrom+iz_start+iz_source~iz_celltype, value.var="repliseqTime_avg")

  if(F) {
    pdf(file="reports/repliseq_termination.pdf", width=80, height=40)
    celltype_pal = c("npc"="#E41A1C", "esc"="#E97F02")
    repliseq_gg = repliseq_df %>% dplyr::filter(repliseq_celltype=="npc", grepl("chr6$", repliseq_chrom) & dplyr::between(repliseq_start, 6e7, 9e7)) %>% dplyr::mutate(chrom=repliseq_chrom)
    repliseqTime_gg = repliseqTime_df %>% dplyr::filter(grepl("chr6$", repliseqTime_chrom) & dplyr::between(repliseqTime_start, 6e7, 9e7)) %>% dplyr::mutate(chrom=repliseqTime_chrom)
    repliseqIZ_gg = repliseqIZ_df %>% dplyr::filter(iz_source=="avg" & grepl("chr6$", iz_chrom) & dplyr::between(iz_start, 6e7, 9e7)) %>% dplyr::mutate(chrom=iz_chrom)
    repliseqIZmerged_gg = repliseqIZmerged_df %>% dplyr::filter(iz_source=="avg" & grepl("chr6$", iz_chrom) & dplyr::between(iz_start, 6e7, 9e7)) %>% dplyr::mutate(chrom=iz_chrom)
    transcription_gg = transcription_df %>% dplyr::filter(grepl("chr6$", tena2020transcription_chrom) & dplyr::between(tena2020transcription_start, 6e7, 9e7)) %>% dplyr::mutate(chrom=tena2020transcription_chrom)

    ggplot() +
      geom_tile(aes(y=repliseq_fraction, x=repliseq_start, fill=repliseq_value), data=repliseq_gg) +
      geom_line(aes(x=repliseqTime_start, y=repliseqTime_avg, color=repliseqTime_celltype), data=repliseqTime_gg, size=1) +
      geom_segment(aes(x=iz_start, xend=iz_start, y=pmin(npc, esc), yend=pmax(npc, esc), color="npc"), size=0.5, linetype="dashed", data=repliseqIZmerged_gg) +
      geom_point(aes(x=iz_start, y=repliseqTime_avg, color=iz_celltype), size=3,  data=repliseqIZ_gg) +
      geom_rect(aes(xmin=tena2020transcription_start, xmax=tena2020transcription_end, ymin=-2*(tena2020transcription_celltype_num-1)-2, ymax=-2*(tena2020transcription_celltype_num-1), alpha=tena2020transcription_rpkm), fill="#FF8C00", data=transcription_gg) +
      scale_fill_gradientn(colours=c("#666666", "#CCCCCC", "#FFFFFF", "#00FFFF", "#000066"), values=c(0, 0.1, 0.3, 0.5, 1)) +
      scale_color_manual(values=celltype_pal)+
      facet_wrap(~chrom, ncol=1, scale="free_x")
    dev.off()
  }

  # Create random regions
  random_lengths = rdc_df %>% dplyr::select(rdc_chrom, rdc_length)
  random_df = random_lengths %>%
    dplyr::inner_join(as.data.frame(genome_info) %>% tibble::rownames_to_column("rdc_chrom"), by=c("rdc_chrom")) %>%
    dplyr::group_by(rdc_chrom) %>%
    dplyr::do((function(z){
      zz<<-z
      n = 5
      z.start = sample(z$seqlengths[1], nrow(z)*n)
      data.frame(rdc_chrom=z$rdc_chrom[1], seqlengths=z$seqlengths[1], start=z.start, rdc_length=rep(z$rdc_length, each=n)) %>%
        dplyr::mutate(rdc_start=ifelse(start+rdc_length>=seqlengths, start-rdc_length, start), rdc_end=start+z$rdc_length, signal="Random positions", rdc_gene="") %>%
        dplyr::mutate(rdc_cluster=paste0(rdc_chrom, "_", 1:n()))
    })(.)) %>%
    dplyr::ungroup()

  random2_df = genes_df %>%
    dplyr::filter(!grepl(re_genes, gene_id) & gene_length>2e5 & gene_chrom %in% paste0("chr", 1:19)) %>%
    dplyr::mutate(rdc_chrom=gene_chrom, rdc_cluster=paste0(rdc_chrom, "_", 1:n()), rdc_start=gene_start, rdc_end=gene_end, rdc_length=gene_length, signal="Non-RDC genes", rdc_gene=gene_id)
  rdc_df = dplyr::bind_rows(rdc_df, random2_df, random_df) %>% dplyr::select(rdc_chrom, rdc_start, rdc_end, rdc_length, rdc_cluster, rdc_gene, signal)
  rdcMargins_df = data.frame()
  for(m in c(0, 5e5, 1e6, 5e6)) {
    rdcMargins_df = dplyr::bind_rows(rdcMargins_df, rdc_df %>% dplyr::mutate(margin=m))
  }

  #
  # Merge different datasets
  #
  repliseqTime_ranges = GenomicRanges::makeGRangesFromDataFrame(repliseqTime_df %>% dplyr::mutate(seqnames=repliseqTime_chrom, start=repliseqTime_start, end=repliseqTime_start), keep.extra.columns=T)
  repliseqIZ_ranges = GenomicRanges::makeGRangesFromDataFrame(repliseqIZ_df %>% dplyr::mutate(seqnames=iz_chrom, start=iz_start, end=iz_start), keep.extra.columns=T)
  repliseqIZmerged_ranges = GenomicRanges::makeGRangesFromDataFrame(repliseqIZmerged_df %>% dplyr::mutate(seqnames=iz_chrom, start=iz_start, end=iz_start), keep.extra.columns=T)
  repliseq_ranges = GenomicRanges::makeGRangesFromDataFrame(repliseq_df %>% dplyr::mutate(seqnames=repliseq_chrom, start=repliseq_start, end=repliseq_start), keep.extra.columns=T)
  rdc_ranges = GenomicRanges::makeGRangesFromDataFrame(rdcMargins_df %>% dplyr::mutate(seqnames=rdc_chrom, start=rdc_start-margin, end=rdc_end+margin), keep.extra.columns=T)
  genes_ranges = GenomicRanges::makeGRangesFromDataFrame(genes_df %>% dplyr::mutate(seqnames=gene_chrom, start=gene_start, end=gene_end), keep.extra.columns=T)
  transcription_ranges = GenomicRanges::makeGRangesFromDataFrame(transcription_df %>% dplyr::mutate(seqnames=tena2020transcription_chrom, start=tena2020transcription_start, end=tena2020transcription_end), keep.extra.columns=T)
  breaks_ranges = GenomicRanges::makeGRangesFromDataFrame(breaks_df %>% dplyr::mutate(seqnames=junction_chrom, start=junction_start, end=junction_end, qvalue=qvalue.sample_vs_max), keep.extra.columns=T)
  ctcf_ranges = GenomicRanges::makeGRangesFromDataFrame(ctcf_df %>% dplyr::mutate(seqnames=ctcfdb_chrom, start=ctcfdb_start, end=ctcfdb_end), keep.extra.columns=T)
  tad_ranges = GenomicRanges::makeGRangesFromDataFrame(tad_df %>% dplyr::mutate(seqnames=tad_chrom, start=tad_start, end=tad_end), keep.extra.columns=T)
  gc_ranges = GenomicRanges::makeGRangesFromDataFrame(gc_df %>% dplyr::mutate(seqnames=gc_chrom, start=gc_start, end=gc_end), keep.extra.columns=T)
  # rpkm_ranges = GenomicRanges::makeGRangesFromDataFrame(rpkm_df %>% dplyr::mutate(seqnames=rpkm_chrom, start=rpkm_start, end=rpkm_end), keep.extra.columns=T)

  #
  # Merge different data frames
  #
  repliseq2rdc_df = as.data.frame(IRanges::mergeByOverlaps(repliseq_ranges, rdc_ranges)) %>% dplyr::select(-dplyr::matches("_ranges\\."))
  repliseqTime2rdc_df = as.data.frame(IRanges::mergeByOverlaps(repliseqTime_ranges, rdc_ranges)) %>% dplyr::select(-dplyr::matches("_ranges\\."))
  genes2rdc_df = as.data.frame(IRanges::mergeByOverlaps(genes_ranges, rdc_ranges)) %>% dplyr::select(-dplyr::matches("_ranges\\."))
  repliseqIZ2rdc_df = as.data.frame(IRanges::mergeByOverlaps(repliseqIZ_ranges, rdc_ranges)) %>% dplyr::select(-dplyr::matches("_ranges\\."))
  repliseqIZmerged2rdc_df = as.data.frame(IRanges::mergeByOverlaps(repliseqIZmerged_ranges, rdc_ranges)) %>% dplyr::select(-dplyr::matches("_ranges\\."))
  transcription2rdc_df = as.data.frame(IRanges::mergeByOverlaps(transcription_ranges, rdc_ranges)) %>% dplyr::select(-dplyr::matches("_ranges\\."))
  # rpkm2rdc_df = as.data.frame(IRanges::mergeByOverlaps(rpkm_ranges, rdc_ranges)) %>% dplyr::select(-dplyr::matches("_ranges\\."))
  breaks2rdc_df = as.data.frame(IRanges::mergeByOverlaps(breaks_ranges, rdc_ranges)) %>% dplyr::select(-dplyr::matches("_ranges\\."))
  ctcf2rdc_df = as.data.frame(IRanges::mergeByOverlaps(ctcf_ranges, rdc_ranges)) %>% dplyr::select(-dplyr::matches("_ranges\\."))
  tad2rdc_df = as.data.frame(IRanges::mergeByOverlaps(tad_ranges, rdc_ranges)) %>% dplyr::select(-dplyr::matches("_ranges\\."))

  repliseqIZ2rdc_ranges = GenomicRanges::makeGRangesFromDataFrame(repliseqIZ2rdc_df %>% dplyr::mutate(seqnames=rdc_chrom, start=rdc_start, end=rdc_end), keep.extra.columns=T)
  repliseqIZ2rdc2transcription_df = as.data.frame(IRanges::mergeByOverlaps(repliseqIZ2rdc_ranges, transcription_ranges)) %>%
    dplyr::select(-dplyr::matches("_ranges\\.")) %>%
    dplyr::filter(iz_celltype==tena2020transcription_celltype) %>%
    dplyr::arrange(dplyr::desc(tena2020transcription_rpkm)) %>%
    dplyr::distinct(rdc_cluster, iz_celltype, iz_start, margin, signal, .keep_all=T)

  breaks2transcription_df = as.data.frame(IRanges::mergeByOverlaps(breaks_ranges, transcription_ranges)) %>%
    dplyr::select(-dplyr::matches("_ranges\\.")) %>%
    dplyr::group_by(tena2020transcription_celltype, tena2020transcription_gene) %>%
    dplyr::summarize(tena2020transcription_rpkm=max(tena2020transcription_rpkm), qvalue=max(qvalue.sample_vs_max))
  breaks2transcription_df %>%
    dplyr::group_by(tena2020transcription_celltype) %>%
    dplyr::summarize(R=cor(qvalue, tena2020transcription_rpkm, method="spearman", use="pairwise.complete.obs"))


  #
  # T-SNE
  #
  long_genes_ranges = GenomicRanges::makeGRangesFromDataFrame(genes_df %>% dplyr::mutate(seqnames=gene_chrom, start=gene_start-2e5, end=gene_end+2e5) %>% dplyr::filter(gene_length>2e5), keep.extra.columns=T)
  long_genes_ranges$gene_rdc = countOverlaps(long_genes_ranges, rdc_ranges[rdc_ranges$signal=="RDC"]) > 1
  long_genes_df = as.data.frame(long_genes_ranges)
  breaks2genes_df = as.data.frame(IRanges::mergeByOverlaps(long_genes_ranges, breaks_ranges)) %>%
    dplyr::group_by(gene_id, gene_chrom, gene_start, gene_end, gene_strand, gene_length, gene_rdc) %>%
    dplyr::summarize(hit=any(qvalue<1e-5))
  breaks2genes_ranges = GenomicRanges::makeGRangesFromDataFrame(breaks2genes_df %>% dplyr::mutate(seqnames=gene_chrom, start=gene_start-2e5, end=gene_end+2e5) %>% dplyr::filter(gene_length>2e5), keep.extra.columns=T)
  repliseq2genes_cluster_df = as.data.frame(IRanges::mergeByOverlaps(breaks2genes_ranges, repliseqTime_ranges)) %>%
    dplyr::select(-dplyr::matches("_ranges\\.")) %>%
    dplyr::filter(repliseqTime_celltype=="npc") %>%
    dplyr::group_by(repliseqTime_celltype, gene_id, hit, gene_rdc) %>%
    dplyr::do((function(z){
      zz <<- z

      n = 50
      as.data.frame(matrix(approx(z$repliseqTime_start, z$repliseqTime_avg, n=n)$y, nrow=1, dimnames=list(row.names=1, col.names=paste0("part_", 1:n)) ))
    })(.)) %>%
    dplyr::ungroup()

  repliseq2genes_cluster_matrix = jitter(as.matrix(repliseq2genes_cluster_df %>% dplyr::select(dplyr::matches("^part_"))))
  repliseq2genes_cluster_tsne = Rtsne(repliseq2genes_cluster_matrix, perplexity=8)
  x = repliseq2genes_cluster_df %>% dplyr::mutate(x=repliseq2genes_cluster_tsne$Y[,1], y=repliseq2genes_cluster_tsne$Y[,2])
  g = ggplot(x) +
    geom_point(aes(x=x, y=y, gene_id=gene_id, color=gene_rdc))
    # ggrepel::geom_text_repel(aes(x=x, y=y, label=gene_id))
  g
  # p = plotly::ggplotly(g)
  # htmlwidgets::saveWidget(plotly::as_widget(p), "index3.html")



  #
  # RANDOM FOREST DATA
  # ==================================
  genes_df$gene_rdc = countOverlaps(GenomicRanges::makeGRangesFromDataFrame(genes_df %>% dplyr::mutate(seqnames=gene_chrom, start=gene_start, end=gene_end)), rdc_ranges[rdc_ranges$signal=="RDC"]) > 1
  genes_long_df = genes_df %>% dplyr::filter(gene_length>2e5 & grepl("chr[0-9]+", genes_chrom))
  genes_promoter_ranges = GenomicRanges::makeGRangesFromDataFrame(rbind(genes_long_df %>% dplyr::mutate(seqnames=gene_chrom, start=gene_start-1e5, end=gene_start), genes_long_df %>% dplyr::mutate(seqnames=gene_chrom, start=gene_end, end=gene_end+1e5)) %>% dplyr::mutate(seqnames=gene_chrom), keep.extra.columns=T)
  genes_exact_ranges = GenomicRanges::makeGRangesFromDataFrame(genes_long_df %>% dplyr::mutate(seqnames=gene_chrom, start=gene_start, end=gene_end), keep.extra.columns=T)
  genes_extended_ranges = GenomicRanges::makeGRangesFromDataFrame(genes_long_df %>% dplyr::mutate(seqnames=gene_chrom, start=gene_start-2e5, end=gene_end+2e5), keep.extra.columns=T)
  genes_extended2_ranges = GenomicRanges::makeGRangesFromDataFrame(genes_long_df %>% dplyr::mutate(seqnames=gene_chrom, start=gene_start-5e5, end=gene_end+5e5), keep.extra.columns=T)

  breaks2genes_df = as.data.frame(IRanges::mergeByOverlaps(genes_exact_ranges, breaks_ranges)) %>%
    dplyr::select(-dplyr::matches("_ranges\\.")) %>%
    dplyr::mutate(qvalue=qvalue.sample_vs_max, qvalue_norm=-log10(qvalue), qvalue_norm=ifelse(qvalue_norm>20, 20, qvalue_norm)) %>%
    dplyr::group_by(gene_id, .drop=F) %>%
    dplyr::summarize(qvalue_norm=ifelse(n()>0, max(qvalue_norm, na.rm=T), 0))

  gc2genes_df = as.data.frame(IRanges::mergeByOverlaps(genes_promoter_ranges, gc_ranges)) %>%
    dplyr::group_by(gene_id, .drop=F) %>%
    dplyr::summarize(gc_freq=ifelse(n()>0, max(gc_freq, na.rm=T), 0.5)) %>%
    dplyr::ungroup()

  repliseq2genes_df = as.data.frame(IRanges::mergeByOverlaps(genes_extended_ranges, repliseqTime_ranges)) %>%
    dplyr::select(-dplyr::matches("_ranges\\.")) %>%
    dplyr::group_by(gene_id, repliseqTime_celltype, .drop=F) %>%
    dplyr::summarize(repliseqTime_slope=ifelse(n()>0, max(diff(abs(repliseqTime_avg)), na.rm=T), 0), repliseqTime_max=ifelse(n()>0, max(repliseqTime_avg), 0), repliseqTime_min=ifelse(n()>0, min(repliseqTime_avg), 0), repliseqTime_diff=ifelse(n()>0, repliseqTime_max-repliseqTime_min, 0))
  # repliseq2genes_ranges = GenomicRanges::makeGRangesFromDataFrame(repliseq2genes_df %>% dplyr::mutate(seqnames=gene_chrom, start=gene_start, end=gene_end), keep.extra.columns=T)

  transcription2genes_df = as.data.frame(IRanges::mergeByOverlaps(genes_exact_ranges, transcription_ranges)) %>%
    dplyr::select(-dplyr::matches("_ranges\\.")) %>%
    dplyr::group_by(gene_id, tena2020transcription_celltype, .drop=F) %>%
    dplyr::summarize(tena2020transcription_rpkm=ifelse(n()>0, max(tena2020transcription_rpkm), 0))

  tad2genes_df = as.data.frame(IRanges::mergeByOverlaps(genes_extended_ranges, tad_ranges)) %>%
    dplyr::select(-dplyr::matches("_ranges\\.")) %>%
    dplyr::group_by(gene_id, tad_celltype, .drop=F) %>%
    dplyr::summarize(A_TAD_all=ifelse(n()>0, all(tad_compartment=="A"), F), B_TAD_all=ifelse(n()>0, all(tad_compartment=="B"), F), A_TAD_any=ifelse(n()>0, any(tad_compartment=="A"), F), B_TAD_any=ifelse(n()>0, any(tad_compartment=="B"), F)) %>%
    dplyr::ungroup()

  IZ2genes_df = as.data.frame(IRanges::mergeByOverlaps(genes_extended2_ranges, repliseqIZ_ranges)) %>%
    dplyr::group_by(gene_id, iz_celltype, .drop=F) %>%
    dplyr::summarize(IZ_count=n(), IZ_dist=ifelse(n()>0, min(abs(gene_start+gene_length/2 - iz_start)/gene_length), 0), IZ_outside=ifelse(n()>0, sum(iz_start>gene_end | iz_start<gene_start), 0), IZ_inside=ifelse(n()>0, sum(iz_start>gene_end | iz_start<gene_start), 0)) %>%
    dplyr::ungroup()

  load("data/methylation/data_all.rda")
  methylation_ranges = GenomicRanges::makeGRangesFromDataFrame(methylation_df %>% dplyr::mutate(seqnames=methylation_chrom, start=methylation_start, end=methylation_end), keep.extra.columns=T)
  meth2genes_df = as.data.frame(IRanges::mergeByOverlaps(genes_promoter_ranges, methylation_ranges)) %>%
    dplyr::group_by(gene_id, methylation_celltype, methylation_type, .drop=F) %>%
    dplyr::summarize(methylation_maxscore=ifelse(n()>0, max(methylation_score, na.rm=T), 0)) %>%
    dplyr::ungroup() %>%
    reshape2::dcast(...~methylation_type, value.var="methylation_maxscore")

    breaks2transcription2repliseq2genes2tad2IZ2gc2meth_df = genes_long_df %>%
      dplyr::left_join(breaks2genes_df, by="gene_id") %>%
      dplyr::left_join(gc2genes_df, by="gene_id") %>%
      dplyr::left_join(meth2genes_df, by="gene_id") %>%
      dplyr::left_join(transcription2genes_df, by=c("gene_id", "methylation_celltype"="tena2020transcription_celltype")) %>%
      dplyr::left_join(repliseq2genes_df, by=c("gene_id", "methylation_celltype"="repliseqTime_celltype")) %>%
      dplyr::left_join(tad2genes_df, by=c("gene_id", "methylation_celltype"="tad_celltype")) %>%
      dplyr::left_join(IZ2genes_df, by=c("gene_id", "methylation_celltype"="iz_celltype")) %>%
      dplyr::mutate(
        gc_freq=tidyr::replace_na(gc_freq, 0.5),
        H3K27ac=tidyr::replace_na(H3K27ac, 0),
        H3K27me3=tidyr::replace_na(H3K27me3, 0),
        H3K4me3=tidyr::replace_na(H3K4me3, 0),
        H3K9me3=tidyr::replace_na(H3K9me3, 0),
        tena2020transcription_rpkm=tidyr::replace_na(tena2020transcription_rpkm, 0),
        repliseqTime_slope=tidyr::replace_na(repliseqTime_slope, 0),
        repliseqTime_max=tidyr::replace_na(repliseqTime_max, 0),
        repliseqTime_min=tidyr::replace_na(repliseqTime_min, 0),
        repliseqTime_diff=tidyr::replace_na(repliseqTime_diff, 0),
        A_TAD_all=tidyr::replace_na(A_TAD_all, F),
        B_TAD_all=tidyr::replace_na(B_TAD_all, F),
        A_TAD_any=tidyr::replace_na(A_TAD_any, F),
        B_TAD_any=tidyr::replace_na(B_TAD_any, F),
        IZ_count=tidyr::replace_na(IZ_count, 0),
        IZ_dist=tidyr::replace_na(IZ_dist, 0),
        IZ_outside=tidyr::replace_na(IZ_outside, 0),
        IZ_inside=tidyr::replace_na(IZ_inside, 0))



  # Find genes different between NPC and ESC
  genes_diff = transcription2genes_df %>%
    dplyr::group_by(tena2020transcription_celltype) %>%
    dplyr::mutate(z_score=(tena2020transcription_rpkm-mean(tena2020transcription_rpkm))/sd(tena2020transcription_rpkm)) %>%
    dplyr::ungroup() %>%
    reshape2::dcast(gene_id ~ tena2020transcription_celltype, value.var="tena2020transcription_rpkm") %>%
    dplyr::filter(xor(npc>=0.2, esc>=0.2) | abs(npc-esc)>=0.5) %>%
    .$gene_id
  # plot(density(with(genes_diff, (npc-esc)[abs(npc-esc)<5])))
  # length(genes_diff)
  # genes_diff = repliseq2genes_df %>%
  #   reshape2::dcast(gene_id ~ repliseqTime_celltype, value.var="repliseqTime_diff") %>%
  #   dplyr::mutate(diff=npc-esc, z_score=(diff-mean(diff))/sd(diff)) %>%
  #   dplyr::filter(abs(z_score)>=1) %>%
  #   .$gene_id

  ggplot(transcription2genes_df %>% dplyr::filter(tena2020transcription_rpkm<=4)) +
    geom_density(aes(x=tena2020transcription_rpkm, color=tena2020transcription_celltype), bw=0.03, n=2^15) +
    geom_vline(xintercept=0.2)


  #
  # Random forest
  #
  # prediction_cols = c("repliseqTime_avg", "A_TAD_all", "B_TAD_all", "A_TAD_any", "B_TAD_any", "tena2020transcription_rpkm")
  # prediction_cols = c("repliseqTime_avg", "tena2020transcription_rpkm", "gene_length")
  prediction_cols = c("repliseqTime_max", "repliseqTime_min", "repliseqTime_diff", "repliseqTime_slope", "tena2020transcription_rpkm", "gene_length", "A_TAD_all", "B_TAD_all", "A_TAD_any", "B_TAD_any", "IZ_count", "IZ_dist", "IZ_outside", "IZ_inside", "gc_freq", "H3K27me3", "H3K4me3", "H3K9me3")
  # prediction_cols = c("repliseqTime_diff", "repliseqTime_slope", "tena2020transcription_rpkm", "gene_length", "IZ_dist", "gc_freq", "H3K27me3", "H3K4me3", "H3K9me3")
  prediction_cols = c("repliseqTime_diff", "repliseqTime_slope", "tena2020transcription_rpkm", "gene_length", "IZ_dist", "A_TAD_all")

  all_df = breaks2transcription2repliseq2genes2tad2IZ2gc2meth_df %>%
    dplyr::filter(!is.na(methylation_celltype))
    # dplyr::filter(gene_id %in% genes_diff)
  all_df = all_df[rowSums(is.na(all_df[,prediction_cols]))==0,]

  oob3 = data.frame()
  perf_results = data.frame()
  varimp = data.frame()
  for(celltype in c("esc", "npc")) {
    all_df_ct = all_df %>% dplyr::filter(methylation_celltype==celltype)
    all_df_ct$y = factor(ifelse(all_df_ct$gene_rdc, "breaks", "no breaks"))
    # all_df_ct$y = factor(ifelse(all_df_ct$qvalue_norm>=3, "breaks", "no breaks"))
    for(i in 1:50) {
      print(paste(celltype, i))

      training_size = round(sum(all_df_ct$y=="breaks")*0.8)
      training_hit_i = sample(which(all_df_ct$y=="breaks"), training_size)
      training_non_i = sample(setdiff(1:nrow(all_df_ct), which(all_df_ct$y=="breaks")), training_size)
      data_training = all_df_ct[c(training_hit_i, training_non_i),]

      testing_hit_i = setdiff(which(all_df_ct$y=="breaks"), training_hit_i)
      testing_non_i = sample(setdiff(which(all_df_ct$y!="breaks"), training_non_i), length(testing_hit_i))
      data_testing = all_df_ct[c(testing_hit_i, testing_non_i),]


      combs = combn(setdiff(prediction_cols, "gene_length"), 2)
      for(i in 1:ncol(combs)) {
        i.cols = c("gene_length", combs[,i])

        input_training = as.matrix(data_training[,i.cols])
        model = randomForest(y=data_training$y, x=input_training, importance=T)

        oob3 = rbind(oob3, data.frame(celltype=celltype, cols=paste(i.cols, collapse="|"), oob=mean(model$err.rate[,1])))
      }

      for(col in c("all")) {
        i.cols = setdiff(prediction_cols, col)

        input_training = as.matrix(data_training[,i.cols])
        model = randomForest(y=data_training$y, x=input_training, importance=T)

        # library(party)
        # model = party::ctree(y ~ repliseqTime_slope, data=data_training, controls=cforest_control(mtry=2, mincriterion=0))
        # plot(model)

        y_ = predict(model, data_testing[,i.cols], type = "prob")
        pred <- ROCR::prediction(as.matrix(y_[,2]), as.matrix(as.numeric(data_testing$y)))
        perf <- ROCR::performance(pred, "tpr", "fpr")
        auc <- ROCR::performance(pred, "auc")@y.values[[1]][1]
        perf_results = dplyr::bind_rows(perf_results, data.frame(celltype=celltype, TPR=perf@y.values[[1]], FPR=perf@x.values[[1]], auc=auc, missing_column=col, i=i))

        if(col=="all") {
          varimp_i = as.data.frame(randomForest::importance(model)) %>%
            tibble::rownames_to_column("variable") %>%
            dplyr::mutate(i=i, celltype=celltype, MeanDecreaseAccuracy=MeanDecreaseAccuracy, MeanDecreaseGini=MeanDecreaseGini) %>%
            dplyr::select(variable, i, celltype, MeanDecreaseAccuracy, MeanDecreaseGini)

          varimp_i$R = NA_real_
          for(col2 in prediction_cols) {
            cor_i = cor(data_testing[[col2]], data_testing$qvalue_norm, method="spearman", use="pairwise.complete.obs")
            varimp_i$R[varimp_i$variable==col2] = cor_i
          }

          varimp = rbind(varimp, varimp_i %>% reshape2::melt(id.vars=c("i", "variable", "celltype"), variable.name="metric"))

        }
      }
    }
  }

  varimp.sorted = as.character(varimp %>%
    dplyr::filter(metric=="MeanDecreaseAccuracy") %>%
    dplyr::group_by(variable, metric) %>%
    dplyr::summarize(value.sd=sd(value), value.mean=mean(value)) %>%
    dplyr::arrange(value.mean) %>%
    dplyr::distinct(dplyr::desc(variable)) %>%
    .$variable)
  varimp = varimp %>%
    dplyr::mutate(variable=factor(variable, varimp.sorted)) %>%
    dplyr::arrange(variable)

  perf_results.sum = perf_results %>%
    dplyr::group_by(celltype, missing_column, FPR) %>%
    dplyr::summarize(TPR=mean(TPR, na.rm=T)) %>%
    dplyr::group_by(celltype, missing_column) %>%
    dplyr::mutate(TPR=smoother::smth.gaussian(TPR, window=3)) %>%
    dplyr::ungroup()

  pdf(file="reports/random_forest_npc2esc_rpkm.pdf", width=12, height=6)
  ggplot(varimp) +
    geom_hline(yintercept=0, color="#000000") +
    geom_boxplot(aes(x=variable, y=value, fill=celltype)) +
    coord_flip() +
    facet_wrap(~metric, scales="free")

  ggplot(oob3) +
    geom_boxplot(aes(x=cols, y=oob, fill=celltype)) +
    coord_flip()


  perf_results.ggplot = perf_results.sum %>%
    dplyr::filter(missing_column=="all")
  ggplot(perf_results.ggplot) +
    geom_line(aes(x=FPR, y=TPR, color=celltype), size=2) +
    geom_abline(slope=1) +
    coord_cartesian(xlim=c(0, 1), ylim=c(0,1))
  dev.off()



  model = randomForest(y=all_df$y, x=all_df[,prediction_cols], importance=T)
  y_ = predict(model, all_df[,prediction_cols], type = "prob")
  table(all_df$y, y_[,1]>y_[,2])

  transcription2repliseqTime_df = as.data.frame(IRanges::mergeByOverlaps(transcription_ranges, repliseqTime_ranges)) %>%
    dplyr::select(-dplyr::matches("_ranges\\.")) %>%
    dplyr::mutate(tena2020transcription_length=tena2020transcription_end-tena2020transcription_start) %>%
    dplyr::left_join(transcription2rdc_df %>% dplyr::filter(margin==0) %>% dplyr::distinct(tena2020transcription_celltype, rdc_cluster, tena2020transcription_gene, signal), by=c("tena2020transcription_celltype", "tena2020transcription_gene")) %>%
    dplyr::inner_join(tad2rdc_df %>% dplyr::filter(margin==0) %>% dplyr::distinct(tad_celltype, tad_compartment, rdc_cluster, signal), by=c("tena2020transcription_celltype"="tad_celltype", "rdc_cluster", "signal")) %>%
    dplyr::filter(repliseqTime_celltype==tena2020transcription_celltype & !is.na(signal) & tena2020transcription_length>1e5) %>%
    dplyr::group_by(tena2020transcription_celltype, signal, tena2020transcription_gene, tena2020transcription_rpkm, tena2020transcription_length) %>%
    dplyr::summarize(repliseqTime_avg=max(repliseqTime_avg), tad_A=ifelse(any(tad_compartment=="A"), "A-TAD", "No A-TAD"), tad_B=ifelse(any(tad_compartment=="B"), "B-TAD", "No B-TAD")) %>%
    dplyr::ungroup()

  transcription2repliseqTime_df %>%
    dplyr::group_by(tena2020transcription_celltype, signal) %>%
    dplyr::summarize(cor=cor(tena2020transcription_rpkm, repliseqTime_avg, method="spearman", use="pairwise.complete.obs"))

  pdf(file="reports/transcription_vs_rpkm.pdf", width=8, height=6)
  transcription2repliseqTime_df %>%
    dplyr::filter(tena2020transcription_celltype=="npc") %>%
    # dplyr::filter(rpkm_value<5) %>%
    dplyr::mutate(tena2020transcription_rpkm=cut(tena2020transcription_rpkm, 10), tena2020transcription_length=cut(tena2020transcription_length, 5)) %>%
    ggplot() +
    geom_boxplot(aes(x=tena2020transcription_rpkm, y=repliseqTime_avg, fill=signal))
    # geom_beeswarm(aes(x=tena2020transcription_rpkm, y=repliseqTime_avg, group=signal), dodge.width=0.8, size=0.5) +
    # facet_wrap(~tad_A)
  dev.off()


  if(F)
  {
    pdf(file="reports/ctcf_counts.pdf", width=10, height=8)
    x = ctcf2rdc_df %>%
      dplyr::group_by(rdc_cluster, margin, signal) %>%
      dplyr::summarize(CTCF=1e6*n()/(rdc_length+margin*2)) %>%
      dplyr::ungroup() %>%
      dplyr::arrange(margin) %>%
      dplyr::mutate(margin_text=paste0(round(margin/1e6, 1), "M")) %>%
      dplyr::mutate(margin_text=factor(margin_text, unique(margin_text)))
    ggplot(x) +
      geom_boxplot(aes(x=margin_text, y=CTCF, fill=signal)) +
      labs(y="Experimenta CTCF count per million", x="") +
      theme_bw(base_size = 30) +
      theme(legend.position="bottom", axis.title.x=element_blank(), legend.text=element_text(size=20),legend.title=element_text(size=20))
    dev.off()

    pdf(file="reports/tad_counts.pdf", width=10, height=8)
    x = tad2rdc_df %>%
      dplyr::group_by(tad_celltype, tad_compartment, rdc_cluster, margin, signal) %>%
      dplyr::summarize(TAD=1e6*n()/(rdc_length+margin*2)) %>%
      dplyr::ungroup() %>%
      dplyr::arrange(margin) %>%
      dplyr::mutate(margin_text=paste0(round(margin/1e6, 1), "M")) %>%
      dplyr::mutate(margin_text=factor(margin_text, unique(margin_text)))
    ggplot(x) +
      geom_boxplot(aes(x=margin_text, y=TAD, fill=signal)) +
      labs(y="Experimenta TAD count per million", x="") +
      theme_bw(base_size = 30) +
      theme(legend.position="bottom", axis.title.x=element_blank(), legend.text=element_text(size=20),legend.title=element_text(size=20)) +
      facet_wrap(~tad_celltype+tad_compartment, scales="free")
    dev.off()

    pdf(file="reports/iz_counts.pdf", width=10, height=8)
    x = repliseqIZ2rdc2transcription_df %>%
      dplyr::filter(signal=="RDC") %>%
      dplyr::group_by(rdc_cluster, signal, margin) %>%
      dplyr::mutate(transcribed=paste0(unique(iz_celltype[tena2020transcription_rpkm>=0.1]), collapse="/"))  %>%
      dplyr::mutate(transcribed=ifelse(transcribed=="", "none", transcribed)) %>%
      dplyr::group_by(rdc_cluster, iz_celltype, signal, margin, transcribed) %>%
      dplyr::summarize(IZ=1e6*n()/rdc_length[1]) %>%
      dplyr::ungroup() %>%
      dplyr::arrange(margin) %>%
      dplyr::mutate(margin_text=paste0(round(margin/1e6, 1), "M")) %>%
      dplyr::mutate(margin_text=factor(margin_text, unique(margin_text)))
    ggplot(x) +
      geom_boxplot(aes(x=margin_text, y=IZ, fill=transcribed), outlier.alpha=0) +
      geom_beeswarm(aes(x=margin_text, y=IZ, group=transcribed), dodge.width=0.8, size=0.5) +
      labs(y="Repliseq peaks (per Mb)", x="") +
      theme_bw(base_size = 30) +
      facet_wrap(~iz_celltype, scales="free_y") +
      theme(legend.position="bottom", axis.title.x=element_blank(), legend.text=element_text(size=20),legend.title=element_text(size=20))


    y = x %>%
      reshape2::dcast(rdc_cluster+signal+margin_text+transcribed~iz_celltype, value.var="IZ")
    ggplot(y) +
      geom_point(aes(x=esc, y=npc, color=transcribed)) +
      geom_abline(slope=1) +
      theme_bw(base_size = 30) +
      labs(x="Repliseq peaks in ESC (per Mb)", y="Repliseq peaks in NPC (per Mb)") +
      facet_wrap(~margin_text)


    x = repliseqIZ2rdc_df %>%
      dplyr::group_by(rdc_cluster, iz_celltype, signal, margin) %>%
      dplyr::summarize(IZ=1e6*n()/rdc_length) %>%
      dplyr::ungroup() %>%
      dplyr::arrange(margin) %>%
      dplyr::mutate(margin_text=paste0(round(margin/1e6, 1), "M")) %>%
      dplyr::mutate(margin_text=factor(margin_text, unique(margin_text)))
    ggplot(x) +
      geom_boxplot(aes(x=margin_text, y=IZ, fill=signal)) +
      labs(y="Repliseq peaks (per Mb)", x="") +
      theme_bw(base_size = 30) +
      facet_wrap(~iz_celltype) +
      theme(legend.position="bottom", axis.title.x=element_blank(), legend.text=element_text(size=20),legend.title=element_text(size=20))
    dev.off()
  }

  pdf(file="reports/repliseq_to_rdc_new4.pdf", width=25, height=60)
  re_genes = "Ctnna"
  # re_genes = "Lsamp|Ctnna|Csmd1"
  # re_genes = "1700048|O20Rik|Ackr2|Astn2|Auts2|Bai3|Ccser1|Cdh13|Celf4|Csmd1|Ctnna2|Ctnnd2|Dcc|Dock1|Fhit|Gpc6|Grid2|Grik2|Hdac9|Lrp1b|Magi2|Naaladl2|Nbea|Nkain2|Nrg3|Nrxn1|Pard3b|Pcdh9|Prkg1|Rbfox1|Rbms3|Sdk1|Tenm4|Tpgs2|Wwox"
  # re_genes = "[A-Za-z]*"

  repliseqIZ2rdc_gg = repliseqIZ2rdc_df %>% dplyr::filter(margin==1e6 & signal=="RDC" & grepl(re_genes, rdc_gene))
  repliseqIZmerged2rdc_gg = repliseqIZmerged2rdc_df %>% dplyr::filter(margin==1e6 & signal=="RDC" & grepl(re_genes, rdc_gene))
  repliseq2rdc_gg = repliseq2rdc_df %>% dplyr::filter(margin==1e6 & signal=="RDC" & grepl(re_genes, rdc_gene)) %>% dplyr::mutate(repliseq_celltype=="npc")
  repliseqTime2rdc_gg = repliseqTime2rdc_df %>% dplyr::filter(margin==1e6 & signal=="RDC" & grepl(re_genes, rdc_gene))
  rdc_gg = rdcMargins_df %>% dplyr::filter(margin==1e6 & signal=="RDC" & grepl(re_genes, rdc_gene))
  genes2rdc_gg = genes2rdc_df %>% dplyr::filter(margin==1e6 & signal=="RDC" & grepl(re_genes, gene_id))
  transcription2rdc_gg = transcription2rdc_df %>% dplyr::filter(margin==1e6 & signal=="RDC" & grepl(re_genes, rdc_gene)) %>%
    dplyr::mutate(tena2020transcription_rpkm=ifelse(tena2020transcription_rpkm>=20, 20, tena2020transcription_rpkm))
  # rpkm2rdc_gg = rpkm2rdc_df %>% dplyr::filter(margin==1e6 & signal=="RDC" & grepl(re_genes, rdc_gene))
  ctcf2rdc_gg = ctcf2rdc_df %>% dplyr::filter(margin==1e6 & signal=="RDC" & grepl(re_genes, rdc_gene))
  tad2rdc_gg = tad2rdc_df %>% dplyr::filter(margin==1e6 & signal=="RDC" & grepl(re_genes, rdc_gene))
  breaks2rdc_gg = breaks2rdc_df %>% dplyr::filter(margin==1e6 & signal=="RDC" & grepl(re_genes, rdc_gene)) %>%
    dplyr::mutate(qvalue_norm = -log10(qvalue), ifelse(!is.finite(qvalue_norm), max(qvalue_norm[is.finite(qvalue_norm)]), qvalue_norm)) %>%
    dplyr::group_by(rdc_cluster) %>%
    dplyr::mutate(qvalue_max=max(qvalue_norm), qvalue_norm=qvalue_norm/qvalue_max) %>%
    dplyr::filter(junction_end-junction_start < 15e3) %>%
    dplyr::ungroup()

  celltype_pal = c("npc"="#E41A1C", "esc"="#E97F02")
  # , alpha=breaks2rdc_gg$qvalue_max/max(breaks2rdc_gg$qvalue_max, na.rm=T)
  ggplot() +
    geom_line(aes(x=junction_start, y=17+qvalue_norm*5), data=breaks2rdc_gg, size=1, alpha=breaks2rdc_gg$qvalue_max)  +
    geom_tile(aes(y=repliseq_fraction, x=repliseq_start, fill=repliseq_value), data=repliseq2rdc_gg) +
    geom_line(aes(x=repliseqTime_start, y=repliseqTime_avg, color=repliseqTime_celltype), data=repliseqTime2rdc_gg, size=1) +
    geom_ribbon(aes(x=repliseqTime_start, ymin=repliseqTime_min, ymax=repliseqTime_max, color=repliseqTime_celltype), data=repliseqTime2rdc_gg, alpha=0.1) +
    geom_segment(aes(x=iz_start, xend=iz_start, y=pmin(npc, esc), yend=pmax(npc, esc), color="npc"), size=0.2, linetype="twodash", data=repliseqIZmerged2rdc_gg) +
    geom_point(aes(x=iz_start, y=repliseqTime_avg, color=iz_celltype), size=3,  data=repliseqIZ2rdc_gg) +
    geom_rect(aes(xmin=rdc_start, xmax=rdc_end, ymin=-2, ymax=0), data=rdc_gg) +
    geom_text(aes(x=(rdc_end+rdc_start)/2, y=-1), label="rdc", data=rdc_gg) +
    # geom_rect(aes(xmin=ctcfdb_start, xmax=ctcfdb_end, ymin=-4, ymax=-2), data=ctcf2rdc_gg, color="#E7298A") +
    geom_rect(aes(xmin=tad_start, xmax=tad_end, ymin=-5+tad_celltype_num, ymax=-4+tad_celltype_num), fill=c("A"="#FF0000","B"="#0000FF")[tad2rdc_gg$tad_compartment], color="#FFFFFF00", data=tad2rdc_gg, color="#E7298A") +
    geom_rect(aes(xmin=gene_start, xmax=gene_end, ymin=-6, ymax=-4), fill="#FF8C00", color="#000000", data=genes2rdc_gg) +
    geom_text(aes(x=(gene_start+gene_end)/2, y=-5, label=gene_id), data=genes2rdc_gg, size=5) +

    geom_rect(aes(xmin=tena2020transcription_start, xmax=tena2020transcription_end, ymin=-9+tena2020transcription_celltype_num, ymax=-8+tena2020transcription_celltype_num, alpha=tena2020transcription_rpkm), fill="#377EB8", data=transcription2rdc_gg) +
    geom_text(aes(x=rdc_start, y=-8.5+tena2020transcription_celltype_num, label=paste(tena2020transcription_celltype)), data=transcription2rdc_gg %>% dplyr::distinct(rdc_cluster, tena2020transcription_celltype, .keep_all=T)) +

    # geom_rect(aes(xmin=rpkm_start, xmax=rpkm_end, ymin=-11+rpkm_celltype_num, ymax=-10+rpkm_celltype_num, alpha=rpkm_value), fill="#377EB8", data=rpkm2rdc_gg) +
    # geom_text(aes(x=rdc_start, y=-10.5+rpkm_celltype_num, label=paste(rpkm_celltype, "(RPKM)")), data=rpkm2rdc_gg %>% dplyr::distinct(rdc_cluster, rpkm_celltype, .keep_all=T)) +

    scale_fill_gradientn(colours=c("#666666", "#CCCCCC", "#FFFFFF", "#00FFFF", "#000066"), values=c(0, 0.1, 0.3, 0.5, 1)) +
    scale_color_manual(values=celltype_pal) +
    scale_x_continuous(breaks=scale_1mb) +
    coord_cartesian(ylim=c(-9, 22)) +
    facet_wrap(rdc_cluster~., scales="free_x", ncol=5, strip.position="left")


  # +
  #     facet_wrap(celltype~., scales="free_x", ncol=5, strip.position="left")
  # p = plotly::ggplotly(g)
  # htmlwidgets::saveWidget(plotly::as_widget(p), "index.html")
  dev.off()
}