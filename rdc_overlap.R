devtools::load_all('breaktools/')
library(ggplot2)
library(GenomicRanges)
library(IRanges)
library(dplyr)
library(readr)
library(rtracklayer)
library(tidyr)
library(ggthemr)
library(ggfortify)

main = function() {
  effective_size = 1.87e9
  sgRNA_length = 19
  extsize = 1e5
  maxgap = extsize*2
  exttype = "symetrical"
  threshold_qvalue = 1e-2
  threshold_pileup = 1
  slocal = 1e7
  llocal = 1e7
  bait_region=6e6

  #
  # Offtargets
  #
  # bowtie2 -k 30 -N 1 -x mm10/mm10 --end-to-end --very-sensitive -c acgagcatttccaaccc

  #
  # Chromosomes sizes
  #
  sizes_df = readr::read_tsv("genomes/mm10/annotation/mm10.chrom.sizes", col_names=c("sizes_chrom", "sizes_length")) %>%
    dplyr::mutate(sizes_effective=sizes_length/sum(sizes_length)*effective_size)

  #
  # Read RDC
  #
  rdc_pnas_df = readr::read_tsv("data/rdc_pnas.tsv")
  rdc_pnas_ranges = GenomicRanges::makeGRangesFromDataFrame(rdc_pnas_df %>% dplyr::mutate(seqnames=rdc_chrom, start=rdc_start, end=rdc_end), keep.extra.columns=T)
  chain_mm9_mm10 = import.chain("genomes/mm9/mm9ToMm10.over.chain")
  rdc_pnas_ranges = unlist(rtracklayer::liftOver(rdc_pnas_ranges, chain_mm9_mm10))
  rdc_pnas_df = as.data.frame(rdc_pnas_ranges) %>%
    dplyr::distinct(rdc_chrom, rdc_start, rdc_end, .keep_all=T) %>%
    dplyr::select(dplyr::matches("rdc_")) %>%
    dplyr::mutate(rdc_name=paste(rdc_chrom, rdc_start, rdc_end))
  rdc_pnas_ranges = GenomicRanges::makeGRangesFromDataFrame(rdc_pnas_df %>% dplyr::mutate(seqnames=rdc_chrom, start=rdc_start, end=rdc_end), keep.extra.columns=T)

  #
  # Read TLX files
  #
  samples_df = readr::read_tsv("data/tlx_samples_extended.tsv")
  tlx_df = tlx_read_many(samples_df)
  tlx_df = tlx_remove_rand_chromosomes(tlx_df)
  tlx_df = tlx_mark_bait_chromosome(tlx_df)
  tlx_df = tlx_mark_bait_junctions(tlx_df, bait_region)
  baits_df = tlx_identify_baits(tlx_df, breaksite_size=sgRNA_length)
  # tlx_df = tlx_mark_repeats(tlx_df, repeatmasker_df)
  # tlx_df = tlx_df %>% dplyr::filter(!tlx_is_bait_junction)
  tlx_df = tlx_df %>%
    dplyr::select(-Seq) %>%
    dplyr::mutate(tlx_id=1:n()) %>%
    dplyr::ungroup()

  #
  # 1. Compare library sizes
  #
  libsizes_df = tlx_df %>%
    dplyr::group_by(tlx_group, tlx_sample, tlx_group_i, tlx_control) %>%
    dplyr::summarize(library_size=n()) %>%
    dplyr::mutate(Treatment=paste0(ifelse(tlx_control, "DMSO (", ""), tlx_group, ifelse(tlx_control, ")", ""))) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(library_factor=max(library_size)/library_size)
  ggplot(libsizes_df) +
    geom_bar(aes(x=tlx_group, y=library_size, fill=Treatment, group=paste0(tlx_group_i, tlx_control)), position="dodge", stat="identity") +
    scale_fill_manual(values=c("APH"="#1F78B4", "DMSO (APH)"="#A6CEE3", "HU"="#E31A1C", "DMSO (HU)"="#FB9A99")) +
    labs(x="", y="Library size", title="Library size comparisons of Hydroxyurea/Aphidicolin treatment and their respective libraries") +
    theme_grey(base_size=18)


  #
  # 2. Compare average distance between junctions
  #
  dist_df = tlx_df %>%
    dplyr::filter(!tlx_is_bait_junction & tlx_is_bait_chromosome) %>%
    dplyr::group_by(tlx_sample, tlx_group, tlx_group_i, tlx_control, Rname) %>%
    dplyr::summarize(dist=diff(sort(Junction))) %>%
    dplyr::ungroup() %>%
    dplyr::filter(dist<5e5) %>%
    dplyr::mutate(Treatment=paste0(ifelse(tlx_control, "DMSO (", ""), tlx_group, ifelse(tlx_control, ")", ""))) %>%
    dplyr::arrange(dplyr::desc(tlx_group), tlx_group_i, dplyr::desc(tlx_control)) %>%
    dplyr::mutate(tlx_sample=factor(tlx_sample, unique(tlx_sample)))

  ggplot(dist_df) +
    ggridges::geom_density_ridges(aes(x=dist, y=tlx_sample, fill=Treatment), bandwidth=1e4) +
    scale_fill_manual(values=c("APH"="#1F78B4", "DMSO (APH)"="#A6CEE3", "HU"="#E31A1C", "DMSO (HU)"="#FB9A99")) +
    labs(y="", title="Average distance between individual junctions") +
    theme_gray(base_size=18)

  #
  # 3. Example of breaksites from APH
  #
  tlxcov_df = tlx_coverage(tlx_df, group="sample", extsize=extsize, exttype="symmetrical") %>%
    dplyr::inner_join(libsizes_df %>% dplyr::select(tlx_sample, library_factor), by="tlx_sample")
  roi_df = readr::read_tsv("data/roi.tsv") %>% dplyr::filter(grepl("Ccser|Ctnna2|Grid2", roi_gene))
  tlxcov_ranges = GenomicRanges::makeGRangesFromDataFrame(tlxcov_df %>% dplyr::mutate(seqnames=tlxcov_chrom, start=tlxcov_start, end=tlxcov_end), keep.extra.columns=T)
  roi_ranges = GenomicRanges::makeGRangesFromDataFrame(roi_df %>% dplyr::mutate(seqnames=roi_chrom, start=roi_start, end=roi_end), keep.extra.columns=T)
  tlxcov_roi_df = as.data.frame(IRanges::mergeByOverlaps(tlxcov_ranges, roi_ranges)) %>%
    dplyr::mutate(sample_name=paste0(ifelse(tlx_control, "DMSO (", ""), tlx_group, ifelse(tlx_control, ")", ""))) %>%
    dplyr::mutate(tlxcov_pileup_norm=tlxcov_pileup*library_factor)
  ggplot() +
    geom_step(aes(x=tlxcov_start, y=tlxcov_pileup, group=paste(tlx_group, tlx_group_i, tlx_control), color=sample_name), data=tlxcov_roi_df %>% dplyr::mutate(Normalization="Unnormalized")) +
    geom_step(aes(x=tlxcov_start, y=tlxcov_pileup_norm, group=paste(tlx_group, tlx_group_i, tlx_control), color=sample_name), data=tlxcov_roi_df %>% dplyr::mutate(Normalization="Library size normalized")) +
    facet_grid(Normalization~roi_gene, scales="free") +
    scale_color_manual(values=c("APH"="#1F78B4", "DMSO (APH)"="#A6CEE3", "HU"="#E31A1C", "DMSO (HU)"="#FB9A99")) +
    theme_grey(base_size=18)

  #
  # 4. PCA samples
  #
  tlx_hist = tlx_df %>%
    dplyr::filter(Rname=="chr6") %>%
    dplyr::inner_join(sizes_df, by=c("Rname"="sizes_chrom")) %>%
    dplyr::group_by(tlx_sample, Rname, sizes_length) %>%
    dplyr::do((function(d){
      dd<<-d
      h = hist(d$Rstart, plot=F, breaks=c(seq(1, d$sizes_length[1], by=extsize), d$sizes_length[1]))
      data.frame(tlx_sample=d$tlx_sample[1], hist_chrom=d$Rname[1], hist_start=h$breaks[-length(h$breaks)], hist_end=h$breaks[-1]-1, hist_count=h$count)
    })(.)) %>%
    dplyr::group_by(tlx_sample) %>%
    dplyr::arrange(dplyr::desc(hist_count)) %>%
    dplyr::mutate(is_top100=1:n()<=50) %>%
    dplyr::group_by(hist_chrom, hist_start, hist_end) %>%
    dplyr::filter(any(is_top100)) %>%
    dplyr::ungroup() %>%
    dplyr::inner_join(libsizes_df %>% dplyr::select(tlx_sample, library_factor), by="tlx_sample") %>%
    dplyr::mutate(hist_count_norm=hist_count*library_factor) %>%
    reshape2::dcast(hist_chrom + hist_start + hist_end ~ tlx_sample, value.var="hist_count_norm")


  tlx_mat = tlx_hist %>%
    dplyr::select(-(hist_chrom:hist_end)) %>%
    t() %>%
    as.data.frame() %>%
    tibble::rownames_to_column("rowname") %>%
    dplyr::mutate(sample=rowname, rowname=sample) %>%
    dplyr::inner_join(samples_df, by="sample") %>%
    tibble::column_to_rownames("rowname") %>%
    dplyr::mutate(group_name=paste0(ifelse(control, "DMSO (", ""), group, ifelse(control, ")", "")))
  pca_res = prcomp(tlx_mat %>% dplyr::select(dplyr::matches("V[0-9]+")), scale.=T)
  ggplot2::autoplot(pca_res, data=tlx_mat, colour="group_name", size=5) +
    ggrepel::geom_text_repel(aes(label=sample), size=8) +
    scale_color_manual(values=c("APH"="#1F78B4", "DMSO (APH)"="#A6CEE3", "HU"="#E31A1C", "DMSO (HU)"="#FB9A99")) +
    labs(title=paste0("PCA of top 50 most abundant bins (", extsize ,") from each sample"), color="Sample") +
    theme_grey(base_size=18)

  #
  # Export data for IGV
  #
  # tlxcov_export_df = tlx_coverage(tlx_df, group="sample", extsize=5e4, exttype="symmetrical")
  # for(s in unique(tlxcov_export_df$tlx_sample)) {
  #   smpl = samples_df %>% dplyr::filter(sample==s)
  #   tlxcov_export_ranges = GenomicRanges::makeGRangesFromDataFrame(tlxcov_export_df %>% dplyr::filter(tlx_sample==smpl$sample) %>% dplyr::select(seqnames=tlxcov_chrom, start=tlxcov_start, end=tlxcov_end, score=tlxcov_pileup), keep.extra.columns=T)
  #   tlxcov_file_basename = paste0(smpl$group, " - ", ifelse(smpl$control, "DMSO", "APH"), "_", smpl$group_i, " - ", smpl$sample, ".bedgraph")
  #   tlxcov_file = paste0("data/pileup/", tlxcov_file_basename)
  #   rtracklayer::export.bedGraph(tlxcov_export_ranges, tlxcov_file)
  # }
  # ggplot()

  #
  # 5. Read MACS2 results
  #
  macs_df = tlx_macs2(tlx_df, grouping="group", effective_size=2*effective_size, extsize=extsize, maxgap=maxgap, exttype=exttype, qvalue=0.01, pileup=threshold_pileup, slocal=slocal, llocal=llocal, exclude_bait_region=T, exclude_repeats=F)
  macs_ranges = GenomicRanges::makeGRangesFromDataFrame(macs_df %>% dplyr::mutate(seqnames=macs_chrom, start=macs_start, end=macs_end), ignore.strand=T)
  macs_reduced_df = as.data.frame(GenomicRanges::reduce(macs_ranges)) %>% dplyr::mutate(macs_chrom=seqnames, macs_start=start, macs_end=end) %>% dplyr::select(-width, -strand)
  macs_reduced_ranges = GenomicRanges::makeGRangesFromDataFrame(macs_reduced_df, ignore.strand=T, keep.extra.columns=T)
  writeLines(macs_df %>% dplyr::mutate(pos=paste0(macs_chrom, ":", macs_start, "-", macs_end)) %>% .$pos)

  tlx_ranges = GenomicRanges::makeGRangesFromDataFrame(tlx_df %>% dplyr::mutate(seqnames=Rname, start=Rstart, end=Rend), ignore.strand=T, keep.extra.columns=T)
  tlx_macs_df = as.data.frame(IRanges::mergeByOverlaps(tlx_ranges, macs_reduced_ranges))
  tlx_macs_df.sum = tlx_macs_df %>%
    dplyr::inner_join(samples_df, by=c("tlx_sample"="sample")) %>%
    dplyr::group_by(tlx_sample, group, group_i, control, macs_chrom, macs_start, macs_end) %>%
    dplyr::summarize(breaks_count=n()) %>%
    dplyr::inner_join(libsizes_df %>% dplyr::select(tlx_sample, library_factor), by="tlx_sample") %>%
    dplyr::mutate(breaks_count_norm=breaks_count*library_factor) %>%
    # dplyr::group_by(tlx_sample, group, control, macs_chrom, macs_start, macs_end) %>%
    # dplyr::summarize(breaks_count_norm.sd=sd(breaks_count_norm), breaks_count_norm=mean(breaks_count_norm)) %>%
    dplyr::mutate(Treatment=paste0(ifelse(control, "DMSO (", ""), group, ifelse(control, ")", "")), macs_pos=paste0(macs_chrom, ":", macs_start, "-", macs_end))
  ggplot(tlx_macs_df.sum) +
    geom_point(aes(x=group, y=breaks_count_norm, color=Treatment),size=2.5, position=position_jitterdodge(jitter.width=0), show.legend=F) +
    scale_fill_manual(values=c("APH"="#1F78B4", "DMSO (APH)"="#A6CEE3", "HU"="#E31A1C", "DMSO (HU)"="#FB9A99")) +
    labs(x="", y="Library size", title="Library size comparisons of Hydroxyurea/Aphidicolin treatment and their respective libraries") +
    theme_grey(base_size=18) +
    facet_wrap(~macs_pos, scales="free")


  #
  # Overlap with RDC
  #
  venn_chromosomes = c("chr6", "All")
  plot.new()
  grid::pushViewport(grid::plotViewport(layout=grid::grid.layout(nrow=2, ncol=1)))
  for(chr in venn_chromosomes) {
    macs_peaks_ranges = GenomicRanges::makeGRangesFromDataFrame(macs_df %>% dplyr::mutate(seqnames=macs_chrom, start=macs_start, end=macs_end), keep.extra.columns=T)
    macs2rdc_df = as.data.frame(IRanges::mergeByOverlaps(macs_peaks_ranges, rdc_pnas_ranges))
    macs2rdc_df = dplyr::bind_rows(macs2rdc_df, rdc_pnas_df %>% dplyr::anti_join(macs2rdc_df %>% dplyr::select(rdc_name)))
    macs2rdc_df = dplyr::bind_rows(macs2rdc_df, macs_df %>% dplyr::anti_join(macs2rdc_df %>% dplyr::select(macs_name)))
    macs2rdc_df = macs2rdc_df %>% dplyr::mutate(common_name=ifelse(!is.na(rdc_name), rdc_name, macs_name))
    if(chr != "All") macs2rdc_df = macs2rdc_df %>% dplyr::filter(rdc_chrom==chr | macs_chrom==chr)

    venn_size = 6
    venn_pallete = "Pastel2"
    venn_list = list(
      RDC=macs2rdc_df %>% dplyr::filter(!is.na(rdc_name)) %>% .$common_name,
      HU=macs2rdc_df %>% dplyr::filter(!is.na(macs_name) & macs_group=="HU") %>% .$common_name,
      APH=macs2rdc_df %>% dplyr::filter(!is.na(macs_name) & macs_group=="APH") %>% .$common_name)

    grid::pushViewport(grid::plotViewport(layout.pos.row=which(chr==venn_chromosomes), layout.pos.col=1))
    p = VennDiagram::venn.diagram(
      x=venn_list, height=venn_size, width=venn_size,
      margin=0.1,
      cat.cex=venn_size/1.5, cex=venn_size, main.cex=venn_size,
      lwd=2, lty='blank', fill=RColorBrewer::brewer.pal(8, venn_pallete)[1:length(venn_list)], cat.fontface="bold", filename=NULL, main=paste0("Considering ", chr, " chromosome(s)"))
    grid::grid.draw(p)
    grid::popViewport()
  }

  "

  #
  # Cluster breaks clusters into larger regions of clusters
  #
  macs_peaks_reduced_df = as.data.frame(GenomicRanges::reduce(GenomicRanges::makeGRangesFromDataFrame(macs_peaks_df %>% dplyr::mutate(macs_start=macs_start-1e5, macs_end=macs_end+1e5)))) %>%
    dplyr::select(reduced_chrom=seqnames, reduced_start=start, reduced_end=end) %>%
    dplyr::mutate(reduced_loci=paste0(reduced_chrom, ":", reduced_start, "-", reduced_end))
  macs_peaks_reduced_ranges = GenomicRanges::makeGRangesFromDataFrame(macs_peaks_reduced_df %>% dplyr::mutate(seqnames=reduced_chrom, start=reduced_start, end=reduced_end), keep.extra.columns=T)
  macs_peaks_ranges = GenomicRanges::makeGRangesFromDataFrame(macs_peaks_df %>% dplyr::mutate(seqnames=macs_chr, start=macs_start, end=macs_end), keep.extra.columns=T)

  tlxcov_ranges = GenomicRanges::makeGRangesFromDataFrame(tlxcov_df %>% dplyr::mutate(seqnames=tlxcov_chrom, start=tlxcov_start, end=tlxcov_end), keep.extra.columns=T)
  tlxcov_reduced_df = as.data.frame(IRanges::mergeByOverlaps(tlxcov_ranges, macs_peaks_reduced_ranges)) %>%
    dplyr::select(-dplyr::matches("macs_peaks_reduced_ranges|tlxcov_ranges")) %>%
    dplyr::mutate(tlxcov_start=ifelse(tlxcov_start<reduced_start, reduced_start-1, tlxcov_start), tlxcov_end=ifelse(tlxcov_end<reduced_start, reduced_start, tlxcov_end)) %>%
    dplyr::mutate(tlxcov_start=ifelse(tlxcov_start>reduced_end, reduced_end, tlxcov_start), tlxcov_end=ifelse(tlxcov_end>reduced_end, reduced_end+1, tlxcov_end))

  #
  # Plot break cluster and pileup
  #
  macs_ggplot_df = as.data.frame(IRanges::mergeByOverlaps(macs_peaks_ranges, macs_peaks_reduced_ranges)) %>%
    dplyr::select(-dplyr::matches("macs_peaks_ranges|macs_peaks_reduced_ranges")) %>%
    dplyr::mutate(Condition=ifelse(macs_control, "DMSO", "HU"))
  macs_ggplot_ymax = 50
  ggplot(macs_ggplot_df) +
    geom_rect(aes(xmin=macs_start, xmax=macs_end, ymin=macs_ggplot_ymax/10*(macs_control+macs_group_i/2), ymax=macs_ggplot_ymax/10*(macs_control+(macs_group_i+1)/2), fill=Condition), alpha=0.5, color="#222222") +
    geom_step(aes(x=tlxcov_start, y=tlxcov_pileup, color=Condition, group=paste0(tlx_sample)), data=tlxcov_reduced_df %>% dplyr::mutate(Condition=ifelse(tlx_control, "DMSO", "HU"))) +
    geom_text(aes(x=macs_end/2+macs_start/2, y=macs_ggplot_ymax/10*(macs_control+macs_group_i/2+1/4), label=macs_pileup)) +
    facet_wrap(~reduced_loci, scales="free_x") +
    scale_color_manual(values=c("DMSO"="#E41A1C", "HU"="#377EB8")) +
    scale_fill_manual(values=c("DMSO"="#E41A1C", "HU"="#377EB8")) +
    guides(color="none") +
    coord_cartesian(ylim=c(0, macs_ggplot_ymax))

  rdc2emily_df = as.data.frame(IRanges::mergeByOverlaps(macs_peaks_ranges, rdc_pnas_ranges)) %>%
    dplyr::select(-dplyr::matches("_ranges\\."))

  venn_size = 6
  venn_pallete = "Pastel2"
  venn_list = list(RDC=c(rdc2emily_df$macs_name, rdc_pnas_df$rdc_name), Emily=macs_peaks_df$macs_name)

  plot.new()
  p = VennDiagram::venn.diagram(
    x=venn_list, height=venn_size, width=venn_size,
    margin=0.1,
    cat.cex=venn_size/1.5, cex=venn_size, main.cex=venn_size,
    lwd=2, lty='blank', fill=RColorBrewer::brewer.pal(8, venn_pallete)[1:length(venn_list)], cat.fontface="bold", filename=NULL)
  grid::grid.draw(p)



  tlxcov_df = tlx_coverage(tlx_df, group="group", extsize=1e6, exttype="symmetrical") %>%
    dplyr::filter(tlx_group=="HU") %>%
    dplyr::mutate(tlxcov_pileup_log10=log10(tlxcov_pileup), tlxcov_pileup_log10=ifelse(tlxcov_pileup_log10>=0, tlxcov_pileup_log10, 0))
  ggplot(tlxcov_df %>% dplyr::filter(tlxcov_chrom=="chr6")) +
    geom_step(aes(tlxcov_start, tlxcov_pileup, color=tlx_control)) +
    geom_rect(aes(ymin=-1, ymax=1, xmin=rdc_start, xmax=rdc_end), data=rdc_pnas_df %>% dplyr::filter(rdc_chrom=="chr6") %>% dplyr::arrange(rdc_start))

  tlx_df.f = tlx_df %>%
    dplyr::filter(tlx_control & Rname=="chr6")
  table(tlx_df$tlx_group, tlx_df$tlx_control)
  libsize.f = tlx_df.f %>%
    dplyr::group_by(tlx_group) %>%
    dplyr::summarize(libsize=n()) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(factor=max(libsize)/libsize)

  tlxcov_df = tlx_coverage(tlx_df.f, group="group", extsize=1.2e6, exttype="symmetrical") %>%
    dplyr::inner_join(libsize.f, by="tlx_group") %>%
    dplyr::mutate(tlxcov_pileup=tlxcov_pileup*factor) %>%
    dplyr::group_by(tlx_group) %>%
    dplyr::mutate(tlxcov_smooth=smoother::smth.gaussian(tlxcov_pileup, window=20)) %>%
    dplyr::ungroup()
  ggplot(tlxcov_df) +
    # geom_step(aes(tlxcov_start, tlxcov_pileup, color=tlx_group)) +
    geom_line(aes(tlxcov_start, tlxcov_smooth, color=tlx_group)) +
    geom_line(aes(x=repliseqTime_start, y=-repliseqTime_avg), data=repliseqTime_df %>% dplyr::filter(repliseqTime_chrom=="chr6")) +
    geom_rect(aes(ymin=-1, ymax=0, xmin=rdc_start, xmax=rdc_end), data=rdc_pnas_df %>% dplyr::filter(rdc_chrom=="chr6") %>% dplyr::arrange(rdc_start))

}