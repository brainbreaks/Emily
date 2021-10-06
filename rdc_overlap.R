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

  pdf("reports/meeting_2021-10-06.pdf", width=11.69, height=8.27, paper="a4r")

  #
  # Offtargets
  #
  # bowtie2 -k 30 -N 1 -x mm10/mm10 --end-to-end --very-sensitive -c acgagcatttccaaccc

  #
  # Chromosomes sizes
  #
  sizes_df = readr::read_tsv("genomes/mm10/annotation/mm10.chrom.sizes", col_names=c("sizes_chrom", "sizes_length")) %>%
    dplyr::mutate(sizes_effective=sizes_length/sum(sizes_length)*effective_size)

  genes_ranges = rtracklayer::import("genomes/mm10/annotation/refGene.bed")
  genes_ranges = genes_ranges[!grepl("_rev",genes_ranges$name)]
  values(genes_ranges) = data.frame(gene_name=genes_ranges$name)

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
    scale_fill_manual(values=color_scheme) +
    labs(x="", y="Library size", title="Library size comparisons of Hydroxyurea/Aphidicolin treatment and their respective libraries") +
    theme_grey(base_size=14) +
    guides(fill=guide_legend(nrow=2, byrow=TRUE)) +
    theme(legend.position="bottom")


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
    scale_fill_manual(values=color_scheme) +
    labs(y="", title="Average distance between individual junctions") +
    theme_gray(base_size=14) +
    guides(fill=guide_legend(nrow=2, byrow=TRUE)) +
    theme(legend.position="bottom")

  #
  # 3. Example of breaksites from APH
  #
  tlxcov_df = tlx_coverage(tlx_df, group="sample", extsize=extsize, exttype="symmetrical") %>%
    dplyr::inner_join(libsizes_df %>% dplyr::select(tlx_sample, library_factor), by="tlx_sample")
  roi_df = readr::read_tsv("data/roi.tsv") %>% dplyr::filter(grepl("Ccser|Ctnna2|Grid2", roi_gene))
  tlxcov_ranges = GenomicRanges::makeGRangesFromDataFrame(tlxcov_df %>% dplyr::mutate(seqnames=tlxcov_chrom, start=tlxcov_start, end=tlxcov_end), keep.extra.columns=T)
  roi_ranges = GenomicRanges::makeGRangesFromDataFrame(roi_df %>% dplyr::mutate(seqnames=roi_chrom, start=roi_start, end=roi_end), keep.extra.columns=T)
  tlxcov_roi_df = as.data.frame(IRanges::mergeByOverlaps(tlxcov_ranges, roi_ranges)) %>%
    dplyr::mutate(Treatment=paste0(ifelse(tlx_control, "DMSO (", ""), tlx_group, ifelse(tlx_control, ")", ""))) %>%
    dplyr::mutate(tlxcov_pileup_norm=tlxcov_pileup*library_factor)
  ggplot() +
    geom_step(aes(x=tlxcov_start/1e6, y=tlxcov_pileup, group=paste(tlx_group, tlx_group_i, tlx_control), color=Treatment), data=tlxcov_roi_df %>% dplyr::mutate(Normalization="Unnormalized")) +
    geom_step(aes(x=tlxcov_start/1e6, y=tlxcov_pileup_norm, group=paste(tlx_group, tlx_group_i, tlx_control), color=Treatment), data=tlxcov_roi_df %>% dplyr::mutate(Normalization="Library size normalized")) +
    facet_grid(Normalization~roi_gene, scales="free") +
    labs(title="Examples of known APH clusters", y="Reads", x="Mbp") +
    scale_color_manual(values=color_scheme) +
    theme_grey(base_size=14)

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
    ggrepel::geom_text_repel(aes(label=sample), size=6) +
    scale_color_manual(values=color_scheme) +
    labs(title=paste0("PCA of top 50 most abundant bins (", extsize ,") from each sample"), color="Sample") +
    theme_grey(base_size=14)

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
  # 5. Compare breaks for MACS2 hits
  #
  macs_df = tlx_macs2(tlx_df, grouping="group", effective_size=effective_size, extsize=extsize, maxgap=maxgap, exttype=exttype, qvalue=0.01, pileup=threshold_pileup, slocal=slocal, llocal=llocal, exclude_bait_region=T, exclude_repeats=F)
  macs_ranges = GenomicRanges::makeGRangesFromDataFrame(macs_df %>% dplyr::mutate(seqnames=macs_chrom, start=macs_start, end=macs_end), ignore.strand=T)
  macs_reduced_df = as.data.frame(GenomicRanges::reduce(macs_ranges)) %>% dplyr::mutate(macs_chrom=seqnames, macs_start=start, macs_end=end) %>% dplyr::select(-width, -strand)
  macs_reduced_ranges = GenomicRanges::makeGRangesFromDataFrame(macs_reduced_df, ignore.strand=T, keep.extra.columns=T)
  writeLines(macs_df %>% dplyr::mutate(pos=paste0(macs_chrom, ":", macs_start, "-", macs_end)) %>% .$pos)
  # table(macs_df$macs_group)

  tlx_ranges = GenomicRanges::makeGRangesFromDataFrame(tlx_df %>% dplyr::mutate(seqnames=Rname, start=Rstart, end=Rend), ignore.strand=T, keep.extra.columns=T)
  tlx_macs_df = as.data.frame(IRanges::mergeByOverlaps(tlx_ranges, macs_reduced_ranges)) %>% dplyr::select(-dplyr::matches("_ranges\\."))
  tlx_macs_ranges = GenomicRanges::makeGRangesFromDataFrame(tlx_macs_df %>%  dplyr::mutate(seqnames=macs_chrom, start=macs_start, end=macs_end), keep.extra.columns=T, ignore.strand=T)
  tlx_named_macs_df = leftJoinByOverlaps(tlx_macs_ranges, genes_ranges) %>%
    dplyr::distinct(tlx_sample, macs_chrom, macs_start, macs_end, Rname, Junction, .keep_all=T) %>%
    dplyr::inner_join(samples_df, by=c("tlx_sample"="sample"))

  tlx_macs_df.sum1 = tlx_named_macs_df %>%
    dplyr::group_by(tlx_sample, macs_chrom, macs_start, macs_end, group, gene_name, control) %>%
    dplyr::summarize(breaks_count=n()) %>%
    dplyr::ungroup()
  tlx_macs_df.sum2 = tlx_macs_df.sum1 %>%
    dplyr::distinct(macs_chrom, macs_start, macs_end, gene_name) %>%
    tidyr::crossing(tlx_macs_df.sum1 %>% dplyr::distinct(tlx_sample, group, control)) %>%
    dplyr::mutate(breaks_count=0) %>%
    dplyr::anti_join(tlx_macs_df.sum1 %>% dplyr::select(-breaks_count))
  tlx_macs_df.sum = dplyr::bind_rows(tlx_macs_df.sum1, tlx_macs_df.sum2) %>%
    dplyr::group_by(macs_chrom, macs_start, macs_end) %>%
    dplyr::filter(any(breaks_count>0)) %>%
    dplyr::ungroup() %>%
    dplyr::inner_join(libsizes_df %>% dplyr::select(tlx_sample, library_factor), by="tlx_sample") %>%
    dplyr::mutate(breaks_count_norm=breaks_count*library_factor) %>%
    # dplyr::group_by(tlx_sample, group, control, macs_chrom, macs_start, macs_end) %>%
    # dplyr::summarize(breaks_count_norm.sd=sd(breaks_count_norm), breaks_count_norm=mean(breaks_count_norm)) %>%
    dplyr::mutate(Treatment=paste0(ifelse(control, "DMSO (", ""), group, ifelse(control, ")", "")), macs_pos=paste0(macs_chrom, ":", macs_start, "-", macs_end)) %>%
    dplyr::mutate(Treatment=factor(Treatment)) %>%
    dplyr::mutate(Facet=paste(gene_name, "|", macs_pos))

  ggplot(tlx_macs_df.sum) +
    geom_boxplot(aes(x=group, y=breaks_count_norm, fill=Treatment), position=position_dodge2(preserve="single")) +
    geom_point(aes(x=group, y=breaks_count_norm, fill=Treatment), size=1.8, pch=21, position=position_jitterdodge(jitter.width=0.2), show.legend=F) +
    scale_color_manual(values=color_scheme, drop=F) +
    scale_fill_manual(values=color_scheme, drop=F) +
    labs(x="", y="Normalized junctions count", title="Number of junctions in macs hits comparison") +
    theme_grey(base_size=14) +
    facet_wrap(~Facet, scales="free_y") +
    scale_x_discrete(drop=F) +
    guides(fill=guide_legend(nrow=2, byrow=TRUE)) +
    theme(legend.position="bottom", strip.text=element_text(size=7), axis.text = element_text(size=7))


  #
  # 6. Overlap with RDC
  #
  venn_chromosomes = c("chr6", "All")
  plot.new()
  grid::pushViewport(grid::plotViewport(layout=grid::grid.layout(nrow=1, ncol=2)))
  for(chr in venn_chromosomes) {
    macs_peaks_ranges = GenomicRanges::makeGRangesFromDataFrame(macs_df %>% dplyr::mutate(seqnames=macs_chrom, start=macs_start, end=macs_end), keep.extra.columns=T)
    macs2rdc_df = as.data.frame(IRanges::mergeByOverlaps(macs_peaks_ranges, rdc_pnas_ranges))
    macs2rdc_df = dplyr::bind_rows(macs2rdc_df, rdc_pnas_df %>% dplyr::anti_join(macs2rdc_df %>% dplyr::select(rdc_name)))
    macs2rdc_df = dplyr::bind_rows(macs2rdc_df, macs_df %>% dplyr::anti_join(macs2rdc_df %>% dplyr::select(macs_name)))
    macs2rdc_df = macs2rdc_df %>% dplyr::mutate(common_name=ifelse(!is.na(rdc_name), rdc_name, macs_name))
    if(chr != "All") macs2rdc_df = macs2rdc_df %>% dplyr::filter(rdc_chrom==chr | macs_chrom==chr)

    venn_size = 2
    venn_pallete = "Pastel2"
    venn_list = list(
      RDC=macs2rdc_df %>% dplyr::filter(!is.na(rdc_name)) %>% .$common_name,
      HU=macs2rdc_df %>% dplyr::filter(!is.na(macs_name) & macs_group=="HU") %>% .$common_name,
      APH=macs2rdc_df %>% dplyr::filter(!is.na(macs_name) & macs_group=="APH") %>% .$common_name)

    grid::pushViewport(grid::plotViewport(layout.pos.col=which(chr==venn_chromosomes), layout.pos.row=1))
    p = VennDiagram::venn.diagram(
      x=venn_list, height=venn_size, width=venn_size,
      margin=0.1,
      cat.cex=venn_size/1.5, cex=venn_size, main.cex=venn_size,
      lwd=2, lty='blank', fill=RColorBrewer::brewer.pal(8, venn_pallete)[1:length(venn_list)], cat.fontface="bold", filename=NULL, main=paste0("Considering ", chr, " chromosome(s)"))
    grid::grid.draw(p)
    grid::popViewport()
  }
  dev.off()
}