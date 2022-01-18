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
source("circos.R")


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
  color_scheme = c("APH"="#1F78B4", "DMSO (APH)"="#A6CEE3", "HU"="#E31A1C", "DMSO (HU)"="#FB9A99")

  #zz
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
  genes_reduced_ranges = GenomicRanges::reduce(genes_ranges)
  genes_reduced_ranges$range_id = 1:length(genes_reduced_ranges)
  genes_ranges = as.data.frame(IRanges::mergeByOverlaps(genes_ranges, genes_reduced_ranges)) %>%
    dplyr::arrange(dplyr::desc(genes_ranges.width)) %>%
    dplyr::distinct(range_id, .keep_all=T) %>%
    dplyr::mutate(seqnames=genes_reduced_ranges.seqnames, start=genes_reduced_ranges.start, end=genes_reduced_ranges.end) %>%
    dplyr::distinct(seqnames, start, end, gene_name) %>%
    GenomicRanges::makeGRangesFromDataFrame(ignore.strand=T, keep.extra.columns=T)


  #
  # Read RDC
  #
  rdc_pnas_df = readr::read_tsv("data/rdc_pnas.tsv")
  rdc_pnas_ranges = GenomicRanges::makeGRangesFromDataFrame(rdc_pnas_df %>% dplyr::mutate(seqnames=rdc_chrom, start=rdc_start, end=rdc_end), keep.extra.columns=T)
  chain_mm9_mm10 = import.chain("genomes/mm9/mm9ToMm10.over.chain")
  rdc_pnas_ranges = unlist(rtracklayer::liftOver(rdc_pnas_ranges, chain_mm9_mm10))
  rdc_pnas_df = as.data.frame(rdc_pnas_ranges) %>%
    dplyr::distinct(rdc_chrom, rdc_start, rdc_end, .keep_all=T) %>%
    dplyr::mutate(rdc_chrom=seqnames, rdc_start=start, rdc_end=end) %>%
    dplyr::select(dplyr::matches("rdc_")) %>%
    dplyr::mutate(rdc_name=paste(rdc_chrom, rdc_start, rdc_end))
  readr::write_tsv(rdc_pnas_df %>% dplyr::select(-rdc_name), file="data/rdc_pnas_mm10.tsv", quote_escape=F, na="")
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
    dplyr::mutate(Junction1=Junction+200) %>%
    dplyr::mutate(tlx_id=1:n()) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(tlx_sample_num=match(tlx_df$tlx_sample, c("VI021", "VI028", "VI048", "VI054", "VI020", "VI047", "VI053", "LC0024", "JF025", "JF026", "JF027", "JF028")))

  libsizes_df = tlx_df %>%
    dplyr::group_by(tlx_group, tlx_sample, tlx_sample_num, tlx_group_i, tlx_control) %>%
    dplyr::summarize(library_size=n()) %>%
    dplyr::mutate(Treatment=paste0(ifelse(tlx_control, "DMSO (", ""), tlx_group, ifelse(tlx_control, ")", ""))) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(library_factor=max(library_size)/library_size)


  pdf("reports/meeting_2021-10-14.pdf", width=11.69, height=8.27, paper="a4r")
  #
  # 1. Compare library sizes
  #
  ggplot(libsizes_df) +
    geom_bar(aes(x=tlx_group, y=library_size, fill=Treatment, group=paste0(tlx_control, tlx_sample_num, tlx_group_i)), position="dodge", color="#EEEEEE", stat="identity") +
    geom_text(aes(x=tlx_group, y=library_size, group=paste0(tlx_control, tlx_sample_num, tlx_group_i, Treatment), label=tlx_sample), vjust=-1, position=position_dodge(width=0.9), size=4) +
    scale_fill_manual(values=color_scheme) +
    scale_y_continuous(labels=scales::label_number(accuracy=1, scale=1e-3, suffix="K")) +
    labs(x="", y="Library size", title="Library size comparisons of Hydroxyurea/Aphidicolin treatment and their respective libraries") +
    theme_grey(base_size=14) +
    guides(fill=guide_legend(nrow=2, byrow=TRUE)) +
    theme(legend.position="bottom")

  #
  # 2. Compare average distance between junctions
  #
  dist_df = tlx_df %>%
    dplyr::filter(!tlx_is_bait_junction & tlx_is_bait_chromosome) %>%
    dplyr::group_by(tlx_sample, tlx_sample_num, tlx_group, tlx_group_i, tlx_control, Rname) %>%
    dplyr::summarize(dist=diff(sort(Junction))) %>%
    dplyr::ungroup() %>%
    dplyr::inner_join(libsizes_df %>% dplyr::select(tlx_sample, library_factor), by="tlx_sample") %>%
    dplyr::mutate(dist=dist/library_factor) %>%
    dplyr::mutate(Treatment=paste0(ifelse(tlx_control, "DMSO (", ""), tlx_group, ifelse(tlx_control, ")", ""))) %>%
    dplyr::arrange(dplyr::desc(tlx_sample_num)) %>%
    dplyr::mutate(tlx_sample=factor(tlx_sample, unique(tlx_sample)), dist_log10=log10(dist))

  ggplot(dist_df) +
    ggridges::stat_density_ridges(aes(x=dist, y=tlx_sample, fill=Treatment, height=..density..), geom="density_ridges", bandwidth=1e3, scale=1, n=1e4, from=0, to=1e5) +
    scale_fill_manual(values=color_scheme) +
    labs(y="", title="Density of distance between junctions", x="Distance between two neibohring junctions (Kb, sample normalized)") +
    theme_gray(base_size=14) +
    guides(fill=guide_legend(nrow=2, byrow=TRUE)) +
    theme(legend.position="bottom") +
    scale_x_continuous(labels=scales::label_number(accuracy=1, scale=1e-3, suffix="Kb"))

  #
  # 3. Example of breaksites from APH
  #
  # pdf("reports/early_replicating_pileup_2021-10-19.pdf", width=11.69, height=8.27, paper="a4r")
  pileup = c("Late"=extsize, "Early"=1e6)
  roi_df = readr::read_tsv("data/roi.tsv") %>%
    dplyr::filter(grepl("Ccser|Ctnna2|Grid2", roi_gene) | grepl("Early", roi_description)) %>%
    dplyr::mutate(roi_description=ifelse(grepl("Early", roi_description), "Early", "Late")) %>%
    dplyr::mutate(roi_start=ifelse(grepl("Early", roi_description), roi_start-1e6, roi_start)) %>%
    dplyr::mutate(roi_end=ifelse(grepl("Early", roi_description), roi_end+1e6, roi_end))
  roi_ranges = GenomicRanges::makeGRangesFromDataFrame(roi_df %>% dplyr::mutate(seqnames=roi_chrom, start=roi_start, end=roi_end), keep.extra.columns=T)
  for(desc in unique(tlxcov_roi_df$roi_description)) {
    tlxcov_df = tlx_coverage(tlx_df, group="sample", extsize=pileup[desc], exttype="symmetrical") %>%
      dplyr::inner_join(libsizes_df %>% dplyr::select(tlx_sample, library_factor), by="tlx_sample")
    tlxcov_ranges = GenomicRanges::makeGRangesFromDataFrame(tlxcov_df %>% dplyr::mutate(seqnames=tlxcov_chrom, start=tlxcov_start, end=tlxcov_end), keep.extra.columns=T)
    tlxcov_roi_df = as.data.frame(IRanges::mergeByOverlaps(tlxcov_ranges, roi_ranges)) %>%
      dplyr::mutate(Treatment=paste0(ifelse(tlx_control, "DMSO (", ""), tlx_group, ifelse(tlx_control, ")", ""))) %>%
      dplyr::mutate(tlxcov_pileup_norm=tlxcov_pileup*library_factor) %>%
      dplyr::arrange(tlxcov_start) %>%
      dplyr::group_by(roi_gene) %>%
      dplyr::filter(1:n()>rle(tlxcov_pileup)$lengths[1]) %>%
      dplyr::ungroup()

    tlxcov_roi_fdf = tlxcov_roi_df %>% dplyr::filter(roi_description==desc)
    p = ggplot() +
      # geom_step(aes(x=tlxcov_start, y=tlxcov_pileup, group=paste(tlx_group, tlx_group_i, tlx_control), color=Treatment), data=tlxcov_roi_fdf %>% dplyr::mutate(Normalization="Unnormalized")) +
      # geom_step(aes(x=tlxcov_start, y=tlxcov_pileup_norm, group=paste(tlx_group, tlx_group_i, tlx_control), color=Treatment), data=tlxcov_roi_fdf %>% dplyr::mutate(Normalization="Sample normalized")) +
      # facet_grid(Normalization~roi_gene, scales="free") +
      geom_step(aes(x=tlxcov_start, y=tlxcov_pileup_norm, group=paste(tlx_group, tlx_group_i, tlx_control), color=Treatment), data=tlxcov_roi_fdf %>% dplyr::mutate(Chromosome=tlxcov_chrom)) +
      facet_wrap(Chromosome~roi_gene, scales="free_x") +
      labs(title="Examples of known APH clusters", y="Junctions", x="Position on chromosome (Mbp)") +
      scale_color_manual(values=color_scheme) +
      theme_grey(base_size=14) +
      scale_x_continuous(labels=scales::label_number(accuracy=0.1, scale=1e-6, suffix="Mb"), breaks=scales::breaks_width(0.5e6))
    print(p)
  }
  # dev.off()

  #
  # 4. PCA samples
  #
  pca_binsize = 1e6
  tlx_hist_df = tlx_df %>%
    dplyr::filter(Rname=="chr6") %>%
    dplyr::inner_join(sizes_df, by=c("Rname"="sizes_chrom")) %>%
    dplyr::group_by(tlx_sample, Rname, sizes_length) %>%
    dplyr::do((function(d){
      dd<<-d
      h = hist(d$Junction, plot=F, breaks=c(seq(1, d$sizes_length[1], by=pca_binsize), d$sizes_length[1]))
      data.frame(tlx_sample=d$tlx_sample[1], hist_chrom=d$Rname[1], hist_start=h$breaks[-length(h$breaks)], hist_end=h$breaks[-1]-1, hist_count=h$count)
    })(.)) %>%
    dplyr::group_by(tlx_sample) %>%
    dplyr::arrange(dplyr::desc(hist_count)) %>%
    dplyr::mutate(is_top100=1:n()<=50) %>%
    dplyr::group_by(hist_chrom, hist_start, hist_end) %>%
    dplyr::filter(any(is_top100)) %>%
    dplyr::ungroup()


  tlx_hist_ranges = GenomicRanges::makeGRangesFromDataFrame(tlx_hist_df %>% dplyr::mutate(seqnames=Rname, start=hist_start, end=hist_end), ignore.strand=T, keep.extra.columns=T)
  tlx_named_hist_df = leftJoinByOverlaps(tlx_hist_ranges, genes_ranges) %>%
    dplyr::distinct(tlx_sample, hist_chrom, hist_start, hist_end, .keep_all=T) %>%
    dplyr::group_by(tlx_sample, hist_chrom, gene_name) %>%
    dplyr::summarize(hist_count=sum(hist_count))
  tlx_named_hist_wide_df = tlx_named_hist_df %>%
    dplyr::inner_join(libsizes_df %>% dplyr::select(tlx_sample, library_factor), by="tlx_sample") %>%
    dplyr::mutate(hist_count_norm=hist_count*library_factor) %>%
    reshape2::dcast(hist_chrom + gene_name ~ tlx_sample, value.var="hist_count_norm")
  tlx_mat = tlx_named_hist_wide_df %>%
    dplyr::mutate(rowname=paste(gene_name)) %>%
    tibble::column_to_rownames("rowname") %>%
    dplyr::select(-(hist_chrom:gene_name)) %>%
    t() %>%
    as.data.frame() %>%
    tibble::rownames_to_column("rowname") %>%
    dplyr::mutate(sample=rowname, rowname=sample) %>%
    dplyr::inner_join(samples_df, by="sample") %>%
    tibble::column_to_rownames("rowname") %>%
    dplyr::mutate(group_name=paste0(ifelse(control, "DMSO (", ""), group, ifelse(control, ")", "")))
  pca_res = prcomp(tlx_mat %>% dplyr::select(-(sample:group_name)), scale.=T)

  ggplot2::autoplot(pca_res, data=tlx_mat, colour="group_name", size=5, loadings=F, loadings.colour='blue', loadings.label=F, loadings.label.size=1) +
    ggrepel::geom_text_repel(aes(label=sample), size=6) +
    scale_color_manual(values=color_scheme) +
    labs(title=paste0("PCA of top 50 most abundant bins (", scales::number(pca_binsize, scale=1e-6,  suffix="Mbp"), ") from each sample"), color="Sample") +
    theme_grey(base_size=14)
  ggplot2::autoplot(pca_res, data=tlx_mat, colour="group_name", size=5, loadings=T, loadings.colour='blue', loadings.label=T, loadings.label.size=1) +
    ggrepel::geom_text_repel(aes(label=sample), size=6) +
    scale_color_manual(values=color_scheme) +
    labs(title=paste0("PCA of top 50 most abundant bins (", scales::number(pca_binsize, scale=1e-6,  suffix="Mbp"), ") from each sample"), color="Sample") +
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
  macs_df = tlx_macs2(tlx_df, grouping="group", effective_size=effective_size, extsize=extsize, maxgap=maxgap, exttype=exttype, qvalue=threshold_qvalue, pileup=threshold_pileup, slocal=slocal, llocal=llocal, exclude_bait_region=T, exclude_repeats=F)
  macs_ranges = GenomicRanges::makeGRangesFromDataFrame(macs_df %>% dplyr::mutate(seqnames=macs_chrom, start=macs_start, end=macs_end), ignore.strand=T)
  macs_reduced_df = as.data.frame(GenomicRanges::reduce(macs_ranges)) %>% dplyr::mutate(macs_chrom=seqnames, macs_start=start, macs_end=end) %>% dplyr::select(-width, -strand)
  macs_reduced_ranges = GenomicRanges::makeGRangesFromDataFrame(macs_reduced_df, ignore.strand=T, keep.extra.columns=T)
  writeLines(macs_df %>% dplyr::mutate(pos=paste0(macs_chrom, ":", macs_start, "-", macs_end)) %>% .$pos)
  # table(macs_df$macs_group)


  #
  # 7. Overlap with RDC
  #
  venn_chromosomes = c("chr6", "All")
  grid::pushViewport(grid::plotViewport(layout=grid::grid.layout(nrow=1, ncol=2)))
  plot.new()
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
      margin=0,
      cat.cex=venn_size/1.5, cex=venn_size, main.cex=venn_size/1.5,
      lwd=2, lty='blank', fill=RColorBrewer::brewer.pal(8, venn_pallete)[1:length(venn_list)], cat.fontface="bold", main.fontfamily="Helvetica", cat.fontfamily="Helvetica", sub.fontfamily="Helvetica", filename=NULL, main=paste0("Overlap with RDC across\n", chr, " chromosome(s)"),
      cat.dist=rep(-0.03, length(venn_list))
    )
    grid::grid.draw(p)
    grid::popViewport()
  }

  #
  # 6. Plot circos
  #
  cytoband_path = file.path("genomes/mm10/annotation/cytoBand.txt")
  # groups_n = length(unique(tlx_df$tlx_group))
  # layout(matrix(1:groups_n, 1, groups_n))
  for(gr in unique(tlx_df$tlx_group)) {
    tlx_df.gr = tlx_df %>% dplyr::filter(tlx_group==gr)
    macs_df.gr = macs_df %>% dplyr::filter(macs_group==gr)
    baits_df.gr = baits_df %>% dplyr::distinct(bait_group, .keep_all=T) %>% dplyr::filter(bait_group==gr)
    hits_df.gr = macs_df.gr %>% dplyr::select(chrom=macs_chrom, start=macs_start, end=macs_end)
    bait_region.gr = baits_df.gr %>% dplyr::mutate(start=bait_start-bait_region/2, end=bait_end+bait_region/2) %>% dplyr::select(chrom=bait_chrom, start, end)

    links_df.gr = macs_df.gr %>%
      dplyr::inner_join(baits_df.gr, by=c("macs_group"="bait_group")) %>%
      dplyr::mutate(color="#74ADD180") %>%
      dplyr::select(chrom1=macs_chrom, start1=macs_start, end1=macs_end, chrom2=bait_chrom, start2=bait_start, end2=bait_end, color)

    plot_circos(
      input=tlx_df.gr %>% dplyr::filter(!tlx_control) %>% dplyr::select(chrom=Rname, start=Junction, end=Junction),
      control=tlx_df.gr %>% dplyr::filter(tlx_control) %>% dplyr::select(chrom=Rname, start=Junction, end=Junction),
      bait_region=bait_region.gr,
      title=gr,
      annotations=list("hits"=hits_df.gr),
      cytoband_path=cytoband_path,
      chromosomes="chr6",
      links=links_df.gr,
      circos_bw=extsize,
      cex=1.3)
  }

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
    dplyr::arrange(group, control) %>%
    dplyr::mutate(Treatment=factor(Treatment, unique(Treatment))) %>%
    dplyr::mutate(Facet=paste(gene_name, "\n", macs_pos)) %>%
    dplyr::group_by(Facet) %>%
    dplyr::mutate(breaks_count_norm.max=max(breaks_count_norm)) %>%
    dplyr::mutate(FacetGroup=dplyr::case_when(
      breaks_count_norm.max<=20~"<=20 breaks",
      breaks_count_norm.max<=50~"<=50 breaks",
      # breaks_count_norm.max<=100~"<=100 breaks",
      breaks_count_norm.max<=200~"<=200 breaks",
      breaks_count_norm.max>200~">200 breaks",
    )) %>%
    dplyr::group_by(FacetGroup) %>%
    dplyr::mutate(FacetGroupCount=n()) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(breaks_count_norm.max) %>%
    dplyr::mutate(FacetGroup=factor(FacetGroup, unique(FacetGroup)))


  plist = lapply(split(tlx_macs_df.sum, f=tlx_macs_df.sum$FacetGroup), FUN=function(df) {
    g = ggplot(df) +
      geom_boxplot(aes(x=group, y=breaks_count_norm, fill=Treatment), position=position_dodge2(preserve="single")) +
      geom_point(aes(x=group, y=breaks_count_norm, fill=Treatment), size=1.8, pch=21, position=position_jitterdodge(jitter.width=0.2), show.legend=F) +
      scale_color_manual(values=color_scheme, drop=F) +
      scale_fill_manual(values=color_scheme, drop=F) +
      labs(x=NULL, y=NULL, title=df$FacetGroup[1]) +
      theme_grey(base_size=10) +
      facet_wrap(~Facet, nrow=5) +
      scale_x_discrete(drop=F) +
      guides(fill=guide_legend(ncol=1, byrow=TRUE)) +
      theme(legend.position="none", strip.text=element_text(size=7), axis.text=element_text(size=7.5))
    g
    # ggplot_gtable(ggplot_build(g))

  })

  plist[[1]] = plist[[1]] +labs(y="Junctions (sample normalized)")
  plist[[length(plist)+1]] = cowplot::get_legend(plist[[1]] + theme(legend.position="left", legend.box.margin = margin(0, 0, 0, 12)))
  plist_widths = c(2*ceiling(sapply(plist[-length(plist)], function(d) length(unique(d$data$Facet)))/5), 1.5)
  cowplot::plot_grid(plotlist=plist, align="v", axis="t", nrow=1, rel_widths=plist_widths)

  dev.off()
}
