darker_colors = function(colors, l=0.3) {
  cols1 = colorspace::readhex(file=textConnection(paste(colors, collapse = "\n")), class="RGB")
  cols1 = as(cols1, "HLS")
  cols1@coords[, "L"] = pmax(0, cols1@coords[, "L"] - l)
  cols1 = colorspace::hex(as(cols1, "RGB"))
  names(cols1) = names(colors)
  cols1
}

plot_circos = function(input, control, title, cytoband_path, chromosomes, bait_region=NULL, annotations=NULL, links=NULL, circos_bw=1e6-1, cex=5, colors=c(neutral="#999999", input="#FFCB00", control="#FF5700")) {
  cytoband = circlize::read.cytoband(cytoband_path)
  cytoband_df = data.frame(cytoband$chr.len) %>% tibble::rownames_to_column("chrom") %>% dplyr::rename(chrom_length="cytoband.chr.len")

  if(!is.null(control)) {
    has_control = nrow(control)>0
  } else {
    has_control = F
  }

  unknown_chroms = list("data (chrom1)"=unique(setdiff(unique(c(input$chrom, control$chrom)), cytoband_df$chrom)))

  if(!is.null(annotations)) {
    if(is.data.frame(annotations)) {
      annotations = list("annotations"=annotations)
    }

    for(n in names(annotations)) {
      unknown_chroms[[paste0("annotations (", n, ")")]] = unique(setdiff(annotations[[n]]$chrom, cytoband_df$chrom))
    }
  }
  if(!is.null(links)) {
    unknown_chroms[["links (chrom1)"]] = unique(setdiff(links$chrom1, cytoband_df$chrom))
    unknown_chroms[["links (chrom2)"]] = unique(setdiff(links$chrom2, cytoband_df$chrom))
  }

  if(sum(sapply(unknown_chroms, length))>0) {
    unknown_chroms = sapply(unknown_chroms[sapply(unknown_chroms, length)>0], paste, collapse=",")
    unknown_chroms_err = paste("Unknown chromosomes provided when ploting circos plot...\n", paste(paste0(names(unknown_chroms), ":", unknown_chroms), collapse="\n"))
    stop(unknown_chroms_err)
  }


  if(has_control) {
    scale = c(input=1, control=nrow(input)/nrow(control))
    data_sum = rbind(input %>% dplyr::mutate(circos_signal="input"), control %>% dplyr::mutate(circos_signal="control"))
  } else {
    scale = c(input=1)
    data_sum = input %>% dplyr::mutate(circos_signal="input")
  }

  data_sum = data_sum %>%
    dplyr::inner_join(cytoband_df, by="chrom") %>%
    dplyr::group_by(chrom) %>%
    dplyr::mutate(chrom_length=max(c(chrom_length, end))) %>%
    dplyr::group_by(circos_signal, chrom_length, chrom) %>%
    dplyr::do((function(d){
      dd<<-d
      h = hist(d$start, plot=F, breaks=c(seq(1, d$chrom_length[1], by=circos_bw), d$chrom_length[1]))
      data.frame(start=h$breaks[-length(h$breaks)], end=h$breaks[-1]-1, count=scale[d$circos_signal[1]]*h$count, count_log10=ifelse(h$count>0, log10(h$count), 0))
    })(.)) %>%
    dplyr::mutate(start=pmax(1, start), end=pmin(end, chrom_length-1)) %>%
    dplyr::mutate(circos_col=paste0("count_log10.", circos_signal), circos_chrom=chrom) %>%
    dplyr::filter(end<=chrom_length) %>%
    reshape2::dcast(chrom + start + end + chrom_length + circos_chrom ~ circos_col, value.var="count_log10")


  if(has_control) {
    circos_ymax = max(c(data_sum$count_log10.input, data_sum$count_log10.control))
  } else {
    circos_ymax = max(data_sum$count_log10.input)
  }

  circos_ylim = c(0, circos_ymax)
  circos_ylabels = sort(expand.grid(power=10^seq(0, ceiling(circos_ylim[2]), 1), prec=c(1, 2, 5)) %>% dplyr::mutate(y=power*prec) %>% .$y)
  circos_ylabels = circos_ylabels[circos_ylabels<=10^circos_ymax]
  circos_yaxis = log10(circos_ylabels)
  circos_yaxis_pal = circlize::colorRamp2(circos_yaxis, colorRampPalette(rev(RColorBrewer::brewer.pal(5, "Blues")))(length(circos_yaxis)), transparency=0.8)

  colors_darker = darker_colors(colors, 0.2)
  # barplot(rep(1,3), col=colors_darker)
  # barplot(rep(1,3), col=colors)

  circlize::circos.par("gap.degree"=c(rep(1, length(chromosomes)-1), 5))
  par(cex=cex, cex.main=cex)
  circlize::circos.initializeWithIdeogram(cytoband=cytoband_path, chromosome.index=chromosomes, plotType=c("axis", "labels"), axis.labels.cex=cex*0.5)
  circlize::circos.genomicTrack(data_sum, bg.border=NA, ylim=circos_ylim,
      panel.fun = function(region, value, ...) {
        rr<<-region
        vv<<-value
        if(circlize::get.current.chromosome() == intersect(cytoband$chromosome, chromosomes)[1]) {
          circlize::circos.yaxis(at=circos_yaxis[circos_ylabels>=2], labels=circos_ylabels[circos_ylabels>=2], labels.cex=cex*0.22)
        }
        if(length(circos_yaxis)>1) {
          for(ax in 2:length(circos_yaxis)) {
            circlize::circos.rect(xleft=0, xright=cytoband$chr.len[circlize::get.current.chromosome()], ybottom=circos_yaxis[ax-1], ytop=circos_yaxis[ax], col=circos_yaxis_pal(circos_yaxis[ax-1]), border="#00000000")
          }
        }

        if(!is.null(bait_region) & bait_region$chrom==value$circos_chrom[1]) {
          bait_breaks = sapply(rr$start, function(z) { any(z>=bait_region$start) }) & sapply(rr$end, function(z) { any(z<=bait_region$end) })
          colors_neutral = sapply(bait_breaks, function(z) { ifelse(z, colors_darker["neutral"], colors["neutral"]) })
          colors_input = sapply(bait_breaks, function(z) { ifelse(z, colors_darker["input"], colors["input"]) })
          colors_control = sapply(bait_breaks, function(z) { ifelse(z, colors_darker["control"], colors["control"]) })
        } else {
          colors_neutral = colors["neutral"]
          colors_input = colors["input"]
          colors_control = colors["control"]
        }

        if(has_control) {
          count_log10.pmin = pmin(value$count_log10.input, value$count_log10.control)
          circlize::circos.rect(xleft=region$start, xright=region$end+1, ybottom=0, ytop=count_log10.pmin, col=colors_neutral, border=NA)

          f = value$count_log10.input>count_log10.pmin
          if(any(f)) circlize::circos.rect(xleft=region$start[f], xright=region$end[f]+1, ybottom=count_log10.pmin[f], ytop=value$count_log10.input[f], col=colors_input[f], border=NA)

          f = value$count_log10.control>count_log10.pmin
          if(any(f)) circlize::circos.rect(xleft=region$start[f], xright=region$end[f]+1, ybottom=count_log10.pmin[f], ytop=value$count_log10.control[f], col=colors_control[f], border=NA)
        } else {
          circlize::circos.rect(xleft=region$start, xright=region$end+1, ybottom=0, ytop=value$count_log10.input, col=colors_neutral, border=NA)
        }
  })

  if(!is.null(annotations) && length(annotations)>0) {
    for(n in names(annotations)) {
      circlize::circos.genomicTrack(annotations[n], bg.border=NA, ylim=c(0,1), track.height=0.04, cell.padding=c(0,0),
          panel.fun = function(region, value, ...) {
            circlize::circos.rect(xleft=region$start, xright=region$end, ybottom=0, ytop=1, col="#330000", border="#330000")
      })
    }
  }

  if(!is.null(links)) {
    if(!("color" %in% colnames(links))) {
      links$color = "#F46D4380"
    }


    links_sum = links %>%
      dplyr::filter(chrom1 %in% chromosomes & chrom2 %in% chromosomes) %>%
      dplyr::inner_join(cytoband_df %>% setNames(., paste0(colnames(.), "1")), by="chrom1")  %>%
      dplyr::inner_join(cytoband_df %>% setNames(., paste0(colnames(.), "2")), by="chrom2")  %>%
      dplyr::rowwise() %>%
      dplyr::do((function(d){
        dd<<-d
        breaks = seq(1, d$chrom_length1[1], by=circos_bw)
        data.frame(
          chrom1=d$chrom1,
          start1=breaks[which(breaks>d$start1)[1]-1],
          end1=breaks[rev(which(breaks<d$end1))[1]+1],
          chrom2=d$chrom2,
          start2=breaks[which(breaks>d$start2)[1]-1],
          end2=breaks[rev(which(breaks<d$end2))[1]+1],
          color=d$color
        )
      })(.))

    if(nrow(links_sum) > 0) {
      circlize::circos.genomicLink(
        region1=links_sum %>% dplyr::select(chr=chrom1, start=start1, end=end1),
        region2=links_sum %>% dplyr::select(chr=chrom2, start=start2, end=end2),
        col=links_sum$color, border=NA)
    }
  }
  # title(title)

  # mtext(title, side=3, adj=0, line=1.2, cex=2, font=2)
  mtext(title, side=2, padj=1, adj=0.5, cex=cex*2)
  # @ TODO change legend size
  legend("right", title="Reads", legend=names(colors), fill=colors, xjust=1, yjust=1, cex=cex/2)

  circlize::circos.clear()
}