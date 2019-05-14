require(cowplot)
require(data.table)
require(ggplot2)

#' ggkaryo: a class for karyotype plotting and overlaying.
#'
#' ggkaryo full description
#'

ggkaryo <- setRefClass("ggkaryo",
  field = list(
    n_chrom = "numeric",
    hetero = "character",
    chrom_width = "numeric",
    chrom_padding = "numeric",
    giemsa_palette = "character",
    giemsa_levels = "character",
    opposite = "logical",
    data = "list"
  ),

  method = list(
    initialize = function(giemsa, ...,
        n_chrom=24, hetero=c("X", "Y"),
        chrom_width=1, chrom_padding=5,
        giemsa_palette=c(
          "#DDDDDD", "#9A9A9A", "#787878", "#555555", "#333333",
          "#FF0000", "#C4FFFC", "#AFE6FF"),
        giemsa_levels=c(
          "gneg", "gpos25", "gpos50", "gpos75", "gpos100",
          "acen", "gvar", "stalk"),
        opposite=F
      ) {
      "@param giemsa (character) path to Giemsa BED5+ file.
      @param giemsa (data.table) data table with Giemsa BED5+ file info."
      stopifnot(length(giemsa_levels) == length(giemsa_palette))
      stopifnot(chrom_width > 0)
      stopifnot(chrom_padding >= chrom_width)

      callSuper(...,
        n_chrom=n_chrom, hetero=hetero,
        chrom_width=chrom_width, chrom_padding=chrom_padding,
        giemsa_palette=giemsa_palette, giemsa_levels=giemsa_levels,
        opposite=opposite, data=list(tracks=list()))
      names(.self$giemsa_palette) = giemsa_levels

      .self$prep4karyo(giemsa)
    },

    chrom2id = function(chrom) {
      "Converts a chromosome signature (seqname) to a numerical id.
      @param chrom (string) chromosome signature (e.g., 'chr1' or '1')
      @returns numeric: chromosome numerical ID"
      if ( grepl("^chr", chrom) ) chrom = gsub("^chr", "", chrom)
      if ( grepl(":", chrom) ) {
        return(floor(as.numeric(gsub(":", ".", chrom))))
      } else {
        if ( chrom %in% hetero ) {
          return(.self$n_chrom-length(.self$hetero)+which(.self$hetero==chrom))
        } else {
          return(as.numeric(chrom))
        }
      }
    },
    chromID2x = function(chromID) {
      "Retrieve the position of a chromosome on the X axis.
      @param chromID (numeric)
      @returns numeric: chromosome position on the X axis"
      return((chromID-1)*(.self$chrom_width + .self$chrom_padding))
    },

    read_giemsa = function(giemsa) {
      "Reads a Giemsa bed file. Adds chromID, bandID, and X columns.
      @param giemsa (character) path to Giemsa BED5 file.
      @param giemsa (data.table) data table with Giemsa BED5 file info."
      if ( is(giemsa, "character") ) {
        stopifnot(file.exists(giemsa))
        giemsa = fread(giemsa)
      }
      stopifnot(is(giemsa, "data.table"))
      stopifnot(ncol(giemsa) >= 5)
      giemsa = giemsa[, 1:5]
      colnames(giemsa) = c("chrom", "start", "end", "name", "value")
      giemsa[, chromID := unlist(lapply(chrom, .self$chrom2id))]
      giemsa[, x := unlist(lapply(chromID, .self$chromID2x))]
      giemsa[, bandID := paste0(chrom, ":", start, "-", end)]
      return(giemsa)
    },
    prep4bands = function() {
      "Prepares data for plotting chromosome bands.
      Builds .self$data[['bands']] object."

      stopifnot("giemsa" %in% names(.self$data))
      stopifnot(is(.self$data[['giemsa']], "data.table"))

      non_acen_bands = data.table(
        chrom = rep(.self$data[["giemsa"]]$chrom, each = 4),
        chromID = rep(.self$data[["giemsa"]]$chromID, each = 4),
        y = c(t(cbind(
          .self$data[["giemsa"]]$start,
          .self$data[["giemsa"]]$start,
          .self$data[["giemsa"]]$end,
          .self$data[["giemsa"]]$end))),
        x = c(t(cbind(
          .self$data[["giemsa"]]$x,
          .self$data[["giemsa"]]$x+.self$chrom_width,
          .self$data[["giemsa"]]$x+.self$chrom_width,
          .self$data[["giemsa"]]$x))),
        value = rep(.self$data[["giemsa"]]$value, each = 4),
        bandID = rep(.self$data[["giemsa"]]$bandID, each = 4)
      )
      non_acen_bands = non_acen_bands[non_acen_bands$value != "acen",]

      acen_data = .self$data[["giemsa"]][value == "acen", .(
          start = min(start), end = max(end), name = NA, value = "acen",
          x = x[1], bandID = bandID[1]
        ), by = c("chrom", "chromID")]
      acen_bands = data.table(
        chrom = rep(acen_data$chrom, each = 4),
        chromID = rep(acen_data$chromID, each = 4),
        y = c(t(cbind(acen_data$start, acen_data$start,
          acen_data$end, acen_data$end))),
        x = c(t(cbind(acen_data$x, acen_data$x+.self$chrom_width,
          acen_data$x, acen_data$x+.self$chrom_width))),
        value = rep(acen_data$value, each = 4),
        bandID = rep(acen_data$bandID, each = 4)
      )

      .self$data[["bands"]] = rbind(non_acen_bands, acen_bands)
      .self$data[["bands"]][, value := factor(value,
        levels = .self$giemsa_levels)]
      NULL
    },
    prep4boxes = function() {
      "Prepares data for plotting chromosome arm boxes. Chromosome arms are
      identified based on the 'acen' bands that are used to divide each
      chromosomes in two arms. Builds .self$data[['boxes']] object."

      stopifnot("giemsa" %in% names(.self$data))
      stopifnot(is(.self$data[['giemsa']], "data.table"))

      select_chrom_arms = function(chrom_data) {
        chrom_x = .self$chromID2x(chrom_data[1, chromID])
        acen_band_ids = which(chrom_data$value == "acen")
        if ( ! 1 %in% acen_band_ids ) {
          p_arm_data = chrom_data[1:(min(acen_band_ids)-1), .(
              x = c(chrom_x, chrom_x,
                chrom_x+.self$chrom_width, chrom_x+.self$chrom_width, chrom_x),
              y = c(min(start), max(end), max(end), min(start), min(start)),
              arm_id = "p"
            )]
        } else {
          p_arm_data = NULL
        }
        if ( ! nrow(chrom_data) %in% acen_band_ids ) {
          q_arm_data = chrom_data[(max(acen_band_ids)+1):nrow(chrom_data), .(
              x = c(chrom_x, chrom_x,
                chrom_x+.self$chrom_width, chrom_x+.self$chrom_width, chrom_x),
              y = c(min(start), max(end), max(end), min(start), min(start)),
              arm_id = "q"
            )]
        } else {
          q_arm_data = NULL
        }
        return(rbind(q_arm_data, p_arm_data))
      }
      .self$data[["boxes"]] = .self$data[["giemsa"]][,
        select_chrom_arms(.SD), by="chrom"]
      NULL
    },
    prep4labels = function() {
      "Prepares data for plotting chromosome labels.
      Builds .self$data[['chrom_labels']] object."

      stopifnot("giemsa" %in% names(.self$data))
      stopifnot(is(.self$data[['giemsa']], "data.table"))

      .self$data[['chrom_labels']] = .self$data[["giemsa"]][, .(
          x = min(x) + .self$chrom_width/2,
          y = -5 * 10**(ceiling(abs(log10(max(.self$data[["giemsa"]]$end))))-3)
        ), by = c("chrom", "chromID")]
      NULL
    },
    prep4plot = function() {
      "Prepares for plotting. Builds .self$data[['plot']] object."

      stopifnot("bands" %in% names(.self$data))
      stopifnot(is(.self$data[['bands']], "data.table"))

      .self$data[['plot']] = ggplot(.self$data[['bands']], aes(x=x, y=-y/10)
        ) + geom_polygon(aes(fill=value, group=bandID)
        ) + scale_fill_manual(values=.self$giemsa_palette
        ) + theme(axis.line = element_blank(), axis.ticks = element_blank(),
          axis.title = element_blank(), axis.text = element_blank()
        ) + guides(fill = F)
      NULL
    },
    prep4karyo = function(giemsa) {
      "Builds a data.table to plot the ideograms and (optionally) boxes around
      each chromosome arm.
      @param giemsa (character) path to Giemsa BED5 file.
      @param giemsa (data.table) data table with Giemsa BED5 file info."
      .self$data[["giemsa"]] = .self$read_giemsa(giemsa)
      .self$prep4bands()
      .self$prep4boxes()
      .self$prep4labels()
      .self$prep4plot()
      NULL
    },

    bin_track = function(track, size, step, method="within",
      fun.aggreg=mean, ...) {
      "Bins a track based on provided bin size and step.
      Regions from the track are assigned to the bins when they are completely
      include ('within' method) or overlap even partially ('overlap' method).
      @param track (data.table) BED5+ track data table
      @param size (numeric) bin size in nt
      @param step (numeric) bin step in nt
      @param method (string) either 'within' or 'overlap'
      @param fun.aggreg (function) how to aggregate values in bins
      @param ... (mixed) additional parameters to pass to fun.aggreg"
      stopifnot(is(track, "data.table"))
      stopifnot(ncol(track) >= 5)
      stopifnot(method %in% c("within", "overlap"))
      track = track[, 1:5]
      colnames(track) = c("chrom", "start", "end", "name", "value")
      track[, chromID := unlist(lapply(chrom, .self$chrom2id))]

      mk_bins = function(chrom_data, size, step) {
        "Generates a list of bins of given size and step."
        end = max(chrom_data$end, na.rm = T)
        starts = seq(0, end-step, by=step)
        data.table(start = starts, end = starts+size, value = 0)
      }
      select_overlap = function(data, start, end)
        data$start > start | data$end <= end
      select_within = function(data, start, end)
        data$start > start & data$end <= end
      bin_chrom = function(chrom_data, size, step, method,
        fun.aggreg=mean, ...) {
        "Bin chromosome data using given size and step."
        select_regions = ifelse("within"==method, select_within, select_overlap)
        chrom_data = chrom_data[order(start)]
        bins = mk_bins(chrom_data, size, step)
        for ( bi in 1:nrow(bins) ) {
          ri = which(select_regions(chrom_data, bins[bi, start], bins[bi, end]))
          bins[bi, value := fun.aggreg(chrom_data[ri, value], ...)]
        }
        return(bins)
      }
      track[, bin_chrom(.SD, size, step, method, fun.aggreg, na.rm=T),
        by=c("chrom", "chromID")][,
        .(chrom, start, end, paste0("bin_", 1:.N), value)]
    },
    add_track = function(track, step,
      position = "auto", color = "auto", alpha = .5) {
      "Adds a profile to the current ggkaryo plot. The input track must have
      already been binned with a consistent step. A consistent step is needed
      to automatically set any gap to 0 in the profile.
      @param track (string) path to BED5+ file.
      @param track (data.table) BED5+ data table.
      @param step (numerical) bin step in nt.
      @param position (string) one of auto|left|right. 'left' can be used only
        if opposite=T was used when initializing the ggkaryo object.
      @param color (string) either 'auto' or a color string.
      @param alpha (numerical) opacity level."
      stopifnot(position %in% c("auto", "right", "left"))
      if ( is(track, "character") ) {
        stopifnot(file.exists(track))
        track = fread(track)
      }
      stopifnot(is(track, "data.table"))
      stopifnot(ncol(track) >= 5)
      track = track[, 1:5]
      colnames(track) = c("chrom", "start", "end", "name", "value")
      track[, chromID := unlist(lapply(chrom, .self$chrom2id))]
      track = track[is.na(value), value := 0]

      track = track[order(start)]
      track = track[order(chromID)]

      track[, norm := value - min(value, na.rm = T)]
      track[, norm := value / max(value, na.rm = T)]
      track[is.na(norm), norm := 0]

      stopifnot(all(track[, .(v=unique(diff(start))), by=chrom]$v == step))

      set_gaps_to_zero = function(chrom_data, step) {
        "Adds 0s where gaps are detected in the track."
        id = which(diff(chrom_data$start) != step)
        if ( 0 == length(id) ) return(chrom_data)
        
        basepoints = chrom_data[c(id[1], id[1]+1),]
        basepoints[, norm := 0]
        if ( 1 == length(id) ) {
          return(do.call(rbind, list(
            chrom_data[1:id[1],],
            basepoints,
            chrom_data[(id[1]+1):nrow(chrom_data),])))
        } else {
          out = do.call(rbind, list(
            chrom_data[1:id[1],],
            basepoints,
            do.call(rbind, lapply(2:length(id), FUN = function(ii) {
              basepoints = chrom_data[c(id[ii], id[ii]+1),]
              basepoints[, norm := 0]
              rbind(chrom_data[(id[ii-1]+1):id[ii],], basepoints)
            })),
            chrom_data[(id[length(id)]+1):nrow(chrom_data),]
          ))
          return(out)
        }
      }
      track = track[, set_gaps_to_zero(.SD, step), by = chrom]

      add_chrom_ends = function(chrom_data) {
        "Sets the chromosome ends to 0."
        pre = chrom_data[1,]
        pre$value = NA; pre$norm = 0
        pos = chrom_data[nrow(chrom_data),]
        pos$value = NA; pos$norm = 0
        do.call(rbind, list(pre, chrom_data, pos))
      }
      track = track[, add_chrom_ends(.SD), by=chrom]

      nTracks = length(.self$data[['tracks']])
      .self$data[['tracks']][[nTracks+1]] = list(
        data = track, position = position, color = color)
    },

    add_lois = function(data, position, colorName, alpha = 1) {
      "docline"
      NULL
    },

    add_arm_boxes = function() {
      "Adds boxes around chromosome arms."
      .self$data[['plot']] = .self$data[['plot']] + geom_path(
        data=ggk$data[['boxes']], aes(group=paste0(chrom, "_", arm_id)),
        color="black")
    },
    add_chrom_labels = function() {
      "Adds chromosome labels."
      .self$data[['plot']] = .self$data[['plot']] + geom_text(
        data = .self$data[['chrom_labels']], aes(label = chrom), size = 5)
    },

    plot_base = function() {
      "Plots the current ggkaryo object (only basic layers)."
      print(.self$data[['plot']])
    },
    plot_full = function() {
      "Plots the current ggkaryo object with tracks and lois."

      p = .self$data[['plot']]

      nTracks = length(.self$data[['tracks']])
      if ( 0 != nTracks ) {
        for ( trackID in 1:nTracks ) {
          track = .self$data[['tracks']][[trackID]]$data
          track[, x := .self$chromID2x(chromID)]
          track[, x := x+.self$chrom_width+norm*.self$chrom_padding]
          track[, y := start+(end-start)/2]
          p = p + geom_polygon(data = as.data.frame(track),
            aes(group = chrom), fill = "red")
        }
      }

      print(p)
    }
  )
)

# Initialize ggkaryo object (only ideograms)
ggk = ggkaryo("giemsa.bed",
  chrom_width = 0.75, chrom_padding = 5)

# Add boxes around chromosome arms
ggk$add_arm_boxes()

# Add chromosome labels
ggk$add_chrom_labels()

# Add profile
trackData = as.data.table(readRDS("track_test.rds"
  ))[!is.na(value), .(chrom, start, end, name=paste0("bin_", 1:.N), value)]
binnedTrack = ggk$bin_track(trackData, 1e6, 1e5, na.rm=T)
ggk$add_track(binnedTrack, 1e5, "right")

# Add loci of interest
# ggk$add_lois(loiData, "right", "sample")

# Show plot
ggk$plot()

#ggk
