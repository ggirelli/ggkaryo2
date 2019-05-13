library(cowplot)
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

      track = track[order(start)]
      track = track[order(chromID)]

      track[, norm := value - min(value, na.rm = T)]
      track[, norm := value / max(value, na.rm = T)]
      track[is.na(norm), norm := 0]

      set_gaps_to_zero = function(chrom_data, step) {
        id = which(diff(chrom_data$start) != step)
        stopifnot( 0 == length(id) || id == 1 )
        
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

      track = do.call(rbind, by(track, track$chrom, FUN = function(ct) {
        pre = ct[1,]
        pre$value = NA
        pre$norm = 0
        pos = ct[nrow(ct),]
        pos$value = NA
        pos$norm = 0
        do.call(rbind, list(pre, ct, pos))
      }))

      .self$data[['tracks']] = c(.self$data[['tracks']],
        list(data = track, position = position, color = color))
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

    plot = function() {
      "Plots the current ggkaryo object."

      # nTracks = length(.self$data[['tracks']])
      # lapply(1:nTracks, function(trackID) {
      #   track = .self$data[['tracks']][[trackID]]

        
      # })

      print(.self$data[['plot']])
    }
  )
)

# # Initialize ggkaryo object (only ideograms)
# ggk = ggkaryo("/mnt/data/Resources/hg19.giemsa_bands.bed",
#   chrom_width = 0.75, chrom_padding = 5)

# # Add boxes around chromosome arms
# ggk$add_arm_boxes()

# # Add chromosome labels
# ggk$add_chrom_labels()

# # Add profile
# trackData = as.data.table(readRDS("/home/gire/Desktop/Code/src/190507_jesko/profile.rds"
#   ))[, .(chrom, start, end, name=paste0("bin_", 1:.N), value)]
# ggk$add_track(trackData, "right")

# # Add loci of interest
# # ggk$add_lois(loiData, "right", "sample")

# # Show plot
# ggk$plot()

# #ggk
