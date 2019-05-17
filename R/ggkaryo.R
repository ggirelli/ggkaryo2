
#' ggkaryo: a class for karyotype plotting and overlaying.
#'
#' The \\code{ggkaryo} class allows to plot (labeled) ideograms and to overlay
#' them with data track profiles and also highlight loci of interes (lois).
#'
#' @import methods
#' @export ggkaryo
#' @exportClass ggkaryo
#'
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 geom_polygon
#' @importFrom ggplot2 scale_fill_manual
#' @importFrom ggplot2 theme
#' @importFrom ggplot2 element_blank
#' @importFrom ggplot2 guides
#' @importFrom ggplot2 geom_segment
#' @importFrom ggplot2 aes_string
#' @importFrom ggplot2 scale_color_brewer
#' @importFrom ggplot2 guide_legend
#' @import cowplot
#' @import data.table
#' @import methods
#' @import RColorBrewer
#'
#' @field n_chrom (numerical) number of chromosomes, default: 24
#' @field hetero (character) heterosome labels (without "chr"),
#'        default: c("X", "Y")
#' @field chrom_width (numerical) width of the ideograms, default: 1
#' @field chrom_padding (numerical) space between ideograms, default: 5
#' @field track_palette_name (character) name of RColorBrewer palette for track
#'        filling
#' @field lois_palette_name (character) name of RColorBrewer palette for lois
#'        color
#' @field giemsa_palette (character) vector of colors for the giemsa bands
#' @field giemsa_levels (character) vector of giemsa band levels
#' @field opposite (logical) to plot profiles on both sides of the ideogram
#' @field data (list) contains ggkaryo data for plotting
#'
#' @examples
#' require(data.table)
#' require(ggkaryo2)
#'
#' # Load example data
#' data('giemsa', package='ggkaryo2')
#' data('track', package='ggkaryo2')
#' data('lois', package='ggkaryo2')
#'
#' # Plot ideogram
#' ggk = ggkaryo(giemsa)
#' ggk$plot_full()
#'
#' # Plot ideogram with boxes around chromosome arms and labels
#' ggk = ggkaryo(giemsa)
#' ggk$add_arm_boxes()
#' ggk$add_chrom_labels()
#' ggk$plot_full()
#'
#' # Plot ideogram with one profile track
#' ggk = ggkaryo(giemsa)
#' binnedTrack = track
#' ggk$add_track(binnedTrack, 1e5)
#' ggk$plot_full()
#'
#' # Plot ideogram with two profile tracks on the same side
#' ggk = ggkaryo(giemsa)
#' binnedTrack2 = copy(binnedTrack)
#' binnedTrack2[, value := value*abs(rnorm(nrow(binnedTrack2)))]
#' ggk$add_track(binnedTrack, 1e5)
#' ggk$add_track(binnedTrack2, 1e5)
#' ggk$plot_full()
#'
#' # Plot ideogram with two profile tracks on opposite sides
#' ggk = ggkaryo(giemsa, opposite=TRUE)
#' binnedTrack2 = copy(binnedTrack)
#' binnedTrack2[, value := value*abs(rnorm(nrow(binnedTrack2)))]
#' ggk$add_track(binnedTrack, 1e5)
#' ggk$add_track(binnedTrack2, 1e5)
#' ggk$plot_full()
#'
#' # Plot ideogram with two profile tracks on opposite sides and central lois
#' ggk = ggkaryo(giemsa)
#' binnedTrack2 = copy(binnedTrack)
#' binnedTrack2[, value := value*abs(rnorm(nrow(binnedTrack2)))]
#' ggk$add_track(binnedTrack, 1e5)
#' ggk$add_track(binnedTrack2, 1e5)
#' loiData = lois
#' ggk$add_lois(loiData, "center", "sample")
#' ggk$plot_full()
#'

ggkaryo <- setRefClass("ggkaryo",
  fields = list(
    n_chrom = "numeric",
    hetero = "character",
    chrom_width = "numeric",
    chrom_padding = "numeric",
    track_palette_name = "character",
    lois_palette_name = "character",
    giemsa_palette = "character",
    giemsa_levels = "character",
    opposite = "logical",
    data = "list"
  ),
  method = list(
    initialize = function(giemsa, ...,
        n_chrom=24, hetero=c("X", "Y"),
        chrom_width=1, chrom_padding=5,
        track_palette_name="Paired", lois_palette_name="Dark2",
        giemsa_palette=c(
          "#DDDDDD", "#9A9A9A", "#787878", "#555555", "#333333",
          "#FF0000", "#C4FFFC", "#AFE6FF"),
        giemsa_levels=c(
          "gneg", "gpos25", "gpos50", "gpos75", "gpos100",
          "acen", "gvar", "stalk"),
        opposite=FALSE
      ) {
      "Initializer method. See \\code{ggkaryo} description for more details"
      stopifnot(length(giemsa_levels) == length(giemsa_palette))
      stopifnot(chrom_width > 0)
      stopifnot(chrom_padding >= chrom_width)

      callSuper(...,
        n_chrom=n_chrom, hetero=hetero,
        chrom_width=chrom_width, chrom_padding=chrom_padding,
        track_palette_name=track_palette_name,
        lois_palette_name=lois_palette_name,
        giemsa_palette=giemsa_palette, giemsa_levels=giemsa_levels,
        opposite=opposite, data=list(tracks=list(), lois=list()))
      names(.self$giemsa_palette) = giemsa_levels

      .self$prep4karyo(giemsa)
    },

    chrom2id = function(chrom) {
      "Converts a chromosome signature (seqname) to a numerical id.
      \\describe{
        \\item{\\code{chrom}}{
          (string) chromosome signature (e.g., 'chr1' or '1')}
      }
      \\describe{\\item{returns}{numeric: chromosome numerical ID}}"
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
      \\describe{
        \\item{\\code{chromID}}{(numeric)}
      }
      \\describe{\\item{returns}{numeric: chromosome position on the X axis}}"
      return((chromID-1)*(.self$chrom_width + .self$chrom_padding))
    },
    norm2x = function(chromID, norm, position) {
      "Converts normalized score to X coordinate in the ggkaryo plot.
      \\describe{
        \\item{\\code{chromID}}{(numeric)}
        \\item{\\code{norm}}{(numeric) normalized score}
        \\item{\\code{position}}{(character) 'left' or 'right'}
      }
      \\describe{\\item{returns}{numeric: normalized score X coordinate}}"
      padding = .self$chrom_padding
      if ( .self$opposite )
        padding = padding / 2
      if ( "right" == position ) {
        return(.self$chromID2x(chromID)+.self$chrom_width+norm*padding)
      } else {
        stopifnot(.self$opposite)
        return(.self$chromID2x(chromID)-norm*padding)
      }
    },

    read_giemsa = function(giemsa) {
      "Reads a Giemsa bed file. Adds chromID, bandID, and X columns.
      \\describe{
        \\item{\\code{giemsa}}{(character) path to giemsa BED5+ file}
        \\item{\\code{giemsa}}{(data.table) data table with giemsa BED5+ data}
      }
      \\describe{\\item{returns}{data.table: adjusted giemsa data.table}}"
      if ( is(giemsa, "character") ) {
        stopifnot(file.exists(giemsa))
        giemsa = fread(giemsa)
      }
      stopifnot(is.data.table(giemsa))
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
      Builds .self\\$data[['bands']] object."

      stopifnot("giemsa" %in% names(.self$data))
      stopifnot(is.data.table(.self$data[['giemsa']]))

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
      chromosomes in two arms. Builds .self\\$data[['boxes']] object."

      stopifnot("giemsa" %in% names(.self$data))
      stopifnot(is.data.table(.self$data[['giemsa']]))

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
      Builds .self\\$data[['chrom_labels']] object."

      stopifnot("giemsa" %in% names(.self$data))
      stopifnot(is.data.table(.self$data[['giemsa']]))

      .self$data[['chrom_labels']] = .self$data[["giemsa"]][, .(
          x = min(x) + .self$chrom_width/2,
          y = -5 * 10**(ceiling(abs(log10(max(.self$data[["giemsa"]]$end))))-3)
        ), by = c("chrom", "chromID")]
      NULL
    },
    prep4plot = function() {
      "Prepares for plotting. Builds .self\\$data[['plot']] object."

      stopifnot("bands" %in% names(.self$data))
      stopifnot(is.data.table(.self$data[['bands']]))

      .self$data[['plot']] = ggplot(.self$data[['bands']], aes(x=x, y=-y)
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
      \\describe{
        \\item{\\code{giemsa}}{(character) path to giemsa BED5+ file}
        \\item{\\code{giemsa}}{(data.table) data table with giemsa BED5+ data}
      }
      \\describe{\\item{returns}{data.table: adjusted giemsa data.table}}"
      .self$data[["giemsa"]] = .self$read_giemsa(giemsa)
      .self$prep4bands()
      .self$prep4boxes()
      .self$prep4labels()
      .self$prep4plot()
      NULL
    },

    get_color = function(color, trackID) {
      "Extracts, in order, track colors from the .self\\$track_palette_name.
      See RColorBrewer for more details.
      \\describe{
        \\item{\\code{color}}{(character) a color or 'auto'}
        \\item{\\code{trackID}}{(numeric) track number}
      }"
      if ( "auto" == color ) {
        nTracks = max(3, length(.self$data[["tracks"]]))
        color = brewer.pal(nTracks, .self$track_palette_name)[trackID]
      }
      color
    },
    get_next_position = function(position) {
      "Selects position for the next track in such a fashion to balance out
      left/right sides of the ideograms.
      \\describe{
        \\item{\\code{position}}{(character) 'right', 'left', or 'auto'}
      }"
      stopifnot(position %in% c("auto", "right", "left"))
      if ( "auto" == position ) {
        if ( 0 == length(.self$data[['tracks']]) )
          return("right")
        if ( !.self$opposite ) {
          position = "right"
        } else {
          position_count = table(unlist(lapply(.self$data[['tracks']],
            function(x) x$position)))
          position_count = data.table(
            label = c("right", "left"),
            n = c(
              max(0, position_count["right"], na.rm = T),
              max(0, position_count["left"], na.rm = T)
            )
          )
          if ( position_count[1, n] == position_count[2, n] ) {
            return("right")
          } else {
            return(position_count[which.min(position_count$n), label])
          }
        }
      }
      position
    },

    bin_track = function(track, size, step, method="within",
      fun.aggreg=mean, ...) {
      "Bins a track based on provided bin size and step.
      Regions from the track are assigned to the bins when they are completely
      include ('within' method) or overlap even partially ('overlap' method).
      \\describe{
        \\item{\\code{track}}{(data.table) BED5+ track data table}
        \\item{\\code{size}}{(numeric) bin size in nt}
        \\item{\\code{step}}{(numeric) bin step in nt}
        \\item{\\code{method}}{(string) either 'within' or 'overlap'}
        \\item{\\code{fun.aggreg}}{(function) how to aggregate values in bins}
        \\item{\\code{...}}{(mixed) additional parameters to pass to fun.aggreg}
      }
      \\describe{\\item{returns}{data.table: binned track}}"
      stopifnot(is.data.table(track))
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
      Builds .self\\$data[['tracks']].
      \\describe{
        \\item{\\code{track}}{(character) path to BED5+ file}
        \\item{\\code{track}}{(data.table) BED5+ data table}
        \\item{\\code{step}}{(numerical) bin step in nt}
        \\item{\\code{position}}{(character) one of auto|left|right. 'left' can be used
          only if opposite=T was used when initializing the ggkaryo object}
        \\item{\\code{color}}{(character) either 'auto' or a color string}
        \\item{\\code{alpha}}{(numerical) opacity level.}
      }"
      position = .self$get_next_position(position)
      stopifnot(alpha <= 1 & alpha > 0)
      if ( is(track, "character") ) {
        stopifnot(file.exists(track))
        track = data.table::fread(track)
      }
      stopifnot(is.data.table(track))
      stopifnot(ncol(track) >= 5)
      track = track[, 1:5]
      colnames(track) = c("chrom", "start", "end", "name", "value")
      track[, chromID := unlist(lapply(chrom, .self$chrom2id))]
      track[is.na(value), value := 0]

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
        data = track, position = position, color = color, alpha = alpha)
    },

    add_lois = function(loiData, position, colorName, alpha = 1) {
      "Adds details on Loci Of Interest (loi) to the current ggkaryo plot.
      Builds .self\\$data[['lois']].
      \\describe{
        \\item{\\code{loiData}}{(character) path to BED5+ loi file}
        \\item{\\code{loiData}}{(data.table) data.table with BED5+ loi data}
        \\item{\\code{position}}{(character) either 'left', 'right' or 'center'}
        \\item{\\code{colorName}}{(character) column with color factors}
        \\item{\\code{alpha}}{(numeric) opacity level}
      }"
      stopifnot(position %in% c("left", "right", "center"))
      if ( !.self$opposite ) stopifnot(position %in% c("right", "center"))
      stopifnot(alpha <= 1 & alpha > 0)
      if ( is(loiData, "character") ) {
        stopifnot(file.exists(loiData))
        loiData = data.table::fread(loiData)
      }
      stopifnot(is.data.table(loiData))
      stopifnot(ncol(loiData) >= 5)
      stopifnot(colorName %in% names(loiData))
      loiData = loiData[, .SD, .SDcols=c(1:5, which(names(loiData)==colorName))]
      colnames(loiData) = c("chrom", "start", "end", "name", "value", colorName)
      loiData[, chromID := unlist(lapply(chrom, .self$chrom2id))]
      loiData[is.na(value), value := 0]
      loiData[, colorName] = factor(unlist(loiData[, .SD, .SDcols=colorName]))

      .self$data[['lois']] = list(data=loiData, position=position,
        color=colorName, alpha=alpha)
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
      require(cowplot)
      print(.self$data[['plot']])
    },
    add_lois_overlay = function(p) {
      "Overlays track profiles to a ggkaryo plot.
      \\describe{
        \\item{\\code{p}}{(ggplot)}
      }"
      if ( 0 != length(.self$data[['lois']]) ) {
        lois = .self$data[['lois']]

        padding = .self$chrom_padding
        if ( .self$opposite ) padding = padding/2

        if ( lois$position == "right" ) {
          loiData[, x := .self$chromID2x(lois$data$chromID)+.self$chrom_width]
          loiData[, xend := .self$chromID2x(lois$data$chromID
            )+.self$chrom_width+padding]
        }
        if ( lois$position == "left" ) {
          loiData[, x := .self$chromID2x(lois$data$chromID)-padding]
          loiData[, xend := .self$chromID2x(lois$data$chromID)]
        }
        if ( lois$position == "center" ) {
          loiData[, x := .self$chromID2x(lois$data$chromID)+.self$chrom_width]
          loiData[, xend := .self$chromID2x(lois$data$chromID)]
        }

        loiData[, y := (start+end)/2]
        `+.uneval` <- function(a,b) {
          `class<-`(modifyList(a,b), "uneval")
        }
        p = p + geom_segment(data=loiData, aes(xend=xend, yend=-y
            ) + aes_string(color=lois$color), alpha=lois$alpha
          ) + scale_color_brewer(palette=.self$lois_palette_name
          ) + guides(color = guide_legend(title=lois$color))
      }
      return(p)
    },
    add_track_overlay = function(p) {
      "Overlays track profiles to a ggkaryo plot.
      \\describe{
        \\item{\\code{p}}{(ggplot)}
      }"
      nTracks = length(.self$data[['tracks']])
      if ( 0 != nTracks ) {
        for ( trackID in 1:nTracks ) {
          track = .self$data[['tracks']][[trackID]]
          track$data[, x := .self$norm2x(chromID,norm,track$position), by=chrom]
          track$data[, y := start+(end-start)/2]
          p = p + geom_polygon(data = as.data.frame(track$data),
            aes(group = chrom), fill = get_color(track$color, trackID),
            alpha = track$alpha)
        }
      }
      return(p)
    },
    plot_full = function() {
      "Plots the current ggkaryo object with tracks and lois."
      p = .self$data[['plot']]
      p = .self$add_lois_overlay(p)
      p = .self$add_track_overlay(p)
      require(cowplot)
      print(p)
    }
  )
)
