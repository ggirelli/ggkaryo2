#' ggkaryo2: A package for ggplot-compatible overlaying of karyotype and tracks
#'
#' ggkaryo2 is fully based on ggplot. It allows to plot a simple karyotype, or a
#' karyotype with one/multiple data-track profile(s). This can be useful when
#' visualizing multiple informations, as it allows to easily identify the
#' localization (e.g., genomic coordinate, location relative to
#' centromere/telomeres, ...) of any region of interest based on a data-track
#' profile. The karyotype is automatically built from a giemsa-staining bed
#' file, and the data-track are added as additional bed files. ggkaryo2 supports
#' simultaneously plotting, onto a karyotype, any number of data-tracks as
#' profiles, and to highlight loci of interest. Being fully based on ggplot, the
#' users can easily customize a ggkaryo2 plot and add their own custom layers to
#' it.
#' 
#' @section ggkaryo2 classes:
#' ggkaryo
#'
#' @docType package
#' @name ggkaryo2
NULL