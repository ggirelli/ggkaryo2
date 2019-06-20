ggkaryo2
===

[![DOI](https://zenodo.org/badge/185577333.svg)](https://zenodo.org/badge/latestdoi/185577333)

An R package to overlay karyotype and data-track profiles in a `ggplot`-compatible manner.

Installation
-------------

To **install**, run the following:

```R
require("devtools")
devtools::install_github("ggirelli/ggkaryo2")
```

To **uninstall** run the following from within the repository folder:

```R
remove.packages("ggkaryo2")
```

Usage
----------

```R
require(ggkaryo2)
require(data.table)

# Load example data
data('giemsa', package='ggkaryo2')
data('track', package='ggkaryo2')
data('lois', package='ggkaryo2')

# Plot ideogram
ggk = ggkaryo(giemsa)
ggk$plot_full()

# Plot ideogram with boxes around chromosome arms and labels
ggk = ggkaryo(giemsa)
ggk$add_arm_boxes()
ggk$add_chrom_labels()
ggk$plot_full()

# Plot ideogram with one profile track
ggk = ggkaryo(giemsa)
binnedTrack = track
ggk$add_track(binnedTrack, 1e5)
ggk$plot_full()

# Plot ideogram with two profile tracks on the same side
ggk = ggkaryo(giemsa)
binnedTrack2 = copy(binnedTrack)
binnedTrack2[, value := value*abs(rnorm(nrow(binnedTrack2)))]
ggk$add_track(binnedTrack, 1e5)
ggk$add_track(binnedTrack2, 1e5)
ggk$plot_full()

# Plot ideogram with two profile tracks on opposite sides
ggk = ggkaryo(giemsa, opposite=T)
binnedTrack2 = copy(binnedTrack)
binnedTrack2[, value := value*abs(rnorm(nrow(binnedTrack2)))]
ggk$add_track(binnedTrack, 1e5)
ggk$add_track(binnedTrack2, 1e5)
ggk$plot_full()

# Plot ideogram with two profile tracks on opposite sides and central lois
ggk = ggkaryo(giemsa)
binnedTrack2 = copy(binnedTrack)
binnedTrack2[, value := value*abs(rnorm(nrow(binnedTrack2)))]
ggk$add_track(binnedTrack, 1e5)
ggk$add_track(binnedTrack2, 1e5)
loiData = lois
ggk$add_lois(loiData, "center", "sample")
ggk$plot_full()
```

Contributing
---

We welcome any contributions to `ggkaryo2`. Please, refer to the [contribution guidelines](https://github.com/ggirelli/ggkaryo2/blob/master/CONTRIBUTING.md) if this is your first time contributing! Also, check out our [code of conduct](https://github.com/ggirelli/ggkaryo2/blob/master/CODE_OF_CONDUCT.md).

License
---

```
MIT License
Copyright (c) 2019 Gabriele Girelli
```
