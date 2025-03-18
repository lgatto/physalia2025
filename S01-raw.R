## BiocManager::install("Spectra")

library("Spectra")



data.frame()
## tibble()

spd <- DataFrame(msLevel = c(1L, 2L),
                 rtime = c(1.1, 1.2))
spd

spd$mz <- list(c(100, 103.2, 104.3, 106.5),
               c(45.6, 120.4, 190.2))

spd$intensity <- list(c(200, 100, 500.1, 200),
                      c(12.3, 130, 34.1))

spd

sp0 <- Spectra(spd)

sp0

spectraVariables(sp0)

spectraData(sp0)

peaksData(sp0)

peaksData(sp0)[[1]]

mz(sp0)
intensity(sp0)


sp <- Spectra(f)
sp

spectraData(sp)$msLevel
sp$msLevel

msLevel(sp)

rtime(sp)

peaksData(sp)

head(peaksData(sp)[[1]])

precursorIntensity(sp)

table(msLevel(sp))

msLevel(sp)[1000:1010]

precursorIntensity(sp)[1000:1010]

sp[1000:1010]

plotSpectra(sp[1000:1003])


## Ex

sp <- Spectra(f)

## https://rformassspectrometry.github.io/Spectra

## The chromatogram at the top displays the total ion current along the
## retention time. The vertical line identifies one scan in particular at
## retention time 1800.68 seconds (the 2807th scan).

## 1. The chromatogram can be created by extracting the `totIonCurrent` and
##    `rtime` variables for all MS1 spectra. Annotate the spectrum of interest.


sel1 <- msLevel(sp) == 1

plot(rtime(sp)[sel1], tic(sp)[sel1], type = "l")

filterMsLevel(sp, 1)

with(spectraData(filterMsLevel(sp, 1)),
     plot(rtime, totIonCurrent, type = "l"))

plot(filterMsLevel(sp, 1)$rtime,
     filterMsLevel(sp, 1)$totIonCurrent,
     type = "l")

abline(v = rtime(sp)[2807], col = "red")

## Plot MS1 spectrum 2807 (tip: use plotSpectra())


plotSpectra(sp[2807])

plotSpectra(sp[2807], xlim = c(400, 1000))


## Plot the MS1 and/or the 10 MS2 scans

sp2 <- filterPrecursorScan(sp, 2807)

plotSpectra(sp2[1], xlim = c(400, 1000))

plotSpectra(sp2)

plotSpectra(sp2[-1])

## Ex: plot the MS1 scan and highlight the precursors that were selected
## hint: see precursorMz()


plotSpectra(sp2[1], xlim = c(400, 1000))

abline(v = precursorMz(sp2)[-1], col = "grey")
abline(v = precursorMz(sp2)[2], col = "red")

plotSpectra(sp2[1], xlim = c(521.1, 522.5), type = "l")
abline(v = precursorMz(sp2)[2], col = "red")

length(mz(sp2[7])[[1]])

mzLabel <- function(z) {
    zpd <- peaksData(z)[[1]]
    lbs <- format(zpd[, "mz"], digits = 4)
    lbs[zpd[, "intensity"] < 1e5] <- ""
    lbs
}


plotSpectra(sp2[7],
            xlim = c(126, 132),
            labels = mzLabel)


mzLabel2 <- function(z) {
    zpd <- peaksData(z)[[1]]
    lbs <- format(zpd[, "intensity"], digits = 4)
    lbs[zpd[, "intensity"] < 1e5] <- ""
    lbs
}

plotSpectra(sp2[7],
            xlim = c(126, 132),
            labels = mzLabel2)


## Ex: find MS2 spectra in sp that have the same precursor M/Z.
## Hint: precursorMz()

sp3 <- filterMsLevel(sp, 2L)

anyDuplicated(precursorMz(sp3))

i <- which(precursorMz(sp3)[37] == precursorMz(sp3))

plotSpectraOverlay(sp3[i],
                   col = c("red", "steelblue"),
                   xlim = c(300, 400),
                   ylim = c(0, 5000))

plotSpectraMirror(sp3[i][1], sp3[i][2],
                  xlim = c(100, 800))



k <- which(duplicated(precursorMz(sp3)))


plotSpectraOverlay(sp3[c(86, 101, 102)],
                   col = c("red", "steelblue", "gree"),
                   ylim = c(0, 1e4))

plotSpectra(sp3[c(86, 101, 102)])

## BiocManager::install("RforMassSpectrometry/SpectraVis")
library(SpectraVis)

browseSpectra(sp)

plotlySpectra(sp2[2])

## more ...

(fls <- dir(system.file("sciex", package = "msdata"), full.names = TRUE))

sp_sciex <- Spectra(fls)

table(dataOrigin(sp_sciex))

sp_sciex

setBackend(sp_sciex, MsBackendMemory())

sp_sciex

filterDataOrigin(sp_sciex, fls[2])

filterRt(sp_sciex, rt = c(175, 189))

sp_sciex |>
    filterDataOrigin(fls[2]) |>
    filterRt(rt = c(175, 189)) |>
    spectraData()

library(BiocParallel)


min(intensity(sp_sciex[1]))

sp_sciex <- filterIntensity(sp_sciex, intensity = c(10, Inf))

sp_sciex

min(intensity(sp_sciex[1]))

sp_sciex@processingQueue

sp_sciex <- reset(sp_sciex)

min(intensity(sp_sciex[1]))

## EX:
##
## 1. Download the 3 first mzML files from the PXD022816 project (using rpx)
## 2. Create a Spectra object containing all these data
## 3. How many scans? How many scans per file? How many MS1 and MS2 scan? How
##    many MS1 and MS2 scan per file?
## 4. Visualise the chromatograms of these 3 files (ideally on one figure).

library(rpx)

px2 <- PXDataset("PXD022816")
pxref(px2)

grep("mzML", pxfiles(px2))[1:3]

f <- pxget(px2, grep("mzML", pxfiles(px2))[1:3])

## f <- dir("data", full.names = TRUE)

library(Spectra)

sp <- Spectra(f)

length(sp)

table(dataOrigin(sp))

table(msLevel(sp))

table(dataOrigin(sp), msLevel(sp))

library(tidyverse)

filterMsLevel(sp, 1L)|>
    spectraData() |>
    as_tibble() |>
    ggplot(aes(x = rtime,
               y = totIonCurrent,
               colour = basename(dataOrigin))) +
    geom_line()


filterMsLevel(sp, 1L) |>
    spectraData() |>
    as_tibble() |>
    ggplot(aes(x = rtime,
               y = totIonCurrent)) +
    geom_line() +
    facet_wrap(~ basename(dataOrigin))
