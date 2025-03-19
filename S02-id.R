library(PSMatch)

idf <- MsDataHub::TMT_Erwinia_1uLSike_Top10HCD_isol2_45stepped_60min_01.20141210.mzid()

## idf <- msdata::ident(full.names = TRUE)

idf

id <- PSM(idf)
id

dim(id)
names(id)

data.frame(id)

## Verify that this table contains 5802 matches for 5343 scans and
## 4938 peptides sequences.

id


nrow(id)

length(unique(id$spectrumID))

length(unique(id$sequence))

unique(id$spectrumFile)

BiocManager::install("cleaver")


calculateFragments("RQCRTDFLNYLR")

summary(id$MS.GF.RawScore)

## Proteins -> peptides -> PSM/scores

## Decoy (reverse proteins) -> peptides -> decoy PSM/scores

table(id$isDecoy)

table(table(id$spectrumID))

i <- which(id$spectrumID == "controllerType=0 controllerNumber=1 scan=1774")

data.frame(id[i, c(1, 30:32)])


data.frame(id[i, ])


id2 <- reducePSMs(id, id$spectrumID)

j <- which(id2$spectrumID == "controllerType=0 controllerNumber=1 scan=1774")

data.frame(id2[j, ])

id2[j, "DatabaseAccess"]

id[!id$isDecoy, ]

id[id$rank == 1, ]

id_filtered <- PSMatch::filterPSMs(id)

describePeptides(id)

describeProteins(id)
describeProteins(id_filtered)

## Compare the distribution of raw identification MS.GF.RawScore of
## the decoy and non-decoy hits. Interpret the figure.

id |>
    ggplot(aes(x = MS.GF.RawScore,
               colour = isDecoy)) +
    geom_density()

## Understanding protein groups with adjacency matrices

id0 <- msdata::ident(full.names = TRUE,
                     pattern = "TMT") |>
    PSM() |>
    filterPsmDecoy() |>
    filterPsmRank()



id0

data.frame(id0[1:10, c("sequence", "DatabaseAccess")])

adj <- makeAdjacencyMatrix(id)
adj[1:10, 1:5]

dim(adj)

cc <- ConnectedComponents(adj)

length(cc)

cc

connectedComponents(cc, 1)




tibble(id = 1:length(cc),
       nr = nrows(cc),
       nc = ncols(cc)) |>
    filter(nc > 1) |>
    filter(nr > 1)

connectedComponents(cc, 2)

connectedComponents(cc, 549)

tibble(id = 1:length(cc),
       nr = nrows(cc),
       nc = ncols(cc)) |>
    filter(nc == 2) |>
    filter(nr == 4)


connectedComponents(cc, 2504)


tibble(id = 1:length(cc),
       nr = nrows(cc),
       nc = ncols(cc)) |>
    filter(nr > 4) |>
    arrange(desc(nc))


cx <- connectedComponents(cc, 1124)

plotAdjacencyMatrix(cx)

plotAdjacencyMatrix(cx, 1)

plotAdjacencyMatrix(cx, 2)

## https://rformassspectrometry.github.io/PSMatch/articles/AdjacencyMatrix.html

cctab <- prioritiseConnectedComponents(cc)

head(cctab)

library(factoextra)

fviz_pca(prcomp(cctab, scale = TRUE, center = TRUE))

## Add id to spectra

## Identify the spectum identifier columns in the sp the id_filtered
## variables needed for integration.


spectraVariables(sp)

head(sp$spectrumId)

head(id_filtered$spectrumID)

table(table(id_filtered$spectrumID))

which(table(id_filtered$spectrumID) == 4)

which(id_filtered$spectrumID == "controllerType=0 controllerNumber=1 scan=5490")


data.frame(id_filtered[14:17, 1:15])

data.frame(id_filtered[14:17, c("modName", "modLocation")])


## px <- PXDataset("PXD000001")
## f <- pxget(px, "TMT_Erwinia_1uLSike_Top10HCD_isol2_45stepped_60min_01-20141210.mzML")
## sp <- Spectra(f)

## idf <- MsDataHub::TMT_Erwinia_1uLSike_Top10HCD_isol2_45stepped_60min_01.20141210.mzid()
## id_filtered <- filterPSMs(id)
## id_filtered <- reducePSMs(id_filtered, id_filtered$spectrumID)

id_filtered <- reducePSMs(id_filtered, id_filtered$spectrumID)

table(table(id_filtered$spectrumID))


sp2 <- joinSpectraData(sp, id_filtered,
                       by.x = "spectrumId",
                       by.y = "spectrumID")

## Verify that the identification data has been added to the correct
## spectra.
## MS1 -> no identication added
## MS2 -> some identication added


msLevel(sp2)
sp2$sequence

all(is.na(filterMsLevel(sp2, 1L)$sequence))

table(is.na(filterMsLevel(sp2, 2L)$sequence))




sp2 <- countIdentifications(sp2)

table(msLevel(sp2), sp2$countIdentifications)

sp2 |>
    filterMsLevel(1L) |>
    spectraData() |>
    as_tibble() |>
    ggplot(aes(x = rtime,
               y = totIonCurrent)) +
    geom_line() +
    geom_point(
        aes(colour = ifelse(countIdentifications == 0,
                            NA, countIdentifications)))


i <- which(sp2$MS.GF.RawScore > 100)[1]



plotSpectra(sp2[i])

plotSpectra(sp2[i], labels = addFragments)

calculateFragments(sp2[i]$sequence)

plotSpectra(sp2[i], labels = addFragments,
            labelCol = "red")

## Compare spectra

## compareSpectra()

sp2

## Create a new Spectra object containing the MS2 spectra with
## sequences "SQILQQAGTSVLSQANQVPQTVLSLLR" and
## "TKGLNVMQNLLTAHPDVQAVFAQNDEMALGALR".

sp_k <- sp2[which(sp2$sequence %in% c("SQILQQAGTSVLSQANQVPQTVLSLLR", "TKGLNVMQNLLTAHPDVQAVFAQNDEMALGALR"))]

sp_k


## Calculate the 5 by 5 similarity matrix between all spectra using
## compareSpectra(). See the ?Spectra man page for details. Draw a
## heatmap of that matrix.

mat <- compareSpectra(sp_k)

mat

colnames(mat) <- rownames(mat) <- strtrim(sp_k$sequence, 3)

mat

pheatmap::pheatmap(mat)


## Visually compare the spectra with the plotting function seen
## previously.

## - filterIntensity()  keep peaks with an intensity > 1e3

filterIntensity(sp_k, 1e3) |>
    plotSpectra()

sp_k2 <- filterIntensity(sp_k, 1e3)

plotSpectraMirror(sp_k2[1], sp_k2[2])


plotSpectraOverlay(sp_k2[3:5],
                   col = c("red", "blue", "green"))

par(mfrow = c(2, 1))

plotSpectraMirror(sp_k2[3], sp_k2[4])

plotSpectraMirror(sp_k2[3], sp_k2[5])


filterIntensity(sp_k, 1e3) |>
    plotSpectra(labels = addFragments,
                labelCol = "red")

## Summary exercise

library(rpx)
px <- PXDataset("PXD022816")

(mzmls <- pxget(px, grep("mzML", pxfiles(px))[1:3]))

library(Spectra)
sp <- Spectra(mzmls)

## See S01-raw.R

(mzids <- pxget(px, grep("mzID", pxfiles(px))[1:3]))

## Check the quality of the identification data by comparing the density of the
## decoy and target PSMs id scores for each file

library(PSMatch)

id <- PSM(mzids)

id

table(basename(dataOrigin(filterMsLevel(sp, 2))))

names(id)

table(id$idFile)

table(basename(dataOrigin(filterMsLevel(sp, 2)))) / table(id$idFile)


ggplot(id,
       aes(x = MetaMorpheus.score,
           colour = isDecoy)) +
    geom_density() +
    facet_wrap(~ spectrumFile)

max(id$PSM.level.q.value)

table(id$isDecoy)


id_filtered <- filterPSMs(id)


i <- grep("scan=12040", sp$spectrumId)

sp$spectrumId[i]
dataOrigin(sp)[i]


paste(sub("^.+_QEP", "QEP", basename(dataOrigin(sp)[i])),
      sub("^.+scan=", "", sp$spectrumId[i]),
      sep = "::")

head(id$spectrumID)
head(id$spectrumFile)