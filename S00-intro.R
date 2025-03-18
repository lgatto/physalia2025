BiocManager::version()

library(Spectra)

## bioconductor.org

install.package("BiocManager")

BiocManager::install("Spectra")

## PRIDE, MassIVE, JProt

library("rpx")

px <- PXDataset("PXD000001")

px
pxfiles(px)
pxtax(px)
pxurl(px)
pxref(px)

f <- pxget(px, "TMT_Erwinia_1uLSike_Top10HCD_isol2_45stepped_60min_01-20141210.mzML")

tail(pxfiles(px))
pxget(px, "erwinia_carotovora.fasta")


px2 <- PXDataset("PXD022816")
pxref(px2)
pxfiles(px2)


library("MsDataHub")
MsDataHub()

ko15.CDF()

library("msdata")

proteomics()
ident()
quant()

f <- proteomics(full.names = TRUE)[4]

## library("pRolocdata")
## library("scpdata")

## Convertion to mzML
## - msconvert (from proteowizard)
## - https://github.com/compomics/ThermoRawFileParser
