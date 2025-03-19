library(QFeatures)

## example(SummarizedExperiment)

data(feat1)

feat1

colData(feat1)

colData(feat1)$X <- c("x1", "x2")

feat1$Y <- c("y1", "y2")

feat1

feat1[[1]]
feat1[["psms"]]

assay(feat1[["psms"]])

colData(feat1[["psms"]])

rowData(feat1[["psms"]])

## Aggregation psms -> peptides

feat1 <- aggregateFeatures(
    feat1,
    i = 1, ## or "psms"
    fcol = "Sequence",
    name = "peptides",
    fun = colMeans,
    na.rm = TRUE
)

feat1

## Ex: verify the calculations

## Ex: aggregate peptides -> proteins

rowData(feat1[[2]])

feat1 <- aggregateFeatures(
    feat1,
    i = 2,
    fcol = "Protein",
    name = "proteins",
    fun = MsCoreUtils::robustSummary
)

## Subsetting and filtering

feat1["ProtA", , ]

filterFeatures(feat1, ~ pval < 0.05)

filterFeatures(feat1, ~ pval < 0.05, keep = TRUE)

filterFeatures(feat1, ~ pval < 0.05, i = 1)

filterFeatures(feat1, ~ location == "Mitochondrion")

## Ex: filter rows that do not localise to the mitochondrion

filterFeatures(feat1, ~ location != "Mitochondrion")

## Creating QFeatures

data(hlpsms)

class(hlpsms)

dim(hlpsms)

hlpsms

hl <- readQFeatures(hlpsms, quantCols = 1:10, name = "psms")

hl

hl[[1]]

assay(hl[[1]])

rowData(hl[[1]])

se <- readSummarizedExperiment(hlpsms, quantCols = 1:10)

se

QFeatures(list(psms = se))


## Ex: create SummarizedExperiment for the CPTAC data

## (f <- msdata::quant(full.names = TRUE))

(f <- MsDataHub::cptac_a_b_peptides.txt())


tmp <- read.delim(f)
class(tmp)
dim(tmp)

names(tmp)

(i <- grep("Intensity\\.", names(tmp)))

cptac_se <- readSummarizedExperiment(tmp,
                                     quantCols = i,
                                     fnames = "Sequence")

cptac_se

colData(cptac_se)

colnames(cptac_se) <- sub("Intensity\\.", "", colnames(cptac_se))

cptac_se$condition <- rep(c("6A", "6B"), each = 3)
cptac_se$id <- rep(7:9, 2)

names(rowData(cptac_se))

keep_var <- c("Sequence", "Proteins",
              "Leading.razor.protein",
              "Score", "Reverse", "PEP",
              "Potential.contaminant")

rowData(cptac_se) <- rowData(cptac_se)[, keep_var]

rowData(cptac_se)

## Processing and analyses

assay(cptac_se)

anyNA(assay(cptac_se))

cptac_se <- zeroIsNA(cptac_se)

anyNA(assay(cptac_se))

assay(cptac_se)

nNA(cptac_se)

barplot(nNA(cptac_se)$nNAcols$nNA)

table(nNA(cptac_se)$nNArow$nNA)

cptac_se <- filterNA(cptac_se, pNA = 4/6)

table(nNA(cptac_se)$nNArow$nNA)

impute(cptac_se, method = "knn")

impute(cptac_se, method = "MinDet")

impute(cptac_se, method = "zero")

MsCoreUtils::imputeMethods()

## QC

table(rowData(cptac_se)$Reverse)

table(rowData(cptac_se)$Potential.contaminant)

## Ex: Visualise the identification score distributions from forward and reverse
## hits and interpret the figure.

library(tidyverse)

rowData(cptac_se) |>
    as_tibble() |>
    ggplot(aes(x = Score,
               colour = Reverse)) +
    geom_density()

rowData(cptac_se) |>
    as_tibble() |>
    ggplot(aes(x = PEP,
               colour = Reverse)) +
    geom_density()

max(rowData(cptac_se)$PEP)

prots <- rowData(cptac_se)$Proteins
names(prots) <- rowData(cptac_se)$Sequence

head(prots)


## Connected components from quant data
library(PSMatch)
adj <- makeAdjacencyMatrix(prots, split = ";")
dim(adj)
cc <- ConnectedComponents(adj)
cc

## Create a QFeatures object

(cptac <- QFeatures(list(peptides = cptac_se)))

colData(cptac) <- colData(cptac_se)


## Ex: Using the filterFeatures() function, filter out the reverse and
## contaminant hits, and also retain those that have a posterior error
## probability smaller than 0.05.

cptac <-
    cptac |>
    filterFeatures(~ Reverse != "+") |>
    filterFeatures(~ Potential.contaminant != "+") |>
    filterFeatures(~ PEP < 0.05)

cptac

cptac <- logTransform(cptac,
                      i = "peptides",
                      name = "log_peptides")

cptac

## Ex: Use the normalize() method to normalise the data. The syntax is the same
## as logTransform(). Use the "center.median" method.

cptac <- cptac |>
    normalize(i = "log_peptides",
              name = "lognorm_peptides",
              method = "center.median")



library(limma)

plotDensities(assay(cptac[[1]]))
plotDensities(assay(cptac[[2]]))
plotDensities(assay(cptac[[3]]))

boxplot(assay(cptac[[2]]))
boxplot(assay(cptac[[3]]))

## Use median to aggregation peptides into protein values

cptac <- aggregateFeatures(
    cptac,
    i = "lognorm_peptides",
    name = "proteins_med",
    fcol = "Leading.razor.protein",
    fun = colMedians,
    na.rm = TRUE)


prcomp(t(assay(cptac[[4]])))

library(factoextra)

cptac[[4]] |>
    filterNA() |>
    assay() |>
    t() |>
    prcomp(scale = TRUE, center = TRUE) |>
    fviz_pca_ind(habillage = cptac$condition)

cptac[[4]] |>
    impute(method = "knn") |>
    assay() |>
    t() |>
    prcomp(scale = TRUE, center = TRUE) |>
    fviz_pca_ind(habillage = cptac$condition)

cptac[[4]] |>
    impute(method = "zero") |>
    assay() |>
    t() |>
    prcomp(scale = TRUE, center = TRUE) |>
    fviz_pca_ind(habillage = cptac$condition)

nNA(cptac[[4]])

longFormat(cptac["P02787ups|TRFE_HUMAN_UPS", ,
                 c("lognorm_peptides", "proteins_med")]) |>
    as_tibble() |>
    mutate(condition = ifelse(grepl("A", colname), "A", "B")) |>
    ggplot(aes(
        x = colname,
        y = value,
        colour = rowname,
        shape = condition)) +
    geom_point(size = 3) +
    geom_line(aes(group = rowname)) +
    facet_grid(~ assay)


cptac["P02787ups|TRFE_HUMAN_UPS", ,
      c("lognorm_peptides", "proteins_med")] |>
    longFormat()


tmp <- impute(cptac, i = 3,
              method = "knn") |>
    aggregateFeatures(
        i = "imputedAssay",
        name = "proteins_med2",
        fcol = "Leading.razor.protein",
        fun = colMedians,
        na.rm = TRUE)

longFormat(tmp["P02787ups|TRFE_HUMAN_UPS", ,
               c("imputedAssay", "proteins_med2")]) |>
    as_tibble() |>
    mutate(condition = ifelse(grepl("A", colname), "A", "B")) |>
    ggplot(aes(
        x = colname,
        y = value,
        colour = rowname,
        shape = condition)) +
    geom_point(size = 3) +
    geom_line(aes(group = rowname)) +
    facet_grid(~ assay)


plot(cptac)


normalize(cptac, "log_peptides",
          name = "logquantiles_peptides",
          method = "quantiles.robust") |>
    aggregateFeatures(
        "logquantiles_peptides",
        name = "proteins_med2",
        fcol = "Leading.razor.protein",
        fun = colMedians,
        na.rm = TRUE) |>
    plot()

## Statistical analysis

prots <- getWithColData(cptac, "proteins_med")
prots

colData(prots)