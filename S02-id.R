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
