load("../Mar.datawrangling/data/RV.GSSPECIES.RData")
rm(definitive_Aphia)
rm(definitive_TSN)
rm(multiMatches)
rm(mysterySpec)
GSSPECIES = GSSPECIES[,c("SPEC", "COMM","CODE")]
test = getTaxaIDs(head(GSSPECIES,25), comm_col = "COMM", sci_col = "SPEC")