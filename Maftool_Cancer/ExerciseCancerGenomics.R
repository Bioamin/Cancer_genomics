#Calling the library
library(maftools)

# lung adenocarcinoma(LUAD)
#loading the mutation file
luad_maf = "/Users/amin/Desktop/Cancer_Maftool/luad_tcga/data_mutations_extended.txt" 
luad = read.maf(maf = luad_maf)

#lung squamous cell carcinoma (LUSC)
#loading the mutation file
lusc_maf = "/Users/amin/Desktop/Cancer_Maftool/lusc_tcga/data_mutations.txt" 
lusc = read.maf(maf = lusc_maf)

#gene summary
print("____________LUAD____Summary______________")
getGeneSummary(luad)
print("____________LUSC____Summary______________")
getGeneSummary(luad)

#A general visualization
#LUAD
plotmafSummary(
  maf = luad, rmOutlier = TRUE, 
  addStat = 'mean', dashboard = TRUE,
  titvRaw = FALSE)
#LUSC
plotmafSummary(
  maf = lusc, rmOutlier = TRUE, 
  addStat = 'mean', dashboard = TRUE,
  titvRaw = FALSE)

#oncoplot for top ten most frequent mutated genes.
#LUAD
oncoplot(maf = luad, top = 10, keepGeneOrder = FALSE, sortByMutation=FALSE)
#LUSC
oncoplot(maf = lusc, top = 10, keepGeneOrder = FALSE, sortByMutation=FALSE)


#LUAD exclusive/co-occurance event analysis on top 10 mutated genes. 
somaticInteractions(maf = luad, top = 10, pvalue = c(0.05, 0.1), 
                    fontSize = 0.9, sigSymbolsSize = 2, sigSymbolsFontSize = 2,
                    nShiftSymbols = 2,showSum = FALSE
                    )
#LUSC exclusive/co-occurance event analysis on top 10 mutated genes. 
somaticInteractions(maf = lusc, top = 10, pvalue = c(0.05, 0.1),
                    fontSize = 0.9, sigSymbolsSize = 2, sigSymbolsFontSize = 2 , 
                    nShiftSymbols = 2, showSum = FALSE)

#LUad exclusive/co-occurance of TP53 and KRAS
somaticInteractions(maf = luad, pvalue = c(0.05, 0.1),
                    genes = c("TP53", "KRAS"))
#LUSC exclusive/co-occurance of TP53 and KRAS
somaticInteractions(maf = lusc, pvalue = c(0.05, 0.1),
                    genes = c("TP53", "KRAS"))
#LUad exclusive/co-occurance of EGFR and KRAS
somaticInteractions(maf = luad, pvalue = c(0.05, 0.1),
                    genes = c("EGFR", "KRAS"))
#LUSC exclusive/co-occurance of EGFR and KRAS
somaticInteractions(maf = lusc, pvalue = c(0.05, 0.1),
                    genes = c("EGFR", "KRAS"))

#LUad exclusive/co-occurance of TP53 and EGFR
somaticInteractions(maf = luad, pvalue = c(0.05, 0.1),
                    genes = c("EGFR", "TP53"))
#LUSC exclusive/co-occurance of TP53 and EGFR
somaticInteractions(maf = lusc, pvalue = c(0.05, 0.1),
                    genes = c("EGFR", "TP53"))
#visualize the distribution of mutant genes in sample
oncoplot(maf = luad,genes = c("EGFR", "TP53","KRAS" ))
oncoplot(maf = lusc,genes = c("EGFR", "TP53","KRAS" ))

#variant allele frequency 
#LUAD
plotVaf(maf = luad)
plotVaf(maf = luad, gene_fs = 1, axis_fs = 1,
         height = 20, width = 20)
plotVaf(maf = luad, genes = c("KRAS","TP53"), height = 15, width = 15  )
#LUSC
plotVaf(maf = lusc, gene_fs = 1, axis_fs = 1)
plotVaf(maf = lusc, genes = c("KRAS","TP53"), height = 15, width = 15  )

#comparing LUAD and LUSC
pt.vs.rt <- mafCompare(m1 = luad, m2 = lusc, m1Name = 'LUAD', m2Name = 'LUSC', 
                       minMut = 5)
print(pt.vs.rt)
#TP53 
lollipopPlot2(m1 = luad, m2 = lusc, 
              gene = "TP53", 
              m1_name = "LUAD", m2_name = "LUSC",              
              domainLabelSize = 1.6, 
              legendTxtSize = 1.3, verbose = FALSE,
              colors = NULL, alpha = 1, 
              axisTextSize = c(1.3,1.3), pointSize = 1.2)

#KRAS 
lollipopPlot2(m1 = luad, m2 = lusc, 
              gene = "KRAS", 
              m1_name = "LUAD", m2_name = "LUSC",              
              domainLabelSize = 0.65, 
              legendTxtSize = 1.3, verbose = FALSE,
              colors = NULL, alpha = 1, 
              axisTextSize = c(1.3,1.3), pointSize = 1.2)


#Appendix
luad.titv = titv(maf = luad, plot = FALSE, useSyn = TRUE)
#plot titv summary
plotTiTv(res = luad.titv)

luad.sig = oncodrive(
  maf = luad,AACol = 'Protein_Change', 
  minMut = 5, pvalMethod = 'zscore')

plotOncodrive(res = luad.sig,
              fdrCutOff = 0.1, 
              useFraction = TRUE, labelSize = 0.5)




