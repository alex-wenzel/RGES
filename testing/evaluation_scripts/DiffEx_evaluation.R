### This produces a plot of the sorted upregulated genes from
### testing/res.df.txt. This plot should show increasing log2 fold change

diffex.unit.upreg = read.csv("DiffEx_unit_upreg.tsv", sep='\t')
plot(diffex.unit.upreg$log2FoldChange)

### This produces a plot of the sorted downregulated genes from
### testing/res.df.txt. This plot should show decreasing log2 fold change

diffex.unit.dnreg = read.csv("DiffEx_unit_dnreg.tsv", sep='\t')
plot(diffex.unit.dnreg$log2FoldChange)

### This code produces a histogram of expression values and expression z-scores

l1k.unit = read.csv("l1k_unit_zscore.tsv", sep='\t')
hist(l1k.unit$LINCSCP_1, breaks=200)
hist(l1k.unit$LINCSCP_1_zscore, breaks=200)
hist(l1k.unit$LINCSCP_1_zscore_norm, breaks=200)
hist(l1k.unit$LINCSCP_1_zscore_norm_pos, breaks=200)

### This code checks the sorting for the LINCS profile, this should plot
### genes by increasing differential expression

plot(l1k.unit$LINCSCP_1)

### DiffEx+profile_upreg.tsv should contain the phenotype expression
### joined with the drug profile, sorted according to ascending abs(diff_ex)

diffex.up.join.prof = read.csv("DiffEx+profile_upreg.tsv", sep='\t')
diffex.dn.join.prof = read.csv("DiffEx+profile_dnreg.tsv", sep='\t')

### This code gets the gene symbols that are present l1k.unit
### but not in diffex.up.join.prof or diffex.dn.join.prof

all.gene.sym.l1k = as.character(l1k.unit$Name_GeneSymbol)
all.gene.sym.up.dn = c(as.character(diffex.up.join.prof$Name_GeneSymbol),
                       as.character(diffex.dn.join.prof$Name_GeneSymbol))
setdiff(all.gene.sym.l1k, all.gene.sym.up.dn)
