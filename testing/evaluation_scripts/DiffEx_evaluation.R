### This produces a plot of the sorted upregulated genes from
### testing/res.df.txt. This plot should show increasing log2 fold change

diffex.unit.upreg = read.csv("DiffEx_unit_upreg.tsv", sep='\t')
x.axis.vals = seq(1, nrow(diffex.unit.upreg))
plot(x.axis.vals, diffex.unit.upreg$log2FoldChange)

### This produces a plot of the sorted downregulated genes from
### testing/res.df.txt. This plot should show decreasing log2 fold change

diffex.unit.dnreg = read.csv("DiffEx_unit_dnreg.tsv", sep='\t')
x.axis.vals = seq(1, nrow(diffex.unit.dnreg))
plot(x.axis.vals, diffex.unit.dnreg$log2FoldChange)

### This code produces a histogram of expression values and expression z-scores

l1k.unit = read.csv("l1k_unit_zscore.tsv", sep='\t')
hist(l1k.unit$LINCSCP_1, breaks=200)
hist(l1k.unit$LINCSCP_1_zscore, breaks=200)
hist(l1k.unit$LINCSCP_1_zscore_norm, breaks=200)
hist(l1k.unit$LINCSCP_1_zscore_norm_pos, breaks=200)
