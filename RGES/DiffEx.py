"""
This script defines a representation for a differential expression result. The
table is required to have the following headers:
    --baseMean (float)
    --log2FoldChange (float)
    --lfcSE (float)
    --stat (float)
    --pvalue (float)
    --padj (float)
    --gene_id (str)
    --symbol (str)
"""

import pandas as pd

class DiffEx:
    """
    Represents a file with differential expression data as described above
    """
    def __init__(self, path, sep='\t'):
        """
        Initializes a DiffEx instance and populates it with data from path

            path (str): Path to differential expression data
            sep (str): The delimiter for the data in path

            returns: None
        """
        self.data = pd.read_csv(path, sep=sep)

    """
    Accession
    """

    def get_up_genes(self):
        """
        Returns a subset of self.data with upregulated genes

            returns (pd.DataFrame): Upregulated genes in self.data
        """
        up_gene_rows = self.data.ix[self.data['log2FoldChange']>0]
        ug_rows_sorted = up_gene_rows.sort_values(by=['log2FoldChange'])
        return ug_rows_sorted['gene_id']  #return id because some symbols missing

    def get_down_genes(self):
        """
        Returns a subset of self.data with downregulated genes

            returns (pd.DataFrame): Downregulated genes in self.data
        """
        dn_gene_rows = self.data.ix[self.data['log2FoldChange']<0]
        ##NOTE unclear if downregulated should be ascending or descending
        dg_rows_sorted = dn_gene_rows.sort_values(by=['log2FoldChange'], ascending=False)
        return dg_rows_sorted[['gene_id', 'log2FoldChange']]  #return id because some symbols missing

if __name__ == "__main__":
    de = DiffEx("testing/res.df.txt")
    #print(de.get_up_genes()[:10])
    print(de.get_down_genes()[:10])
