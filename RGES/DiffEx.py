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

from RGES.L1KGCT import L1KGCTX

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
        #self.data = self.data.sort_values(by=['log2FoldChange'], ascending=False)
        #self.data['phen_rank'] = list(range(1, len(self.data)+1))
        self.data['entrezgene'] = self.data['entrezgene'].astype(str)
        self.data['entrezgene'] = self.data['entrezgene'].apply(lambda x: x[:-2])
        #print(self.data[['log2FoldChange', 'phen_rank']])  #Debug

    """
    Accession
    """

    def get_up_genes(self):
        """
        Returns a subset of self.data with upregulated genes

            returns (pd.DataFrame): Upregulated genes in self.data
        """
        up_gene_rows = self.data.ix[self.data['log2FoldChange']>0]
        ug_rows_sorted = up_gene_rows.sort_values(by=['log2FoldChange'], ascending=False)
        return ug_rows_sorted[['entrezgene', 'log2FoldChange']]  #return id because some symbols missing

    def get_down_genes(self):
        """
        Returns a subset of self.data with downregulated genes

            returns (pd.DataFrame): Downregulated genes in self.data
        """
        dn_gene_rows = self.data.ix[self.data['log2FoldChange']<0]
        dg_rows_sorted = dn_gene_rows.sort_values(by=['log2FoldChange'], ascending=False)
        return dg_rows_sorted[['entrezgene', 'log2FoldChange']]  #return id because some symbols missing

    def get_profile_order(self, l1k_prof, signame):
        """
        Uses an L1000 drug profile to re-order the upregulated and downregulated genes

            l1k_prof (pd.DataFrame): Columns ID_geneid and profile data
            signame (str): Name of the signature in l1k_prof

            returns ((pd.DataFrame, pd.DataFrame)): 
        """
        cols2save = ['entrezgene', 'log2FoldChange', signame]

        up = self.get_up_genes()
        dn = self.get_down_genes()

        up_prof = pd.merge(up, l1k_prof, how='left', left_on='entrezgene', right_index=True)
        up_prof = up_prof[up_prof[signame].notnull()][cols2save]
        up_prof = up_prof.sort_values(by=signame, ascending=True)
        up_prof = up_prof.drop_duplicates(subset=["entrezgene"])
        up_prof.index = list(range(1, len(up_prof.index)+1))

        dn_prof = pd.merge(dn, l1k_prof, how='left', left_on='entrezgene', right_index=True)
        dn_prof = dn_prof[dn_prof[signame].notnull()][cols2save]
        dn_prof = dn_prof.sort_values(by=signame, ascending=True)
        dn_prof = dn_prof.drop_duplicates(subset=["entrezgene"])
        dn_prof.index = list(range(1, len(dn_prof.index)+1))

        return up_prof, dn_prof

    """
    Operators
    """

    def __len__(self):
        return len(self.data.index)

if __name__ == "__main__":
    DIR = "/scratch/alexw/l1k/"

    de = DiffEx("testing/res.df.entrez.txt")
    gctx = L1KGCTX("/scratch/alexw/l1k/LINCS_FULL_GEO/GSE70138_2017-03-06_landmarks_ranked_n118050x972.gctx")
    sig1 = list(gctx.data)[0]
    up_prof, dn_prof = de.get_profile_order(gctx.data, sig1)
    print(up_prof)
    print(dn_prof)
