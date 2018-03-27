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

from L1KGCT import L1KGCT

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
        self.data['entrezgene'] = self.data['entrezgene'].astype(str)
        self.data['entrezgene'] = self.data['entrezgene'].apply(lambda x: x[:-2])

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
        return ug_rows_sorted[['entrezgene', 'log2FoldChange']]  #return id because some symbols missing

    def get_down_genes(self):
        """
        Returns a subset of self.data with downregulated genes

            returns (pd.DataFrame): Downregulated genes in self.data
        """
        dn_gene_rows = self.data.ix[self.data['log2FoldChange']<0]
        dg_rows_sorted = dn_gene_rows.sort_values(by=['log2FoldChange'], ascending=False)
        return dg_rows_sorted[['entrezgene', 'log2FoldChange']]  #return id because some symbols missing

    def get_profile_order(self, l1k_prof):
        """
        Uses an L1000 drug profile to re-order the upregulated and downregulated genes

            l1k_prof (L1KGCT): A drug signature

            returns ((pd.DataFrame, pd.DataFrame)): 
        """
        cols2save = ['entrezgene', 'log2FoldChange', 'Name_GeneSymbol', l1k_prof.name, 'drug_rank']

        up = self.get_up_genes()
        dn = self.get_down_genes()
        up_prof = pd.merge(up, l1k_prof.data, how='left', left_on='entrezgene', right_on='ID_geneid')
        up_prof = up_prof[up_prof['drug_rank'].notnull()][cols2save]

        dn_prof = pd.merge(dn, l1k_prof.data, how='left', left_on='entrezgene', right_on='ID_geneid')
        dn_prof = dn_prof[dn_prof['drug_rank'].notnull()][cols2save]

        return up_prof, dn_prof

if __name__ == "__main__":
    DIR = "/scratch/alexw/l1k/"

    de = DiffEx("testing/res.df.entrez.txt")
    
    up = de.get_up_genes()
    dn = de.get_down_genes()
    up.to_csv(DIR+"DiffEx_unit_upreg.tsv", sep='\t', index=False, na_rep="NA")
    dn.to_csv(DIR+"DiffEx_unit_dnreg.tsv", sep='\t', index=False, na_rep="NA")

    l = L1KGCT("testing/LINCSCP_1.gct", normalized=True)
    up_prof, dn_prof = de.get_profile_order(l)
    up_prof.to_csv(DIR+"DiffEx+profile_upreg.tsv", sep='\t', index=False, na_rep="NA")
    dn_prof.to_csv(DIR+"DiffEx+profile_dnreg.tsv", sep='\t', index=False, na_rep="NA")
