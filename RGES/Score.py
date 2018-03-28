"""
This script implements methods for scoring a differential
expression result's enrichment for a LINCS drug profile
"""

import numpy as np
import pandas as pd

from RGES.DiffEx import DiffEx
from RGES.L1KGCT import L1KGCT, MultiL1KGCT

def get_a(de_join_prof, total_de_genes, signame):
    """
    From the up or down regulated differential expression matrix
    joined with the drug profile, calculate the 'a' term according
    to Lamb et al 2006

        de_join_prof (pd.DataFrame): The differential expression combined
                                    with profile drug rank
        total_de_genes (int): The total number of genes in input differential
                                expression
        signame (str): The name of the drug signature input

        returns (float): The term for 'a' according to the above reference
    """
    terms = []
    total_de_genes = float(total_de_genes)
    t = float(len(de_join_prof.index))
    for j, row in de_join_prof.iterrows():
        terms.append((j/t) - (row[signame+'_drug_rank']/total_de_genes))
    return max(terms)

def get_b(de_join_prof, total_de_genes, signame):
    """
    From the up or down regulated differential expression matrix joined
    with the drug profile, calculate the 'b' term according to 
    Lambe et al 2006

        de_join_prof (pd.DataFrame): The differential expression combined
                                    with drug profile rank
        total_de_genes (int): The total number of genes in input differential
                                expression
        signame (str): The name of the drug signature input

        returns (float): The term for 'b' according to the above reference
    """
    terms = []
    total_de_genes = float(total_de_genes)
    t = float(len(de_join_prof.index))
    for j,row in de_join_prof.iterrows():
        terms.append((row[signame+'_drug_rank']/total_de_genes) - ((j-1)/t))
    return max(terms)

def score(de, lincs_sigs, signame):
    """
    Computes the RGES score for the differential expression stored
    at de_path and the lincs drug profile stored at lincs_path

        de (DiffEx): Location of a differential expression file
        lincs_sigs (MultiL1KGCT): Signature data

        returns (float): The RGES score
    """
    total_genes = len(de)
    sig = lincs_sigs.data[['Name_GeneSymbol', 'ID_geneid', 
                            signame, signame+'_drug_rank']]
    up, dn = de.get_profile_order(sig, signame)

    a_up = get_a(up, total_genes, signame)
    b_up = get_b(up, total_genes, signame)
    a_dn = get_a(dn, total_genes, signame)
    b_dn = get_b(dn, total_genes, signame)

    print("====================================")  #Debug
    print(a_up, b_up, a_dn, b_dn)  #Debug

    es_up = a_up if a_up > b_up else -1*b_up
    es_dn = a_dn if a_dn > b_dn else -1*b_dn
    
    print(es_up, es_dn)  #Debug

    return es_up - es_dn

if __name__ == "__main__":
    print("Score.py")
    de = DiffEx("testing/res.df.entrez.txt")
    lincs_path = "testing/CTPRES_100_concordant_sigs.gct"
    lincs_sigs = MultiL1KGCT(lincs_path, normalized=True)
    sig_name = list(lincs_sigs.metadata.keys())[0]
    print(score(de, lincs_sigs, sig_name))
