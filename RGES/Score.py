"""
This script implements methods for scoring a differential
expression result's enrichment for a LINCS drug profile
"""

import numpy as np
import pandas as pd

from DiffEx import DiffEx
from L1KGCT import L1KGCT

def get_a(de_join_prof, total_de_genes):
    """
    From the up or down regulated differential expression matrix
    joined with the drug profile, calculate the 'a' term according
    to Lamb et al 2006

        de_join_prof (pd.DataFrame): The differential expression combined
                                    with profile drug rank
        total_de_genes (int): The total number of genes in input differential
                                expression

        returns (float): The term for 'a' according to the above reference
    """
    terms = []
    total_de_genes = float(total_de_genes)
    t = float(len(de_join_prof.index))
    for j, row in de_join_prof.iterrows():
        terms.append((j/t) - (row['drug_rank']/total_de_genes))
    return max(terms)

def get_b(de_join_prof, total_de_genes):
    """
    From the up or down regulated differential expression matrix joined
    with the drug profile, calculate the 'b' term according to 
    Lambe et al 2006

        de_join_prof (pd.DataFrame): The differential expression combined
                                    with drug profile rank
        total_de_genes (int): The total number of genes in input differential
                                expression

        returns (float): The term for 'b' according to the above reference
    """
    terms = []
    total_de_genes = float(total_de_genes)
    t = float(len(de_join_prof.index))
    for j,row in de_join_prof.iterrows():
        terms.append((row['drug_rank']/total_de_genes) - ((j-1)/t))
    return max(terms)

def score(de_path, lincs_path):
    """
    Computes the RGES score for the differential expression stored
    at de_path and the lincs drug profile stored at lincs_path

        de_path (str): Location of a differential expression file
        lincs_path (str): Location of a LINCS drug profile

        returns (float): The RGES score
    """
    de = DiffEx(de_path)
    total_genes = len(de)
    prof = L1KGCT(lincs_path, normalized=True)
    up, dn = de.get_profile_order(prof)

    a_up = get_a(up, total_genes)
    b_up = get_b(up, total_genes)
    a_dn = get_a(dn, total_genes)
    b_dn = get_b(dn, total_genes)

    es_up = a_up if a_up > b_up else -1*a_dn
    es_dn = a_dn if a_dn > b_dn else -1*b_dn

    return es_up - es_dn

if __name__ == "__main__":
    print("Score.py")
    de_path = "testing/res.df.entrez.txt"
    lincs_path = "testing/LINCSCP_1.gct"
    print(score(de_path, lincs_path))
