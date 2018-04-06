"""
This script implements methods for scoring a differential
expression result's enrichment for a LINCS drug profile
"""

#import humanfriendly as hf
import numpy as np
import pandas as pd
import time

from RGES.DiffEx import DiffEx
from RGES.L1KGCT import L1KGCTX, MultiL1KGCT

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
    #t = 1417.0  #Debug
    for j, row in de_join_prof.iterrows():
        terms.append((j/t) - (row[signame]/total_de_genes))
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
    #t = 1417.0  #Debug
    for j, row in de_join_prof.iterrows():
        terms.append((row[signame]/total_de_genes) - ((j-1)/t))
    return max(terms)

def score(de, lincs_sigs, signame):
    """
    Computes the RGES score for the differential expression stored
    at de_path and the lincs drug profile stored at lincs_path

        de (DiffEx): Location of a differential expression file
        lincs_sigs (MultiL1KGCT): Signature data

        returns (float): The RGES score
    """
    total_genes = len(lincs_sigs.data.index)
    sig = lincs_sigs.data[[signame]]
    sig = sig.sort_values(by=signame, ascending=False)
    sig = sig.rank(method='first', ascending=False)
    up, dn = de.get_profile_order(sig, signame)
    
    a_up = get_a(up, total_genes, signame)
    b_up = get_b(up, total_genes, signame)
    a_dn = get_a(dn, total_genes, signame)
    b_dn = get_b(dn, total_genes, signame)

    es_up = a_up if a_up > b_up else -1*b_up
    es_dn = a_dn if a_dn > b_dn else -1*b_dn

    return es_up - es_dn

if __name__ == "__main__":
    print("Score.py")
    de = DiffEx("testing/res.df.entrez.txt")
    #lincs_path = "testing/CTPRES_100_concordant_sigs.gct"
    #lincs_sigs = MultiL1KGCT(lincs_path, normalized=True)
    lincs_path = "/scratch/alexw/l1k/LINCS_FULL_GEO/GSE70138_2017-03-06_landmarks_ranked_n118050x972.gctx"
    lincs_sigs = L1KGCTX(lincs_path)
    #times = []
    #for signame in lincs_sigs.metadata.keys():
    for signame in list(lincs_sigs.data):
        t0 = time.time()
        s = score(de, lincs_sigs, signame)
        t1 = time.time()
        #times.append(t1-t0)
        #print(t1-t0)
        #print(s)
        #break
    #LINCS_L = 68960.0
    #print("Mean time per score: "+hf.format_timespan(np.mean(times)))
    #print("Max score time: "+hf.format_timespan(max(times)))
    #print("Min score time: "+hf.format_timespan(min(times)))
    #time_for_lincs = np.mean(times)*LINCS_L
    #print("Would score all of iLINCS ("+str(LINCS_L)+") on this signature in "+str(hf.format_timespan(time_for_lincs)))
