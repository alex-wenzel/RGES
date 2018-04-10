"""
This is a command-line only version of RGES_perm_testing.ipynb
"""

import numpy as np
import pandas as pd

import json
from multiprocessing import Pool
import sklearn
import sys
sys.path.insert(0, '../')

from RGES.DiffEx import DiffEx
from RGES.L1KGCT import L1KGCTX
from RGES.Score import score

PHENOTYPE_SIGNATURE_PATH = "/mnt/oncogxA/Alex/l1k/DEG_SC_5um_entrezgene.txt"

#DRUG_PROFILE_PATH = "/mnt/oncogxA/Alex/l1k/10x_ilincs_sigs_top500_ranked_n500x978.gctx"
#DRUG_PROFILE_PATH = "/mnt/oncogxA/Alex/l1k/LINCS_FULL_GEO_RANKED/GSE70138_2017-03-06_landmarks_ranked_n118050x972.gctx"
DRUG_PROFILE_PATH = "/mnt/oncogxA/Alex/l1k/LINCS_CGS/GSE106127_CGS_ranked_n33839x978.gctx"

#OUTPATH = "ilincs_top500_permutations.json"
#OUTPATH = "GSE70138_perms.json"
OUTPATH = "GSE106127_perms.json"

DE = DiffEx(PHENOTYPE_SIGNATURE_PATH)
LINCS = L1KGCTX(DRUG_PROFILE_PATH)

merge_l2fc = lambda x: -1.0*x['log2fc.y'] if not np.isnan(x['log2fc.y']) else x['log2FoldChange']
DE.data['log2FoldChange'] = DE.data.apply(merge_l2fc, axis=1)

def mt_score_CHILD(signame):
    """Returns the RGES score for signame based on DE and LINCS"""
    return ((signame, score(DE, LINCS, signame)))

def mt_score(procs):
    """Returns a dictionary of {drug_profile: score}"""
    p = Pool(processes=procs)
    res = p.map(mt_score_CHILD, list(LINCS.data))
    p.close()
    p.join()
    return {r[0]: r[1] for r in res}

def shuffle_sigs():
    LINCS.data = LINCS.data.apply(sklearn.utils.shuffle, axis=0)
    LINCS.data.index = map(str, sorted(map(int, list(LINCS.data.index))))

PERMUTATIONS = 100
PROCESSES = 32

res_d = {signame: [] for signame in list(LINCS.data)}

for i in range(PERMUTATIONS):
    print("Starting permutation "+str(i))
    print("\tShuffling drug signature rankings...")
    shuffle_sigs()
    print("\tCalculating RGES for all profiles...")
    p_res = mt_score(PROCESSES)
    for signame in p_res.keys():
        res_d[signame].append(p_res[signame])
    open(OUTPATH, 'w').write(json.dumps(res_d))
open(OUTPATH, 'w').write(json.dumps(res_d))
