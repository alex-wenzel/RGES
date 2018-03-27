"""
This script defines a representation for an input .gct file containing
L1000-assay mRNA expression data. This representation assumes the format
as provided by the iLINCS browser, which can be downloaded here:
http://www.ilincs.org/ilincs/signaturesL1000/LDG-1188/search/
"""

import pandas as pd

class L1KGCT:
    """
    This class represents a .gct file for L1000 mRNA expression as
    described above. This class should only be used to load .gct files that
    contain a single LINCS signature. For multi-signature files, use
    MultiL1KGCT below
    """
    def __init__(self, path, normalized=False):
        """
        Initializes an L1KGCT instance and calls the load function on the
        input path

            path (str): A path to a GCT with L1000 data
            normalize (bool): If True, apply normalization as described
                                below in self.normalize()

            returns: None
        """
        self.name = ""
        self.metadata = {}
        self.data = None

        self.path = path
        self.normalized = normalized

        self.load()
        if self.normalized:
            self.normalize()

    """
    Loading Data
    """

    def load(self):
        """
        Reads data from self.path and populates self.data and metadata fields

            returns: None
        """
        self.data = pd.read_csv(self.path, sep='\t', skiprows=[0,1,3,4,5,6,7,8,9])
        self.data['ID_geneid'] = self.data['ID_geneid'].astype(str)
        for i, line in enumerate(open(self.path)):
            if i in [0,1]:
                continue
            if i > 9:
                break
            lv = line.strip('\n').split('\t')
            self.metadata[lv[0]] = lv[-1]
        self.name = self.metadata['id']
        self.data = self.data.sort_values(by=self.name)
        self.data['drug_rank'] = range(1, len(self.data)+1) 

    """
    Normalize Data
    """

    def normalize(self):
        """
        This function transforms self.data such that each value
        becomes a z-score constrained between [-3, 3] and then shifted by 6
        such that all values are in the interval [0, 6]

            returns: None
        """
        pname = self.metadata['id']
        raw = self.data[pname]
        self.data[pname+'_zscore'] = (raw - raw.mean())/raw.std(ddof=0)
        self.data[pname+'_zscore_norm'] = self.data[pname+'_zscore']
        self.data.loc[self.data[pname+'_zscore_norm'] > 3, pname+'_zscore_norm'] = 3
        self.data.loc[self.data[pname+'_zscore_norm'] < -3, pname+'_zscore_norm'] = -3
        self.data[pname+'_zscore_norm_pos'] = self.data[pname+'_zscore_norm']+3

    """
    Accession
    """

class MultiL1KGCT:
    """
    This class represents a .gct file for 1 or more signatures from the
    L1000 assay drug signature database. 
    """
    def __init__(self, path, normalized=False):
        """
        Initalizes a MultiL1KGCT instance and call the load function
        on path

            path (str): A path to a file with one more L1000 signatures
            normalize (bool): If true, apply normalization as described below
                            in self.normalize()

            returns: None
        """
        self.names = []
        self.metadata = {}  #{sig: {metadata}}
        self.data = None  #pd.DataFrame

        self.path = path
        self.normalized = normalized

        self.load()
        if self.normalized:
            self.normalize()

    """
    Laoding Data
    """

    def load(self):
        """
        Reads data from self.path and populates self.data and metadata fields
        for each LINCS signature in the file

            returns: None
        """
        sig_md = {}
        self.data = pd.read_csv(self.path, sep='\t', skiprows=[0,1,3,4,5,6,7,8,9])
        self.data['ID_geneid'] = self.data['ID_geneid'].astype(str)
        for i, line in enumerate(open(self.path)):
            if i in [0,1]:
                continue
            if i > 9:
                break
            lv = line.strip('\n').split('\t')
            sig_md[lv[0]] = lv[4:]
        for i, sig in enumerate(sig_md['id']):
            self.metadata[sig] = {key: sig_md[key][i] for key in sig_md.keys()
                                        if key != "id"}

    """
    Normalize Data
    """

    def normalize(self):
        """
        This function transforms self.data such that each value
        becomes a z-score constrained between [-3, 3] and then shifted by 6
        such that all values are in the interval [0, 6]

            returns: None
        """
        zscore = lambda x: (x-x.mean())/x.std(ddof=0)
        col_d = {sig: zscore for sig in self.metadata.keys()}
        self.data = self.data.transform(col_d)

        norm = lambda x: 6 if x>3 else (0 if x<-3 else x+3)
        col_d = {sig: norm for sig in self.metadata.keys()}
        self.data = self.data.transform(col_d)

if __name__ == "__main__":
    DIR = "/scratch/alexw/l1k/"    

    l = L1KGCT("testing/LINCSCP_1.gct", normalized=True)
    l.data.to_csv(DIR+"l1k_unit_zscore.tsv", sep='\t', index=False, na_rep="NA")

    ml = MultiL1KGCT("testing/CTPRES_100_concordant_sigs.gct", normalized=True)
    ml.data.to_csv(DIR+"mul1k_unit_zscore.tsv", sep='\t', index=False, na_rep="NA")
