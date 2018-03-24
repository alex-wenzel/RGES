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
    described above. 
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
        self.perturbagen_id = None
        self.compound = None
        self.treatment = None
        self.concentration = None
        self.cell_line = None
        self.time = None
        self.factor = None
        self.data = pd.Series()

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
        """

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

    """
    Accession
    """

if __name__ == "__main__":
    print("L1KGCT.py")
