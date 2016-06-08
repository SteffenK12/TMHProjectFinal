import os
import sys

from Helix import Helix
from ModuleMethods import is_transmembrane_helix
from PDBHelixParser import PDBHelixParser

__author__ = "Kevin Menden"
__date__ = '08.06.2016'
"""
Validation of the model
"""


def isTransmembraneProtein(file):
    # parser = argparse.ArgumentParser(description="Membrane Plane Finder")
    # parser.add_argument('pdb')
    pdbParser = PDBHelixParser(file)
    pdbParser.parse_pdb_file()
    structure = pdbParser.structure         # The whole structure of the PDB file
    print("Parsed PDB file succesfully")
    raw_helices = pdbParser.proteinHelixSequences

    # Convert raw helices into Helix objects
    helix_set = []
    for h in raw_helices:
        helix_set.append(Helix(h))
    print("Found " + str(len(helix_set)) + " helices")

    # Predict transmembrane helices
    tmh_set = []
    for h in helix_set:
        if is_transmembrane_helix(h) == 'tm':
            tmh_set.append(h)
    print(len(tmh_set) / len(helix_set))
    if len(tmh_set) < 3 or len(tmh_set)/len(helix_set) < 0.1:
        return False
    else:
        return True


if __name__ == "__main__":

    tm = 0
    glob = 0

    directory = "validation_set_tm"
    for subdir, dirs, files in os.walk(directory):
        for file in files:
            if file.endswith(".pdb"):
                f = directory + "/" + file
                if isTransmembraneProtein(f):
                    tm += 1
                else:
                    glob += 1

    print("Tramsmembrane: " + str(tm))
    print("Globular: " + str(glob))