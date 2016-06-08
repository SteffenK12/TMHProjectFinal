from Helix import Helix
from ModuleMethods import is_transmembrane_helix
from PDBHelixParser import PDBHelixParser
from PDBextractor import SecondaryStructureCounter

__author__ = "Kevin Menden"
__date__ = '24.05.2016'

from sklearn import svm, cross_validation
import numpy as np
import matplotlib.pyplot as plt
from Bio.PDB.PDBParser import PDBParser
from XMLParser import parseXML
from sklearn.externals import joblib
from matplotlib import cm as cm



hydrophobic_residues = ["F", "G", "I", "L", "M", "V", "W", "Y"]
hydrophilic_residues = ["A", "C", "D", "E", "H", "K", "N", "P", "Q", "R", "S", "T"]

def calculate_hydrophobicity(seq):
    """
    Calculate the hydrophobicity factor for an amino acid sequence
    :param seq: the input amino acid sequence
    :return: the hydrophobicity factor
    """
    phobic = 0
    seq = seq.upper()
    for elem in seq:
        if elem in hydrophobic_residues:
            phobic += 1

    return phobic




def calculate_features(sequence):
    """
    Calculates all features of a given amino acid sequence
    :param sequence: amino acid sequence
    :return: array of features
    """
    hp_factor = calculate_hydrophobicity(sequence)
    length_factor = len(sequence)

    return [hp_factor, length_factor]

def make_training_set(tm_set, glob_set):
    """
    Make a training set in correct format for the SVM
    :param tm_set:
    :param glob_set:
    :return: the training set, the set of class labels
    """
    training = []
    lables = []

    # Append transmembrane set
    for elem in tm_set:
        training.append(calculate_features(elem))
        lables.append("tm")
    # Appen globular set
    for elem in glob_set:
        training.append(calculate_features(elem))
        lables.append("glob")

    return training, lables



tm_set = []

glob_set = []

# Parse PDB files and extract helices
structureCounter = SecondaryStructureCounter("globular_set")
structureCounter.parse_all_pdb_files()

new_helices = structureCounter.proteinHelixSequences
for prot in new_helices:
    for res in prot:
        glob_set.append(structureCounter.convertResiduesToLetters(res))
        # print(structureCounter.convertResiduesToLetters(res))

# Parse XML file containing helices
alpha_helices = parseXML("pdbtmall.xml")
print(len(alpha_helices))
for i in range(2500):
    tm_set.append(alpha_helices[i])


training, class_labels = make_training_set(tm_set, glob_set)

# Training
clf = svm.SVC(kernel='rbf', class_weight={"glob": 10})
clf.fit(X=training, y=class_labels)

# Cross-validation
scores = cross_validation.cross_val_score(clf, X=training, y=class_labels)
print("Accuracy: " + str(np.mean(cross_validation.cross_val_score(clf, X=training, y=class_labels))))


# # Make plot
# nf_tm = open("tm_file.txt", "w")
# nf_glob = open("glob_file.txt", "w")
# x1 = []
# y1 = []
# x2 = []
# y2 = []
# for elem in tm_set:
#     tmp = calculate_features(elem)
#     x1.append(tmp[0])
#     y1.append(tmp[1])
#     nf_tm.write(str(tmp[0]) + "\t" + str(tmp[1]) + "\n")
#
# for elem in glob_set:
#     tmp = calculate_features(elem)
#     x2.append(tmp[0])
#     y2.append(tmp[1])
#     nf_glob.write(str(tmp[0]) + "\t" + str(tmp[1]) + "\n")
#
# nf_tm.close()
# nf_glob.close()
# fig = plt.figure()
# tm_array = np.array([y1, x1])
# glob_array = np.array([y2, x2])
# max_length = max(max(y1), max(y2))
# H, xedges, yedges = np.histogram2d(tm_array[0], tm_array[1], bins=range(max_length))
# plt.pcolor(xedges, yedges, H, cmap="Blues")
#
# # ax2.colorbar()
# plt.show()



# plt.xlabel("Hydrophobicity factor")
# plt.ylabel("Length")
# plt.show()

# Save the model
joblib.dump(clf, "tmh_predictor_weighted.pkl")


