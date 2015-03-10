"""
This module computate BDT Score
"""

import sys
import os
import logging
import ROOT
import pickle
import numpy as np

class BDTScore(object):

    def __init__(self, year):

        self.BDT_DIR = str(os.getenv('HIGGSTAUTAU_DIR'))
        self.BDT_DIR += '/bdts/'
        self.clf_suffix = '_125_NONISOL_ebz_%s_' % ( year % 1000 )

        self.clfs_vbf = [None,None]
        self.clfs_boosted = [None,None]


        ## VBF Classifier
        clf_vbf_0_name = self.BDT_DIR + 'clf_vbf' + self.clf_suffix + '0.pickle'
        if os.path.isfile(clf_vbf_0_name):
            with open(clf_vbf_0_name, 'rb') as f:
                self.clfs_vbf[0] = pickle.load(f)
        else :
            raise ValueError("could not open {0}".format(clf_vbf_0_name))

        clf_vbf_1_name = self.BDT_DIR + 'clf_vbf' + self.clf_suffix + '1.pickle'
        if os.path.isfile(clf_vbf_1_name):
            with open(clf_vbf_1_name, 'rb') as f:
                self.clfs_vbf[1] = pickle.load(f)
        else :
            raise ValueError("could not open {0}".format(clf_vbf_1_name))

         ## Boosted Classifier
        clf_boosted_0_name = self.BDT_DIR + 'clf_boosted' + self.clf_suffix + '0.pickle'
        if os.path.isfile(clf_boosted_0_name):
            with open(clf_boosted_0_name, 'rb') as f:
                self.clfs_boosted[0] = pickle.load(f)
        else :
            raise ValueError("could not open {0}".format(clf_boosted_0_name))

        clf_boosted_1_name = self.BDT_DIR + 'clf_boosted' + self.clf_suffix + '1.pickle'
        if os.path.isfile(clf_boosted_1_name):
            with open(clf_boosted_1_name, 'rb') as f:
                self.clfs_boosted[1] = pickle.load(f)
        else :
             raise ValueError("could not open {0}".format(clf_boosted_1_name))

    def get_score(self, EventNumber, partition, category, inputs):

        if category == 'VBF':
            clf = self.clfs_vbf[(EventNumber + partition) % 2]
        elif category == 'Boosted':
            clf = self.clfs_boosted[(EventNumber + partition) % 2]
        else :
            raise ValueError("No classifier for {0} cateogry".format(category))

        # logistic tranformation used by TMVA (MethodBDT.cxx)
        score = -1 + 2.0 / (1.0 +
                            np.exp(-clf.n_estimators *
                                    clf.learning_rate *
                                    clf.decision_function(inputs) / 1.5))

        return score



