# -*- coding: utf-8 -*-
"""
@author: Yann Jobic
"""

import numpy as np
import scipy
from scipy.sparse import coo_matrix,coo_array

#Quand on ajoute, c'est dans l'ordre I,J,val
#et dans data c'est val,I,J.
#A changer

class Triplets:
    def __init__(self):
        self.data = (np.array([], dtype=float), (np.array([], dtype=int), np.array([], dtype=int)))
    def __str__(self):
        return str(self.data)
    def append(self, I, J, val):
        self.data = (np.append(self.data[0],val),
                 (np.append(self.data[1][0],I),np.append(self.data[1][1],J)))
    def to_coo_array(self):
        return coo_array(self.data)

