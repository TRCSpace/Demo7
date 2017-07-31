# NodeGen.py
# This script takes the inputs specified above and generates nodes to solve the problem

import numpy

class NodeGen:

    def __init__(self):

        self.fuel = 'C3H8'
        self.oxidiser = 'O2'

        # i1:
        self.T_f = 100
        # i2:
        self.T_ox = 100
        # j
        self.OFR = 2
        # i3:
        self.T_r = (1/(self.OFR+1))*(self.T_f + self.OFR * self.T_ox)
        # k:
        self.P_c = 3.447e6
        # l:
        self.mdot_t = 1.25 # adjusted for approx 750 lbf thrust.

nodes = NodeGen()
