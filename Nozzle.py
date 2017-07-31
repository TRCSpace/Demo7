from Combustor import Combustor
from Component import Component
import Nodes

import math
from solid import *
from solid.utils import *
from pprint import pprint

class Nozzle(Component):

    def __init__(self, nodes, combThermo):

        self.P_c = nodes.P_c
        # First, need procedure to fix thr.P (by fixing gamma!):
        self.gm10 = combThermo.gamma
        self.gm15 = self.gm10  # Initial guess
        self.err15 = 1  # Relative error
        self.gas = combThermo.gas

        # Initial guess for pressure:
        self.P20g = self.P_c * (1 + 0.5 * (self.gm10 - 1)) ** (-1 * self.gm10 / (self.gm10 - 1))
        self.T20g = combThermo.T0 / (1 + 0.5 * (self.gm10 - 1))

        # 10 is combustor 20 is throat 30 is exit

        while self.err15 > 0.001:
            # Initial calculation for throat equilibrium point: Let the initial conditions be those in the combustor
            self.gas = combThermo.gas

            # Expand the gas to estimated throat T & P. Then observe the change in properties as the gas releases more
            # heat as it reacts
            self.gas.TP = self.T20g, self.P20g
            self.gas.equilibrate('HP')

            self.gm20 = self.gas.cp_mass / self.gas.cv_mass
            self.gm15n = 0.5 * (self.gm10 + self.gm20)
            # P20g = P_c*(1+0.5*(gm10-1))**(-gm10/(gm10-1))
            self.P20g = self.P_c * (1 + 0.5 * (self.gm15n - 1)) ** (-self.gm15n / (self.gm15n - 1))
            self.T20g = combThermo.T0 / (1 + 0.5 * (self.gm15n - 1))
            self.err15 = abs(self.gm15n - self.gm15) / self.gm15n
            self.gm15 = self.gm15n

        self.P20 = self.gas.P
        self.T20 = self.gas.T
        self.P0 = self.P_c
        self.T0 = combThermo.T0
        self.cp20 = self.gas.cp_mass
        self.gm20 = self.gas.cp_mass / self.gas.cv_mass
        self.Mr20 = self.gas.mean_molecular_weight
        self.rho20 = self.gas.density
        self.c020 = math.sqrt(self.gm20 * 8314 * self.T0 / self.Mr20)
        self.v20 = self.c020
        self.mu20 = self.gas.viscosity
        self.k_th20 = self.gas.thermal_conductivity
        self.Pr20 = self.cp20 * self.mu20 / self.k_th20

        self.Y20 = self.gas.mass_fraction_dict(1e-6)
        self.X20 = self.gas.mole_fraction_dict(1e-6)

        self.mbar_thr20 = (self.gm20 / (math.sqrt(self.gm20 - 1))) * \
                        ((0.5 * (self.gm20 + 1)) ** (-0.5 * (self.gm20 + 1) / (self.gm20 - 1)))

        ## Exit Equilibrium
        # First, need procedure to fix gamma!:
        self.gm30 = self.gm20
        self.gm25 = 0.5 * (self.gm20 + self.gm30)
        self.err25 = 1

        self.T20eq = self.T20
        self.P20eq = self.P20
        self.gm20eq = self.gm20
        self.T020eq = self.T20eq * (1 + 0.5 * (self.gm20 - 1))
        self.P020eq = self.P20eq * ((1 + 0.5 * (self.gm20 - 1)) ** (self.gm20 / (self.gm20 - 1)))

        # self.P0 = T20eq*(1+0.5*(gm20-1))
        # self.T0 = P20eq*((1+0.5*(gm20-1))**(gm20/(gm20-1)))

        # Initial guess for gamma & Mach number:
        # P_e > 0.4 atm to avoid flow separation (Summerfield Criterion)
        self.P30 = (5 / 14.7) * 101325  # Set arbitrarily at 5 psi
        self.M30g = math.sqrt((((self.P020eq / self.P30) ** ((self.gm25 - 1) / self.gm25)) - 1) * (2 / (self.gm25 - 1)))
        self.T30g = self.T0 / (1 + (0.5 * (self.M30g ** 2) * (self.gm25 - 1)))
        # self.gas = thrThermo.gas

        while self.err25 > 0.001:
            # Initial calculation for throat equilibrium point:
            self.gas.TP = [self.T30g, self.P30]
            self.gas.equilibrate('HP')

            self.gm30 = self.gas.cp_mass / self.gas.cv_mass
            self.gm25n = 0.5 * (self.gm20 + self.gm30)

            self.M30g = math.sqrt((((self.P020eq / self.P30) ** ((self.gm25n - 1) / self.gm25n)) - 1) *
                                  (2 / (self.gm25n - 1)))
            self.T30g = self.T0 / (1 + (0.5 * (self.M30g ** 2) * (self.gm25 - 1)))

            self.err25 = abs(self.gm25 - self.gm25n) / self.gm25
            self.gm25 = self.gm25n

        self.P30 = self.P30
        self.T30 = self.gas.T
        self.M30 = self.M30g
        self.c030 = math.sqrt(self.gm30 * 8314 * self.T0 / self.gas.mean_molecular_weight)
        self.v30 = self.M30g * self.c030
        self.T0 = self.T0
        self.cp30 = self.gas.cp_mass
        self.gm30 = self.cp30 / self.gas.cv_mass
        self.Mr30 = self.gas.mean_molecular_weight
        self.rho30 = self.gas.density
        self.mu30 = self.gas.viscosity
        self.k_th30 = self.gas.thermal_conductivity
        self.Pr30 = self.cp30 * self.mu30 / self.k_th30

        self.mbar30 = (self.gm30 * self.M30 / (math.sqrt(self.gm30 - 1))) * \
                         ((0.5 * (self.gm30 + 1) * self.M30 ** 2) ** (-0.5 * (self.gm30 + 1) / (self.gm30 - 1)))
        self.Y30 = self.gas.mass_fraction_dict(1e-6)
        self.X30 = self.gas.mole_fraction_dict(1e-6)


# Change to make nozzle!!
def RenderToModel(combLength, combDiam, combThickness):
    # Your code here!
    # converting from metres to mm:
    combLength = combLength * 1000
    combDiam = combDiam * 1000
    combThickness = combThickness * 1000

    combustor = cylinder(r=combDiam/2, h=combLength)
    combustor -= translate([0, 0, -1])(cylinder(r=combDiam/2, h=combLength+2))
    combustor = translate([0, 0, -combLength])(combustor)

    return combustor


def print_r(the_object):
    print("CLASS: ", the_object.__class__.__name__, " (BASE CLASS: ", the_object.__class__.__bases__, ")")
    pprint(vars(the_object))


# Start of script:

combustor = Combustor(Nodes.nodes)
nozzle = Nozzle(Nodes.nodes, combustor)


if __name__ == '__main__':
    combCAD = RenderToModel(combustor.L_tot_combustor, combustor.d_comb, 0.02)
    SEGMENTS = 48
    scad_render_to_file(combCAD, file_header='$fn = %s;' % SEGMENTS, include_orig_code=True)
