# ! /usr/bin/env python
# -*- coding: UTF-8 -*-
from __future__ import division
import os
import sys
import re
from solid import *
from solid.utils import *

from pprint import pprint
import cantera as ct
import CoolProp as cp

from Component import Component
import Nodes


class Combustor(Component):

    def __init__(self, nodes):

        nodes.T_r = nodes.T_r
        self.OFR = nodes.OFR
        self.P_c = nodes.P_c
        self.mdot_t = nodes.mdot_t
        self.fuel = nodes.fuel
        self.oxidiser = nodes.oxidiser
        # Beginning Cantera analysis (equilibrium combustion analysis):
        self.gas = ct.Solution('gri30.xml')
        self.gas.transport_model = 'Multi'
        self.ct_OFR_string = self.fuel + ':1, ' + self.oxidiser + ':' + str(self.OFR)
        self.gas.Y = self.ct_OFR_string
        self.gas.TP = nodes.T_r, self.P_c
        self.gas.equilibrate('HP')
        # Extracting gas properties from the gas Solution object - storing them permanently, so that ct doesn't alter
        # them when we call it in the future
        self.gamma = self.gas.cp_mass / self.gas.cv_mass
        self.cp = self.gas.cp_mass
        self.mu = self.gas.viscosity
        nodes.T0 = self.gas.T
        self.P0 = self.gas.P
        self.Mr = self.gas.mean_molecular_weight
        self.rho = self.gas.density_mass
        self.c0 = math.sqrt(self.gamma * 8314 * nodes.T0 / self.Mr)
        self.k_th = self.gas.thermal_conductivity
        self.Pr = self.cp * self.mu / self.k_th
        self.Y = self.gas.mass_fraction_dict(1e-6)
        self.X = self.gas.mole_fraction_dict(1e-6)

        # Placeholder variables - will be operated on by other functions:
        nodes.T = nodes.T0
        self.P = self.P0
        self.M = 0
        self.v = self.M * self.c0

        self.d_inj_hole = 1e-3
        self.rho_amb = 2.5  # estimated value for combustion gases
        # Half angle between like impinging jets:
        self.theta_jet = 20

        # Determining Required velocity for jet atomisation:
        # Wrap into a function
        self.mdot_f = nodes.mdot_t * (1 / (nodes.OFR + 1))
        # Propellant Density:
        self.rho_f = cp.CoolProp.PropsSI('D', 'P', self.P_c, 'T', nodes.T_r, nodes.fuel)
        # Propellant viscosity:
        self.mu_f = cp.CoolProp.PropsSI('viscosity', 'P', self.P_c, 'T', nodes.T_r, nodes.fuel)
        # Propellant Surface tension:
        self.sigma_f = cp.CoolProp.PropsSI('surface_tension', 'P', self.P_c, 'T', nodes.T_r, nodes.fuel)

        # These constants are the slope and y-intercept of the log-linear fit for the Reynolds - Ohnesorge fit for
        # instantaneous droplet atomisation conditions (bit of a mouthful!)
        self.m = -1.3315
        self.c = 3.2621
        # m and c will need to be updated if we need to look for wave breakup instead of instantaneous atomisation, as
        # the relative velocities for this breakup regime may result in unstable combustion (risk of blowout)

        # Fuel Ohnesorge Number, based on droplet diameter:
        self.Oh_f = self.mu_f / (math.sqrt(self.rho_f * self.sigma_f * self.d_inj_hole))
        # log base 10 of fuel Reynolds Number required for instant atomisation:
        self.log10Re_f = (self.c - math.log10(self.Oh_f)) / (-self.m)
        # Fuel Reynolds Number:
        self.Re_f = 10 ** self.log10Re_f
        # Critical velocity (required) for instant atomisation of droplet
        self.u_cr = self.Re_f * self.mu_f / (self.rho_f * self.d_inj_hole)

        # Cd_sphere determination: wrap into a function
        self.Cd_sphere = 0.14  # good approximation for turbulent flow (Re > 10**6)

        # Temporary measure: taken value manually from main script. Will need to be
        # linked to Combustor.rho
        # rho_amb = 2.4831

        # Inviscid Weber Number:
        self.We_crit_inviscid = 8 / self.Cd_sphere
        self.d_fmax = 8 * self.sigma_f / (self.Cd_sphere * self.rho_amb * self.u_cr ** 2)
        self.Oh_f = self.mu_f / math.sqrt(self.rho_f * self.sigma_f * self.d_fmax)

        # Correlation below has a 20# error - find a better one if possible
        self.We_jet_crit = self.We_crit_inviscid + (14 * self.Oh_f)
        # Accounting for the error that may exist:
        self.We_jet_crit_max = 1.2 * self.We_jet_crit
        self.We_jet_crit_min = 0.8 * self.We_jet_crit

        # # Conduct secondary droplet breakup calculations:
        # Droplet secondary breakup: has been wrapped into a function - new
        # function requires testing before formal integration

        self.u_rel = self.u_cr * 1.0
        self.u_jet = self.u_rel / (2 * math.sin(math.radians(self.theta_jet)))

        # gas-based Weber number - (dynamic pressure/surface tension) force ratio
        # We_g  = rho_g *(u_jf **2)*d_jf /sigma_f
        # Calculating critical Weber number:
        # Source: https://www.researchgate.net/publication/257803315

        # Empirical constants:
        self.Cd = 5  # Note: NOT drag coefficient!!
        self.Ck = 8
        self.d_pri = self.d_inj_hole
        self.r_pri = self.d_pri / 2  # primary droplet radius

        self.tau = 2 * self.rho_f * (self.r_pri ** 2) / (self.Cd * self.mu_f)
        self.omega = math.sqrt((self.Ck * self.sigma_f / (self.rho_f * self.r_pri ** 3)) - (1 / (self.tau ** 2)))

        # Working out velocity components:
        # Need to develop a way of capturing the fan spread of the droplets - if we
        # can do so, then it will mean that we can really conserve momentum, and
        # work out more accurate values for droplet axial and radial valocities
        # (assuming no azimuthal component)

        # Currently, if momentum conservation is applied, and coalescence is
        # assumed, then the result a straight axial jet.

        # Known that a straight turbulent axial jet will diverge with half angle of
        # about 12 degrees. Could use this as an initial approximation of the fan
        # half angle. Still probably a better method out there than this crude
        # approximation. Makes conceptual sense that the fan half angle would be
        # larger...
        self.theta_fan = 12
        # Conceptually, fan half angle will be somewhere between 12 degrees and the
        # the impinging jet half angle.

        # Maximum radial component - at the edge of the fan

        # Based on assumption that velocity magnitude remains constant for all
        # droplets (i.e. only change in direction)
        self.u_r2max = math.sin(math.radians(self.theta_jet)) * self.u_jet
        # At the fan centre, the radial component is assumed to be zero
        self.u_r2min = 0
        self.u_z2min = math.cos(math.radians(self.theta_jet)) * self.u_jet
        self.u_z2max = self.u_jet
        self.u_r2_avg = 0.5 * (self.u_r2max + self.u_r2min)
        self.u_z2_avg = 0.5 * (self.u_z2max + self.u_z2min)

        # Developing a (Hybrid) Timescale:
        self.Lscale = math.pi * self.d_pri  # relevant turbulent length scale
        # relative turbulence intensity (true, but add source to backup)
        self.I_r = 0.25
        self.t_sbu = math.pi / self.omega  # breakup timescale
        self.u_rel = self.u_cr

        # turbulent flow timescale
        self.t_turb = (1 / self.I_r) * (self.Lscale / self.u_rel)
        # hybrid (combined) timescale
        self.t_hyb = ((1 / self.t_sbu) + (1 / self.t_turb)) ** -1

        # A,B,C,D,E are constituent variables created to make the code more
        # readable
        self.A = 1 * (1 + math.exp(-math.pi / (self.omega * self.tau)))
        self.B = math.exp(-self.t_hyb / self.tau)
        self.C = self.omega * self.t_hyb
        self.D = math.cos(self.C) + (math.sin(self.C) / (self.omega * self.tau))
        self.E = (4 * math.pi / self.A) * (1 - (self.D * self.B))

        # Critical Weber number (accounting for turbulence and density effects)
        self.We_crit = self.E * (1 + (self.rho_amb / self.rho_f))
        # Max stable droplet diameter (theoretical)
        self.d_sbu = self.We_crit * self.sigma_f / (self.rho_amb * (self.u_jet ** 2))
        # Breakup length, based on average velocity (u_r and u_d)
        self.L_sbu = self.u_jet * self.t_sbu


        ## Function to output combustion times
        # Input structs: Nodes, Reactant.fuel,Pintle.output,DBU.fuel

        # function[output] = ArrayFuelDropletComb(Nodes, ReactantSpecies, Combustor, PintleOutput, DBUProp)

        # u_f = PintleOutput.vel_f     # Unused at present
        # u_ox = PintleOutput.vel_ox   # Unused at present

        nodes.T_prop = nodes.T_f
        self.P_c = nodes.P_c

        ## Declare Initial Conditions:
        nodes.T_b = cp.CoolProp.PropsSI('T','P',self.P_c,'Q',0,nodes.fuel)
        self.rho_f = cp.CoolProp.PropsSI('D', 'P', self.P_c, 'T', nodes.T_prop, nodes.fuel)
        self.cp_f = cp.CoolProp.PropsSI('Cpmass', 'P', self.P_c, 'T', nodes.T_prop, nodes.fuel)
        self.h_fg = abs(cp.CoolProp.PropsSI('H','P',self.P_c,'Q',1,nodes.fuel) -
                        cp.CoolProp.PropsSI('H','P',self.P_c,'Q',0,nodes.fuel))
        nodes.T_i = nodes.T_prop  # FUEL temperature, NOT flame temperature

        nodes.T_af = nodes.T0
        self.k_g = self.k_th
        self.rho_g = self.rho
        self.mu_g = self.mu
        self.cp_g = self.cp

        # u_rel = DBUProp.u_rel    # Unused at present
        self.u_avg = self.u_rel
        self.Pr = self.cp_g * self.mu_g / self.k_g

        self.d_sbu = self.d_sbu

        ## Calculating Pre-heat Time:
        self.Re_d = self.rho_g * self.u_avg * self.d_sbu / self.mu_g
        self.Nu = 2 + (0.4 * (self.Re_d ** 0.5) * (self.Pr ** (1 / 3)))
        self.A = self.rho_f * self.cp_f * self.d_sbu ** 2
        self.B= 6 * self.k_g * self.Nu

        # Find a way to account for non-adiabatic effects!
        nodes.T_amb = nodes.T_af * 1
        # Non-dimensional temperature
        nodes.T_nd = (nodes.T_amb - nodes.T_i) / (nodes.T_amb - nodes.T_b)

        self.Dsq = self.d_sbu ** 2
        self.dt = 1e-6
        self.C1 = 0.5
        self.D = self.d_sbu
        self.C2 = 0.5

        # math.log used here was the reallog function in MATLAB - ensure that results obtained here do not conflict
        # with the MATLAB ones
        self.t_ph = (self.A / self.B) * math.log(nodes.T_nd)

        # Calculate Vaporisation Time:
        # Method can be found in Fundamentals of Combustion Processes, McAllister, Chapter 8: Droplet Combustion
        self.t_vap = 0  # Initialise counter variable for droplet vaporisation
        while self.Dsq >= 0:
            # Instantaneous droplet Reynolds number:
            self.Re_d = self.rho_g * self.u_avg * self.D / self.mu_g
            # Terms in theoretical formula to calculate droplet vapourisation time
            self.beta0 = 4 * self.k_g * (nodes.T_amb - nodes.T_b) / (self.rho_f * self.h_fg * self.C1)
            self.beta = (1.6 * self.k_g / (self.rho_f * self.h_fg)) * (nodes.T_amb - nodes.T_b) * (self.Re_d ** 0.5) * (self.Pr ** (1 / 3))
            # d(D^2)/dt - rate of change of the droplet's diameter squared with time
            self.dDsq_dt = -2 * (self.C1 * self.beta0) - self.beta
            # Instantaneous square of droplet diameter
            self.Dsq = self.Dsq + (self.dDsq_dt * self.dt)
            # Break the while loop if there Dsq < 0, in order to stop it from throwing a math domain error
            if self.Dsq < 0:
                break
            # Instantaneous droplet diameter
            self.D = math.sqrt(self.Dsq)
            # Counter to monitor vaporisation time
            self.t_vap = self.t_vap + self.dt


        # Calculate Burn Time:
        self.beta0pr = 4 * self.k_g * (nodes.T_af - nodes.T_b) / (self.rho_f * self.h_fg * self.C2)
        # Formula for quiescent combustion:
        self.t_burn = (self.d_sbu ** 2) / self.beta0pr

        # Calculating relevant variables
        self.t_tot = self.t_ph + self.t_vap + self.t_burn
        self.L_comb = self.t_tot * self.u_avg


def RenderToModel():
    # Your code here!
    combustor = cylinder(r=100, h=100)
    combustor -= translate([0,0,-1])(cylinder(r=80, h=102))
    combustor = translate([0, 0, -100])(combustor)

    return combustor


def print_r(the_object):
    print("CLASS: ", the_object.__class__.__name__, " (BASE CLASS: ", the_object.__class__.__bases__, ")")
    pprint(vars(the_object))


# Start of script:

combustor = Combustor(Nodes.nodes)

if __name__ == '__main__':
    combCAD = RenderToModel()
    SEGMENTS = 48
    scad_render_to_file(combCAD, file_header='$fn = %s;' % SEGMENTS, include_orig_code=True)
