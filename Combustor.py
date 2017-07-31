# ! /usr/bin/env python
# -*- coding: UTF-8 -*-
from __future__ import division
import os
import sys
import re
from solid import *
from solid.utils import *

from Component import Component
# Assumes SolidPython is in site-packages or elsewhwere in sys.path

class Combustor(Component):
  pass

SEGMENTS = 48

def RenderToModel():
    # Your code here!
    combustor = cylinder(r=100, h=100)
    combustor -= translate([0,0,-1])(cylinder(r=80, h=102))
    combustor = translate([0, 0, -100])(combustor)

    return combustor

if __name__ == '__main__':
    combCAD = RenderToModel()
    scad_render_to_file(combCAD, file_header='$fn = %s;' % SEGMENTS, include_orig_code=True)
