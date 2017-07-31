from Component import Component
import solid

class Combustor(Component):
  pass


def RenderToModel(geoms):

  c = solid.cylinder(r=10, h=100)
  c -= solid.cylinder(r=2, h=100)