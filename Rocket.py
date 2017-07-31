from Component import *
from Combustor import *
from Fuel import *
from Injector import *
from Nozzle import *
from O2 import *


# Represents an entire rocket
class Rocket(Component):
  def __init__(self):
    Component.__init__(self)

    self.subcomponents['fuel'] = Fuel();
    self.subcomponents['o2'] = O2();
    self.subcomponents['nozzle'] = Nozzle();
    self.subcomponents['injector'] = Injector();
    self.subcomponents['combustor'] = Combustor();

