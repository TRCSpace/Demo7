from Component import Component

class Tank(Component):
  def __init__(self):
    Component.__init__(self)
    
    # Placeholder parameters for testing.
    # All diameters in mm.
    self.parameters['diameter'] = 77
    # diameter = combustor diameter
    self.parameters['length'] = 77
    # length defined by volume of propellant required
    self.parameters['wall_thickness'] = 77
