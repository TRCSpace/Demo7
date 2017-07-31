from Component import Component

class Tank(Component):
  def __init__(self):
    Component.__init__(self)
    
    # Placeholder paremeters for testing.
    self.parameters['diameter'] = 77;
    self.parameters['length'] = 77;
    self.parameters['wall_thickness'] = 77;
    
  def Forward(self):
    # This math is all made up...   Take with a large grain of salt...
    
    output = {}
    output['volume'] = self.parameters['diameter'] ** 2 / 2 * self.parameters['length']
    output['mass'] = output['volume'] * self.parameters['fuel_density']
    
    return output;
