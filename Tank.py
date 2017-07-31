from Component import Component

class Tank(Component):
  def __init__(self):
    Component.__init__(self)
    
    # Placeholder paremeters for testing.
    self.parameters['diameter'] = 77;
    self.parameters['length'] = 77;
    self.parameters['wall_thickness'] = 77;
    
