from Tank import Tank

class O2(Tank):
  def __init__(self):
    Tank.__init__(self)
    
    # Made up...
    self.parameters['fuel_density'] = 1.03;
    
