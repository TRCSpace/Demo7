from Tank import Tank

class Fuel(Tank):
  def __init__(self):
    Tank.__init__(self)
    
    # Made up...
    self.parameters['fuel_density'] = 2.5;

