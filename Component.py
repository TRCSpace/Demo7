# Basic component.  Implements nothing, subclasses do all the work.

class Component:
  def __init__(self):
    # Stores the current value of all design parameters, and estimates of all calculated values.
    self.parameters = {};
    self.subcomponents = {};

  def Forward(self):
    pass
  
  def Backward(self):
    pass
  
  def RenderToModel(self):
    pass
  
  
  # Extra "nice to have" features
  def __str__():
     return str(self.parameters)
