# Basic component.  Implements nothing, subclasses do all the work.

class Component(object):
  def __init__(self):
    # Stores the current value of all design parameters, and estimates of all calculated values.
    self.parameters = new dict();

  def Forward():
    pass
  
  def Backward():
    pass
  
  def RenderToModel():
    pass
  
