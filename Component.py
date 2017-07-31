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

  #  Pretty print component, all parameters, and all subcomponents
  def __repr__(self):
     def indent(text, amount, ch=' '):
       padding = amount * ch
       return ''.join(padding+line for line in text.splitlines(True))

     def render(i):
       return indent('\n'.join([k + ': ' + repr(v) for k,v in i.iteritems()]), 2)

     result = self.__class__.__name__

     if len(self.parameters):
       result = result + '\n' + render(self.parameters)

     if len(self.subcomponents):
       result = result + '\n' + render(self.subcomponents)

     return result
