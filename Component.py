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
    # If not overridden, simply render all the subcomponents as a rendering of this component.
    sub_component_models = [x.RenderToModel() for x in self.subcomponents.values()]
    
    #Filter places we don't have models
    sub_component_models = [x for x in sub_component_models if x is not None]
    
    if len(sub_component_models):
      return reduce((lambda x, y: x + y), sub_component_models)
    else:
      return None
  
  
  # Extra "nice to have" features

  #  Pretty print component, all parameters, and all subcomponents
  def __repr__(self):
     def indent(text, amount, ch=' '):
       padding = amount * ch
       return ''.join(padding + line for line in text.splitlines(True))

     def render(i):
       return indent('\n'.join([k + ': ' + repr(v) for k,v in i.iteritems()]), 2)

     result = self.__class__.__name__

     if len(self.parameters):
       result = result + '\n' + render(self.parameters)

     if len(self.subcomponents):
       result = result + '\n' + render(self.subcomponents)

     return result
