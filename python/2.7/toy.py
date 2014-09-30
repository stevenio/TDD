__author__ = 'stevenw'
import numpy
class Toy(object):
  def __init__(self):
    self.a = numpy.array([1])
  def get_prop(self):
    return self.a

toy1 = Toy()
tmp = toy1.get_prop()

print("tmp = ", tmp)
toy1.a = '1'
print("tmp = ", tmp)


