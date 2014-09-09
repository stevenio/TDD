"""
fun with "with/as" clause
author: YHW
date: Sep-09-2014
"""

class TraceBlock:
  def message(self, arg):
    print("running "+arg)

  def __enter__(self):
    print("starting with block")
    return self

  def __exit__(self, exc_type, exc_value, exc_traceback):
    if exc_type is None:
      print("Existed normally!\n")
    else:
      print('raise an exception! '+str(exc_type))
      return False # Propagate
if __name__ == '__main__':
  a = TraceBlock()
  print(a.message('hello'))
  with TraceBlock() as action:
    action.message('test 1')
    print('reached')

  with TraceBlock() as action:
    action.message('test 2')
    #raise TypeError
    print('not reached')
  
  print("Multiple context manager in 2.7\n")
  with TraceBlock() as a, TraceBlock() as b:
    print(a.message('1'), b.message('2'))
