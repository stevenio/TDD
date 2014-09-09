"""
fun with python try/except, try/finally, with/as statements
author: YHW
date: Sep-09-2014
"""
#########################################
# Compatibility with python 3
#########################################
from __future__ import print_function, division
#########################################

def fetcher(obj, index):
  return obj[index]

x = 'spam'
try:
  _tmp = fetcher(x,4)
  print('result:', _tmp)
except IndexError:
  print('Error Hint: index out of bound!!!')


class NewException(Exception):
  def __str__(self):
    return  "str ..."
def grail():
  raise NewException()

try:
  grail()
except NewException:
  print('new exception')
finally:
  print('Done!')
  

print(1/2)
