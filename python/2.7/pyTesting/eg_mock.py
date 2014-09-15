"""
fun with mock
author: YHW
date: sep-09-2014
"""

from mock import MagicMock

## usage 1
class ProductionClass(object): pass
thing = ProductionClass()
thing.method = MagicMock(return_value=3)
thing.method(3,4,5, key='value')
thing.method.assert_called_with(3,4,5,key='value')

## useage 2
def mock_search(self):
  class MockSearchQuerySet(SearchQuerySet):
    def __iter__(self):
      return iter(["foo", "bar", "baz"])
  return MockSearchQuerySet()
      

