import nose
from umath import multiply

def setup_func():
  "set up test fixtures"

def teardown_func():
  "tear down test fixtures"

@nose.with_setup(setup_func, teardown_func)
## Alternative way to add setup|teardown
#test.setup = setup_func
#test.teardown = teardown_func
def test_numbers_3_4():
    assert multiply(3,4) == 12 

