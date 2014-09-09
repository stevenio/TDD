import nose


def setup_func():
  "set up test fixtures"

def teardown_func():
  "tear down test fixtures"

@nose.with_setup(setup_func, teardown_func)
def test():
  "test ..."

## Alternative way to add setup|teardown
#test.setup = setup_func
#test.teardown = teardown_func
