from moha import *  


def test_recursive_timer():
    @timer.with_section('Foo')
    def factorial(n):
        if n <= 1:
            return 1
        else:
            return factorial(n-1)*n
    assert factorial(4) == 24
