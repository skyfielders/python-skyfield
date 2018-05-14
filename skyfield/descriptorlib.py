from functools import update_wrapper

class reify(object):
    """Adapted from Pyramid's `reify()` memoizing decorator."""
    def __init__(self, method):
        self.method = method
        update_wrapper(self, method)

    def __get__(self, instance, objtype=None):
        if instance is None:
            return self
        value = self.method(instance)
        instance.__dict__[self.__name__] = value
        return value
