# class generators with base
# symbols cls argument!

from symengine import symbols


class GeneratingSet:
    def __init__(self, name):
        self._symbols_arr = symbols(f"{name}(0:100)")
        for i in range(len(self._symbols_arr)):
            self._symbols_arr[i].base = name
            self._symbols_arr[i].index = i

    def __getitem__(self, i):
        return self._symbols_arr[i]

def is_indexed(x):
    return hasattr(x, "index")

