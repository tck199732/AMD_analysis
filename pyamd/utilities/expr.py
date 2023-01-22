
import numpy as np
import sympy
import inspect
from functools import wraps

class expression:
    """ A class for mathematical calculation using sympy.
    """
    
    @staticmethod
    def make_function(vars, expr, **kwargs):
        """ Convert a string to python function
        Parameter
        ---------
        vars : collection of string 
            order is unimportant, a non-repeating sequence of string represening variables of the function. Note that the order of args followed the order of occurence in `vars`.
        expr : string
            the string of the mathematic expression which is to be accepted by sympy
        """
        kw = {
            'modules': 'numpy',
        }
        kw.update(kwargs)
        symbols = list(set(vars))
        f = sympy.lambdify(symbols, expr, **kw)
        fargs = [p.name for p in inspect.signature(f).parameters.values()]
    
        @wraps(f)
        def wrapped(*args_, **kwargs_):
            mapping = {v : a for v, a in zip(vars, args_)}
            print(mapping)

            args = [args_[vars.index(str(f))] for f in fargs]
            print(fargs)
            print(args)
            return f(*args, **kwargs_)
        return wrapped

    @staticmethod
    def error_propagation(vars, expr, error_prefix='d', **kwargs):
        """ Construct the function for calculating propagated error of a function represented by `expr` with respect to the variables `vars`
        Parameters
        ----------
        vars : collection of str
            a collection of variable name to be converted to `sympy.symbols`
        expr : str
            to be converted to `sympy` expression
        Return
        ------
        python function 
        """
        kw = {
            'modules': 'numpy',
        }
        kw.update(kwargs)
        vars = list(set(vars))
        dexpr = sympy.sqrt(sum([
            sympy.diff(expr, var) ** 2 * sympy.symbols(f'{error_prefix}{var}') ** 2 for var in vars
        ]))
        return sympy.lambdify(vars, dexpr, **kw)
        
        
if __name__ == '__main__':
    f = expression.make_function(['a', 'x','b'], 'x**2 + b/a', modules= 'numpy')
    