
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
            args = [args_[vars.index(str(f))] for f in fargs]
            return f(*args, **kwargs_)
        return wrapped

    @staticmethod
    def error_propagation(vars, pars, expr, error_prefix='d', **kwargs):
        """ Construct the function for calculating propagated error of a function represented by `expr`. `vars` specified the full list of variable that enters the expression. `pars` indicates the elements in `vars` which will not enter error calculation. 
        Parameters
        ----------
        vars : collection of str
            a collection of variable name to be converted to `sympy.symbols`
        pars : collection of str
            a collection of parameter name to be converted to `sympy.symbols`. It must be a subset of `vars`, symbols in this collection does not enter error calculation.
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

        pars = list(set(pars))
        vars = list(set(vars))
        dvars = [f'{error_prefix}{var}' for var in vars if not var in pars]

        dexpr = sympy.sqrt(sum([
            sympy.diff(expr, var) ** 2 * sympy.symbols(f'{error_prefix}{var}') ** 2 for var in vars if not var in pars
        ]))

        return expression.make_function(vars + dvars, dexpr, **kw)
        
        
if __name__ == '__main__':
    
    expr = 'c * x / (1 + exp((x - b0) / db))'
    vars = ['x', 'c', 'b0', 'db']
    f = expression.make_function(vars, expr)

    errs = [0.7, 0.9, 0.1]
    x = np.linspace(0.,10, 21)
    x = 0.5 * (x[1:] + x[:-1])
    y = f(x, 7.803, 9.328, 1.) * 10


    f = expression.error_propagation(vars, pars=['x'], expr=expr)
    yerr = f(x, 7.803, 9.328, 1., *errs) * 10.

    import matplotlib.pyplot as plt
    fig, ax = plt.subplots()
    ax.errorbar(x, y, yerr=yerr, fmt='--')
    plt.show()

    
