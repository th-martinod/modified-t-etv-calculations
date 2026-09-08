import sympy as smp
from sympy.printing.latex import LatexPrinter
from IPython.display import Math, display


class ContinuumMechanicsPrinter(LatexPrinter):

    def _print_AppliedUndef(self, expr):
        # v_1(x,y,z,t) -> v_1
        # c_11(x,y,z,t) -> c_{11}
        return self._print(smp.Symbol(expr.func.__name__))

    def _print_Derivative(self, expr):

        vars_ = expr.variables

        # first derivatives
        if len(vars_) == 1:
            var = vars_[0]
            f = expr.expr

            return (
                rf"\partial_{{{self._print(var)}}}"
                rf" {self._print(f)}"
            )

        # higher-order derivatives
        idx = "".join(self._print(v) for v in vars_)

        return (
            rf"\partial_{{{idx}}}"
            rf" {self._print(expr.expr)}"
        )


def latex_cm(expr):
    return ContinuumMechanicsPrinter().doprint(expr)


def show(expr, label=None):

    latex = latex_cm(expr)

    if label is not None:
        latex = rf"{label} = {latex}"

    display(Math(latex))