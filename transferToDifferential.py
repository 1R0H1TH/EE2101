from sympy.integrals.transforms import inverse_laplace_transform, InverseLaplaceTransform
from sympy.abc import s, t
from sympy import pprint, fraction, Function, Wild, Derivative, Eq, Pow

def transfer_to_differential(tf, fun_X = Function('X'), fun_F = Function('F')):
    tf = fraction(tf)
    res = Eq(inverse_laplace_transform(tf[1] * fun_X(s), s, t), inverse_laplace_transform(tf[0] * fun_F(s), s, t))

    wf = Wild('w')
    ilw = InverseLaplaceTransform(wf, s, t, None)
    
    for exp in res.find(ilw):
        e = exp.match(ilw)[wf]
        args = e.args
        if len(args) == 2:
            p = 1 if not isinstance(args[0], Pow) else args[0].args[1]
            newexp = Derivative(Function(args[1].name.lower())(t), t, p)
            res = res.replace(exp, newexp)
        elif len(args) == 1:
            newexp = Function(e.name.lower())(t)
            res = res.replace(exp, newexp)

    return res