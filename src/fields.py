# fields.py

import sympy as smp


def scalar_field(name, deps):
    return smp.Function(name)(*deps)


def vector_field(name, deps, dim=3):
    return smp.Matrix(
        dim,
        1,
        lambda i, j: smp.Function(f"{name}_{i+1}")(*deps)
    )


def tensor_field(name, deps, rows=3, cols=3):
    return smp.Matrix(
        rows,
        cols,
        lambda i, j: smp.Function(f"{name}_{i+1}{j+1}")(*deps)
    )

def double_contract(A, B):
    return smp.trace(A.T*B)

def material_derivative(phi, v, coords, t):

    result = smp.diff(phi, t)

    for i in range(3):
        result += v[i]*smp.diff(phi, coords[i])

    return smp.simplify(result)


def upper_convected(A, grad_v, t):

    return (
        smp.diff(A, t)
        - grad_v.T*A
        - A*grad_v
    )


def tensor_derivative(f, A):

    rows, cols = A.shape

    return smp.Matrix(
        rows,
        cols,
        lambda i,j: smp.diff(f, A[i,j])
    )

