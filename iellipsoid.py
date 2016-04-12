import math
import numpy as np
from numpy import linalg

from pymol.cgo import BEGIN, COLOR, TRIANGLES, VERTEX, NORMAL, END
from pymol import cmd


def vertex(a1, a2, a3, u, v, M, r0):
    vrtx = M.dot(np.array([
        a1 * math.cos(u) * math.cos(v),
        a2 * math.cos(u) * math.sin(v),
        a3 * math.sin(u)
        ]))
    nrml = M.dot(np.array([
        math.cos(u) * math.cos(v) / a1,
        math.cos(u) * math.sin(v) / a2,
        math.sin(u) / a3
        ]))
    return vrtx + r0, nrml


def iellipsoid(sele, name='iellipsoid', col=[0.5, 0.5, 0.5], u_segs=12, v_segs=12, scale=0.0004):
    data = cmd.get_coords(sele)
    
    r0 = data.mean(axis=0)
    x, y, z = (data - r0).transpose()
    Jxx = sum(y ** 2 + z ** 2)
    Jyy = sum(x ** 2 + z ** 2)
    Jzz = sum(x ** 2 + y ** 2)
    Jxy, Jxz, Jyz = sum(x * y), sum(x * z), sum(y * z)
    
    ws, vs = linalg.eig(np.array([
        [Jxx, -Jxy, -Jxz],
        [-Jxy, Jyy, -Jyz],
        [-Jxz, -Jyz, Jzz]
        ]))
    M = linalg.inv(vs)
    a1, a2, a3 = ws * scale

    mesh = [BEGIN, TRIANGLES, COLOR]
    mesh.extend(col)
    dU = math.pi / u_segs
    dV = 2 * math.pi / v_segs
    U = -math.pi / 2
    for Y in range(0, u_segs):
        V = math.pi
        for X in range(0, v_segs):

            (x1, y1, z1), (n1x, n1y, n1z) = vertex(a1, a2, a3, U, V, M, r0)
            (x2, y2, z2), (n2x, n2y, n2z) = vertex(a1, a2, a3, U + dU, V, M, r0)
            (x3, y3, z3), (n3x, n3y, n3z) = vertex(a1, a2, a3, U + dU, V + dV, M, r0)
            (x4, y4, z4), (n4x, n4y, n4z) = vertex(a1, a2, a3, U, V + dV, M, r0)

            mesh.extend([NORMAL, n1x, n1y, n1z, VERTEX, x1, y1, z1])
            mesh.extend([NORMAL, n2x, n2y, n2z, VERTEX, x2, y2, z2])
            mesh.extend([NORMAL, n4x, n4y, n4z, VERTEX, x4, y4, z4])
            mesh.extend([NORMAL, n2x, n2y, n2z, VERTEX, x2, y2, z2])
            mesh.extend([NORMAL, n3x, n3y, n3z, VERTEX, x3, y3, z3])
            mesh.extend([NORMAL, n4x, n4y, n4z, VERTEX, x4, y4, z4])

            V += dV
        U += dU
    mesh.append(END)

    cmd.load_cgo(mesh, name)


cmd.extend('iellipsoid', iellipsoid)
