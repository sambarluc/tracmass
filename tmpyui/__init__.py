__all__ = ['tmsim', 'tmdata']

from tmpyui.tmdata import tmdata
from tmpyui.tmsim import tmsim

def _get_geometry(grid, geometry):
    """
    Convenience function returning XG and YG of mitgcm
    """
    from numpy import zeros

    if geometry=="curvilinear":
        ni = grid.i_g.size + 1
        nj = grid.j_g.size + 1
        xG = zeros((nj, ni))
        yG = zeros((nj, ni))
        # copy to array, we need numpy fancy indexing
        xG[:-1, :-1] = grid.XG.to_masked_array()
        yG[:-1, :-1] = grid.YG.to_masked_array()
        cs = grid.CS.to_masked_array()
        sn = grid.SN.to_masked_array()
        dxG = grid.dxG.to_masked_array()
        dyG = grid.dyG.to_masked_array()

        # Fill (approximate) end points of the grid
        xG[:-1, -1] = xG[:-1, -2] + dxG[:, -1] * cs[:, -1]
        xG[-1, :-1] = xG[-2, :-1] - dyG[-1, :] * sn[-1, :]
        # we lack the last metric at the NE corner, so we use the
        # nearby metric
        xG[-1, -1] = xG[-1, -2] + dxG[-1, -1] * cs[-1, -1]

        yG[-1, :-1] = yG[-2, :-1] + dyG[-1, :] * cs[-1, :]
        yG[:-1, -1] = yG[:-1, -2] + dxG[:, -1] * sn[:, -1]
        yG[-1, -1] = yG[-2, -1] + dyG[-1, -1] * cs[-1, -1]
    elif geometry=="cartesian":
        ni = grid.i.size + 1
        nj = grid.j.size + 1
        xG = zeros((nj, ni))
        yG = zeros((nj, ni))
        # copy to array, we need numpy fancy indexing
        xG[:-1, :-1] = grid.XG.to_masked_array()
        yG[:-1, :-1] = grid.YG.to_masked_array()
        dxG = grid.dxG.to_masked_array()
        dyG = grid.dyG.to_masked_array()

        # Fill (approximate) end points of the grid
        xG[:-1, -1] = xG[:-1, -2] + dxG[:, -1]
        xG[-1, :] = xG[-2, :]

        yG[-1, :-1] = yG[-2, :-1] + dyG[-1, :]
        yG[:, -1] = yG[:, -2]
    else:
        raise ValueError("Grid geometry not recognised.")
    return xG, yG

def _xy2grid(x, y, Ax, Ay, Bx, By, Cx, Cy, Dx, Dy):
    from numpy import sqrt
    # form Mathematica
    gx = (1/(2*((Ay-By)*(Cx-Dx)-(Ax-Bx)*(Cy-Dy)))) * \
         (Ay*Cx-2*Ay*Dx+By*Dx-Bx*Dy+Ay*x-By*x+Cy*x-Dy*x+Bx*y-Cx*y+Dx*y-Ax*(Cy-2*Dy+y) + \
         sqrt((By*(Dx-x)+Cy*x-Dy*(Bx+x)+Ay*(Cx-2*Dx+x)+(Bx-Cx+Dx)*y-Ax*(Cy-2*Dy+y))**2 + \
         4*((Ay-By)*(Cx-Dx)-(Ax-Bx)*(Cy-Dy))*(Ay*(Dx-x)+Dy*x-Dx*y+Ax*(-Dy + y))))
         
    gy = (1/(2*((-(By-Cy))*(Ax-Dx)+(Bx-Cx)*(Ay-Dy)))) * \
         (2*Ay*Bx-2*Ax*By-Ay*Cx+Ax*Cy+By*Dx-Bx*Dy-Ay*x+By*x-Cy*x+Dy*x+Ax*y-Bx*y+Cx*y-Dx*y + \
         sqrt((By*(Dx-x)+Cy*x-Dy*(Bx+x)+Ay*(Cx-2*Dx+x)+(Bx-Cx+Dx)*y-Ax*(Cy-2*Dy+y))**2 + \
         4*((Ay-By)*(Cx-Dx)-(Ax-Bx)*(Cy-Dy))*(Ay*(Dx-x)+Dy*x-Dx*y+Ax*(-Dy+y))))
    return gx, gy
