from cython cimport wraparound

@wraparound(False)
def transform(float[:] ii,
              float[:] jj,
              float[:] kk,
              double[:, ::1] X,
              double[:, ::1] Y,
              double[:, ::1] dX,
              double[:, ::1] dY,
              double[:, ::1] CS,
              double[:, ::1] SN,
              double[:, :, ::1] Z,
              double[:, ::1] stx,
              double[:, ::1] sty,
              double[:, ::1] stz):

    cdef int m, ii_int, jj_int, kk_int
    cdef double nx, ny, px, py, pxpy, pxny, nxpy, nxny, nz
    cdef double i, j ,k, cosij, sinij, dxij, dyij
    cdef int N = ii.size

    for m in range(N):
        i = ii[m]
        j = jj[m]
        k = kk[m]

        ii_int = int(i)
        jj_int = int(j)
        kk_int = int(k)

        nx = i % 1
        ny = j % 1
        nz = k % 1

        # inside the cell
        if (nx >= 1e-9) & (ny >= 1e-9):
            cosij = CS[jj_int, ii_int]
            sinij = SN[jj_int, ii_int]
            dxij = dX[jj_int, ii_int]
            dyij = dY[jj_int, ii_int]

            # x
            stx[m, 0] = X[jj_int, ii_int] + \
                        cosij * dxij * nx - sinij * dyij * ny

            # y
            sty[m, 0] = Y[jj_int, ii_int] + \
                        sinij * dxij * nx + cosij * dyij * ny

            # z
            if nz >= 1e-9:
                stz[m, 0] = (1.0 - nz) * Z[kk_int, jj_int, ii_int] + \
                            nz * Z[kk_int+1, jj_int, ii_int]
            else:
                stz[m, 0] = Z[kk_int, jj_int, ii_int]
        # at the SW corner
        elif (nx < 1e-9) & (ny < 1e-9):
            # x
            stx[m, 0] = X[jj_int, ii_int]

            # y
            sty[m, 0] = Y[jj_int, ii_int]

            # z
            if nz >= 1e-9:
                stz[m, 0] = (1.0 - nz) * Z[kk_int, jj_int, ii_int] + \
                            nz * Z[kk_int+1, jj_int, ii_int]
            else:
                stz[m, 0] = Z[kk_int, jj_int, ii_int]
        # on W face
        elif nx < 1e-9:
            cosij = CS[jj_int, ii_int]
            sinij = SN[jj_int, ii_int]
            dyij = dY[jj_int, ii_int]

            # x
            stx[m, 0] = X[jj_int, ii_int] - sinij * dyij * ny

            # y
            sty[m, 0] = Y[jj_int, ii_int] + cosij * dyij * ny


            # z
            if nz >= 1e-9:
                stz[m, 0] = (1.0 - nz) * Z[kk_int, jj_int, ii_int] + \
                            nz * Z[kk_int+1, jj_int, ii_int]
            else:
                stz[m, 0] = Z[kk_int, jj_int, ii_int]
        # on S face
        else:
            cosij = CS[jj_int, ii_int]
            sinij = SN[jj_int, ii_int]
            dxij = dX[jj_int, ii_int]

            # x
            stx[m, 0] = X[jj_int, ii_int] + cosij * dxij * nx

            # y
            sty[m, 0] = Y[jj_int, ii_int] + sinij * dxij * nx

            # z
            if nz >= 1e-9:
                stz[m, 0] = (1.0 - nz) * Z[kk_int, jj_int, ii_int] + \
                            nz * Z[kk_int+1, jj_int, ii_int]
            else:
                stz[m, 0] = Z[kk_int, jj_int, ii_int]
