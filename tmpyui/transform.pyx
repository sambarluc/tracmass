from cython cimport wraparound

@wraparound(False)
def transform(float[:] ii,
              float[:] jj,
              float[:] kk,
              double[:, ::1] X,
              double[:, ::1] Y,
              double[:, :, :] Z,
              double[:, ::1] stx,
              double[:, ::1] sty,
              double[:, ::1] stz):

    cdef int m, ii_int, jj_int, kk_int
    cdef double nx, ny, px, py, pxpy, pxny, nxpy, nxny, nz
    cdef double i, j ,k
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
            px = 1.0 - nx
            py = 1.0 - ny
            pxpy = px * py
            pxny = px * ny
            nxpy = nx * py
            nxny = nx * ny

            # x
            stx[m, 0] = pxpy * X[jj_int, ii_int] + pxny * X[jj_int+1, ii_int] + \
                        nxpy * X[jj_int, ii_int+1] + nxny * X[jj_int+1, ii_int+1]

            # y
            sty[m, 0] = pxpy * Y[jj_int, ii_int] + pxny * Y[jj_int+1, ii_int] + \
                        nxpy * Y[jj_int, ii_int+1] + nxny * Y[jj_int+1, ii_int+1]

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
            py = 1.0 - ny

            # x
            stx[m, 0] = py * X[jj_int, ii_int] + ny * X[jj_int+1, ii_int]

            # y
            sty[m, 0] = py * Y[jj_int, ii_int] + ny * Y[jj_int+1, ii_int]


            # z
            if nz >= 1e-9:
                stz[m, 0] = (1.0 - nz) * Z[kk_int, jj_int, ii_int] + \
                            nz * Z[kk_int+1, jj_int, ii_int]
            else:
                stz[m, 0] = Z[kk_int, jj_int, ii_int]
        # on S face
        else:
            px = 1.0 - nx

            # x
            stx[m, 0] = px * X[jj_int, ii_int] + nx * X[jj_int, ii_int+1]

            # y
            sty[m, 0] = px * Y[jj_int, ii_int] + nx * Y[jj_int, ii_int+1]

            # z
            if nz >= 1e-9:
                stz[m, 0] = (1.0 - nz) * Z[kk_int, jj_int, ii_int] + \
                            nz * Z[kk_int+1, jj_int, ii_int]
            else:
                stz[m, 0] = Z[kk_int, jj_int, ii_int]
