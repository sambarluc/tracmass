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
    cdef int N = ii.size

    for m in range(N):

        ii_int = int(ii[m])
        jj_int = int(jj[m])
        kk_int = int(kk[m])

        nx = ii[m] - ii_int
        ny = jj[m] - jj_int
        nz = kk[m] - kk_int
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
        stz[m, 0] = (1.0 - nz) * Z[kk_int, jj_int, ii_int] + \
                    nz * Z[kk_int+1, jj_int, ii_int]
