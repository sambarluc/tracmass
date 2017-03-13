from __future__ import division, print_function


def tmdata(mitdir, tmtracks, tstart, **xmitgcm):
    """
    mitdir:    Path to the MITgcm simulations results, used to load the grid information.
    tmdir:     File containing the Lagrangian tracks computed by tracmass.
    tstart:    Beginning time of the simulation, as a string with format "2010-01-24 12:20"
    **xmitgcm: Keywords arguments passed to xmitgcm.mitgcmds to load MITgcm grid files.
    """

    import numpy as np
    import xarray as xr
    from xmitgcm import open_mdsdataset as mitgcmds

    grid = mitgcmds(mitdir, read_grid=True, **xmitgcm)

    ni = grid.i_g.size + 1
    nj = grid.j_g.size + 1
    xG = np.zeros((nj, ni))
    yG = np.zeros((nj, ni))
    # copy to array, we need numpy fancy indexing
    xG[:-1, :-1] = grid.XG.to_masked_array()
    yG[:-1, :-1] = grid.YG.to_masked_array()
    # We need the last lines
    dxG = grid.dxG.isel(i=-1)
    dyG = grid.dyG.isel(j=-1)
    # Fill (approximate) end points of the grid
    xG[-1, :] = xG[-2, :]
    xG[:-1, -1] = xG[:-1, -1] + dxG
    xG[-1, -1] = xG[-2, -1]
    yG[-1, :-1] = yG[-1, :-1] + dyG
    yG[:, -1] = yG[:, -2]
    yG[-1, -1] = yG[-1, -2]
    # tracmass has opposite Z order, Zu is the lower interface
    Z = grid.Zu[::-1].to_masked_array()
    dZ = grid.drF[::-1].to_masked_array()

    tmbin = np.fromfile(tmtracks, '>f4')
    if (tmbin.size % 5) != 0:
        raise ValueError("Something is wrong in the size of the tracmass file.")
    tmbin = np.reshape(tmbin, (tmbin.size//5, 5))
    ids = np.unique(np.int32(tmbin[:, 0]))
    ntracks = ids.size
    tsteps = np.datetime64(tstart) + \
             np.asarray(tmbin[:, 1]*86400, dtype="timedelta64[s]")
    tcoord = np.unique(tsteps)

    tracks = xr.Dataset()
    tracks.attrs["MITgcm_dir"] = mitdir
    tracks.attrs["tracmass_file"] = tmtracks
    tracks.coords["p_id"] = ("p_id", ids)
    tracks.coords["time"] = ("time", np.unique(tsteps))
    tracks["xtrack"] = xr.DataArray(np.empty((ntracks, tcoord.size))*np.nan,
                                         coords={"p_id": ids,
                                                 "time": tcoord},
                                         dims=["p_id", "time"])
    tracks["ytrack"] = xr.DataArray(np.empty((ntracks, tcoord.size))*np.nan,
                                         coords={"p_id": ids,
                                                 "time": tcoord},
                                         dims=["p_id", "time"])
    tracks["ztrack"] = xr.DataArray(np.empty((ntracks, tcoord.size))*np.nan,
                                         coords={"p_id": ids,
                                                 "time": tcoord},
                                         dims=["p_id", "time"])

    for thisid in ids:
        thisind = np.where(tmbin[:, 0]==thisid)[0]
        ii = tmbin[thisind, 2]
        ii_int = np.int32(ii)
        jj = tmbin[thisind, 3]
        jj_int = np.int32(jj)
        kk = tmbin[thisind, 4]
        kk_int = np.int32(kk)

        nx = ii - ii_int
        px = 1.0 - nx
        ny = jj - jj_int
        py = 1.0 - ny

        # x
        st = px * (py * xG[jj_int, ii_int] + ny * xG[jj_int+1, ii_int]) + \
             nx * (py * xG[jj_int, ii_int+1] + ny * xG[jj_int+1, ii_int+1])
        tracks["xtrack"].loc[{"p_id": [thisid], "time": tsteps[thisind]}] = \
                      np.atleast_2d(st)

        # y
        st = px * (py * yG[jj_int, ii_int] + ny * yG[jj_int+1, ii_int]) + \
             nx * (py * yG[jj_int, ii_int+1] + ny * yG[jj_int+1, ii_int+1])
        tracks["ytrack"].loc[{"p_id": [thisid], "time": tsteps[thisind]}] = \
                      np.atleast_2d(st)

        # z
        st = Z[kk_int] + (kk - kk_int) * dZ[kk_int]
        tracks["ztrack"].loc[{"p_id": [thisid], "time": tsteps[thisind]}] = \
                      np.atleast_2d(st)

    return tracks
