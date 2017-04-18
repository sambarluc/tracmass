from __future__ import division, print_function


def tmdata(mitdir, tmtracks, tstart, ids=None, inds=False, **xmgcm):
    """
    mitdir:    Path to the MITgcm simulations results, used to load the grid information.
    tmdir:     File containing the Lagrangian tracks computed by tracmass.
    tstart:    Beginning time of the MITgcm simulation, as a string with
               format "2010-01-24 12:20"
    ids:       Particle id numbers to load. If None (default), all particles are loaded
    inds:      Return tracks as indices too, default: False.
    **xmgcm:   Keywords arguments passed to xmgcm.mitgcmds to load MITgcm grid files.

    Returns:
        tracks: xarray Dataset containing time, id and position information of each particle.
                If inds=True, return also the tracks as indices, not as physical coordinates.
    """

    import numpy as np
    import xarray as xr
    from xmitgcm import open_mdsdataset as mitgcmds
    from . import _get_geometry

    xmgcm["swap_dims"] = False
    grid = mitgcmds(mitdir, read_grid=True, iters=[], **xmgcm)
    
    xG, yG = _get_geometry(grid, xmgcm["geometry"])

    # We have to take into account the shaved bottom cells
    # tracmass returns a "normalised" vertical coordinate
    # which has to be translated into real space by taking
    # into account hFacC
    dZ = (grid.drF * grid.hFacC).to_masked_array()
    Z = np.zeros((dZ.shape[0] + 1, dZ.shape[1], dZ.shape[2]))
    Z[1:, ...] = np.cumsum(dZ, axis=0).filled(0)
    # tracmass has opposite Z order
    Z = Z[::-1, ...]

    tmbin = np.fromfile(tmtracks, '>f4')
    if (tmbin.size % 5) != 0:
        print("Something is wrong in the size of the tracmass file.")
        tmbin = tmbin[:(tmbin.size//5)*5]
    tmbin = np.reshape(tmbin, (tmbin.size//5, 5))
    all_ids = np.unique(np.int32(tmbin[:, 0]))
    if ids is None:
        ntracks = all_ids.size
        ids = all_ids
    else:
        ids = np.atleast_1d(np.squeeze(ids))
        ntracks = np.size(ids)
    tsteps = np.datetime64(tstart) + \
             np.asarray(tmbin[:, 1]*86400, dtype="timedelta64[s]").astype("timedelta64[ns]")
    tcoord = np.unique(tsteps)

    tracks = grid.drop(grid.data_vars.keys())
    tracks.attrs["MITgcm_dir"] = mitdir
    tracks.attrs["tracmass_file"] = tmtracks
    tracks.coords["id"] = ("id", ids)
    tracks.coords["time"] = ("time", np.unique(tsteps))
    # we also store the full mesh with edges coordinates
    tracks.coords["i_p1"] = ("i_p1", np.arange(xG.shape[1]))
    tracks.coords["j_p1"] = ("j_p1", np.arange(xG.shape[0]))
    tracks.coords["XG_p1"] = (("j_p1", "i_p1"), xG)
    tracks.coords["YG_p1"] = (("j_p1", "i_p1"), yG)

    tracks["xtrack"] = xr.DataArray(np.empty((ntracks, tcoord.size))*np.nan,
                                         coords={"id": ids,
                                                 "time": tcoord},
                                         dims=["id", "time"])
    tracks["ytrack"] = xr.DataArray(np.empty((ntracks, tcoord.size))*np.nan,
                                         coords={"id": ids,
                                                 "time": tcoord},
                                         dims=["id", "time"])
    tracks["ztrack"] = xr.DataArray(np.empty((ntracks, tcoord.size))*np.nan,
                                         coords={"id": ids,
                                                 "time": tcoord},
                                         dims=["id", "time"])
    if inds:
        tracks["itrack"] = xr.DataArray(np.empty((ntracks, tcoord.size))*np.nan,
                                             coords={"id": ids,
                                                     "time": tcoord},
                                             dims=["id", "time"])
        tracks["jtrack"] = xr.DataArray(np.empty((ntracks, tcoord.size))*np.nan,
                                             coords={"id": ids,
                                                     "time": tcoord},
                                             dims=["id", "time"])
        tracks["ktrack"] = xr.DataArray(np.empty((ntracks, tcoord.size))*np.nan,
                                             coords={"id": ids,
                                                     "time": tcoord},
                                             dims=["id", "time"])

    for thisid in ids:
        print("Process id: %d (total: %d)" % (thisid, ids.size))
        thisind = np.where(tmbin[:, 0]==thisid)[0]
        # NOTE: we can use these indices directly, because the grid in tracmass
        # has been defined starting from zero, similarly to python's indexing.
        ii = tmbin[thisind, 2]
        ii_int = np.int32(ii)
        jj = tmbin[thisind, 3]
        jj_int = np.int32(jj)
        kk = tmbin[thisind, 4]
        kk_int = np.int32(kk)

        if inds:
            tracks["itrack"].loc[{"id": [thisid], "time": tsteps[thisind]}] = \
                          np.atleast_2d(ii)
            tracks["jtrack"].loc[{"id": [thisid], "time": tsteps[thisind]}] = \
                          np.atleast_2d(jj)
            tracks["ktrack"].loc[{"id": [thisid], "time": tsteps[thisind]}] = \
                          np.atleast_2d(kk)

        nx = ii - ii_int
        px = 1.0 - nx
        ny = jj - jj_int
        py = 1.0 - ny
        nz = kk - kk_int
        pz = 1.0 - nz

        # x
        st = px * (py * xG[jj_int, ii_int] + ny * xG[jj_int+1, ii_int]) + \
             nx * (py * xG[jj_int, ii_int+1] + ny * xG[jj_int+1, ii_int+1])
        tracks["xtrack"].loc[{"id": [thisid], "time": tsteps[thisind]}] = \
                      np.atleast_2d(st)

        # y
        st = px * (py * yG[jj_int, ii_int] + ny * yG[jj_int+1, ii_int]) + \
             nx * (py * yG[jj_int, ii_int+1] + ny * yG[jj_int+1, ii_int+1])
        tracks["ytrack"].loc[{"id": [thisid], "time": tsteps[thisind]}] = \
                      np.atleast_2d(st)

        # z
        st = pz * Z[kk_int, jj_int, ii_int] + nz * Z[kk_int+1, jj_int, ii_int]
        tracks["ztrack"].loc[{"id": [thisid], "time": tsteps[thisind]}] = \
                      np.atleast_2d(st)

    return tracks.where(np.isfinite(tracks.xtrack))
