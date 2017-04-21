from __future__ import division, print_function


def tmdata(mitdir, tmtracks, tstart, ids=None, **xmgcm):
    """
    mitdir:    Path to the MITgcm simulations results, used to load the grid information.
    tmdir:     File containing the Lagrangian tracks computed by tracmass.
    tstart:    Beginning time of the MITgcm simulation, as a string with
               format "2010-01-24 12:20"
    ids:       Particle id numbers to load. If None (default), all particles are loaded
    **xmgcm:   Keywords arguments passed to xmgcm.mitgcmds to load MITgcm grid files.

    Returns:
        tracks: xarray Dataset containing time, id and position information of each particle.
                If inds=True, return also the tracks as indices, not as physical coordinates.
    """

    import numpy as np
    import xarray as xr
    from xmitgcm import open_mdsdataset as mitgcmds
    from . import _get_geometry
    from transform import transform

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

    tracks = xr.open_dataset(tmtracks, mask_and_scale=False)
    # we load it into memory now, because while we run the file may change
    # if the simulation is still running, and this causes NETcdf issues
    # In general, very little memory is required to load the tracks.
    tracks.load()
    all_ids = tracks.id[(tracks.itrack != 0).any(dim="time")]
    if ids is None:
        ntracks = all_ids.size
        ids = all_ids
    else:
        ids = np.atleast_1d(np.squeeze(ids))
        ntracks = np.size(ids)
    tcoord = np.datetime64(tstart) + tracks.time

    for k,v in grid.drop(grid.data_vars.keys()).items():
        tracks.coords[k] = v
    tracks.attrs["MITgcm_dir"] = mitdir
    tracks.attrs["tracmass_file"] = tmtracks
    tracks.coords["time"] = ("time", tcoord)
    # we also store the full mesh with edges coordinates
    tracks.coords["i_p1"] = ("i_p1", np.arange(xG.shape[1]))
    tracks.coords["j_p1"] = ("j_p1", np.arange(xG.shape[0]))
    tracks.coords["XG_p1"] = (("j_p1", "i_p1"), xG)
    tracks.coords["YG_p1"] = (("j_p1", "i_p1"), yG)

    tracks["xtrack"] = xr.DataArray(np.empty((tcoord.size, ntracks))*np.nan,
                                         coords={"time": tcoord,
                                                 "id": ids},
                                         dims=["time", "id"])
    tracks["ytrack"] = xr.DataArray(np.empty((tcoord.size, ntracks))*np.nan,
                                         coords={"time": tcoord,
                                                 "id": ids},
                                         dims=["time", "id"])
    tracks["ztrack"] = xr.DataArray(np.empty((tcoord.size, ntracks))*np.nan,
                                         coords={"time": tcoord,
                                                 "id": ids},
                                         dims=["time", "id"])

    print("Converting particle trajectories...")
    stx = np.empty((tracks.itrack.shape[0], 1))
    sty = np.empty_like(stx)
    stz = np.empty_like(stx)
    for thisid in ids:
        trid = tracks.sel(id=thisid, drop=True)
        # NOTE: we can use these indices directly, because the grid in tracmass
        # has been defined starting from zero, similarly to python's indexing.
        ii = trid.itrack.to_masked_array().filled(0)
        jj = trid.jtrack.to_masked_array().filled(0)
        kk = trid.ktrack.to_masked_array().filled(0)

        transform(ii, jj, kk, xG, yG, Z, stx, sty, stz)
        
        tracks["xtrack"].loc[{"id": [thisid], "time": tcoord}] = \
                      stx
        tracks["ytrack"].loc[{"id": [thisid], "time": tcoord}] = \
                      sty
        tracks["ztrack"].loc[{"id": [thisid], "time": tcoord}] = \
                      stz

    print("Done.")

    return tracks.where(tracks.itrack != 0)
