from __future__ import division, print_function
import numpy as np
from xmitgcm import open_mdsdataset as mitgcmds
from os.path import join, isfile
from time import sleep


class tmsim(object):
    """
    A class that deals with tracmass.
    Tracmass must be compiled manually (at least for now),
    then tmsim can configure and run a simulation.
    """
    _empty = np.array([], int)

    def __init__(self, basedir, projdir, ptemplate, inname, mitgcmdir):
        """
        Define a tmsim object.
        basedir:      base directory, where the "runtrm" executable
                      and the "projects" directory are found.
        projdir:      directory of tracmass project template,
                      containing namelist, radgrid, etc.
                      (to be found inside the "projects" subdirectory)
        ptemplate:    template file name with namelist for tracmass,
                      this will be edited to split the computation
                      over multiple cpus.
        inname:       filename of the tracmass namelist, which will be
                      fed to the executable
        mitgcmdir:    Path to MITgcm data
        """
        
        tmexe = join(basedir, "runtrm")
        if not isfile(tmexe):
            raise ValueError("Cannot find 'runtrm' executable.")
        self.exepath = tmexe

        with open(join(basedir, "projects", projdir, ptemplate), 'r') as f:
            self._infile = f.readlines()
        self.project = projdir
        self.projdir = join(basedir, "projects", projdir)

        inname = inname.rstrip()
        if inname.endswith(".in"):
            self.nmlname = inname
        else:
            self.nmlname = inname + ".in"
            
        if not mitgcmdir.endswith("/"):
            mitgcmdir += "/"
        self.mitgcmdir = mitgcmdir

        self.ii = None
        self.jj = None
        self.kk = None
        self.indlist = None

    def ijklims(self, ii=_empty, jj=_empty, kk=_empty,
                inds=True, geometry="curvilinear"):
        """
        Define the ii, jj, kk indices from which particles will be
        released.
        ii, jj, kk: list, array, or scalar with indices (i, j, k) from which
                    particles will be released (must be integers).
        inds:       if True (default), ii, jj and kk are interpreted
                    as indices, otherwise they are interpreted as lists of
                    exact release positions
        geometry:   only used when inds=False, it describes the MITgcm
                    grid geometry. At the moment, only "curvilinear" and
                    "cartesian" have been implemented.
        """
        self.seed_inds = inds

        if inds:
            if not np.issubdtype(self.ii.dtype, np.integer):
                raise TypeError("Indices must be integers (i-indices).")
            self.ii = np.atleast_1d(np.squeeze(ii))
            if not np.issubdtype(self.jj.dtype, np.integer):
                raise TypeError("Indices must be integers (j-indices).")
            self.jj = np.atleast_1d(np.squeeze(jj))
            if not np.issubdtype(self.kk.dtype, np.integer):
                raise TypeError("Indices must be integers (k-indices).")
            self.kk = np.atleast_1d(np.squeeze(kk))
        else:
            # if MITgcm coordinates are passed, we have to load the model grid
            # and then translate into the "normalised" index coordinates of
            # tracmass
            from . import _get_geometry, _xy2grid
            from matplotlib.path import Path
            from itertools import product

            ii = np.atleast_1d(np.squeeze(ii))
            jj = np.atleast_1d(np.squeeze(jj))
            kk = np.atleast_1d(np.squeeze(kk))

            if (ii.size != jj.size) or (ii.size != kk.size):
                raise ValueError("If inds=False, ii, jj and kk must have "
                                 "all the same dimension.")
            
            grid = mitgcmds(self.mitgcmdir, read_grid=True,
                            iters=[], prefix=["UVEL"], swap_dims=False,
                            geometry=geometry)
            xG, yG = _get_geometry(grid, geometry)
            dZ = (grid.drF * grid.hFacC).to_masked_array()
            zG = np.zeros((dZ.shape[0] + 1, dZ.shape[1], dZ.shape[2]))
            zG[1:, ...] = np.cumsum(dZ, axis=0).filled(0)
            # tracmass has opposite Z order
            zG = zG[::-1, ...]
            self.ii = np.zeros(ii.size) * np.nan
            self.jj = np.zeros(ii.size) * np.nan
            self.kk = np.zeros(ii.size) * np.nan
            for nn, (xx, yy, zz) in enumerate(zip(ii, jj, kk)):
                for jj, ii in product(range(xG.shape[0]-1), range(xG.shape[1]-1)):
                    bbPath = Path([[xG[jj, ii], yG[jj, ii]],
                                   [xG[jj, ii+1], yG[jj, ii+1]],
                                   [xG[jj+1, ii+1], yG[jj+1, ii+1]],
                                   [xG[jj+1, ii], yG[jj+1, ii]]])
                    if bbPath.contains_point((xx, yy)):
                        nx, ny = _xy2grid(xx, yy,
                                          xG[jj, ii], #Ax
                                          yG[jj, ii], #Ay
                                          xG[jj, ii+1], #Bx
                                          yG[jj, ii+1], #By
                                          xG[jj+1, ii+1], #Cx
                                          yG[jj+1, ii+1], #Cy
                                          xG[jj+1, ii], #Dx
                                          yG[jj+1, ii]) #Dy
                        z_here = zG[:, jj, ii]
                        if (zz > z_here.max()) or (zz <= z_here.min()):
                            print("Point outside vertical bounds at x,y,z=%.2f,%.2f,%.2f" %
                                  (xx, yy, zz))
                            break
                        kk = np.where(z_here > zz)[0][-1]
                        nz = (zz - z_here[kk]) / (z_here[kk+1] - z_here[kk])
                        self.ii[nn] = ii + nx
                        self.jj[nn] = jj + ny
                        self.kk[nn] = kk + nz
            self.ii = self.ii[np.isfinite(self.ii)]
            self.jj = self.jj[np.isfinite(self.jj)]
            self.kk = self.kk[np.isfinite(self.kk)]


    def _decompose(self, ncpu):
        if (self.ii is None) or (self.jj is None) or (self.kk is None):
            raise ValueError("To perform a domain decomposition, first "
                             " define the indices from where the particles"
                             " will be released.")
        if self.seed_inds:
            ii, ji, ki = np.meshgrid(self.ii, self.jj, self.kk)
            self.indlist = zip(ii.ravel(), ji.ravel(), ki.ravel())
        else:
            self.indlist = zip(self.ii, self.jj, self.kk)
        #nsplit = len(indlist) // ncpu
        #self.mpinds = []
        #for nc in range(ncpu-1):
        #    self.mpinds.append(indlist[nc*nsplit:(nc+1)*nsplit])
        #self.mpinds.append(indlist[(nc+1)*nsplit:])

    def run(self, ncpu=1):
        """
        Run simulation on ncpu cores. Each realese grid cell is run separately.
        """
        if not np.issubdtype(type(ncpu), np.integer):
            raise TypeError("Number of cpus must be an integer.")

        from subprocess32 import Popen

        def check_active(procs):
            # check which processes are running
            active = []
            for p in procs:
                if p.poll() is None:
                    active.append(p)
            return active, len(active)

        self._decompose(ncpu)

        procs = []
        ijkproc = []
        outf = [open(join(self.projdir, ("output_%.6d.txt" % nproc)), 'w')
                for nproc, _ in enumerate(self.indlist)]
        for nproc, ijk in enumerate(self.indlist):
            # generate namelist for the grid point
            if self.seed_inds:
                nml = self._make_namelist(i1=ijk[0], i2=ijk[0],
                                          j1=ijk[1], j2=ijk[1],
                                          k1=ijk[2], k2=ijk[2],
                                          cn=nproc,
                                          ddir=self.mitgcmdir)
            else:
                nml = self._make_namelist(cn=nproc,
                                          ddir=self.mitgcmdir)
                f = open(join(self.projdir, "seed.txt"), 'w')
                f.write((3*"%10.2f" + 2*"%6i" + "%12i\n") %
                        (ijk[0], ijk[1], ijk[2], 5, 0, 1))
                f.close()


            # write nml to file
            f = open(join(self.projdir, self.nmlname), 'w')
            f.writelines(nml)
            f.close()
            # start tracmass
            procs.append(Popen([self.exepath, self.project,
                                self.nmlname[:-3]],
                               stdout=outf[nproc], stderr=outf[nproc]))
            ijkproc.append(ijk)
            # sleep to let the process start and read the namelist
            # this is a poor trick, but as long as we dont want to use
            # a direct interface between python and fortran, it does the job
            if ncpu > 1:
                sleep(0.5)
            # check which processes are running
            active, na = check_active(procs)
            while na >= ncpu:
                # if we wait for a specific process to finish,
                # we may leave one or more CPUs empty even for
                # long time. This is certainly not professional
                # but works.
                sleep(5)
                active, na = check_active(procs)

        # wait for all to finish and store exit code
        self.exit = {}
        for n, (ijk, p) in enumerate(zip(ijkproc, procs)):
            self.exit[ijk] = p.wait()
            outf[n].close()

    def _make_namelist(self, **kwargs):

        from string import Template
        outfile = []
        for l in self._infile:
            t = Template(l)
            outfile.append(t.safe_substitute(kwargs))
        return outfile
