from __future__ import division, print_function
import numpy as np
from xmitgcm import open_mdsdataset as mitgcmds
from os.path import join, isfile
from os import access, X_OK
from time import sleep


class tmsim(object):
    """
    A class that deals with tracmass.
    Tracmass must be compiled manually (at least for now),
    then tmsim can configure and run a simulation.
    """
    _empty = np.array([], int)

    def __init__(self, tmexe, projdir, ptemplate, inname):
        """
        Define a tmsim object.
        tmexe:        full path to the tracmass "runtrm" executable
        projdir:      directory of tracmass project template,
                      containing namelist, radgrid, etc.
        ptemplate:    template file with namelist for tracmass,
                      this will be edited to split the computation
                      over multiple cpus.
        inname:       filename of the tracmass namelist
        """
        
        if not tmexe.endswith("runtrm"):
            tmexe = join(tmexe, "runtrm")
        if not isfile(tmexe):
            raise ValueError("Cannot find 'runtrm' executable.")
        if not access(tmexe, X_OK):
            raise ValueError("The given 'runtrm' is not executable.")
        self.exepath = tmexe

        with open(join(projdir, ptemplate), 'r') as f:
            self._infile = f.readlines()
        self.projdir = projdir

        inname = inname.rstrip()
        if inname.endswith(".in"):
            self.nmlname = inname
        else:
            self.nmlname = inname + ".in"
            
        self.ii = None
        self.jj = None
        self.kk = None
        self.indlist = None

    def ijklims(self, ii=_empty, jj=_empty, kk=_empty):
        """
        Define the ii, jj, kk indices from which particles will be
        released.
        ii, jj, kk: list, array, or scalar with indices (i, j, k) from which
                    particles will be released (must be integers).
        """
        self.ii = np.atleast_1d(ii)
        if not np.issubdtype(self.ii.dtype, np.integer):
            raise TypeError("Indices must be integers (i-indices).")
        self.jj = np.atleast_1d(jj)
        if not np.issubdtype(self.jj.dtype, np.integer):
            raise TypeError("Indices must be integers (j-indices).")
        self.kk = np.atleast_1d(kk)
        if not np.issubdtype(self.kk.dtype, np.integer):
            raise TypeError("Indices must be integers (k-indices).")

    def _decompose(self, ncpu):
        if (self.ii is None) or (self.jj is None) or (self.kk is None):
            raise ValueError("To perform a domain decomposition, first "
                             " define the indices from where the particles"
                             " will be released.")
        ii, ji, ki = np.meshgrid(self.ii, self.jj, self.kk)
        self.indlist = zip(ii.ravel(), ji.ravel(), ki.ravel())
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

        self._decompose(ncpu)

        procs = []
        ijkproc = []
        outf = [open(join(self.projdir, ("output_%.6d.txt" % nproc)), 'w')
                for nproc, _ in enumerate(self.indlist)]
        for nproc, ijk in enumerate(self.indlist):
            # generate namelist for the grid point
            nml = self._make_namelist(i1=ijk[0], i2=ijk[0],
                                      j1=ijk[1], j2=ijk[1],
                                      k1=ijk[2], k2=ijk[2],
                                      cn=nproc)
            # write nml to file
            f = open(join(self.projdir, self.nmlname), 'w')
            f.writelines(nml)
            f.close()
            # start tracmass
            procs.append(Popen([self.exepath],
                               stdout=outf[nproc], stderr=outf[nproc]))
            ijkproc.append(ijk)
            # sleep to let the process start and read the namelist
            # this is a poor trick, but as long as we dont want to use
            # a direct interface between python and fortran, it does the job
            if ncpu > 1:
                sleep(0.5)
            # check which processes are running
            active = []
            for n,p in zip(ijkproc, procs):
                if p.poll() is None:
                    active.append(p)
            # wait if there are too many active
            if len(active) >= ncpu:
                active[0].wait()

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
