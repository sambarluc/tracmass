from __future__ import division, print_function
import numpy as np
from xmitgcm import open_mdsdataset as mitgcmds


class tmsim(object):
    """
    A class that deals with tracmass.
    Tracmass must be compiled manually (at least for now),
    then tmsim can configure and run a simulation.
    """

    from string import Template
    input_file_template=Template(
"""
&INIT_NAMELIST_VERSION
      gridvernum = $gridvernum,
/
&INIT_GRID_DESCRIPTION
         gcmname = '$gcmname',
       gcmsource = '$gcmsource',
        gridname = '$gridname',
      gridsource = '$gridsource',
        griddesc = '$griddesc',
       indatadir = '$indatadir',
/
&INIT_CASE_DESCRIPTION
        casename = '',
        casedesc = '',
/
&INIT_GRID_SIZE
             imt =   $imt,
             jmt =   $jmt,
              km =   $kmt,
             nst =      2,
         subgrid =      0,
     subgridimin =     -1,
     subgridimax =     -1,
     subgridjmin =     -1,
     subgridjmax =     -1,
     subgridfile =     '',
       subgridid =      1,
/
&INIT_TIME
            ngcm    =    $ngcm,
            dtgcm   =   $dtgcm,
            intmin  =  $intmin,
            intspin = $intspin,
            endJD   =   $endJD,
            iter    =    $iter,
/
&INIT_WRITE_TRAJS
      twritetype =      2,
           kriva = $kriva,
      outdatadir = '',
     outdatafile = '$outdatafile',
 intmininoutfile =  0,
/
&INIT_SEEDING
! nff:            1 = Follow trajectories forward
!                 2 = Follow trajectories backward
!                 3 = Follow trajectories both ways.
             nff =  $nff,
! isec:           1 = Seed particles meridional from the east (y-z)
!                 2 = Seed particles zonal from the north (x-z)
!                 3 = Seed particles horiz from top (x-y)
!                 4 = Seed particles in the middle of T-box
            isec = $isec,
! idir:           1 = follow positive direction (eastward/northward)  
!                -1 = follow negative direction (westward/southward)
!                 0 = both directions
            idir = $idir,
! === Number of trajectories can be set by
! nqua:           1 = constant number of particles in all boxes
!                     (partQuant in # particles / gridcell)
!                 2 = Each particle reflects flux at seeding. 
!                     (partQuant in m3s-1. per particle)
!                 3 = Each particle reflects volume at seeding.
!                     (partQuant in m3 per particle)
!                 4 = Set number of particles in seed file.
!                     (partQuant is not used)
            nqua = $nqua,
! === Particles seeded (particles/cell or m3s-1/particle or m3/particle)
       partquant = $partquant,
        ntracmax =  $ntracmax,
    loneparticle =     -1,
! === initial directions all in MODEL COORDINATES ===
! Method for seeding particles.
! seedtype:       1 = Seed an area defined by ist, jst, and kst.
!                 2 = Use a list to define which cells to seed.
!                 3 = Use a 2-D mask file.
        seedtype =      1,
            ist1 =  $ist1,
            ist2 =  $ist2,
            jst1 =  $jst1,
            jst2 =  $jst2,
! k goes from 1 (bottom) upwards
            kst1 =  $kst1,
            kst2 =  $kst2,
            tst1 =  $tst1,
            tst2 =  $tst2,
         seeddir =     '',
        seedfile =     '',
     varseedfile =      0,
! 1: seed between tst1 and tst2
! 2: seed in file-defined time definition
        seedtime =      0,
! 1: All particles are seeded at each seeding time step
! 2: Each seed position is used once
         seedall =      0,
! 1: seed at i,j,k positions
! 2: seed at x,y,z positions (?)
         seedpos =      1,
       seedparts =      0,
     seedpart_id =      0,
/
&INIT_KILLZONES
            nend =      2,
            ienw =      1,
            iene =      1,
            jens =      1,
            jenn =      1,
           timax = 36500.0,
/
&INIT_TEMP_SALT
           tmin0 =  -50.0,
           tmax0 =  400.0,
           smin0 = -500.0,
           smax0 =  400.0,
           rmin0 = -100.0,
           rmax0 =  500.0,
           tmine =  -50.0,
           tmaxe =  400.0,
           smine = -150.0,
           smaxe =  500.0,
           rmine = -100.0,
           rmaxe =  500.0,
/
&INIT_DIFFUSION
              ah = 2000.0,
              av =    0.0,
/
&INIT_SEDIMENT
        partdiam =  0.001,
            rhos = 2620.0,
           cwamp =   20.0,
           twave =    8.0,
         critvel =    0.1,
/
""")

    def __init__(self, inoutdir):
        self.inoutdir = inoutdir
      self.gridvernum = 0
         self.gcmname = 'none'
       self.gcmsource = 'none'
        self.gridname = 'none'
      self.gridsource = 'none'
        self.griddesc = 'none'
       self.indatadir = 'none'
# grid definition
             self.imt = 0
             self.jmt = 0
              self.km = 0
# Define the time parameters:
# ngcm   : number of hours between GCM output
# dtgcm  : time step in seconds of the GCM simulation
# intmin : initial time step number of GCM, the simulation will start from this
# intspin : seeding time in steps (of output of GCM simulation)
# endJD  : number of days of the particle tracking simulation
# iter   : number of time steps between GCM output
            self.ngcm = 0
           self.dtgcm = 0.0
          self.intmin = 0
         self.intspin = 0
           self.endJD = 0
            self.iter = 0
# kriva:          0 = no writing
#                 1 = write at time intervals of gcm datasets (each ints)
#                 2 = write at each time iteration
#                 3 = write all the time
#                 4 = write only start and end positions
#                 5 = write at chosen intervals
#                 6 = write each spatial grid-crossing 
           self.kriva = 0
     self.outdatafile = ''
# nff:            1 = Follow trajectories forward
#                 2 = Follow trajectories backward
#                 3 = Follow trajectories both ways.
             self.nff = 1
# isec:           1 = Seed particles meridional from the east (y-z)
#                 2 = Seed particles zonal from the north (x-z)
#                 3 = Seed particles horiz from top (x-y)
#                 4 = Seed particles in the middle of T-box
            self.isec = 4
# idir:           1 = follow positive direction (eastward/northward)  
#                -1 = follow negative direction (westward/southward)
#                 0 = both directions
            self.idir = 0
# === Number of trajectories can be set by
# nqua:           1 = constant number of particles in all boxes
#                     (partQuant in # particles / gridcell)
#                 2 = Each particle reflects flux at seeding. 
#                     (partQuant in m3s-1. per particle)
#                 3 = Each particle reflects volume at seeding.
#                     (partQuant in m3 per particle)
#                 4 = Set number of particles in seed file.
#                     (partQuant is not used)
            self.nqua = 0
# === Particles seeded (particles/cell or m3s-1/particle or m3/particle)
       self.partquant = 0
# Maximum number of particles
        self.ntracmax = 1e6
# Seeding area
            self.ist1 = 0
            self.ist2 = 0
            self.jst1 = 0
            self.jst2 = 0
# k goes from 1 (bottom) upwards
            self.kst1 = 0
            self.kst2 = 0
# seed between tst1 and tst2
            self.tst1 = 0
            self.tst2 = 0
