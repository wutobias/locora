
This is *LoCorA*, a package that provides the Local Correlation Analysis method, both as a python-library and a command line tool. The method is a post-processing tool for molecular dynamics trajectories and is used for investigating temporal properties of water molecules on molecular surfaces. Its main application is the calculation of residence times and rotational time constants of water molecules that occupy the first hydration layer of a protein surface. The typical application of *LoCorA* is the analysis of water molecules in a protein binding pocket. The following pamphlet should give a technical introduction to the program and will be complemented by an upcoming paper. However, some more information on the background can be found in the sources provided in the *Background* section.


Requirements
============
* Linux Operating System
* GNU C compiler
* python 2.7 (Anaconda2 recommended)
* python packages: mdtraj, numpy, scipy, matplotlib, multiprocessing
* python packages used for installation only: pip, setuptools


Installation
============

First, download the *LoCorA* package from GitHub
    
```
git clone https://github.com/wutobias/locora
```

Second, go into the locora directory and install the program

```
python setup.py install
```

That's all!


Background
==========

The background on the method and a case study is presented in our recent work [*(Schiebel et al.)*](https://www.nature.com/articles/s41467-018-05769-2) about a experimental and computational investigation on some trypsin structures, featuring several high-resolution neutron structures. The calculations in the paper were carried out with an earlier version of *LoCorA*, which can be obtained upon request. Also, more background on the method and some prelimanary data from my study about thrombin-trypsin selectivity, can be found on my poster [link].

Usage
=====

This is a brief introduction into the functionalities of the program. The main executable of *LoCorA* is called run_locora. Type run_locora --help into your shell in order to get an overview about all options and a brief description of their individual meaning.

Run mode
--------
*LoCorA* has two possible modes of action, either "grid_solvent" or "process_data". The mode must be initialized with the *--mode* option and is mandatory for each run of the run_locora program. The "grid_solvent" mode builds the non-fixed coordinate system and performs all coordinate transformations on the water molecules. Internally, the program uses fractional coordinate representations of the water molecules. The water molecules are treated as if they were on a three dimensional grid (therefore it is called "grid_solvent"). Although the command line version does not support grid output, the coordinates of the water molecules can be easily made descrete through the python API (i.e. write them on a grid) and saved in Data Explorer (DX) file format. Every coordinate axis in the non-fixed coordinate system is defined by a group of atoms and specific operations performed to build the coordinate system axis. The group of atoms must be specified by the user (more on atom defintions in the following paragraphs) and ideally should constitute a rigid substructure, e.g. a phenyl ring or a carboxylate group. 
Note, that the program will write several files on your disk, which can be become rather large (~100-500MB each) depending on the size of the coordinate system and number of water molecules in the coordinate system. These files contain all information about the water molecules and the translation and rotation coordinates of the water molecules.
The other run mode, "process_data", uses the data from the previous "grid_solvent" run (hence, one must always run "grid_solvent" before "process_data") and it is used to calculate translational as well as rotational time constants. This run mode comes with several different optional arguments, which will be introduced in the following paragraphs.

Input file
----------
For every run of *LoCorA*, it is required to provide an input file by using the *--input* option (if not, an error will be raised). The input file contains all definitions of the non-fixed coordinate system as well as the locations of the files containing information about your system. Please note, that it is highly recommended to use the exact same input file for both run modes. However, in some cases it can be useful to make some changes to the input file between running the grid_solvent and process_data commands (e.g. if the grid_solvent mode is run for several individual chunks of the trajectory, but the process_data mode must be performed on the full trajectory.)
In the following, every single parameter for the input file is explained, at the end you will find a working example of an input file.

**trajin**
The path to the trajectory file. It must be in a file format which can be read by mdtraj, the netcdf, dcd or xtc file formats usually work quite well.

**parm**
The path to the parameter file. It must be in a file format which can be read by mdtraj, however currently only amber parameter files (prmtop) are tested. Therefore, it cannot be guaranteed that other file formats will work as well.

**start**
Start the trajectory analysis at this frame in the MD trajecory.

**stop**
Perform the analysis only up to this frame in the MD trajectory. When set to -1, it is performed up to the last frame of the trajectory.

**xx**
This is a mdtraj-style selection mask, which specifies the atoms (from the trajectory and parameter files above) used for construction of the x-axis in the non-fixed coordinate system. The axis is defined as the average orientation of the connecting vectors between all possible combinations of atoms in the mask. Usually, it is a good idea to use a group of atoms which have a fixed orientation with respect to the other atoms defining the z-axis. For instance, for a tyrosine amino acid side chain, one would use two atoms on opposite positions of the phenyl ring (i.e. 1,4 neigboring). Note, that for mdtraj-style selection masks, residue numbering starts with 0.

**zz**
This is a mdtraj-style selection mask, which specifies the atoms (from the trajectory and parameter files above) used for construction of the x-axis in the non-fixed coordinate system. The axis is defined as the average orientation of the normal vectors of the planes constructed from all possible combinations of the atoms in this mask. As already mentioned for the defintion of the other axis, it is wise to use a group of atoms which constitutes a substructure that is rigid with respect to the atoms defining the x-axis. Note, that for mdtraj-style selection masks, residue numbering starts with 0.

**center**
This is a mdtraj-style selection mask, which specifies the atoms used for construction of the center of the fractional coordinate system (i.e. (0.5|0.5|0.5) in fractional coordinates). The center is calculated from the geometrical mean of the atom postions. If *center* is not specified, the atoms used to define the x and z axis are used. Note, that for mdtraj-style selection masks, residue numbering starts with 0.

**xx_ref**
This is a mdtraj-style selection mask, which specifies the atoms for construction of an internal reference orientation vector. In some cases, it can be useful to have an internal point of refernce for the construction of the coordinate system. If the group of atoms that constitutes the non-fixed coordinate system tends to be rather mobile, the non-fixed coordiante system could undergo rotations about one or more of the coordinate axis. If the rotation of one of the axis around +/-180° would lead to an arrangement of the atoms, which cannot be distinguised from the arrangement at 0°, the non-fixed coordinate system would be rotated although the arrangement of the atoms remains (physically) unchanged. This would lead to an apparent rotation of the coordinate system. In order to curcumvent this, it is possible construct a reference vector spanned by the geometric center of the set of reference atoms *xx_ref* and the center of the non-fixed coordinate system. The dot product between this reference vector and the x-axis must satisfy the condition $ x \cdot x_{ref} > 0 $, if not, then the vector defining the x-axis is multiplied by -1.

**zz_ref**
This is a mdtraj-sytle selection mask, and basically has the same meaning as *xx_ref* but for the z-axis. The vector defining the z-axis is multiplied by -1, if the inequality $ z \cdot z_{ref} > 0 $ is not fulfilled.

**water**
This is a mdtraj-style selection mask, which specifies the water molecules which will be included into the analysis.

**dims**
This is a triple of floating point numbers, which specifies the length of the box edges. Only the water molecules which are within this box (constructed using the non-fixed coordinate system) are recorded and saved to disk. Therefore, it is important to have dimensions that are big enough to encompass the full first hydration layer of the atoms of interest. However, the box edges should not be too big, because then the program will consume more computational resources (RAM memory as well as disk space) and needs more time to finish. Usually, (10,10,10) should be sufficient to capture the first hydration layer of most amino acid sidechains.

**unitcell**
In this file, information about the orientation of the non-fixed coordinate system is stored as a (3x3) matrix for each frame.

**pop**
This file contains the number of water molecules that is present in the box for each frame

**frames**
This file contains the frame number (i.e. the time stamp), in which a water molecule was found in the box (consecutively written for each frame).

**center**
The Cartesian coordinates of the center of the non-fixed coordinates system for each frame.

**origin**
The Cartesian coordinates of the origin of the non-fixed coordinates system for each frame

**O_idxs**
The oxygen atom indices for all water molecules, which occupy the box (consecutively written for each frame).

**theta, phi, pis**
The Euler angles of each water molecule occupying the box (consecutively written for each frame).

**xx1_wat, xx2_wat, yy_wat, zz_wat**
The water coordinate system axis (see Background and sources therein for more details) for each water molecule occupying the box (consecutively written for each frame).

**O_frac, H1_frac, H2_frac**
Fractional coordinates of oxygen and hydrogen atoms of water molecules that occupy the box (consecutively written for each frame).

**Example input file**

 
```
### Trajectory and Parameter files ###
### ------------------------------ ###
trajin my_traj.nc
parm   my_parameter.prmtop
start  0
stop   -1

### Coordinate system definitions (selection masks) ###
### ----------------------------------------------- ###
xx          resid 205 and (name CZ  or name CG)
zz          resid 205 and (name CD1 or name CD2 or name CE1 or name CE2)
center_sele resid 205 and (name CD1 or name CD2 or name CE1 or name CE2)
xx_ref      None
zz_ref      resid 192 and (name CG or name CD1 or name NE1 or name CE2 or name CZ2 or name CH2 or name CZ3 or name CE3 or name CD2)
water       water
dims        10 10 10

### File locations for file writing during grid_solvent and reading during process_data ###
### ----------------------------------------------------------------------------------- ###
unitcell   unitcell.dat
pop        pop.dat
frames     frames.dat
center     center.dat
origin     origin.dat
O_idxs     O_idxs.dat
theta      theta.dat
phi        phi.dat
psi        psi.dat
xx1_wat    xx1_wat.dat
xx2_wat    xx2_wat.dat
yy_wat     yy_wat.dat
zz_wat     zz_wat.dat
O_frac     O_frac.dat
H1_frac    H1_frac.dat
H2_frac    H2_frac.dat
```

Command Line Options
====================
In order to see all possible command line options, type
 

```
run_locora --help
```

Also note, that not all options are being used for every run mode.
 

```
usage: run_locora [-h] -i INPUT -m {grid_solvent,process_data} [-noim]
                  [-c CUTOFF] [-ts TIMESTEP] [-ms MINSTEP]
                  [-tr TRANSIENT TRANSIENT] [-np NPROC] [-pre PREFIX]
                  [-w WINDOW] [-b BOOTSTRAP] [-pl] [-ai] [-lp LEGENDRE] [-de]
                  [-v]

Executable for solute-solvent local correlation analysis.

required arguments:
  -i INPUT, --input INPUT
                        Input file.

optional arguments:
  -h, --help            show this help message and exit
  -m {grid_solvent,process_data}, --mode {grid_solvent,process_data}
                        Run mode.
  -noim, --noimage      Swtich off re-imaging during calculations.
  -c CUTOFF, --cutoff CUTOFF
                        Use only water molecules that within cutoff (in Ang.),
                        around the center of the fractional coordinate system.
                        Default is 3.5 Ang. Only used in mode=process_data.
  -ts TIMESTEP, --timestep TIMESTEP
                        Timestep per MD frame in ps. Default is 1 ps. Only
                        used in mode=process_data.
  -ms MINSTEP, --minstep MINSTEP
                        Minium number of occurancies per window, in order to
                        analyse lifetime distribution. Default is 1. Only used
                        in mode=process_data.
  -tr TRANSIENT TRANSIENT, --transient TRANSIENT TRANSIENT
                        Scan transient re-entering time of water molecules in
                        multiples of timestep. Start with -tr[0] and end with
                        -tr[1]. Default is [0, 1]. Only used in
                        mode=process_data.
  -np NPROC, --nproc NPROC
                        Total number of procesess. Default is 1(=no
                        multiprocessing).
  -pre PREFIX, --prefix PREFIX
                        Output prefix. Default is ''. Only used in
                        mode=process_data.
  -w WINDOW, --window WINDOW
                        window size used for calculation of time correlation.
                        Default is 1000. Only used in mode=process_data.
  -b BOOTSTRAP, --bootstrap BOOTSTRAP
                        Number of bootstrapping resample steps. If
                        bootstrap=0, no bootstrapping will be carried out and
                        averages/standard deviation are calculated from
                        subsequent windows of length --window. If bootstrap>0,
                        --boostrap bootstrapping iterations are performed for
                        windows of legnth --window. Default is 0. Only used in
                        mode=process_data.
  -pl, --plot           Turn on plotting. Default is off (=no plotting). Only
                        used in mode=process_data.
  -ai, --anisotropic    Perform anisotropic analysis. Default is off (=no
                        anisotropic analysis). Only used in mode=process_data.
  -lp LEGENDRE, --legendre LEGENDRE
                        Legendre polynomial used for calculation of
                        orientational lifetimes. Can be either 1 or 2. Default
                        is 1(=First order Legendre polynomial). Only used in
                        mode=process_data.
  -de, --double-exp     Perform double exponential fitting. Only used in
                        mode=process_data.
  -v, --verbose         Verbosity output. Default is off.
```

**-i, --input**
The location of the input file as outlined above.

**-m, --mode**
The run mode, can be either "grid_solvent" (the first step) or "process_data" (the second step).

**-noim, --noimage**
Switches off wrapping/imaging of coordinates back into the simulation cell. Usually it is not healthy to switch it off.

**-c, --cutoff**
A water molecule must be less than *--cutoff* away from the center of the non-fixed coordinate system. If this condition is fulfilled, then it populates the solvent layer.

**-ts, --timestep**
This is the physical timestep between MD frames in picoseconds.

**-ms, --minstep**
Minimum number of MD steps that a water molecule must occupy the first hydration layer in order to populate it.

**-tr, --transient***
The maximum number of transient MD steps in between leaving and re-entering the solvent layer. If the water molecule re-enters the solvent layer in less than *--transient* number of steps, it is treated as not having left the solvent layer. The value of this option must be an integer tuple *(value1 value2)*, since the data processing is carried out for all transient steps between *value1* and *value2*.

**-np, --nproc**
The number of processes started during multiprocessing.

**-pre, --prefix**
Prefix used for names of the files written to disk.

**-w, --window**
All time constants are calculated as averages over windows taken from the full trajectory. The length of these windows is given by this option. The averaging is carried out over sequentiel, non-overlapping windows, if the *--bootstrap* option is set to 0 (see below).

**-b, --bootstrap**
If the value of this option is >0, the windows will be selected using bootstrapping. This effectively deactivates the use of sequentiel, non-overlapping windows. The length of the windows is specified with the *--windows* option. This approach is more robust then the use of sequentiel windows and is recommended.

**-pl, --plot**
Generate some plots showing distributions of water density and lifetimes.

**-ai, --anisotropic**
Perform anisotropic analysis. If this option is activated, the non-fixed coordinate system is subdevided into its octants and the analysis is performed individually for each one of them.

**-lp, --legendre**
This option specifies the order of the Legendre polynom used for calculation of orientational lifetimes.

**-de, --double-exp**
Use a weighted double exponential fit instead of only a single exponential.

**-v, --verbose**
Verbose output or not. Note, that in *--mode=grid_solvent*, PDB files with the coordinates of the water molecules in the box (as defined by the non-fixed coordinate system) are written to disk, if the verbose mode is active. This might generate a lot of data...

Contact
=======

Main Author of the Program
--------------------------
* Tobias Hüfner-Wulsdorf, Philipps-Universität Marburg, tobias.wulsdorf@gmail.com 

Scientific Supervision
----------------------
* Gerhard Klebe, Philipps-Universität Marburg, klebe@staff.uni-marburg.de
