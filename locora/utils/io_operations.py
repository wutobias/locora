import os
import numpy as np
from string import ascii_uppercase
from locora.utils.misc import are_you_numpy
from collections import OrderedDict
import gzip

class write_files(object):

    def __init__(self, Delta=None, Frac2Real=None, Bins=None, Origin=None, \
                 Value=None, XYZ=None, X=None, Y=None, Z=None, Format='PDB', \
                 Filename=None, Nan_fill=-1.0):

        """
        This class can write different file types.
        currently only dx and pdb are supported.
        """

        self._delta     = Delta
        self._frac2real = Frac2Real
        self._bins      = Bins
        self._origin    = Origin
        self._value     = Value
        self._x         = X
        self._y         = Y
        self._z         = Z
        self._format    = Format
        self._filename  = Filename
        self._xyz       = XYZ
        self._nan_fill  = Nan_fill

        if type(self._filename) != str:

            self._filename  = 'output.'
            self._filename  += self._format

        self._writers = {
                        'PDB'  : self._write_PDB,
                        'DX'   : self._write_DX,
                        }

        data = self._writers[self._format]()

        o = open(self._filename, "w")
        o.write(data)
        o.close()

    def _merge_x_y_z(self):

        return np.stack( ( self._x, self._y, self._z ), axis=1 )


    def _write_PDB(self):

        """
        Write a PDB file.
        This is intended for debugging. It writes all atoms
        as HETATM of element X with resname MAP.
        """

        if are_you_numpy(self._xyz):

            if self._xyz.shape[-1] != 3:

                raise TypeError(
                    "XYZ array has wrong shape.")

        else:

            if not ( are_you_numpy(self._x) or are_you_numpy(self._y) or are_you_numpy(self._z) ):

                raise TypeError(
                    "If XYZ is not given, x,y and z coordinates\
                     must be given in separate arrays.")

            else:

                self._xyz = self._merge_x_y_z()

        if self._value == None:

            self._value = np.zeros( len(self._xyz), dtype=float )

        data = 'REMARK File written by write_files.py\n'

        for xyz_i, xyz in enumerate(self._xyz):

            #iterate over uppercase letters
            chain_id    = ascii_uppercase[( len(str(xyz_i+1)) / 5 )]

            atom_counts = xyz_i - ( len(str(xyz_i+1)) / 6 ) * 100000
            resi_counts = xyz_i - ( len(str(xyz_i+1)) / 5 ) * 10000
            data += \
            '%-6s%5d %4s%1s%3s %1s%4d%1s   %8.3f%8.3f%8.3f%6.2f%6.2f          \n' \
            %('HETATM',atom_counts+1,'X','', 'MAP', chain_id, resi_counts+1, '', xyz[0], xyz[1], xyz[2], 0.00, float( self._value[xyz_i] ) )

        data += 'END\n'

        return data


    def _write_DX(self):

        """
        Writes DX files according to openDX standard.
        """

        if not ( are_you_numpy(self._origin) or are_you_numpy(self._bins) ):

            raise TypeError(
            "Origin and bins must be given.")

        #This means not (a XOR b) or not (a or b)
        if are_you_numpy(self._delta) == are_you_numpy(self._frac2real) :

            raise TypeError(
            "Either delta or frac2real must be given.")

        if are_you_numpy(self._delta):

            self._frac2real = np.zeros((3,3), dtype=float)

            np.fill_diagonal(self._frac2real, self._delta)

        data = '''object 1 class gridpositions counts %d %d %d
origin %8.4f %8.4f %8.4f
delta %8.4f %8.4f %8.4f
delta %8.4f %8.4f %8.4f
delta %8.4f %8.4f %8.4f
object 2 class gridconnections counts %d %d %d
object 3 class array type float rank 0 items %d data follows
''' %(self._bins[0], self._bins[1], self._bins[2],\
      self._origin[0], self._origin[1], self._origin[2],\
      self._frac2real[0][0], self._frac2real[0][1], self._frac2real[0][2],\
      self._frac2real[1][0], self._frac2real[1][1], self._frac2real[1][2],\
      self._frac2real[2][0], self._frac2real[2][1], self._frac2real[2][2],\
      self._bins[0],   self._bins[1],   self._bins[2],\
      self._bins[2] * self._bins[1] * self._bins[0])

        i = 0
        for x_i in range(0, self._bins[0]):

            for y_i in range(0, self._bins[1]):

                for z_i in range(0, self._bins[2]):

                    ### writing an integer instead of float
                    ### saves us some disk space
                    if np.isnan(self._value[x_i][y_i][z_i]):

                        data += str(self._nan_fill) + " "

                    else:

                        if self._value[x_i][y_i][z_i] == 0.0:

                            data += "0 " 

                        else:

                            data += str(self._value[x_i][y_i][z_i]) + ' '
                
                    i += 1

                    if i == 3:

                        data += '\n'
                        i = 0
        return data


class PDB(object):

  """
  Class that reads a pdb file and provides pdb type data structure.
  """

  def __init__(self, Path):

    self.path = Path

    self.crd  = list()
    self.B    = list()

    with open(self.path, "r") as PDB_file:

      for i, line in enumerate(PDB_file):

        if not (line[0:6].rstrip() == 'ATOM' or line[0:6].rstrip() == 'HETATM'):

          continue

        if i <= 9999:

          #Coordinates
          self.crd.append(list())
          self.crd[-1].append(float(line.rstrip()[30:38]))
          self.crd[-1].append(float(line.rstrip()[38:46]))
          self.crd[-1].append(float(line.rstrip()[46:54]))

          #B-Factors
          self.B.append(line.rstrip()[54:59])

        if 9999 < i <= 99999:

          #Coordinates
          self.crd.append(list())
          self.crd[-1].append(float(line.rstrip()[31:39]))
          self.crd[-1].append(float(line.rstrip()[39:47]))
          self.crd[-1].append(float(line.rstrip()[47:55]))

          #B-Factors
          self.B.append(line.rstrip()[55:60])

        if i > 99999:

          #Coordinates
          self.crd.append(list())
          self.crd[-1].append(float(line.rstrip()[33:41]))
          self.crd[-1].append(float(line.rstrip()[41:49]))
          self.crd[-1].append(float(line.rstrip()[49:57]))

          #B-Factors
          self.B.append(line.rstrip()[57:62])

    self.crd  = np.array(self.crd)
    self.B    = np.array(self.B)