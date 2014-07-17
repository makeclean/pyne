#!/usr/bin/env python
""" Module for the production of PartiSn input decks. PartiSn is a discrete
ordinates code produced by Los Almos National Laboratory (LANL). Can be used
to produce neutron, photon, or coupled neutron photon prblems, adjoint or
forward or time dependent problems can be run.

The module is designed to operate on either 2D or 3D meshes, and produce the 
appropriate input. It would be lovely if we eventually manage to get it working
with 1D as well as this appears to be a common mode of operation for PartiSn.

The Input class is the on being worked on currently and should need the least work
to improve. Fundamental inputs to the PartiSn class are: 
    cell_fracs, a list of cell fractions with the number of materials 
                as produced by DG
    mesh, a PyNE mesh instance including materials
    bxslib, the filename of the PartiSn cross section file 

Next should be the Output class to read the output file and rtflux file

If PyTaps not installed then this module will not work.
"""

from __future__ import print_function, division
import sys
import collections
import string
import struct
import math
import os
import linecache
import datetime
from warnings import warn
from pyne.utils import VnVWarning
import itertools

import numpy as np
import tables

from pyne.material import Material
from pyne.material import MultiMaterial
from pyne import nucname
from pyne.binaryreader import _BinaryReader, _FortranRecord

warn(__name__ + " is not yet V&V compliant.", VnVWarning)

# Mesh specific imports
try:
    from itaps import iMesh
    HAVE_PYTAPS = True
except ImportError:
    warn("the PyTAPS optional dependency could not be imported. "
                  "All aspects of the PartiSn module are not imported.",
                  VnVWarning)
    HAVE_PYTAPS = False

if HAVE_PYTAPS:
    from pyne.mesh import Mesh, StatMesh, MeshError, IMeshTag

def _string_width(string,char_len):
    """This functions takes as argument an arbitrarly long string
    and will delimit it by newlines every n characters
    
    """
    
    # 

    # x is now array of strings
    x = [string[i:i+char_len] for i in range(0, len(string), char_len)]
    # join the strings by new lines rather than array of strings
    return_string = "".join(part+"\n" for part in x)
    return return_string


class PartisnInput():

    def _block1(mesh):
        """This function reads a structured  mesh object and returns the 1st data
        block of input, the 1st data block contains the definition of the 
        problem initialisation, specifically, the number of energy groups, the
        Sn order 
        
        Parameters
        ----------
        mesh : PyNE Mesh object
        
        Returns
        -------
        block1_str : str
        A string containing the full block 1 definition for a PartiSn input
        """
        
        block1_str = "/A# block 1 \n"
        block1_str += "t \n"
        return block1_str


    def _block2(self, mesh, bounds):
        """This function reads a structured  mesh object and returns the 2nd data
        block of input, the 2nd data block contains the full definition for the 
        problem geometry and material assignments

        This may seem bizarre, it is however correct that the PartiSn material
        assignments are done here, by assigning an integer correspondence 
        to the number of unique material combinations (which themselves are
        not defined until block4)
        
        Parameters
        ----------
        mesh : PyNE Mesh object
        
        Returns
        -------
        block2_str : str
        A string containing the full block 2 definition for a PartiSn input
        """
        
        block2_str = "/A# block 2 \n"
        block2_str += PartisnInput._partisn_geom(self, mesh, bounds)
        block2_str += PartisnInput._partisn_material(self,mesh)
        block2_str += "t \n"

        return block2_str

    def _block3(self, bxslib):
        """This function reads a structured  mesh object and returns the 3rd data
        block of input, the 3rd data block specifies the cross section library
        and edits that you would like. This will be a wall of text that is valid
        PartiSn input
        
        We define the nuclides that make up the "Permenant Materials" here along
        with the edits that were found in the cross section file.

        Parameters
        ----------
        mesh : PyNE Mesh object
        
        Returns
        -------
        block3_str : str
        A string containing the full block 3 definition for a PartiSn input
        """
        
        block3_str = "/A# block 3 \n"
        print (PartisnInput._read_bxslib(self, bxslib))
        block3_str += "lib="+bxslib
        # the next line is the number of neutron energy groups in the file bxslib
        # need to find a way to query the file to get the number of photon and 
        # neutron energy groups
        block3_str += "lng=175"
#        block3_str += PartisnInput._read_bxslib(self, bxslib)
        block3_str += "t \n"
        return block3_str

    def _block4(self):
        """This function reads a structured  mesh object and returns the 4th data
        block of input, the 4th data block specifies the material definitions
        from the problem. The "pure" materials which make up the problem are defined
        and then the mixtures that make up the zones are required, i.e. each unique
        mixture needs to be represented
        
        First we create the matls array, which contains the materials that will later be 
        mixed by volume fraction from Discretise Geom step. These materials represent
        pure (The PartiSn manual refers to them as "Permenant Materials")

        Parameters
        ----------
        mesh : PyNE Mesh object
        
        Returns
        -------
        block4_str : str
        A string containing the full block 4 definition for a PartiSn input
        """
        
        block4_str = "/A# block 4 \n"
        block4_str += "t \n"
        return block4_str


    def _block5(self):
        """This function reads a structured  mesh object and returns the 5th data
        block of input, the 5th data block specifies the source behaviour and 
        normalisation
        
        Parameters
        ----------
        mesh : PyNE Mesh object
        
        Returns
        -------
        block5_str : str
        A string containing the full block 5 definition for a PartiSn input
        """
        
        block5_str = "/A# block 5 \n"
        block5_str += "t \n"
        return block5_str
    
    def _partisn_geom(self, mesh, bounds):
        """This function reads a structured  mesh object and returns the mesh 
        portion of a PartiSn input deck
        
        Parameters
        ----------
        mesh : PyNE Mesh object
        
        Returns
        -------
        geom : str
        A string containing the PartiSn mesh boundaries and geometry that 
        can be written directly as valid syntax 
        """

        x_coords = mesh.structured_get_divisions("x")
        y_coords = mesh.structured_get_divisions("y")
        z_coords = mesh.structured_get_divisions("z")        

        # coarse bin boundaries
        geom  = "/ coarse bins \n"
        xmesh =  "xmesh =" + " ".join(format(x, "f") for x in x_coords) 
        geom  += _string_width(xmesh,50)
        ymesh = "ymesh =" + " ".join(format(x, "f") for x in y_coords) 
        geom  += _string_width(ymesh,50)
        zmesh = "zmesh =" + " ".join(format(x, "f") for x in z_coords) 
        geom  += _string_width(zmesh,50)


        # number of fine mesh intervals between each coarse bin
        geom += "/ fine bins \n"
        xints = "xints =" + " ".join(format(x, "d") for x in bounds[0])
        geom  += _string_width(xints,50)
        yints = "yints =" + " ".join(format(x, "d") for x in bounds[1])
        geom  += _string_width(yints,50)
        zints = "zints =" + " ".join(format(x, "d") for x in bounds[2])
        geom  += _string_width(zints,50)

        # now write out the material assignments, integer material reference
        # number from 0 (void) to N

        # for now write out all zones 0
        geom += "/ material assignments \n"
        zones = "zones= " + " ".join(format(0,"d") for x in mesh)
        geom += _string_width(zones,50)
        
        geom += "* temporary placeholder for geom \n"

        return geom

    def _read_bxslib(self,filename):
        """ This function reads a supplied binary Partisn cross section file
        (bxslib) and provides a list of the possible materials and prepared edits
        from the cross section file. This is particularly hackish, and someone should
        look into this with much more time than me.

        Parameters
        ----------
        filename : the filename of the bxsfile 
        
        Returns
        -------
        edits : a list of the possible edit and material names
        """
        
        bxslib = open(filename,'rb')

        string = ""
        edits = ""

        xs_names=[]

        # 181st byte is the start of xsnames
        bxslib.seek(180)
        done = False
        while not done:
            for i in range(0,8):
                bytes = bxslib.read(1)
                pad1=struct.unpack('s',bytes)[0]
                if '\x00' in pad1:
                    done = True
                    return xs_names
                string += pad1
            xs_names.append(string)
            string=""

    def _partisn_material(self,mesh):
        """This function reads a structured  mesh object and returns the material
        assignnent portion of a PartiSn input deck
        
        Parameters
        ----------
        mesh : PyNE Mesh object
        
        Returns
        -------
        material : str
        A string containing the PartiSn material assignments and is valid
        PartiSn syntax
        """
                
        material = "* temporary placeholder for material \n"
        
        return material

    def write_partisn(self, mesh, bounds, filename, bxslib):
        """This function reads a structured mesh object and returns the 
        the complete PartiSn input deck and writes it out to file

        Parameters
        ----------
        mesh : PyNE Mesh object
        
        Returns
        -------
        filename : str
        A file to write the output to, note this will always overwrite the
        file if it exists
        """
        
        output_data = ("   1    0    0 \n" 
                       "* PartiSn input deck produced automatically by PyNE \n"
                       "* pyne.partisn module \n")
        output_data += PartisnInput._block1(self)
        output_data += PartisnInput._block2(self, mesh, bounds)
        output_data += PartisnInput._block3(self, bxslib)
        output_data += PartisnInput._block4(self)
        output_data += PartisnInput._block5(self)
        output_data += ("/ ********************* \n"
                        "* You must produce your own edits\n")

        f = open(filename,'w')
        f.write(output_data) # python will convert \n to os.linesep
        f.close()
