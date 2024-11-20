import sys, yaml

import openmm as mm
from openmm import app
from openmm import unit
from openmm.app.pdbfile import PDBFile
import numpy as np


def read_input(fileinp):
    """ Reads input from file 

    Parameters
    ----------
    fileinp : str
        Input file in yml format

    """
    try:
        with open(fileinp) as stream:
            try:
                print(yaml.safe_load(stream))
            except yaml.YAMLError as e:
                print (e)
    except FileNotFoundError as e:
        print ("\nFile %s not found"%fileinp)
        print (e)


