#!/usr/bin/env python2.7

#(c) 2015 by Authors
#This file is a part of SpectroGene program.
#Released under the BSD license (see LICENSE file)

"""
This script does all the necessary preparations
and invokes SpectroGene
"""

import os
import sys

LIB_DIR = "lib"

#Check Python version
if sys.version_info[:2] != (2, 7):
    print("Error: SpectroGene requires Python version 2.7 ({0}.{1} detected)."
          .format(sys.version_info[0], sys.version_info[1]))
    sys.exit(-1)

#Setting executable paths
spectrogene_root = os.path.dirname(os.path.realpath(__file__))
lib_absolute = os.path.join(spectrogene_root, LIB_DIR)
os.environ["SPECTROGENE_LIB"] = lib_absolute
sys.path.insert(0, lib_absolute)
sys.path.insert(0, spectrogene_root)

#Spectrogene entry point
from spectrogene.main import main
sys.exit(main())
