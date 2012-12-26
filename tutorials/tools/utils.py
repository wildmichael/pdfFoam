#------------------------------------------------------------------------------
#                pdfFoam: General Purpose PDF Solution Algorithm
#                   for Reactive Flow Simulations in OpenFOAM
#
# Copyright (C) 2012 Michael Wild, Heng Xiao, Patrick Jenny,
#                    Institute of Fluid Dynamics, ETH Zurich
#------------------------------------------------------------------------------
# License
#    This file is part of pdfFoam.
#
#    pdfFoam is free software: you can redistribute it and/or modify it
#    under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 2 of the License, or
#    (at your option) version 3 of the same License.
#
#    pdfFoam is distributed in the hope that it will be useful, but WITHOUT
#    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
#    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
#    for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with pdfFoam.  If not, see <http://www.gnu.org/licenses/>.
#
# Script
#    utils.py
#
# Description
#    Common utility functions.
#
#------------------------------------------------------------------------------
import os.path as _op
import matplotlib as _mp

def parseMargins(strmargins, defaults=None):
  "Parse margins specification"
  marginNames = 'left right bottom top hspace wspace'.split()
  if not defaults:
    defaults = {}
    for n in marginNames:
      defaults[n] = _mp.rcParams['figure.subplot.'+n]
  margins = dict(defaults)
  for i, m in enumerate(strmargins.split(',')):
    if len(m):
      margins[marginNames[i]] = float(m)
  return margins

def getConfig(case=_op.curdir):
  """Sets defaults and reads case.ini

  Reads the following files, in the order specified:
   - `defaults.ini` located next to this script
   - Searches upwards in the directory tree, reading all defaults.ini files
     that are found. The files are read in reverse order, i.e. files more
     "distant" from the case directory are read before the ones "close" to it.
   - If present, the file `<case>/system/case.ini` is read.
  """
  import ConfigParser, io
  defaults = _op.join(_op.dirname(__file__), 'defaults.ini')
  assert _op.isfile(defaults), "Can't find defaults.ini"
  config = ConfigParser.SafeConfigParser()
  # assemble list of config files to read
  files = []
  d = _op.normpath(_op.abspath(case))
  f = _op.join(case, 'system', 'case.ini')
  if _op.isfile(f):
    files.append(f)
  while d != _op.sep:
    f = _op.join(d, 'defaults.ini')
    if _op.isfile(f):
      files.append(f)
    d = _op.dirname(d)
  files.append(defaults)
  # read in reverse order
  for f in reversed(files):
    config.read(f)
  return config
