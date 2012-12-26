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
#    FlameMasterParser.py
#
# Description
#    Functions to parse the FlameMaster output.
#
#------------------------------------------------------------------------------
import re
_gridPointsRe = re.compile(r'^gridPoints = (?P<nData>\d+)')
_bodyRe = re.compile(r'^body$')
_trailerRe = re.compile(r'^trailer$')
_blockStartRe = re.compile(r'^(?P<name>\S+)')

def parse(file):
  """Parse a FlameMaster file. Returns the header, footer, a lists of names and
  variable definitions and a dictionary with variable values."""
  import sys
  import re
  # initialize
  if type(file) == str:
    file = open(file, 'rt')
  header = []
  footer = []
  values = {}
  names = []
  definitions = []
  nData = -1
  inBody = False
  inHead = True
  i = 0
  for l in file:
    i += 1
    # match 'gridPoints = <nData>'
    m = _gridPointsRe.match(l)
    if m:
      nData = int(m.group('nData'))
      header.append(l)
      continue
    # match 'body'
    m = _bodyRe.match(l)
    if m:
      inBody = True
      inHead = False
      header.append(l)
      continue
    # match 'trailer'
    m = _trailerRe.match(l)
    if m:
      inBody = False
      footer.append(l)
      continue
    # match the start of a data block in the body
    m = _blockStartRe.match(l)
    if m and inBody:
      if nData < 0:
        raise AssertionError('Hit body section without encountering gridPoints')
      name = m.group('name').replace('-','_')
      names.append(name)
      definitions.append(l.strip())
      values[name] = []
      continue
    # ok, still here and in the body section -> read values
    if inBody:
      values[names[-1]].extend(
          map(float,l.split()))
      continue
    # if we are still here, print the line
    if inHead:
      header.append(l)
    elif not inBody:
      footer.append(l)
    else:
      raise AssertionError('Corrupt data file (line %d)',i)
  return ''.join(header), ''.join(footer), names, definitions, values
