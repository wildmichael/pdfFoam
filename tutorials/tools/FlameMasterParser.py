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
