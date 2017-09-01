#!env/bin/python

from optparse import OptionParser
from gui import main

# Version
_version = '1.3.0'

# Command line argument parsing
descr = 'CoverView GUI v'+_version
parser = OptionParser(usage='CoverView/env/bin/gui <options>', version=_version, description=descr)
parser.add_option('-i', default=None, dest='input', action='store', help="Input data file names prefix")
parser.add_option('-r', default=None, dest='ref', action='store', help="Reference genome file")
(options, args) = parser.parse_args()

print '\n--------------------------------------------------'
print 'Welcome to CoverView GUI.'
print 'Dataset opened: {}_*'.format(options.input)
print 'Go to http://127.0.0.1:5000/ in your web browser.'
print 'Press CTRL+C to quit.'
print '--------------------------------------------------\n'

main.run(options.input, options.ref)