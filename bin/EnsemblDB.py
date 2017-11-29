#!env/bin/python

from optparse import OptionParser
from ensembldb import main

# Version
_version = '1.4.1'

# Command line argument parsing
descr = 'CoverView ensembl_db v'+_version
parser = OptionParser(usage='CoverView-1.4.1/ensembl_db <options>', version=_version, description=descr)
parser.add_option('-i', "--input", default=None, dest='input', action='store',help="Input filename (list of ENST IDs)")
parser.add_option('-o', "--output", default=None, dest='output', action='store', help="Output filename prefix")
parser.add_option('-e', "--ens", default=None, dest='ensembl', action='store', help="Ensembl release version")
(options, args) = parser.parse_args()

main.run(options, _version)

