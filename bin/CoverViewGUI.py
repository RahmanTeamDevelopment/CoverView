#!env/bin/python

from optparse import OptionParser
from gui_ import main
from gui_ import parsers

# Version
_version = '1.4.2'

# Command line argument parsing
descr = 'CoverView GUI v'+_version
parser = OptionParser(usage='CoverView-{}/gui <options>'.format(_version), version=_version, description=descr)
parser.add_option('-i', '--input', default=None, dest='input', action='store', help="Input data file names prefix")
parser.add_option('-r', '--reference', default=None, dest='ref', action='store', help="Reference genome file")
(options, args) = parser.parse_args()

meta = parsers.read_metadata(options.input)
if meta['config_opts']['outputs']['regions'] == False or meta['config_opts']['outputs']['profiles'] == False or meta['config_opts']['only_flagged_profiles'] == True:
    print '\nCoverView GUI cannot open the input data specified.'
    print 'The following configuration settings are required for generating CoverView output to be accepted by the GUI:'
    print '\n[outputs]'
    print 'regions_file = true'
    print 'profiles_file = true'
    print 'only_flagged_profiles = false\n'
    quit()

print '\n----------------------------------------------------------------'
print 'Welcome to CoverView GUI {}'.format(_version)
print 'Dataset read: {}_*'.format(options.input)
print 'To open the GUI, go to http://127.0.0.1:5000/ in a web browser'
print 'Press CTRL+C to quit'
print '----------------------------------------------------------------\n'

main.run(options.input, options.ref)