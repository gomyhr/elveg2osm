#! /usr/bin/env python2

'''elveg_all Elveg_archive.zip [XXXX [YYYY [...]]]'''

import sys
import os

filename = sys.argv[1]

# Unzip archive if necessary
if filename[-4:] == '.zip':
    # Assume that it is a zip file
    dirname = filename[:-4]
    if not os.path.isdir(dirname):
        os.mkdir(dirname)
        os.system('unzip -d {0} {1}'.format(dirname, filename))
else:
    dirname = filename

# Decide which kommunes to work on    
if len(sys.argv) > 2:
    kommune_numbers = sys.argv[2:]
else:
    allfiles = os.listdir(dirname)
    kommune_numbers = [fn[0:4] for fn in allfiles if fn [4:] == 'Elveg.zip']
    kommune_numbers.sort()

# Iterate over municipalities
for kn in kommune_numbers:
    sys.stdout.write("Processing municipality: {0}\n".format(kn))
    sys.stdout.flush()
    # Unzip municipality files (if directory does not exist)
    kommune_dir = os.path.join(dirname, kn)
    if not os.path.isdir(kommune_dir):
        #os.mkdir(kommune_dir)
        zipfile = os.path.join(dirname, kn + "Elveg.zip")
        os.system('unzip -o -d {0} {1} >/dev/null'.format(kommune_dir, zipfile))
    # Convert SOSI file to OSM using sosi2osm
    sosifile = os.path.join(kommune_dir, kn + 'Elveg.SOS')
    osmfile = os.path.join(kommune_dir, kn + 'Elveg_default.osm')
    fartfile = os.path.join(kommune_dir, kn + 'Fart.txt')
    hoydefile = os.path.join(kommune_dir, kn + 'Hoyde.txt')
    osmoutput = os.path.join(kommune_dir, kn + 'Elveg.osm')
    logfile = os.path.join(kommune_dir, kn + 'elveg2osm.log')
    os.system('sosi2osm {0} default.lua >{1}'.format(sosifile, osmfile))
    os.system('./elveg2osm.py {0} {1} >{2} 2>&1'.format(kommune_dir, kn, logfile))
