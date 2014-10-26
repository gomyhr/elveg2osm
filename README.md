elveg2osm
=========

Conversion from Elveg data to openstreetmap

Usage:
elveg2osm.py XXXXElveg_default.osm XXXXFart.txt XXXXHoyde.txt output.osm

XXXXElveg_default.osm is output from "sosi2osm XXXXElveg.SOS default.lua >XXXXElveg_default.osm"

Requirements:
sosi2osm:      For converting the initial SOSI file to osm-format (without changing the tags)
   - Source code at https://github.com/Gnonthgol/sosi2osm
   - Ubuntu PPA at http://ppa.launchpad.net/saltmakrell/osm/ubuntu/
   - Available in Debian unstable
osmapis:       Python module used for reading/writing/abstracting OSM-files
   - Source code at https://github.com/xificurk/osmapis
geographiclib: Python module used for computing distances between nodes
   - Home page at http://geographiclib.sourceforge.net/
   - Souce code at https://pypi.python.org/pypi/geographiclib
