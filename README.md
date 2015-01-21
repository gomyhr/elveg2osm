elveg2osm
=========

Conversion from Elveg data to openstreetmap

#Usage:
`elveg2osm.py dir XXXX`

The directory dir contains a file called XXXXElveg\_default, which is the
output of the command
`sosi2osm XXXXElveg.SOS default.lua >XXXXElveg_default.osm`.
and at least the file `XXXXFart.txt`.

XXXX is the 4-digit number representing the municipality listed at http://www.statkart.no/Kunnskap/Fakta-om-Norge/Fylker-og-kommuner/Tabell/

#Requirements:
- sosi2osm:      For converting the initial SOSI file to osm-format (without changing the tags)
   - Source code at https://github.com/Gnonthgol/sosi2osm
   - Ubuntu PPA at http://ppa.launchpad.net/saltmakrell/osm/ubuntu/
   - Available in Debian unstable
- osmapis:       Python module used for reading/writing/abstracting OSM-files
   - Source code at https://github.com/xificurk/osmapis
- geographiclib: Python module used for computing distances between nodes
   - Home page at http://geographiclib.sourceforge.net/
   - Source code at https://pypi.python.org/pypi/geographiclib

