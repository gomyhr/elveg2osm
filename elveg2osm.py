#! /usr/bin/env python2
import sys
import os
import osmapis
import csv
import numpy as np
import geographiclib.geodesic as gg

# Output have the following temporary features:
# - Some ways have tags DEBUG=*. Those have Elveg tags I'm unsure about how to tag in OSM.


# Add useful (for our purpose) methods to the osmapis.OSM class
class ElvegOSM(osmapis.OSM):

    def __init__(self, items=()):
        # First call the parent's __init__
        super(ElvegOSM, self).__init__(items)

        # Generate dict with TRANSID as key and is as value
        self.wayid_dict = {}
        for wayid,way in self.ways.iteritems():
            transid = way.tags['TRANSID']
            self.wayid_dict[transid] = wayid

    def way_nodes_from_transid(self, transid):
        wayid = self.wayid_dict[transid]
        way = self.ways[wayid]
        node_ids = way.nds
        nodes = [osmobj.nodes[nid] for nid in node_ids]
        return nodes

    def distances_from_transid(self, transid):
        global g
        nodes = self.way_nodes_from_transid(transid)
        node_distances = []
        distance_so_far = 0.
        prev_lon = nodes[0].lon
        prev_lat = nodes[0].lat
        
        for i,nd in enumerate(nodes):
            #az1,az2,d_from_previous = g.inv(prev_lon, prev_lat, nd.lon, nd.lat)
            ggresults = gg.Geodesic.WGS84.Inverse(prev_lat, prev_lon, nd.lat, nd.lon)
            d_from_previous = ggresults['s12']
            if i != 0 and d_from_previous < 0.5:
                # Report if very short distance
                warn(u"Short distance ({2}) for transid {0} to node No. {1}".format(transid, i,d_from_previous))
            distance_so_far += d_from_previous
            node_distances.append(distance_so_far)
            # Prepare previous coordinates for next round
            prev_lon = nd.lon
            prev_lat = nd.lat

        return node_distances

class ElvegNode(osmapis.Node):

    def __init__(self, attribs={}, tags={}):
        osmapis.Node.__init__(self, attribs, tags)
        # Make sure the class counter is as low as the lowest existing ID
        # This should probably have been done in osmapis.Node
        if self.id is not None:
            self.__class__._counter = min(self.__class__._counter, self.id)
            
class ElvegWay(osmapis.Way):

    def __init__(self, attribs={}, tags={}, nds=()):
        osmapis.Way.__init__(self, attribs, tags, nds)
        # Make sure the class counter is as low as the lowest existing ID
        # This should probably have been done in osmapis.Way
        if self.id is not None:
            self.__class__._counter = min(self.__class__._counter, self.id)


# Override default classes in osmapis.py
osmapis.wrappers["osm"]  = ElvegOSM
osmapis.wrappers["node"] = ElvegNode
osmapis.wrappers["way"]  = ElvegWay

def warn(warning):
    warning = warning.encode('utf-8')
    sys.stderr.write(warning + '\n')

def waynode_from_coord(coord):
    # This assumes that there is only one node for a given
    # coordinates that is part of a way.
    global way_node_ids
    global node_lookup
    way_nodes = [nid for nid in node_lookup[coord] if nid in way_node_ids]
    if len(way_nodes) > 1:
        sys.stderr.write('More than one way nodes at coordinate:\n')
        sys.stderr.write(str(coord) + '\n')
    elif len(way_nodes) == 0:
        #sys.stderr.write('No way nodes at coordinate:\n')
        #sys.stderr.write(str(coord) + '\n')
        return None
    return way_nodes[0]

def merge_nodes(node_id_list):
    global osmobj
    # Record the attributes of the first node
    first_attr = osmobj.nodes[node_id_list[0]].attribs
    del first_attr['id']
    # Join the way_ids lists of the nodes
    way_ids = set()
    for node_id in node_id_list:
        way_ids.update(osmobj.nodes[node_id].way_ids)
    # Join the tags
    merged_tags = {}
    for node_id in node_id_list:
        for key,value in osmobj.nodes[node_id].tags.iteritems():
            if merged_tags.has_key(key):
                # Potential conflict, but only if value is different
                if merged_tags[key] != value:
                    # A conflict for real
                    msg = u"Conflict values when merging tag {0} from node {1}: {2} and {3}".format(
                                                                            key, node_id, merged_tags[key], value)
                    warn(msg)
            # No conflict, so copy tag
            merged_tags[key] = value
        # Delete the node
        osmobj.nodes.pop(node_id)
    # Create a new node
    merged_node = ElvegNode(attribs=first_attr, tags=merged_tags)
    merged_node.way_ids = way_ids
    osmobj.add(merged_node)
    # Replace deleted node_ids with new node_id in all affected ways
    for way_id in way_ids:
        way = osmobj.ways[way_id]
        for i,way_node_id in enumerate(way.nds):
            if way_node_id in node_id_list:
                way.nds[i] = merged_node.id
    
def create_osmtags(elveg_tags):
    '''Create tags based on standard tags in ????Elveg_default.osm'''

    category2highwayclass = {'E': 'trunk',     # Europaveg
                             'R': 'trunk',     # Riksveg
                             'F': 'secondary', # Fylkesveg, could also be primary
                             'K': 'road',      # Kommunal veg
                             'P': 'road',      # Privat veg
                             'S': 'unclassified'}     # Skogsbilveg, possibly more info in the LBVKLASSE tag

    roadOBJTYPEs = set([u'VegSenterlinje', 
                        u'Svingekonnekteringslenke',
                        u'Kj\xf8refelt',
                        u'Kj\xf8rebane'])

    osmtags = dict()

    # Roads and ferry routes share many tags, and are therefore
    # treated together
    if elveg_tags['OBJTYPE'] in roadOBJTYPEs.union([u'Bilferjestrekning']) :

        # Split VNR tag
        # The "vegnummer" tag is optional, but let's assume it is always present for now
        # (i.e. fix it if it causes problems)
        if elveg_tags.has_key('VNR'):
            vegkategori,vegstatus,vegnummer = [s.strip(':;') for s in elveg_tags['VNR'].split()]
        else:
            warn(u"VNR missing for OBJTYPE {OBJTYPE} with TRANSID {TRANSID}".format(**elveg_tags))
            return osmtags

        # There are more vegstatus values than listed in https://wiki.openstreetmap.org/w/images/c/cc/Elveg_SOSI_4.0_2008.pdf
        # There is a more complete list in chapter 7.3.11 in 
        # http://www.statkart.no/Documents/Standard/SOSI-standarden%20del%201%20og%202/SOSI%20standarden/Vegnett.pdf

        if elveg_tags['OBJTYPE'] in roadOBJTYPEs:
            # Set the road category
            if vegstatus in ['V','T','W']: # Eksisterende veg, Veg med midlertidig status, Midlertidig veg mer enn et aar
                osmtags['highway'] = category2highwayclass[vegkategori]
            elif vegstatus == 'A':
                osmtags['highway'] = 'construction'
                osmtags['construction'] = category2highwayclass[vegkategori]
            elif vegstatus == 'G':
                osmtags['FIXME'] = u'Veggrunn, ikke trafikkform\xe5l. Select appropriate road type.'
                osmtags['highway'] = 'road'
            elif vegstatus == 'M':
                osmtags['DEBUG'] = u'M\xf8te- og rasteplasser'
            elif vegstatus in ['P','Q']: # Vedtatt veg, planlagt veg
                osmtags['DEBUG'] = 'Vedtatt (P) eller planglagt (Q): ' + vegstatus
                osmtags['action'] = 'delete'
            else:
                warn(u"Unknown vegstatus {0} for {2} with TRANSID {1}".format(vegstatus,elveg_tags['TRANSID'],elveg_tags['OBJTYPE']))
        elif elveg_tags['OBJTYPE'] == u'Bilferjestrekning':
            # Set the class for the ferry route
            if vegstatus == 'S':
                osmtags['route'] = 'ferry'
                osmtags['class'] = category2highwayclass[vegkategori]
            elif vegstatus in ['E','F']: # Vedtatt fergestrekning, planlagt fergestrekning
                osmtags['DEBUG'] = 'Vedtatt fergestrekning, planlagt fergestrekning ' + vegstatus
                osmtags['action'] = 'delete'
            else:
                warn(u"Ferry route with TRANSID {0} has unknown vegstatus {1}".format(elveg_tags['TRANSID'],vegstatus))

        # Add ref to road kategories Europaveg, Riksveg and Fylkesveg
        if vegkategori == 'E':
            osmtags['ref'] = 'E ' + vegnummer
        elif vegkategori in ['R', 'F']:
            osmtags['ref'] = vegnummer

    # Gang- og sykkelveg. Only a fraction of all of those are in the data. 
    # Nevertheless, include those that are.
    elif elveg_tags['OBJTYPE'] == 'GangSykkelVegSenterlinje':
        osmtags['highway'] = 'cycleway'
        osmtags['foot'] = 'yes'

    # OBJTYPE=Fortau is sometimes used when a Gang- og sykkelveg goes over 
    # in a sidewalk for a while
    # A sidewalk is usually best represented as a sidewalk=* on a road,
    # but at least in the conversion we let it be a separate way.
    elif elveg_tags['OBJTYPE'] == 'Fortau':
        osmtags['highway'] = 'footway'
        osmtags['footway'] = 'sidewalk' 
        osmtags['note'] = 'Consider adding sidewalk as a tag on the road'
    
    # Import OBJTYPE=u'Frittst\xe5endeTrapp'
    # There are many objects in Bergen (1201) and quite a few in 
    # Stavanger (1103) and Sandnes (1102) as well.
    # They are often integrated with the network of footways.
    # There seems to be no consistent direction of the ways, 
    # so do not set incline=up/down
    elif elveg_tags['OBJTYPE'] == u'Frittst\xe5endeTrapp':
        osmtags['highway'] = 'steps'

    # OBJTYPE not handled - add deletion tag and return
    else:
        warn(u"Deleting unimplemented OBJTYPE {OBJTYPE} with TRANSID {TRANSID}".format(**elveg_tags))
        osmtags['DEBUG'] = 'OBJTYPE not handled: ' + elveg_tags['OBJTYPE']
        osmtags['action'] = 'delete'
        return osmtags

    ### Finished switching between OBJTYPEs
    ### From now on we have one of the known OBJTYPEs above

    # Add information about lanes from the VKJORFLT tag (oneway=*, lanes=*)
    if elveg_tags.has_key('VKJORFLT'):
        # This probably only applies to roads and ferry routes - verify that
        if elveg_tags['OBJTYPE'] not in roadOBJTYPEs and elveg_tags['VKJORFLT'] != '1#2'    :
            #print elveg_tags
            warn(u"Processing VKJORFLT tag for OBJTYPE {OBJTYPE} for TRANSID {TRANSID}: {VKJORFLT}".format(**elveg_tags))
        # Use the parse_lanes() function find the correct tags
        lane_tags = parse_lanes(elveg_tags['VKJORFLT'])
        osmtags.update(lane_tags)

    # Import GATENAVN for any type of way, although it would probably only exist for road objects
    # There are some empty GATENAVN values in the data set - do not set a name for those
    if elveg_tags.has_key('GATENAVN') and elveg_tags['GATENAVN'] != '':
        osmtags['name'] = elveg_tags['GATENAVN']

    # Add information about tunnels and bridges from MEDIUM tag
    if elveg_tags.has_key('MEDIUM'):
        # Give a warning if this tag is on a non-road object
        if elveg_tags['OBJTYPE'] not in roadOBJTYPEs.union([u'GangSykkelVegSenterlinje']):
            warn(u"Processing MEDIUM tag for OBJTYPE {OBJTYPE} for TRANSID {TRANSID}: {MEDIUM}".format(**elveg_tags))
        if elveg_tags['MEDIUM'] == 'L':
            osmtags['bridge'] = 'yes'
            osmtags['layer'] = '1'
        elif elveg_tags['MEDIUM'] == 'U':
            osmtags['tunnel'] = 'yes'
            osmtags['layer'] = '-1'
        elif elveg_tags['MEDIUM'] == 'B':
            # B means "through a building".
            # This could be tagged with covered=yes (current tagging
            # for Perleporten in Trondheim), but tunnel=building_passage
            # seems to be preferred.
            warn(u"Processing MEDIUM tag 'B' for OBJTYPE {OBJTYPE} for TRANSID {TRANSID}".format(**elveg_tags))
            osmtags['tunnel'] = 'building_passage'
        else:
            # There should be no other possible values for MEDIUM
            warn(u"Unknown MEDIUM value '{MEDIUM}' for OBJTYPE {OBJTYPE} for TRANSID {TRANSID}".format(**elveg_tags))

    # Add the nvdb:id tag from the TRANSID tag
    # All ways should have a TRANSID (will change to LOKALID with SOSI 4.5)
    osmtags['nvdb:id'] = elveg_tags['TRANSID']


    return osmtags

def parse_lanes(lane_string):
    lane_tags = dict()

    # Strip whitespace from lane_string
    lane_string = lane_string.strip()

    # Early exit for the most coommone values
    if lane_string == '1#2':
        # Most common case - one lane in each direction - no special tags
        pass
    elif lane_string in ('1', '3'):
        # One-way street along way direction
        lane_tags['oneway'] = 'yes'
    elif lane_string in ('2', '4'):
        # One-way street opposite to way direction
        lane_tags['oneway'] = '-1'
    elif lane_string == '1#3':
        # One-way street along way direction
        lane_tags['oneway'] = 'yes'
        lane_tags['lanes'] = '2'
    elif lane_string == '2#4':
        # One-way street along way direction
        lane_tags['oneway'] = '-1'
        lane_tags['lanes'] = '2'
    elif lane_string == '1#3#5':
        # One-way street along way direction
        lane_tags['oneway'] = 'yes'
        lane_tags['lanes'] = '3'
    elif lane_string == '2#4#6':
        # One-way street along way direction
        lane_tags['oneway'] = '-1'
        lane_tags['lanes'] = '3'
    elif lane_string == '':
        # Sometimes this tag is empty -- assume that this means nothing special
        pass
    # TURN LANES
    elif lane_string in ('1V1', '1H1', '1V2','1H2'):
        # Left/right turning lane - mark as one-way
        lane_tags['oneway'] = 'yes'
    elif lane_string in ('2V1','2H1', '2V2', '2H2'):
        # Left/right turning lane - mark as one-way
        lane_tags['oneway'] = '-1'
    elif lane_string in ('1#1V1', '1#1H1', '1#1H1#1V1'):
        # One lane and additional turn lane - ignore the turn lane
        lane_tags['oneway'] = 'yes'
    elif lane_string in ('2#2V1', '2#2H1', '2#2H1#2V1'):
        # One lane and additional turn lane - ignore the turn lane
        lane_tags['oneway'] = '-1'
    elif lane_string in ('1#2#2H1', '1#1H1#2', '1#2#2V1', '1#1V1#2', '1#2#1H1', '1#2#1V1'):
        # One lane each direction and additional turn lane - ignore the turn lane
        pass
    # BIKE LANES
    elif lane_string == '1#2#3S#4S':
        # Two lanes with bike lanes in both directions
        lane_tags['cycleway'] = 'lane'
    elif lane_string == '1#3S':
        # One lane plus bike lane
        lane_tags['oneway'] = 'yes'
        lane_tags['cycleway'] = 'lane'
    elif lane_string == '2#4S':
        # One lane plus bike lane
        lane_tags['oneway'] = '-1'
        lane_tags['cycleway'] = 'lane'
    # BUS LANES
    elif lane_string == '1#3K':
        # One lane plus bus lane
        lane_tags['oneway'] = 'yes'
        lane_tags['lanes'] = '2'
        lane_tags['lanes:psv'] = '1'
    elif lane_string == '2#4K':
        # One lane plus bus lane
        lane_tags['oneway'] = '-1'
        lane_tags['lanes'] = '2'
        lane_tags['lanes:psv'] = '1'
    elif lane_string == '2#4#6K':
        # One lane plus bus lane
        lane_tags['oneway'] = '-1'
        lane_tags['lanes'] = '3'
        lane_tags['lanes:psv'] = '1'
    elif lane_string in ('1K', '3K'):
        # Single bus lane
        lane_tags['oneway'] = 'yes'
        lane_tags['psv'] = 'designated'
    else:
        # TODO: Split lane string into individual lanes
        # Postfix H1, H2, V1, V2 are for turning lanes, 
        # postfix K is for public service vehicles (PSV)
        # postfix O is for "waiting lanes", e.g. at ferry terminals.
        lane_tags['FIXME'] = "Consider adding lane tags based on Elveg data: {0}".format(lane_string)
        #warn("Unhandled VKJORFLT: " + lane_string)
    return lane_tags

def split_way(osmobj, way_id, split_points):
    '''Split way at split points.

    Return list of way ids for the split way. The first id is of the
    original way.

    '''
 
   # Do not go through the hassle, if the way needs no splitting
    if len(split_points) == 0:
        return [way_id]

    # Initialize a list of way id's of the new ways (to be returned)
    # Since the last way is always split off first, the list will be
    # in reverse order, and is turned around at the end.
    splitway_id_list = []

    # Get the way that is to be split
    way = osmobj.ways[way_id]
    transid = way.elveg_tags['TRANSID']

    # Compute distances from start to each node of way
    node_distances = osmobj.distances_from_transid(transid)
    geo_length = node_distances[-1]

    # Compute VPA length and normalize split_points to geographic length
    if way.elveg_tags.has_key("VPA"):
        vpa = [int(n.strip(':;')) for n in way.elveg_tags["VPA"].split()]
    else:
        # These roads are probably not split, so 1.0 is fine, but raise Exception for now
        #corrction_factor = 1.0
        raise KeyError("VPA Elveg tag not present")
    vpa_length = vpa[2] - vpa[1]
    normalization_factor = geo_length / float(vpa_length)
    split_points_normalized = [normalization_factor * sp for sp in split_points]

    # Make sure the normalized split points are sorted
    # (so that we can split off ways from the end of the list)
    split_points_normalized.sort()

    # Loop over the split points, splitting off the last way each time
    while len(split_points_normalized) > 0:
        current_split_point = split_points_normalized.pop()
        upper_split_index = np.searchsorted(node_distances, current_split_point)

        # Find the distance to the nearest nodes
        # (for checking if a new node should be created)
        distance_to_upper = node_distances[upper_split_index] - current_split_point
        distance_to_lower = current_split_point - node_distances[upper_split_index - 1]
        
        # Decide if a new node should be created
        # Reuse node if closer than 0.5 m
        if distance_to_upper < 0.5 or distance_to_lower < 0.5:
            # Verify that we have no negative distances (which is a bug)
            if distance_to_upper < 0. or distance_to_lower < 0.:
                warn(u"Negative distances for TRANSID {0}".format(transid))
            # Reuse closest node
            if distance_to_upper < distance_to_lower:
                split_index = upper_split_index
            else:
                split_index = upper_split_index - 1
            # Create a new way from the split node to the end of the way
            newway_nodes = way.nds[split_index:]
            newway = ElvegWay(tags=way.tags, nds=newway_nodes)
            splitway_id_list.append(newway.id)
            osmobj.ways[newway.id] = newway
            
            # Remove the new way from the old way
            # (the split_index should be included in both ways)
            way.nds = way.nds[:split_index + 1]
            
        else:
            # Find the coordinates for the new split node
            from_node_id = way.nds[upper_split_index - 1]
            to_node_id = way.nds[upper_split_index]
            from_node = osmobj.nodes[from_node_id]
            to_node = osmobj.nodes[to_node_id]
            ggresults = gg.Geodesic.WGS84.Inverse(from_node.lat, from_node.lon, to_node.lat, to_node.lon)
            distance = ggresults['s12']
            azi1 = ggresults['azi1']
            dist_from_last_node = current_split_point - node_distances[upper_split_index - 1]
            ggresults = gg.Geodesic.WGS84.Direct(from_node.lat, from_node.lon, azi1, dist_from_last_node)
            newlon = ggresults['lon2']
            newlat = ggresults['lat2']

            # Create the new node
            split_node = ElvegNode(attribs={"lon": newlon, "lat": newlat})
            if osmobj.nodes.has_key(split_node.id):
                # This should not happen if ElvegNode.__init__() does the right thing
                raise Exception(u"Almost overwrote node {0}\n".format(split_node.id).encode('utf-8'))
            osmobj.nodes[split_node.id] = split_node

            # FOR DEBUGGING WAY SPLITTING
            #osmobj.nodes[split_node.id].tags['newnode'] = 'yes'

            # Create a new way from the split_point to the end of the way
            newway_nodes = [split_node.id] + way.nds[upper_split_index:]
            newway = ElvegWay(tags=way.tags, nds=newway_nodes)
            splitway_id_list.append(newway.id)
            osmobj.ways[newway.id] = newway

            # Remove nodes for the new way from the old way
            way.nds = way.nds[:upper_split_index] + [split_node.id]

    # Finally, add the original way, which is the first segment of the
    # newly split way.
    splitway_id_list.append(way_id)

    # Reverse direction so that first way segment comes first
    return splitway_id_list[::-1]


###########################################################
#           main                                          #
###########################################################
        
# Read input arguments
directory = sys.argv[1]
if len(sys.argv) >= 3:
    kommune_number = sys.argv[2]
else:
    kommune_number = directory.strip('/')[-4:]
    # Check that it is really a number
    kummune_int = int(kommune_number)


# Find the names of the other files
osm_input = os.path.join(directory, kommune_number + 'Elveg_default.osm')
osm_output = os.path.join(directory, kommune_number + 'Elveg.osm')
elveg_fart = os.path.join(directory, kommune_number + 'Fart.txt')
elveg_hoyde = os.path.join(directory, kommune_number + 'Hoyde.txt')
osm_barrier_output = os.path.join(directory, kommune_number + 'detatched_barriers.osm')
osm_deleted_output = os.path.join(directory, kommune_number + 'deleted_elements.osm')

# Loop over speed limits and tags where the whole 
# way where possible. Other places, add to split list
roaddata = {}
with open(elveg_fart, 'rb') as ef:
    # Read first four header lines
    ef_header = ef.next()
    ef_export_line = ef.next()
    ef_some_number = ef.next()
    ef_empty_line = ef.next()
    
    # Then use csv module for reading data
    reader = csv.DictReader(ef, delimiter=';')
    for row in reader:
        transid = row[' TransID']

        fart_start = int(row['Fra'])
        fart_stop =  int(row['   Til'])
        fart_length = fart_stop - fart_start
        fart_limit = row[' Fart']

        if not roaddata.has_key(transid):
            roaddata[transid] = {}
        if not roaddata[transid].has_key('maxspeed'):
            roaddata[transid]['maxspeed'] = []
        roaddata[transid]['maxspeed'].append({'maxspeed': fart_limit,
                                              'start': fart_start,
                                              'stop': fart_stop})
                                              
# Add height limits to roaddata (if the file exists)
if not os.path.isfile(elveg_hoyde):
    warn(u"File {0} does not exist and is not used".format(elveg_hoyde))
else:
    with open(elveg_hoyde, 'rb') as eh:
        # Read first four header lines
        eh_header = eh.next()
        eh_export_line = eh.next()
        eh_empty_line1 = eh.next()
        eh_empty_line2 = eh.next()

        # Then use csv module for reading data
        reader = csv.DictReader(eh, delimiter=';')
        for row in reader:
            transid = row[' TransID']

            height_start = int(row['Fra'])
            height_stop =  int(row['   Til'])
            height_length = height_stop - height_start
            height_limit = row['H\xf8yde']

            if not roaddata.has_key(transid):
                roaddata[transid] = {}
            if not roaddata[transid].has_key('maxheight'):
                roaddata[transid]['maxheight'] = []
            roaddata[transid]['maxheight'].append({'maxheight': height_limit,
                                                   'start': height_start,
                                                   'stop': height_stop})

# TODO: Add information from XXXXAksel.txt to roadddata,
# and add relevant tagging.

# Read OSM file
osmobj = ElvegOSM.load(osm_input)

# Loop through all nodes and move tags to elveg_tags
for nid,node in osmobj.nodes.items():
    node.elveg_tags = node.tags
    node.tags = {}

# Loop through all ways in osmobj and 
# - swap original tags with OSM tags.
# - extract the way length from the Elveg VPA tag and
#   store in roaddata structure
# Important to use items() instead of iteritems() here as we are adding
# items to the obmobj.ways dictionary during the loop.
for wid,w in osmobj.ways.items():
    # Add new tags (using the create_osmtags function)
    w.elveg_tags = w.tags
    osm_tags = create_osmtags(w.elveg_tags)
    w.tags = osm_tags

    # Check that way has VPA Elveg tag
    if not w.elveg_tags.has_key('VPA'):
        warn(u"VPA missing for OBJTYPE {OBJTYPE} with TRANSID {TRANSID}".format(**w.elveg_tags))
        continue

    # Add way length as given by VPA to the roadddata structure
    transid = w.elveg_tags['TRANSID']
    vpa = [int(n.strip(':;')) for n in w.elveg_tags["VPA"].split()]
    # We do not care about those ways where we have no data to add,
    # so move to next if this is the case.
    if not roaddata.has_key(transid):
        continue
    roaddata[transid]['length'] = vpa[2] - vpa[1]
    
    # make a sorted list of meter values, including end
    # points, where some roaddata may change
    end_points = [0, roaddata[transid]['length']]
    for restriction_type in ['maxspeed', 'maxheight']: # Add any new restrictions here
        for endpoint_type in ['start', 'stop']:
            end_points.extend([d[endpoint_type] for d in roaddata[transid].get(restriction_type, [])])
    end_points = list(set(end_points))
    end_points.sort()
    
    # Test endpoints from .txt files against VPA lengths
    # There are several ways where the end point is outside the VPA meter range
    # Remove those TRANSIDs from the roaddata structure and move on to next way
    if end_points[-1] > roaddata[transid]['length']:
        warntemplate = u"Warning: End point {0} m outside of VPA length of road ({1} m) for TRANSID {2}"
        warnstring = warntemplate.format(end_points[-1], roaddata[transid]['length'], transid)
        warn(warnstring)
        del roaddata[transid]
        continue

    # Make a list of intervals, representing the new ways after a split
    # For most ways, there will be only one interval, but whenever
    # the speed limit changes on a way or a height restriction
    # does not apply to the whole way, there will be more than one interval
    interval_list = zip(end_points[:-1],end_points[1:])

    # Make a list of tags (maxheight=*, maxspeed=*)
    # with one list entry per new way interval
    newway_tags = [{} for i in interval_list] # I.e. a list of empty dicts
    for i,interval in enumerate(interval_list):
        for restriction_type in ['maxspeed', 'maxheight']: # Add any new restrictions here
            for j,restr in enumerate(roaddata[transid].get(restriction_type, [])):
                if restr['start'] <= interval[0] and interval[1] <= restr['stop']:
                    newway_tags[i][restriction_type] = restr[restriction_type]

    # DEBUG: Remove later
    #print newway_tags

    # Split the way in osmobj into the right number of segments
    split_points = end_points[1:-1]
    segment_ids = split_way(osmobj, w.id, split_points)

    # Add nvdb:id:part subkey to each part if the Elveg segment has been split
    if len(segment_ids) > 1:
        for i,segment_id in enumerate(segment_ids):
            osmobj.ways[segment_id].tags['nvdb:id:part'] = str(i)
    
    # Add maxheight and maxspeed restrictions
    for i,segment_id in enumerate(segment_ids):
        osmobj.ways[segment_id].tags.update(newway_tags[i])

# Loop through all ways 
# - make a set of those nodes that are part of a way
way_node_ids = set()
for way in osmobj.ways.values():
    way_node_ids.update(way.nds)
# ... and those that are not part of a way
noway_node_ids = set(osmobj.nodes).difference(way_node_ids)

# Add way_id variable to every node which holds the way_ids of all ways
# it is part of.
for node in osmobj.nodes.itervalues():
    node.way_ids = set()
for way in osmobj.ways.itervalues():
    for node_id in way.nds:
        node = osmobj.nodes[node_id]
        node.way_ids.add(way.id)
        


# DATA CHECKING: Check if any way nodes also have tags, or if all tags
# are on duplicate nodes
#for waynode_id in way_nodes:
#    waynode = osmobj.nodes[waynode_id]
#    if len(waynode.tags) > 0:
#        print waynode.tags


# Create OSM object for manual merging of off-way barriers
osmobj_barriers = ElvegOSM()

# Loop through and process all single nodes
for nid in noway_node_ids:
    noway_node = osmobj.nodes[nid]
    coord = (noway_node.lat, noway_node.lon)
    if noway_node.elveg_tags['OBJTYPE'] == 'Vegsperring':
        # Tag the barrier with OSM tags
        vegsperringtype = noway_node.elveg_tags['VEGSPERRINGTYPE']
        if vegsperringtype == 'Betongkjegle':
            noway_node.tags['barrier'] = 'block'
        elif vegsperringtype == u'Bilsperre':
            # This seems to be any type of barrier that has wide enough
            # openings to only stop cars.
            noway_node.tags['barrier'] = 'yes'
        elif vegsperringtype == u'Bussluse':
            noway_node.tags['barrier'] = 'bus_trap'
        elif vegsperringtype == u'L\xe5st bom':
            noway_node.tags['barrier'] = 'gate'
        elif vegsperringtype == u'New Jersey':
            noway_node.tags['barrier'] = 'jersey_barrier'
        elif vegsperringtype == u'R\xf8rgelender':
            # This describes the material more than the actual barrier
            # Similar to barrier=fence, but usually it is possible to
            # walk or bike around
            noway_node.tags['barrier'] = 'yes'
        elif vegsperringtype == u'Steinblokk':
            noway_node.tags['barrier'] = 'block'
        elif vegsperringtype == u'Trafikkavviser':
            # It seems that roads with this kind of barrier are 
            # best tagged as footways in OSM.
            # I suppose the barrier itself could be anything.
            noway_node.tags['barrier'] = 'yes'
        elif vegsperringtype == u'Ukjent':
            noway_node.tags['barrier'] = 'yes'
        else:
            warn(u"Unknown barrier: {0}".format(vegsperringtype))
            noway_node.tags['barrier'] = 'yes'
    elif noway_node.elveg_tags['OBJTYPE'] == 'Kommunedele':
        # We do not use this tag, mark this node for deletion
        noway_node.tags['DEBUG'] = 'Kommunedele'
        noway_node.tags['action'] = 'delete'
    elif noway_node.elveg_tags['OBJTYPE'] == 'Ferjekai':
        # These nodes are not connected to the road network
        # In OSM, they should ideally be on the node between the road and the ferry route.
        noway_node.tags['amenity'] = 'ferry terminal'

# TODO: Add amenity="ferry terminal" on nodes with OBJTYPE=Ferjekai

# Remove all ways and non-way nodes with action=delete and delete unused nodes

# Loop through ways, collect ways with action=delete and
# id of nodes in ways
to_delete = set()
nodes_used = set()
for way in osmobj.ways.itervalues():
    if "action" in way.tags and way.tags['action'] == 'delete':
        to_delete.add(way)
    elif "nvdb:id" in way.tags and len(way.tags) == 1:
        to_delete.add(way)
    elif len(way.tags) == 0:
        to_delete.add(way)
    else:
        for n in way.nds:
            nodes_used.add(n)

# Collects nodes which should be deleted
for node in osmobj.nodes.itervalues():
    if "action" in node.tags and node.tags['action'] == 'delete':
        to_delete.add(node)
    elif (node.id not in nodes_used) and (len(node.tags)) == 0:
        to_delete.add(node)

# Delete objects from output and add them to a separate file
osmobj_deleted = ElvegOSM()
for element in to_delete:
    osmobj.discard(element)
    osmobj_deleted.add(element)
    if hasattr(element, 'elveg_tags'):
        element.tags.update({'Elveg:' + k:v for k,v in element.elveg_tags.iteritems()})

# Copy nodes needed by ways in osmobj_deleted
for delway in osmobj_deleted.ways.itervalues():
    for nid in delway.nds:
        if not osmobj_deleted.nodes.has_key(nid):
            osmobj_deleted.add(osmobj.nodes[nid])

# Make a table with hash of indices of the nodes, in order to identify
# nodes with (exactly) the same coordinates
node_lookup = dict()
for id,node in osmobj.nodes.iteritems():
    key = (node.lat, node.lon)
    if node_lookup.has_key(key):
        node_lookup[key].append(id)
    else:
        node_lookup[key] = [id]

# Merge nodes in same location
for coord, node_id_list in node_lookup.items():
    if len(node_id_list) == 1:
        continue
    else:
        merge_nodes(node_id_list)

# Speed limit cleanup
for id,way in osmobj.ways.iteritems():
    # Remove the default speed limit of 50, since that may
    # be only due to missing reporting
    if (way.tags.get('maxspeed', None) == '50' and 
            way.tags.get('highway',None) not in ('trunk', 'secondary')):
        del way.tags['maxspeed']
    # Remove speed limits for non-roads (footway, cycleway, etc.)
    if (way.tags.has_key('maxspeed') and
            way.tags.get('highway', None) not in ('trunk', 'secondary', 'road', 'unclassified')):
        del way.tags['maxspeed']

# Save barriers that are not merged to other nodes to a separate file
for id,node in osmobj.nodes.items():
    if not node.way_ids:
        osmobj_barriers.nodes[id] = node
        del osmobj.nodes[id]


# TODO: Add turn restrictions from XXXXSving.txt

osmobj.save(osm_output)
osmobj_barriers.save(osm_barrier_output)
osmobj_deleted.save(osm_deleted_output)



