import random
import math
import numpy as np
import os.path
import re
import time
import datetime
import networkx as nx
import pickle

start = time.time()

"""
Header: This script is used to convert a SPRNT input file into a SWMM input file

Required Inputs: 
    SPRNT.spt   - SPRNT input file to parse network data
    Sample.inp  - from which to build out the SWMM input file
    
Optional Inputs:
    IC_H_file               - text file that contains the node initial condition "InitDepth"
    IC_Q_file               - text file that contains the link initial conditions "InitFlow"
    path_to_external_files  - FULLPATH to external files, so the SWMM program can find them
    
Manual Settings:
    ExtFiles_Condition  - Boolean, determines where the forcing data is written (True = External, False = Internal)
                                Writing to External Files is necessary for continuous/complicated forcing data
    IC_condition        - Boolean, determines if the Junction/Conduit initial conditions will be written from files
                                (True = use IC_H/Q files, False = set ICs to 0)
    Cross_section_ID    - Integer, determines which cross-sectional transformation to use
                                (0 = Trapezoidal, 1 = AY-HEC2 Transform)
                                
Outputs:
    SWMM_Final.inp - SWMM input file that can be opened with the GUI or executed with a SWMM wrapper
    ExtFiles/*.dat - Sub-directory containing forcing data in .dat text files (named for the SWMM object they apply to)
"""

def settings():
    """
    This function contains the required input information and manual settings for the SPRNTtoSWMM code
    """
    global start, sprnt_file, swmm_example, swmm_final, path_to_external_files, IC_condition, IC_H_file, IC_Q_file, \
        ExtFiles_Condition, cross_section_ID, cross_section_transformation_tolerance

    # Input Files with Relative Paths
    sprnt_file = r"./SPT_Files/Trinity_1hr_interval.spt"
    swmm_example = r"Lavaca_Sample.inp"
    # IC_H_file = r"./Initial_Condition_Q_H_Data/San_Jacinto_H.txt"
    # IC_Q_file = r"./Initial_Condition_Q_H_Data/San_Jacinto_Q.txt"
    path_to_external_files = 'C:/Users/edt489/Documents/Store_SPRNTtoSWMM/ExtFiles_CorrectedEndTime/Trinity_ExtFiles'

    #
    swmm_final = r"C:/Users/edt489/Documents/Store_SPRNTtoSWMM/ExtFiles_CorrectedEndTime/TRB_extra.inp"
    ExtFiles_Condition = True

    IC_condition = False

    max_segments_per_comid_digits = 4
    cross_section_ID = 0
    cross_section_transformation_tolerance = 0.001
    max_inflows = 500

    return

def read_sprnt_to_contents():
    global sprnt_contents

    f_sprnt = open(sprnt_file, "r")
    sprnt_contents = f_sprnt.readlines()
    read_sprnt_contents = time.time() - start
    print("Read the SPRNT file time: ", read_sprnt_contents)

    return


def global_definitions():
    """This function mostly defines the globals that are used as indices for each of the arrays"""

    global node_id, node_sR, node_n, node_zR, node_hR, node_H_IC, \
            geo_bottomwidth, geo_slope, \
            segment_up, segment_down, segment_length, segment_Q_IC, \
            junction_down, junction_up1, junction_coeff1, junction_up2, junction_coeff2, \
            lateral_source, lateral_timeunit, q_source, q_timeunit, \
            root, root_bc

    node_id = 0
    node_sR = 1
    node_n = 2
    node_zR = 3
    node_hR = 4
    node_H_IC = 5

    geo_bottomwidth = 0
    geo_slope = 1

    segment_up = 0
    segment_down = 1
    segment_length = 2
    segment_Q_IC = 3

    junction_down = 0
    junction_up1 = 1
    junction_coeff1 = 2
    junction_up2 = 3
    junction_coeff2 = 4

    lateral_source = 0
    lateral_timeunit = 1

    q_source = 0
    q_timeunit = 1

    # root = ''
    # root_bc = ''
    return


def find_boundarycondition():
    """
    This function serves two purposes
        1) Find the boundary node - root_location
        2) Determine if the root boundary condition is changing.  If it is I'll need to use a timeseries and I haven't
            put in that functionality yet
    """
    root_bc_0 = 0.0
    root_bc_1 = 0.1
    kk = 0
    for ii in range(len(sprnt_contents)):
        if sprnt_contents[ii].find("boundarycondition") != -1:
            jj = ii + 1
            templine = re.split(' location=| type=|\n', sprnt_contents[jj])
            root_location = templine[1]
            if templine[2] != 'depth':
                print("Different root boundary condition: ", templine[2])
                quit()
            jj = jj + 2
            while sprnt_contents[jj].find('end') == -1:
                try:
                    splitline = sprnt_contents[jj].split()[1].split('=')
                    if kk == 0:
                        root_bc_0 = splitline[1]
                    if kk > 0:
                        root_bc_1 = splitline[1]
                    kk = kk + 1
                except:
                    nothing = ii
                jj = jj + 1
            if root_bc_0 != root_bc_1:
                print("The root boundary condition is changing: ", str(root_bc_0), str(root_bc_1))
                quit()
            else:
                root_bc = root_bc_0

    return root_location, root_bc


def read_sprnt_for_counters():
    """
    This function inspects the contents of the SPRNT input file and creates correctly sized arrays for each of the data
    types.

    This allows us to use np.arrays rather than 2D lists for nodes_name, nodes_geo, segments, and junctions.
    lateralsources and qsources need to be 2D lists because of the inconsistent timeseries entries

    The second dimension of the arrays is hard-coded to reflect that each data type in SPRNT contains a different number
     of parameters.

    In the San Antonio - Guadalupe River basin the SPRNT input file is 3.8 millions lines.
    """
    global total_count, null_count, node_count, segment_count, lateral_count, junction_count, qsource_count, \
        nodes_name, nodes_geo, segments, lateralsources, junctions, qsources, boundaryconditions

    total_count = len(sprnt_contents)
    null_count, node_count, segment_count, lateral_count, junction_count, qsource_count = 0, 0, 0, 0, 0, 0

    for ii in range(total_count):
        if sprnt_contents[ii].find("#") != -1:
            null_count = null_count + 1
        elif sprnt_contents[ii].find('node') != -1:
            node_count = node_count + 1
        elif sprnt_contents[ii].find("segment") != -1:
            segment_count = segment_count + 1
        elif sprnt_contents[ii].find("lateralsource") != -1:
            lateral_count = lateral_count + 1
        elif sprnt_contents[ii].find("junction") != -1:
            junction_count = junction_count + 1
        elif sprnt_contents[ii].find("qsource") != -1:
            qsource_count = qsource_count + 1

    # The sizes of these arrays is network dependent, may need to write some functions to generalize them
    # sys.setrecursionlimit(node_count)
    nodes_geo = np.empty((node_count, 2), dtype=object)
    nodes_name = np.empty((node_count, 6), dtype=object)
    segments = np.empty((segment_count, 4), dtype=object)
    #lateralsources = np.empty((lateral_count, max_time), dtype=object)
    lateralsources = []
    lateralsources.append([])
    junctions = np.empty((junction_count, 5), dtype=object)
    #qsources = np.empty((qsource_count, max_time), dtype=object)
    qsources = []
    qsources.append([])
    # The boundary condition is hard-coded with the value obtained from the SPRNT input file
    # boundaryconditions = (root_raw, root_bc)

    return


def timeseries_length():
    """
    This function looks at the lateralsources and qsources to determine the maximum timeseries value present.

    The max_time value can be used to set the simulation time in the SWMM.inp
    """
    global max_time
    max_time = 0
    for ii in range(len(sprnt_contents)):
        if sprnt_contents[ii].find("lateralsource") != -1:
            # starting the iteration 2 lines down skips the location
            jj = ii + 2
            while sprnt_contents[jj].find('end') == -1:
                try:
                    splitline = sprnt_contents[jj].split()[0].split('=')
                    #print(splitline[1])
                    if int(splitline[1]) > max_time:
                        #print(int(splitline[1]))
                        max_time = int(splitline[1])
                except:
                    nothing = ii
                jj = jj + 1
        if sprnt_contents[ii].find("qsource") != -1:
            # starting the iteration 2 lines down skips the location
            jj = ii + 2
            while sprnt_contents[jj].find('end') == -1:
                try:
                    splitline = sprnt_contents[jj].split()[0].split('=')
                    #print(splitline[1])
                    if int(splitline[1]) > max_time:
                        #print(int(splitline[1]))
                        max_time = int(splitline[1])
                except:
                    nothing = ii
                jj = jj + 1
    return


def read_sprnt_contents():
    """
    This function populates the various data type arrays with the values parsed from the SPRNT input file database.

    The function loops through each line in the database; when it finds the keyword that indicates a certain
    data-type, the corresponding values are added to the np.array.  This function is sort of messily split between
    keyword-value searches and assumptions of .spt line organization.

    lateralsources and qsources are added differently.  Each time a new time-value pair is found it is appended to
    the 2nd dimension of the 2D list.
    """
    global start_date, start_time

    node_number, segment_number, lateralsource_number, junction_number, qsource_number = 0, 0, 0, 0, 0
    for line in range(len(sprnt_contents) - 1):

        # If the line is commented out, ignore it
        if sprnt_contents[line].find("#") != -1:
            continue

        # When you find a node definition - do this
        # Need to replace the nodes population with the AY-transform
        elif sprnt_contents[line].find("node") != -1:
            #print("Found Node")
            splitline = sprnt_contents[line].split()
            for elem in range(len(splitline)):
                try:
                    keyword = splitline[elem].split('=')
                    if keyword[0] == 'id':
                        nodes_name[node_number, node_id] = keyword[1]
                    if keyword[0] == 'sR':
                        nodes_name[node_number, node_sR] = keyword[1]
                    if keyword[0] == 'n':
                        nodes_name[node_number, node_n] = keyword[1]
                    if keyword[0] == 'zR':
                        nodes_name[node_number, node_zR] = keyword[1]
                    if keyword[0] == 'hR':
                        nodes_name[node_number, node_hR] = keyword[1]
                except:
                    continue
            node_number = node_number + 1
            continue

        # When a segment definition is found - do this
        elif sprnt_contents[line].find("segment") != -1:
            #print("Found Segment")
            spline = re.split(' up=| down=| length=| ', sprnt_contents[line])
            segments[segment_number, 0:3] = spline[2:5]
            segment_number = segment_number + 1
            continue

        # When a lateralsource is defined - add the time series
        elif sprnt_contents[line].find("lateralsource") != -1:
            #print("Found Lateral Source")
            if lateralsource_number > 0:
                lateralsources.append([])
            spline = re.split("=|\n", sprnt_contents[line + 1])
            lateralsources[lateralsource_number].append(spline[1])
            nextspline = re.split("=|\n", sprnt_contents[line + 3])
            lateralsources[lateralsource_number].append(nextspline[1])
            temp_line = line + 4
            # Each time series has the source, timeunit, time_number series.  So the time_number series starts at 2
            while sprnt_contents[temp_line].find('t=') != -1:
                timespline = re.split("=|\n| ", sprnt_contents[temp_line])
                for ii in range(len(timespline)):
                    if timespline[ii] == 't':
                        lateralsources[lateralsource_number].append(timespline[ii+1])
                    if timespline[ii] == 'v':
                        lateralsources[lateralsource_number].append(timespline[ii+1])
                temp_line = temp_line + 1
            lateralsource_number = lateralsource_number + 1
            continue


        # When a junction is found - add the connecting nodes and coefficients
        elif sprnt_contents[line].find("junction") != -1:
            #print("Found Junction")
            spline = re.split("down=| up1=| coeff1=| up2=| coeff2=|\n", sprnt_contents[line + 1])
            spline = spline[1:-1]
            for ii in range(len(spline)):
                if spline[ii].find(',') != -1:
                    spline[ii] = spline[ii].strip(',')
            try:
                junctions[junction_number] = spline
            except ValueError:
                junctions[junction_number, junction_down:junction_coeff1 + 1] = spline
            junction_number = junction_number + 1
            continue

        # When a qsource is found - add the time series
        elif sprnt_contents[line].find("qsource") != -1:
            #print("Found Qsource")
            if qsource_number > 0:
                qsources.append([])
            spline = re.split("=|\n", sprnt_contents[line + 1])
            qsources[-1].append(spline[1])
            nextspline = re.split("=|\n", sprnt_contents[line + 3])
            qsources[-1].append(nextspline[1])
            temp_line = line + 4
            # Similarly to lateral sources, the qsource row has source, timeunit, time_number series data
            while sprnt_contents[temp_line].find('t=') != -1:
                timespline = re.split("=|\n| ", sprnt_contents[temp_line])
                #qsources[qsource_number, time_number] = float(timespline[-2])
                for ii in range(len(timespline)):
                    if timespline[ii] == 't':
                        qsources[-1].append(timespline[ii+1])
                    if timespline[ii] == 'v':
                        qsources[-1].append(timespline[ii+1])
                temp_line = temp_line + 1

            qsource_number = qsource_number + 1
            continue

        sum_objects = (node_number + segment_number + lateralsource_number + junction_number + qsource_number)
        if (sum_objects) % 10000 == 0:
            print(sum_objects, "found in read_sprnt_contents")
    return


def read_IC_H_contents():
    with open(IC_H_file, 'r') as H_file:
        next(H_file)
        data = H_file.readlines()
    for line in range(len(data)):
        templine = data[line].split()
        ic_loc = np.where(nodes_name[:, node_id] == templine[0])
        # nodes_name[ic_loc[0][0], node_H_IC] = templine[1]
        nodes_name[ic_loc[0][0], node_H_IC] = 0.0
    return


def read_IC_Q_contents():
    with open(IC_Q_file, 'r') as Q_file:
        next(Q_file)
        data = Q_file.readlines()
    for line in range(len(data)):
        templine = data[line].split()
        if templine[0].find('J') != -1:
            # print("Gotta J", templine[0])
            continue
        ic_loc = np.where(segments[:, segment_up] == templine[0])
        try:
            # print(templine[0], ic_loc[0][0])
            # segments[ic_loc[0][0], segment_Q_IC] = templine[1]
            segments[ic_loc[0][0], segment_Q_IC] = 0.0
        except:
            print(templine[0])
    return


def junction_collapsing():
    """
    This function relates to the SPRNT data structure of multiple nodes existing in the same physical location.  In
    SPRNT these nodes are connected relationally using a junction object.  However, in SWMM a junction and a node are
    synonymous, and only one node/junction can exist in one point in space.

    The junction's downstream node (i.e. the upstream node for the segment LEAVING the junction) is the node that is
    chosen to persist.  The junction's upstream nodes (i.e. the downstream nodes for the segments ENTERING the junction)
    are removed and the downstream nodes for the entering segments are reassigned.
    """
    swmm_junctions = []
    for segment in range(len(segments)):
        from_node_id = segments[segment, segment_up].split('_')[-1]
        to_node_id = segments[segment, segment_down].split('_')[-1]
        if len(from_node_id) > 6:
            segments[segment, segment_up] = from_node_id
        if len(to_node_id) > 6:
            segments[segment, segment_down] = to_node_id
        if to_node_id == 'J':
            for junctionnode in range(len(junctions)):
                if (junctions[junctionnode, junction_up1] == segments[segment, segment_down]) or \
                        (junctions[junctionnode, junction_up2] == segments[segment, segment_down]):
                    segments[segment, segment_down] = junctions[junctionnode, junction_down]
                    swmm_junctions.append(segments[segment,segment_down])
    return swmm_junctions


def junction_adjacent(swmm_junctions):
    junction_adjacent_nodes = []
    for ss in range(len(segments)):
        for jj in range(len(swmm_junctions)):
            if segments[ss, segment_up] == swmm_junctions[jj]:
                junction_adjacent_nodes.append(segments[ss, segment_down])
    return junction_adjacent_nodes


def COMID_Link_Map():
    global node_id_list

    nodes_idx = 1
    comid_link_map = np.empty((len(node_id_list), 2), dtype=object)

    for ii in range(len(node_id_list)):
        if node_id_list[ii].find('_') == -1:
            continue
        elif node_id_list[ii] == '-998877':
            continue
        else:
            spline = re.split("_", node_id_list[ii, 0])
            comid = spline[0]
            comid_link_map[ii] = comid, str(ii+1)
    with open('comid_node_map.pickle', 'wb') as f:
        pickle.dump(comid_link_map, f)
    quit()
    return


def initialization():
    global root, root_bc
    global_definitions()

    read_sprnt_for_counters()
    create_empty_arrays_time = time.time() - start
    print("Created empty arrays/lists for nodes, segments, junctions, lateralsources, qsources:",
          create_empty_arrays_time)

    read_sprnt_contents()
    arrays_populated_time = time.time() - start
    print("The arrays have been populated:", arrays_populated_time)

    # read_IC_H_contents()
    # depth_IC_time = time.time() - start
    # print("The depth initial condition has been read:", depth_IC_time)

    # read_IC_Q_contents()
    # flow_IC_time = time.time() - start
    # print("The flowrate initial condition has been read:", flow_IC_time)

    # COMID_Link_Map()

    timeseries_length()
    read_timeseries = time.time() - start
    print("Found max timeseries length time: ", read_timeseries)

    swmm_junctions = junction_collapsing()
    junction_collapsed_time = time.time() - start
    print("The junctions have been collapsed:", junction_collapsed_time)

    # junction_adjacent_nodes = junction_adjacent(swmm_junctions)
    # print(junction_adjacent_nodes)

    root, root_bc = find_boundarycondition()
    find_root_time = time.time() - start
    print(root, "found as root in:", find_root_time)

    return


def create_graph():
    """
    This function creates a NetworkX object whose edges are the segments of the TRB.  Calling this function after
    the junction_collapse() function ensures that the entire TRB is contained within a contiguous graph.

    NetworkX is strange, the segments need to be added backwards (To_node, From_node) in order for the
    breadth-first-search function to work.
    :return: G = Graph of TRB
    """
    G = nx.DiGraph()

    edges_list = []
    for segment in segments:
        # Need to be added in terms of (To_node, From_node) because of the inexplicable way that bfs_predecessor works
        edges_list.append((segment[segment_down], segment[segment_up]))

    # print(edges_list)
    G.add_edges_from(edges_list)
    return G


def node_search_and_sort(G, root):
    """
    This function contains the breadth-first-search routine to order the TRB nodes in such a way as to easily
    facilitate the elevation calculation.
    :param G: Directed Graph of the TRB
    :param root: The outfall of the TRB
    :return ordered_nodes_list: A list of the nodes of the TRB such that each node is downstream of the ones after it
    """
    bfs_output = list(nx.bfs_predecessors(G, source=root))
    ordered_nodes_list = [root]
    for nn in bfs_output:
        ordered_nodes_list.append(nn[0])
    return ordered_nodes_list


def height_calc(slope, length):
    """
    Calculates the elevation difference between two points based on a trigonomic expression of the diagonal length
    and slope between the two points.
    :param slope: The slope of the triangle's hypotenuse, determined by the listed slope at the upstream node
    :param length: The diagonal length of the triangle's hypotenuse, determined by the length of the segment between
                    the two points
    :return height: The vertical distance between the two points
    """
    angle = math.atan(slope)
    height = length * math.sin(angle)
    return height


def calc_elevation(ordered_nodes_list, root):
    """
    This function uses a sorted list to calculate the elevation at every node in the TRB.
    A for-loop is used through the nodes, the location in the nodes_array, the location of the segment that has that
    node as an upstream boundary, and the downstream node (previous node) of that segment are all determined.
    The height_calc function is then used to determine the increase in elevation for the current node above the previous
    node.
    **The sorted list is organized in a way to guarantee that every time an elevation is needed, it has already been
    calculated by a previous iteration.
    :param ordered_nodes_list: The sorted list of nodes that was returned by the node_search_and_sort() func
    :param root: The outfall of the system - also the "initial condition" for the elevation calculation

    The end result of this function is a nodes_name array that has useful values in the nodes_name[:, node_zR] col
    """
    elevation_counter = 0
    nodes_name[:, node_zR] = 0.0
    for node in ordered_nodes_list:
        # The first node in the ordered list should be the root, so just skip that
        if node == root:
            print(node, root, "found")
            continue
        # Find the index location where the current node exists in the nodes_array
        node_loc = np.where(nodes_name[:, node_id] == node)
        # Find the index location where the current node exists as the upstream boundary in the segments_array
            # Guaranteed to be a unique location for a tree
        segment_loc = np.where(segments[:, segment_up] == node)
        # Find the downstream node of the same segment. This is used to calculate the base elevation for the height_calc
        try:
            previous_node = segments[segment_loc[0][0], segment_down]
        except:
            print(segment_loc)
            quit()
        # Find the index location where the previous node exists in the nodes_array
        previous_node_loc = np.where(nodes_name[:, node_id] == previous_node)
        # The slope is the slope at the upstream node
        slope = float(nodes_name[node_loc[0][0], node_sR])
        # The length is the segment_length between the current node and the previous_node
        length = float(segments[segment_loc[0][0], segment_length])
        # Calculate the height difference between the current node and the previous_node
        height_diff = height_calc(slope, length)
        # Grab the previous_node height
        previous_node_height = nodes_name[previous_node_loc[0][0], node_zR]
        # The elevation of the new node is the sum of the previous_node_height and the height difference
        nodes_name[node_loc[0][0], node_zR] = previous_node_height + height_diff

        elevation_counter = elevation_counter + 1
        if elevation_counter % 1000 == 0:
            print("Calculated the elevation of", elevation_counter, "nodes of", str(len(ordered_nodes_list)))
    return


def node_nullifying(root):
    """
    This array is part of a filtering process that removes the nodes that exist in the SPRNT system but not the SWMM
    system.  Also, the node that is a boundary condition in SRPNT gets a different "OUTFALL" designation in SWMM, so it
    also needs to be filtered.
    """
    for node in range(len(nodes_name)):
        if node % 1000 == 0:
            print("Nullifying Node", node)
        # if (nodes_name[node, node_id] not in segments[:, segment_up]):
        #     nodes_name[node, 0] = '-998877'
        suffix = nodes_name[node, node_id].split('_')[-1]
        if (nodes_name[node, node_id] == root) or (suffix == 'J'):
            nodes_name[node, node_id] = '-998877'
            continue

        if (nodes_name[node, node_zR] == None) or (float(nodes_name[node, node_zR]) == 0.0):
            print(nodes_name[node, :], " not given an elevation")
    return


def elevation_main_func():
    global nodes_list, segment_count
    """
    This function compiles and profiles all of the above functions.

    The scope of the elevation_main() portion of the SPRNTtoSWMM workflow is from the beginning opening of the SPRNT
    .spt file, through calculating the elevations of the nodes, and finally checking if there are any weird nodes
    that evaded elevation calculation.
    """

    TRB_Graph = create_graph()
    graph_time = time.time() - start
    print("The graph has been constructed: ", graph_time)
    # print(TRB_Graph.size())
    # print(segment_count)
    # TRB_Graph_Two = TRB_Graph.to_undirected()
    # print(nx.is_connected(TRB_Graph_Two))

    nodes_list = node_search_and_sort(TRB_Graph, root)
    bfs_time = time.time() - start
    print(nodes_list)
    print("The nodes search and ordering took:", bfs_time)

    calc_elevation(nodes_list, root)
    elevation_calc_time = time.time() - start
    print(nodes_name[:, node_zR])
    print("Calculating the elevation took:", elevation_calc_time)

    node_nullifying(root)
    node_nullifying_time = time.time() - start
    print("The _J nodes and the root have been nullified:", node_nullifying_time)

    for nn in range(len(nodes_name)):
        if (nodes_name[nn, node_zR] < 0.1) and (nodes_name[nn, node_id] != '-998877'):
            print(nodes_name[nn])
    # quit()
    return


def ay_hec2_transform():
    """
    This function contains the mathematical conversion from AYWP intrinsic descriptions of irregular channels to the
    HEC2 standard of Station-Elevation pairs.

    This function is needed when the HAND-derived Area-Depth-TopWidth-WettedPerimeter relationships are used to describe
    the geometries of the nodes in SPRNT.  SWMM uses the HEC2 standard of station-elevation pairs when describing
    irregular channel geometries.
    """
    global station, elev
    station = [0.0]
    elev = []
    for yy in Y:
        elev.append(yy)
    elev.append(0.0)
    for yy in range(1, len(Y)):
        elev.append(Y[-1 - yy])
    for nn in range(1, len(W)):
        offset = (W[nn - 1] - W[nn]) / 2
        station.append(station[nn - 1] + offset)
    reversed_W = list(reversed(W))
    station.append(station[-1] + W[-1])
    for ww in range(1, len(reversed_W)):
        offset = (reversed_W[ww] - reversed_W[ww - 1]) / 2
        station.append(station[-1] + offset)

    for stat in range(1, len(station)):
        if station[stat] < station[stat - 1]:
            print("Station needs to be monotonically increasing")
            print(station[stat], station[stat - 1])
            quit()

    return elev, station


def transform_check(elev, station, A_max, P_max):
    """
    This function is used to ensure that during the AYWP transform to the HEC2 standard that the hydrologically relevant
    properties like area and wetted perimeter have been preserved.
    """
    global minimum_neg, node_counter
    station_cap = [station[0], station[-1]]
    elev_cap = [elev[0], elev[-1]]
    area_channel = np.trapz(elev, station)
    area_cap = np.trapz(elev_cap, station_cap)
    area_difference = abs((A_max - area_cap + area_channel)) / A_max
    if abs(area_difference) > cross_section_transformation_tolerance:
        print("The overall area has not been conserved by ", area_difference, ' percent')
        quit()

    wetted_perimeter = 0
    for ii in range(len(station) - 1):
        x_diff = abs(station[ii + 1] - station[ii])
        y_diff = abs(elev[ii + 1] - elev[ii])
        inc_length = np.sqrt(x_diff ** 2 + y_diff ** 2)
        wetted_perimeter = wetted_perimeter + inc_length

    perimeter_diff = abs((wetted_perimeter - P_max)) / P_max
    if abs(perimeter_diff) > cross_section_transformation_tolerance:
        print("The wetted perimeter has not been conserved by ", perimeter_diff, ' percent')
        quit()

    return


def ay_main():
    """
    This function contains the main workflow of the AY-HEC2 conversion.
    """
    global total_count, A, P, Y, W
    node_counter = 0

    for ii in range(total_count):
        if sprnt_contents[ii].find("intrinsic") != -1:
            A = []
            P = []
            Y = []
            W = []
            jj = ii + 1
            while sprnt_contents[jj].find("end") == -1:
                line = sprnt_contents[jj].split()
                for elem in range(len(line)):
                    line[elem] = line[elem].split('=')[1]
                # print(line)
                A.insert(0, float('{:.3f}'.format(float(line[0]))))
                P.insert(0, float('{:.3f}'.format(float(line[1]))))
                Y.insert(0, float('{:.3f}'.format(float(line[2]))))
                W.insert(0, float('{:.3f}'.format(float(line[3]))))
                jj = jj + 1

            for ww in range(1, len(W)):
                if W[ww - 1] < W[ww]:
                    W[ww] = W[ww - 1]

            nodes_geo[node_counter] = ay_hec2_transform()
            transform_check(nodes_geo[node_counter, 0], nodes_geo[node_counter, 1], A[0], P[0])
            node_counter = node_counter + 1

    return


def trapezoidal_main():
    """
    This function contains the main workflow of the trapezoidal geometric handling
    """
    global total_count, bottomwidth, slope
    node_counter = 0

    for ii in range(total_count):
        if sprnt_contents[ii].find("trapezoidal") != -1:
            line = sprnt_contents[ii].split()
            for elem in range(len(line)):
                try:
                    keyword = line[elem].split('=')
                    if keyword[0] == 'bottomwidth':
                        bottomwidth = float('{:.3f}'.format(float(keyword[1])))
                    if keyword[0] == 'slope':
                        slope = float('{:.3f}'.format(float(keyword[1])))
                except:
                    continue

            nodes_geo[node_counter] = bottomwidth, slope
            #print(nodes_geo[node_counter])
            node_counter = node_counter + 1

    return


def geometry_main():
    """
    This function compiles and profiles the handling of the geometry

    NEED TO ADD A TAG FOR DOING EITHER TRAPEZOIDAL OR AY
    """
    trapezoidal_main()
    #ay_main()
    geometry_time = time.time() - start
    print("The geometries have been found:", geometry_time)
    return


def read_SWMM_contents(inputfilename):
    """
    This function reads the SWMM template as a list of lines

    The necessary data will be inserted into the correct place in this list, and then the list will be re-written
    into another SWMM file.
    """
    global swmm_contents, count

    with open(inputfilename, 'r') as swmmput:
        swmm_contents = swmmput.readlines()

    count = len(swmm_contents)
    return


def calc_end_time():
    global end_date_index, end_time_index, end_date, end_time
    for ii in range(count):
        if swmm_contents[ii].find("START_DATE") != -1:
            start_date_index = ii
            start_date = swmm_contents[ii].split()[-1]
            continue
        if swmm_contents[ii].find("START_TIME") != -1:
            start_time_index = ii
            start_time = swmm_contents[ii].split()[-1]
            continue
        if swmm_contents[ii].find("REPORT_START_DATE") != -1:
            report_start_date_index = ii
            continue
        if swmm_contents[ii].find("REPORT_START_TIME") != -1:
            report_start_time_index = ii
            continue
        if swmm_contents[ii].find("END_DATE") != -1:
            end_date_index = ii
            continue
        if swmm_contents[ii].find("END_TIME") != -1:
            end_time_index = ii
            continue
        if swmm_contents[ii].find("[EVAPORATION]") != -1:
            break

    time_added = datetime.timedelta(hours=max_time)

    templine_date = start_date.split("/")
    templine_time = start_time.split(":")
    start_month, start_day, start_year = templine_date
    start_hours, start_minutes, start_seconds = templine_time
    dd = datetime.datetime(int(start_year), int(start_month), int(start_day), int(start_hours), int(start_minutes),
                           int(start_seconds))
    print(dd)
    end_dd = dd + time_added
    print(end_dd)
    end_year = end_dd.strftime("%Y")
    end_month = end_dd.strftime("%m")
    end_day = end_dd.strftime("%d")
    end_hour = end_dd.strftime("%H")
    end_minute = end_dd.strftime("%M")
    end_second = end_dd.strftime("%S")

    end_date = str(end_month) + "/" + str(end_day) + "/" + str(end_year)
    end_time = str(end_hour) + ":" + str(end_minute) + ":" + str(end_second)

    print(end_date)
    print(end_time)
    return


def hardcode_timeseries_test(series_names, index):
    for name in series_names:
        templine = [name, "0", "2.0", "3500", "10.0"]
        splitline = "      ".join(templine) + "\n"
        swmm_contents.insert(index + 3, splitline)
    return


def hardcode_inflows_test(series_names, index):
    for name in series_names:
        templine = [name, "FLOW", name, "FLOW", "1.0"]
        splitline = "      ".join(templine) + "\n"
        swmm_contents.insert(index + 3, splitline)
    return


def create_inputfile(trialfile, root):
    """
    This function combines the processed data from the np.arrays (organized into the data types from SPRNT) and inserts
    the data line by line into the SWMM inputfile template

    When a header is found in the swmm_contents, the lines beneath that header are populated with the appropriate
    information from the np.arrays. All of the tokens for which no SPRNT corollary exists are be populated with a
    default value.

    **Where ambiguous, values for conduit/xsection dimensions are defined by the values at the upstream SPRNT node.
    """
    global lateralseries_names, qseries_names, nodes_visited, qsource_count, node_id_list

    # # Add the end_date and end_time
    # swmm_contents[end_date_index] = "END_DATE       " + end_date + '\n'
    # swmm_contents[end_time_index] = "END_TIME       " + end_time + '\n'

    # node_id_list = []
    # # SWMM Junctions = SPRNT Nodes
    # junction_index = swmm_contents.index('[JUNCTIONS]\n')
    # print("Creating SWMM Junctions", time.time() - start)
    # # Iterates backwards so they show up in a familiar order in the .inp file
    # for node in nodes_name[::-1]:
    #     # This logical checks the node name, to ensure each node is only added once (i.e. the comID endpoints)
    #     if (node[node_id] not in node_id_list) and (node[node_id] not in root) and (node[node_id] != '-998877'):
    #         try:
    #             # templine = [node[node_id], '{:.3f}'.format(float(node[node_zR])), '0', '{:.3f}'.format(float(node[node_H_IC])), '0', '0']
    #             templine = [node[node_id], '{:.3f}'.format(float(node[node_zR])), '0', '0', '0', '0']
    #         except TypeError:
    #             print(node, "triggered a TypeError")
    #             continue
    #         splitline = "      ".join(templine) + "\n"
    #         swmm_contents.insert(junction_index+3, splitline)
    #     node_id_list.append(node[0])

    #     """
    #     COORDINATES
    #     SWMM requires coordinates for each node bc it functions through a GUI.  SPRNT, as an api, doesn't have the same
    #     requirement, so coordinate data does not exist for the Texas river basins
        
    #     -for each node, randomly select values between 0-1000 for x and y
        
    #     **The index of each keyword must be updated before insertion to reflect the size change due to the previous
    #     section's insertions
    #     """
    #     coordinates_index = swmm_contents.index('[COORDINATES]\n')
    #     x_coor = random.randrange(0, 1000)
    #     y_coor = random.randrange(0, 1000)
    #     templine = [node[0], str(x_coor), str(y_coor)]
    #     splitline = "      ".join(templine) + "\n"
    #     swmm_contents.insert(coordinates_index + 3, splitline)

    # # print(node_id_list)
    # # COMID_Link_Map()
    # # quit()
    # """
    # OUTFALL
    # - name is root
    # - elevation is 0.0
    # - type is FIXED
    # - stage data is root_bc
    # - gated is NO
    # """
    # outfall_index = swmm_contents.index('[OUTFALLS]\n')
    # templine = [root, '0.0', 'FIXED', str(root_bc), 'NO']
    # splitline = "      ".join(templine) + "\n"
    # swmm_contents.insert(outfall_index + 3, splitline)

    # """
    # CONDUITS
    # -The name of the conduit can be from an arbitrary counter
    # -The from_node, to_node, and length are from the segments array
    # -The roughness will be assigned by the roughness in the nodes array with the same id as the segment's upstream node
    # """
    # conduit_name = 1
    # conduits_index = swmm_contents.index('[CONDUITS]\n')
    # print("Creating SWMM Conduits", time.time() - start)

    # for segment in segments[::-1]:
    #     from_node = segment[segment_up]
    #     to_node = segment[segment_down]
    #     length = segment[segment_length]
    #     ic_flow = segment[segment_Q_IC]
    #     # This variable designates the row of the nodes np.array where the from_node can be found
    #     node_loc = np.where(nodes_name[:, node_id] == from_node)
    #     mannings_n = nodes_name[node_loc[0][0], node_n]
    #     # templine = [str(conduit_name), from_node, to_node, length, str(mannings_n), '0', '0', '{:.3f}'.format(float(ic_flow)), '']
    #     templine = [str(conduit_name), from_node, to_node, length, str(mannings_n), '0', '0', '0', '']
    #     splitline = "      ".join(templine) + "\n"
    #     swmm_contents.insert(conduits_index + 3, splitline)

    #     """    
    #     XSECTIONS
    #     -The name of the XSECTION will be the same as the arbitrary counter for conduits.  Meaning they're populated at
    #     the same time
        
    #     For Lavaca, the shape is the transformed HEC2 irregular, and the Tsect name needs to be unique to the conduit
    #     - Tsect name will be the upstream_node
    #     # -The shape is TRAPEZOIDAL
    #     # -Geom1 is the max depth, also arbitrarily defined (probably 100)
    #     # -Geom2 is the bottom width, which is defined at the upstream node for each segment
    #     # -Geom3 is the side slope, also defined at the upstream node for each segment
    #     """
    #     xsections_index = swmm_contents.index('[XSECTIONS]\n')

    #     # Use this templine for when the conduits are irregularly defined
    #     #templine = [str(conduit_name), 'IRREGULAR', from_node]

    #     # Use this templine for when the conduits are trapezoidal
    #     bottom_width = nodes_geo[node_loc[0][0], geo_bottomwidth]
    #     side_slope = nodes_geo[node_loc[0][0], geo_slope]
    #     templine = [str(conduit_name), 'TRAPEZOIDAL', '100', str(bottom_width), str(side_slope), str(side_slope)]
    #     # templine = [str(conduit_name), 'RECT_OPEN', '100', '20', '0', '0', '1']
    #     splitline = "      ".join(templine) + "\n"
    #     swmm_contents.insert(xsections_index + 3, splitline)

    #     conduit_name = conduit_name + 1

    """
    INFLOWS
    -The node name will be the node where the qsource/lateralsource is located
    -Constituent is 'FLOW'
    -Time Series is variable bc there are multiple time series
    -Type 'FLOW'
    -Mfactor, Sfactor = 1.0
    """
    inflows_index = swmm_contents.index('[INFLOWS]\n')
    print("Creating SWMM Inflows", time.time() - start)

    qseries_names = []
    for source in qsources:
        try:
            name = source[0]
        except:
            print("Something messed up with qseries_names")
            print(qsources)

        qseries_names.append(name)

    # hardcode_inflows_test(qseries_names, inflows_index)
        templine = [name, 'FLOW', name, 'FLOW', '1.0', '1.0', '', '']
        splitline = "      ".join(templine) + "\n"
        swmm_contents.insert(inflows_index + 3, splitline)

    # print(len(qseries_names))
    # print(qsource_count)

    # quit()

    lateralseries_names = []
    for source in lateralsources:
        name = source[0]

        templine = [name, 'FLOW', name, 'FLOW', '1.0', '1.0', '', '']
        lateralseries_names.append(name)
        splitline = "      ".join(templine) + "\n"
        swmm_contents.insert(inflows_index + 3, splitline)

    # hardcode_inflows_test(lateralseries_names, inflows_index)

    """
    [TIMESERIES]
    - need to include flow timeseries from lateral sources and qsources
    - we have checked that there is no overlap
    - each timeseries is named after the node at which it occurs
    - on a given line: the name, (time since simulation began, flow value) pairs
    - The lateral sources make the size of the input file explode
    - We also include the T1 series which describes a subset of qsource and lateralsource data
    """
    timeseries_index = swmm_contents.index('[TIMESERIES]\n')
    print("Creating SWMM Timeseries", time.time() - start)


    # hardcode_timeseries_test(qseries_names, timeseries_index)
    # hardcode_timeseries_test(lateralseries_names, timeseries_index)


    # First I should check if the lateralsources and qsources have any overlap

    for name in lateralsources[:]:
        lateralsource_name = name[0]
        if lateralsource_name in qsources[:][:]:
            print("There are multiple inflows to node " + name)
            quit()

    # This is where I need to put a filter for the max_inflows
    for name in lateralseries_names:
        file = '"' + str(name) + '.dat"'
        completeName = os.path.join(path_to_external_files, file)
        templine = [name, 'FILE', file]
        splitline = "      ".join(templine) + "\n"
        swmm_contents.insert(timeseries_index + 3, splitline)

    for name in qseries_names:
        file = '"' + str(name) + '.dat"'
        completeName = os.path.join(path_to_external_files, file)
        templine = [name, 'FILE', file]
        splitline = "      ".join(templine) + "\n"
        swmm_contents.insert(timeseries_index + 3, splitline)

    """
    [TRANSECTS] - this section is needed when the geometry is irregular.  Need to functionalize it an tag it so that
    it's triggered when the "intrinsic" geometry is found. 
    
    """
    # transect_index = swmm_contents.index('[TRANSECTS]\n')
    #
    # error_count = 0
    # for segment in segments[::-1]:
    #     from_node = segment[0]
    #     node_loc = np.where(nodes_name[:, 0] == from_node)
    #     mannings_n = nodes_name[node_loc[0][0], 2]
    #     templine_NC = ['NC', mannings_n, mannings_n, mannings_n]
    #     station = nodes_geo[node_loc[0][0], 1]
    #     elevation = nodes_geo[node_loc[0][0], 0]
    #     stat_past = 0
    #     for stat in station:
    #         if stat >= stat_past:
    #             stat_past = stat
    #         else:
    #             error_count = error_count + 1
    #             break
    #     templine_X1 = ['X1', from_node, str(len(station)), str(station[0]), str(station[-1]), '0.0', '0.0', '0.0', '0.0', '0.0']
    #     templine_GR = ['GR']
    #     for index in range(len(station)):
    #         templine_GR.append(str(elevation[index]))
    #         templine_GR.append(str(station[index]))
    #     splitline_NC = "      ".join(templine_NC) + "\n"
    #     splitline_X1 = "      ".join(templine_X1) + "\n"
    #     splitline_GR = "      ".join(templine_GR) + "\n"
    #     swmm_contents.insert(transect_index + 3, ';\n')
    #     swmm_contents.insert(transect_index + 3, splitline_GR)
    #     swmm_contents.insert(transect_index + 3, splitline_X1)
    #     swmm_contents.insert(transect_index + 3, splitline_NC)
    # print(error_count)

    # Once the SWMM .inp template has been populated, open a new file and write the contents to it
    print("Opening SWMM file and writting data", time.time() - start)
    with open(trialfile, 'w') as newfile:

        for i in range(len(swmm_contents)):
            newfile.write(swmm_contents[i])

    return


def create_external_lateralfiles(extfile, toFile):
    file = str(extfile) + ".dat"
    completeName = os.path.join(path_to_external_files, file)
    with open(completeName, "w") as file1:
        for i in range(len(toFile)):
            file1.write(toFile[i])
        file1.write("\n")
        # lastval = toFile[-1].split()[-1]
        lastval = 0.0
        # print(file, lastval)
        file1.write("     " + str(max_time) + "   " + str(lastval))
    return


def create_external_qsourcefiles(extfile, toFile):
    file = str(extfile) + ".dat"
    completeName = os.path.join(path_to_external_files, file)
    with open(completeName, "w") as file1:
        for i in range(len(toFile)):
            file1.write(toFile[i])
        file1.write("\n")
        # lastval = toFile[-1].split()[-1]
        lastval = 0.1
        # print(file, lastval)
        file1.write("     " + str(max_time) + "   " + str(lastval))
    return


def swmm_create_main():

    read_SWMM_contents(swmm_example)

    calc_end_time()

    # The function argument specifies the name of the SWMM input file that is created
    create_inputfile(swmm_final, root)

    # print(qseries_names)

    for name in lateralseries_names:
        toFile = []
        templine = []
        timeline = []
        valueline = []
        for ii in range(len(lateralsources)):
            row = lateralsources[ii]
            if (row[0] == name):
                if (ii % 100 == 0):
                    print(ii, " of ", len(lateralsources), " lateralsource files written")
                for time in lateralsources[ii][2:-2:2]:
                    timeline.append(str(time))
                for value in lateralsources[ii][3:-1:2]:
                    valueline.append(str(value))

                for timevalue in range(len(timeline)):
                    templine.append(timeline[timevalue])
                    templine.append(valueline[timevalue])
                    templine.append("\n")

                splitline = "      ".join(templine) + "\n"
                splitline = splitline.strip()
                toFile.append(splitline)
        # print(toFile[-1::])
        create_external_lateralfiles(name, toFile)

    for name in qseries_names:
        toFile = []
        templine = []
        timeline = []
        valueline = []
        for ii in range(len(qsources)):
            row = qsources[ii]
            if (row[0] == name):
                if (ii % 100 == 0):
                    print(ii, " of ", len(qsources), " qsource files written")
                for time in qsources[ii][2:-2:2]:
                    timeline.append(str(time))
                for value in qsources[ii][3:-1:2]:
                    valueline.append(str(value))

                for timevalue in range(len(timeline)):
                    templine.append(timeline[timevalue])
                    templine.append(valueline[timevalue])
                    templine.append("\n")

                splitline = "      ".join(templine) + "\n"
                splitline = splitline.strip("   ")
                toFile.append(splitline)
        create_external_qsourcefiles(name, toFile)
    return


def SPRNTtoSWMM_main():
    """
    Main workflow.  Next thing I need to do is have each of these functions take the TRB as arguments and then
    repeat this process for each of the TRBs
    """
    settings()
    read_sprnt_to_contents()
    initialization()
    print("Initialization is over")
    # elevation_main_func()
    # geometry_main()
    swmm_create_main()
    total_time = time.time() - start
    print("Creating the", swmm_final, "took", total_time)

    return

SPRNTtoSWMM_main()



