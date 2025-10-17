import numpy as np
import scipy.io

def AssembleData():

    ##########################
    ### Load Data from GUI ###
    ##########################

    # Load handoff file
    data = scipy.io.loadmat('../Data/ModelHandoff.mat') 
    
    # For debugging purposes (prints intermittently if True)
    printon = False
    # Note that the -1 and many of the +1's in the printon's are to match the
    # indexing of the interface's language, MATLAB

    # Number of nodes in network
    nn = len(data['NetworkData']['X'][0][0][0])
    if printon:
        print('nn = ',nn)

    # Outlet types 'DRT...'
    # Outlet type: (D)erivation box, (R)eceiver, (T)ransmitter, (S)eries outlet, (P)arallel outlet
    OT = 'D'*nn # Initialize every node as a derivation box
    for ind in range(nn):
        if data['NetworkData']['NodeType'][0][0][ind][0] == 2:
            pass                             # If derivation box, do nothing
        elif data['NetworkData']['NodeType'][0][0][ind][0] == 3:
            OT = OT[:ind] + 'R' + OT[ind+1:] # If receiver, replace with 'R'
        elif data['NetworkData']['NodeType'][0][0][ind][0] == 4:
            OT = OT[:ind] + 'T' + OT[ind+1:] # If transmitter, replace with 'T'
        elif data['NetworkData']['NodeType'][0][0][ind][0] == 5:
            OT = OT[:ind] + 'S' + OT[ind+1:] # If series load, replace with 'S'
        elif data['NetworkData']['NodeType'][0][0][ind][0] == 6:
            OT = OT[:ind] + 'P' + OT[ind+1:] # If parallel load, replace with 'P'
        else:
            print('Error: unknown outlet type.')
            exit()
        if printon:
            print('NodeType[',ind+1,'] = ',OT[ind])

    # Load types '0C...0'
    # Load type: -1 = None, 1 = Constant, 2 = Motor, 3 = Series RLC,
    #             4 = Parallel RLC, 5 = Double RLC 6 = Off/Open Circuit, 
    #             7 = Custom
    LT = '0'*nn
    for ind in range(nn):
        if data['NetworkData']['LoadTypenum'][0][0][ind][0] == -1:  # '0'
            pass
        elif data['NetworkData']['LoadTypenum'][0][0][ind][0] == 1: # 'C'
            LT = LT[:ind] + 'C' + LT[ind+1:]                  
        elif data['NetworkData']['LoadTypenum'][0][0][ind][0] == 2: # 'M'
            LT = LT[:ind] + 'M' + LT[ind+1:]
        elif data['NetworkData']['LoadTypenum'][0][0][ind][0] == 3: # 'S'
            LT = LT[:ind] + 'S' + LT[ind+1:]
        elif data['NetworkData']['LoadTypenum'][0][0][ind][0] == 4: # 'P'
            LT = LT[:ind] + 'P' + LT[ind+1:]
        elif data['NetworkData']['LoadTypenum'][0][0][ind][0] == 5: # 'D'
            LT = LT[:ind] + 'D' + LT[ind+1:]
        elif data['NetworkData']['LoadTypenum'][0][0][ind][0] == 6: # 'O'
            LT = LT[:ind] + 'O' + LT[ind+1:]
        elif data['NetworkData']['LoadTypenum'][0][0][ind][0] == 7: # 'X'
            LT = LT[:ind] + 'X' + LT[ind+1:]
        else:
            print('Error: unknown load type.')
            exit()
        if printon:
            print('LoadType[',ind+1,'] = ',LT[ind])

    # Load parameters (LP = [[...],[...],...,[...]])
    # There are a varying number of load parameters depending on the load type,
    # see the documentation for more details
    LP = [[] for ind in range(nn)]
    for ind in range(nn):
        if LT[ind] != '0':
            LP[ind] = list(data['NetworkData']['LoadParam'][0][0][ind][0])
        if printon:
            print('LoadParam[',ind+1,'] = ',LP[ind])

    # Wire lengths
    WireLengths = data['NetworkData']['WireLength'][0][0][0]
    if printon:
        print('Wire lengths = ',WireLengths)

    # Frequencies
    freq = data['NetworkData']['Freq'][0][0][0]
    if printon:
        print('Frequencies = ',freq)

    # Input type voltage/current
    inpVItype = data['NetworkData']['InputVIType'][0][0][0][0]
    if printon:
        print('inpVItype {-1,1} = ',inpVItype)

    # Input value voltage/current
    inpVI = data['NetworkData']['InputVI'][0][0][0]
    if printon:
        print('Input VI = ',inpVI)

    # Type of Rec Model
    whichRec = data['NetworkData']['WhichRec'][0][0][0][0]
    if printon:
        print('Receiver data type = ',whichRec)

    # Rec Model Val
    RecVal = data['NetworkData']['RecVal'][0][0][0]
    if printon:
        print('Receiver data value = ',RecVal)

    # Verify against Spice?
    Spice = data['NetworkData']['Spice'][0][0][0][0]
    if printon:
        print('Spice verification {0,1} = ',Spice)

    # Conductivity of conductors
    sigmavals = data['NetworkData']['Sigma'][0][0]
    if printon:
        print('Sigma = ',sigmavals)

    # Permittivity of insulator
    epsvals = data['NetworkData']['Epsilon'][0][0]
    if printon:
        print('Epsilon = ',epsvals)

    # Number of conductors
    ncond = data['NetworkData']['Ncond'][0][0][0][0]
    if printon:
        print('nconductors = ',ncond)

    # Resolution of points along wire
    wres = data['NetworkData']['WireRes'][0][0][0][0]
    if printon:
        print('wire res = ',wres,' (m)')

    # Wire radii (for each of the conductors)
    wrad = data['NetworkData']['Radius'][0][0]
    if printon:
        print('Wire radii = ',wrad)

    # Separation distance
    dc = data['NetworkData']['DC'][0][0]
    if printon:
        print('DC = ',dc)

    # Hardcode wire formation
    # This was deemed a less impactful detail, yet a burdensome COMSOL conversion
    # factor; the formations used are simply the "equilateral" n-sided shapes.
    if ncond == 4:
        wt = ['square']*nn
    elif ncond == 3:
        wt = ['triangle']*nn
    elif ncond == 2:
        wt = ['line']*nn
    else:
        print('Error: unsupported wire formation detected.')
    if printon:
        print('Wire formation = ',wt[0])

    # Nodes ([[X_1,Y_1,Z_1],[X_2,Y_2,Z_2],...,[X_n,Y_n,Z_n]])    
    node = np.zeros((len(data['NetworkData']['X'][0][0][0]),3))
    for nind in range(len(data['NetworkData']['X'][0][0][0])):
        node[nind] = [data['NetworkData']['X'][0][0][0][nind], data['NetworkData']['Y'][0][0][0][nind], data['NetworkData']['Z'][0][0][0][nind]]
        if printon:
            print('Node[',nind+1,'] = [',node[nind][0],', ',node[nind][1],', ',node[nind][2],']')

    # Wires ([[InitialNodeID_1,FinalNodeID_1],...,[InitialNodeID_n,FinalNodeID_n]])
    wires = data['NetworkData']['Wires'][0][0]-1
    nw = nn-1 # number of wires
    if printon:
        for ind in range(nn-1):
            print('Wire[',ind+1,'] = [',wires[ind][0],', ',wires[ind][1],']')

    ############################################################
    ##### Loading GUI Data Complete... Begin Processing... #####
    ############################################################
    
    # Prepare input type (voltage/current) for usage by code
    if inpVItype == 1:
        UseI = True
    elif inpVItype == -1:
        UseI = False
    else:
        print('Error: Injection type (voltage/current) improperly defined.')
        exit()

    # Identify type of receiver data
    if whichRec == -1:
        whichRec = 'V'
    elif whichRec == 0:
        whichRec = 'I'
    elif whichRec == 1:
        whichRec = 'Z'
    else:
        print('Error: unknown receiver data type.')
        exit()

    # Correct reversed wires (orient wires in direction of power)
    # Method: loads and receiver must be final node, transmitter must be initial
    # node, work backwards and forwards from there
    uncheckedwires = list(range(nw))                # Track unchecked wires
    front = -1                                      # Initialize to non-node ID
    for wirind in range(len(uncheckedwires)):       # Iterate through unchecked wires
        if (OT[wires[wirind][0]] == 'T'):           # If initial node is transmitter
            front = wires[wirind][1]                # Front = corresponding final node
            uncheckedwires.remove(wirind)           # Remove wire from unchecked list
        elif (OT[wires[wirind][1]] == 'T'):         # If final node is transmitter
            tmpwire = wires[wirind][0]              
            wires[wirind][0] = wires[wirind][1]     # Swap wire definition
            wires[wirind][1] = tmpwire
            front = wires[wirind][1]                # Redefine front
            uncheckedwires.remove(wirind)           # Remove
        elif (OT[wires[wirind][0]] == 'R'):         # If initial node is receiver
            tmpwire = wires[wirind][0]
            wires[wirind][0] = wires[wirind][1]     # Swap
            wires[wirind][1] = tmpwire
            uncheckedwires.remove(wirind)           # Remove     
        elif (OT[wires[wirind][1]] == 'R'):         # If final node is receiver
            uncheckedwires.remove(wirind)           # Remove
        elif (LT[wires[wirind][0]] != '0') and (OT[wires[wirind][0]] == 'S'):
            tmpwire = wires[wirind][0]              # If initial node is load
            wires[wirind][0] = wires[wirind][1]     # Swap
            wires[wirind][1] = tmpwire
            uncheckedwires.remove(wirind)           # Remove
        elif (LT[wires[wirind][1]] != '0') and (OT[wires[wirind][1]] == 'S'):
            uncheckedwires.remove(wirind)           # If final node is load, remove

    # Transmitter, receiver, and load wires are removed
    while len(uncheckedwires) > 0:                  # While unchecked wires remain
        countwire = np.zeros((np.max(wires)+1))     # Create list of node involvement
        for wirind in uncheckedwires:               # Iterate through remaining wires
            countwire[wires[wirind][0]] += 1        # Increment node list to track occurrences
            countwire[wires[wirind][1]] += 1           
        for cwind in range(len(countwire)):         # Iterate through node list
            if countwire[cwind] == 1:               # If a node is only present in one remaining wire
                for wind2 in uncheckedwires:        # Iterate through and find the one wire
                    if (wires[wind2][1] == cwind) and (front != cwind):   # If it's final node and =/= front
                        uncheckedwires.remove(wind2)      # Remove
                        break                            
                    elif (wires[wind2][0] == cwind) and (front != cwind): # If it's initial node and =/= front
                        tmpwire = wires[wind2][0] 
                        wires[wind2][0] = wires[wind2][1] # Swap
                        wires[wind2][1] = tmpwire
                        uncheckedwires.remove(wind2)      # Remove
                        break
                    elif (wires[wind2][0] == cwind) and (front == cwind): # If it's initial node and = front
                        uncheckedwires.remove(wind2)      # Remove
                        front = wires[wind2][1]           # Redefine front
                        break
                    elif (wires[wind2][1] == cwind) and (front == cwind): # If it's final node and = front
                        tmpwire = wires[wind2][0]
                        wires[wind2][0] = wires[wind2][1] # Swap
                        wires[wind2][1] = tmpwire
                        front = wires[wind2][1]           # Redefine front
                        uncheckedwires.remove(wind2)      # Remove
                        break

    # Create wire list with node IDs replaced by XYZ nodal locations
    nodes = np.zeros((nw,6))
    for nind in range(nw):
        nodes[nind][range(3)] = node[wires[nind][0]]
        nodes[nind][range(3,6)] = node[wires[nind][1]]


    # Define graph elements (this is the data structure that contains most 
    # of the information used by the code)

    # Define forward and backward paths between nodes
    forward = {}
    backward = {}
    for find in range(nn):                                     # Initialize nn empty slots
        forward.update({find: []})
        backward.update({find: []})
    for find in range(nw):                                     # Fill appropriately
        if len(forward[wires[find][0]]) > 0:                   # If nonempty list
            forward[wires[find][0]].append(wires[find][1])     # Append final node to initial node's forward list
        else:                                                  # If empty list
            forward.update({wires[find][0]: [wires[find][1]]}) # Create forward list for initial node containing final node
        backward.update({wires[find][1]: [wires[find][0]]})    # Define backward list for final node as initial node

    # Define termination and outlet type
    termination = {}
    outlettype = {}
    for otind in range(nn):
        if OT[otind] == 'D':               
            outlettype.update({otind: 'Derivation box'})
            termination.update({otind: False})
        elif OT[otind] == 'R':             
            outlettype.update({otind: 'Receiver'})
            termination.update({otind: True})
            rec = otind                    # Store node ID of receiver
        elif OT[otind] == 'T':             
            outlettype.update({otind: 'Transmitter'})
            termination.update({otind: False})
            trans = otind                  # Store node ID of transmitter
            backward.update({otind: None}) # This finalizes the backwards list as nothing precedes the transmitter
        elif OT[otind] == 'S':             
            outlettype.update({otind: 'Series outlet'})
            termination.update({otind: True})
        elif OT[otind] == 'P':             
            outlettype.update({otind: 'Parallel outlet'})
            termination.update({otind: False})
        else:
            print('Error: outlet type incorrectly defined.')
            exit()

    # Define load type and load parameters
    loadtype = {}
    loadpar = {}
    for ltind in range(nn):
        if LT[ltind] == '0':
            loadtype.update({ltind: None})
            loadpar.update({ltind: None})
        elif LT[ltind] == 'C':
            loadtype.update({ltind: 'Constant'})
            loadpar.update({ltind: list(LP[ltind][0])}) 
        elif LT[ltind] == 'M':
            loadtype.update({ltind: 'Motor'})
            loadpar.update({ltind: list(LP[ltind][0])}) 
        elif LT[ltind] == 'S':
            loadtype.update({ltind: 'Series RLC'})
            loadpar.update({ltind: list(LP[ltind][0])}) 
        elif LT[ltind] == 'P':
            loadtype.update({ltind: 'Parallel RLC'})
            loadpar.update({ltind: list(LP[ltind][0])}) 
        elif LT[ltind] == 'D':
            loadtype.update({ltind: 'Double RLC'})
            loadpar.update({ltind: list(LP[ltind][0])}) 
        elif LT[ltind] == 'O':
            if OT[ltind] == 'R':
                loadtype.update({ltind: None})
            else:
                loadtype.update({ltind: 'Off'})
            loadpar.update({ltind: list(LP[ltind][0])}) 
        elif LT[ltind] == 'X':
            loadtype.update({ltind: 'Custom'})
            loadpar.update({ltind: LP[ltind]}) 
        else:
            print('Error: unexpected load type.')
            exit()

    # Wire radius, length, type, sig, eps, and distance between conductors
    wirerad = {}
    wirelen = {}
    wiretype = {}
    wiredist = {}
    sigmaNew = {}
    epsilonNew = {}
    for wind in range(nn):
        if wind != trans:
            for wind2 in range(nw):
                if wires[wind2][1] == wind:
                    wirerad.update({wind: list(wrad[wind2])})
                    wiretype.update({wind: wt[wind2]})
                    wiredist.update({wind: dc[wind2][0]})
                    wirelen.update({wind: WireLengths[wind2]})
                    #wirelen.update({wind: np.linalg.norm(np.array(nodes[wind2][0:3])-np.array(nodes[wind2][3:6]))}) # Computes physical wire distance
                    sigmaNew.update({wind: list(sigmavals[wind2])})
                    epsilonNew.update({wind: np.double(epsvals[wind2][0])})
        else: # Transmitter has no incoming wire and thus no wire properties
            wirerad.update({wind: None})
            wiretype.update({wind: None})
            wiredist.update({wind: None})
            wirelen.update({wind: None})
            sigmaNew.update({wind: None})
            epsilonNew.update({wind: None})
   
    # Define main path
    path = [rec]
    while trans not in path:
        path.insert(0,backward[path[0]][0])

    # Maintain backward compatibility with earlier model
    M = 1

    # Create graph structure
    graph = {'forward': forward,
             'backward': backward,
             'termination': termination,
             'outlet_type': outlettype,
             'load_type': loadtype,
             'load_parameters': loadpar,
             'wire_radius': wirerad,
             'wire_length': wirelen,
             'wire_type': wiretype,
             'wire_distance': wiredist,
             'transmitter': trans,
             'receiver': rec,
             'path': path,
             'M': M,
             'mean_wire_length': wirelen,
             'Spice': Spice,
             'Sigma': sigmaNew,
             'Epsilon': epsilonNew,
             'Count': 0}

    # Nphase is the number of conductors - 1
    nphase = ncond-1
    
    # Print all results
    if printon:
        print('\n\n',graph,'\n\n')
        exit()

    return graph, freq, nphase, ncond, wres, inpVI, UseI, whichRec, RecVal, node
