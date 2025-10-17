import numpy as np
from scipy.sparse import csc_matrix
from scipy.sparse.linalg import spsolve
from numpy.linalg import inv
import copy
from compVInodes_Functions import compABCDwireto, compLoadImp, compRec, calculate_ABCD_parallel_outlet_load

def compVInodes_Nwire(graph, freq, nwires, UseI, VIin, whichRec, RecVal):
   
    # This code uses the matrix method to solve for the VI at the nodes

    # Set number of wires and phase wires (assume one ground/neutral) 
    nwire = nwires                            # Number of wires
    nphase = int(nwire-1)                     # Number of phase wires
    Cp = 0                                    # Outlet Parasitic capacitance 
    
    # Extracting required information from graph
    mainpath = graph['path']
    forward = graph['forward']
    backward = graph['backward']
    outlet_type = graph['outlet_type']

    # Make definition of backward consistent
    for backind in range(len(backward)):
        if backward[backind] is None:
            backward[backind] = []
    
    # Determine the size of the matrix and vectors
    WireID = np.zeros((len(forward)-1,2))           # Wire IDs for visualization
    WIDC = 0                                        # Count for tracking Wire ID
    WCI = np.zeros((len(forward)-1,1),dtype = int)  # Track which child for indexing
    NumNodes = len(forward)                         # How many nodes in network
    NumVar = 3*nphase*NumNodes                      # Add two unknowns (V/Iin/Iout) for each node, extras will be removed (trans, rec, loads)

    # Finalize expected number of variables/equations for tracking
    for i in range(NumNodes):
            # Assumes node numbering is 0 to n
            NumVar += nphase*len(backward[i])-nphase    # Removes ins current associated with transmitter... could just have -nphase
            NumVar += nphase*len(forward[i])-nphase     # Add additional variables based on splits; load or receiver removes nphase

    VarList = -np.ones((NumVar))              # Will list node associated with V11,V12,V13,I11,I12,I13,I11b,I12b,I13b,...

    # Initialize indexing variables and locate transmitter/receiver
    # Importantly, NodeOrder dictates ordering of solution vector
    Rowcount = 0                              # Track equation number
    NodeOrder = -np.ones((NumNodes))          # Track order of encountered nodes
    for rtind in range(len(mainpath)):
        if outlet_type[mainpath[rtind]] == 'Receiver':
            receiver = mainpath[rtind]
        elif outlet_type[mainpath[rtind]] == 'Transmitter':
            transmitter = mainpath[rtind]

    # Initialize matrix and right hand side vector
    VIMat = np.zeros((NumVar,NumVar),dtype = np.cdouble)  # A in Ax=b, initialized to zeros
    b = np.zeros((NumVar,1),dtype = np.cdouble)
    if len(VIin) != nphase:
        print('Error: for nwire case, we require nwire-1 input voltages or currents.')
    for bind in range(nphase):
        b[bind] = VIin[bind]                              # b in Ax=b, defined as [V/I11in,V/I12in,V/I13in, 0, 0, 0, ..., 0] 

    # Add equations to matrix for initial and final voltages/currents/impedance relations
    if UseI:
        for equind in range(nphase,2*nphase):
            VIMat[Rowcount,equind] = 1                    # Iii = Iii in
            Rowcount += 1
    else:
        for equind in range(nphase):
            VIMat[Rowcount,equind] = 1                    # Vii = Vii in
            Rowcount += 1

    # Add receiver model to matrix depending on user selection
    if whichRec == 'Z':
        Zr = compRec(graph,receiver,freq,nphase,whichRec,RecVal)         # Receiver impedance (nphase x nphase) 
        for equind in range(nphase):                      # Vr = Zr*Ir
            for equind2 in range(nphase):
                VIMat[Rowcount,-nphase+equind2] = -Zr[equind,equind2]
            VIMat[Rowcount,-2*nphase+equind] = 1        
            Rowcount += 1
    elif whichRec == 'V':
        Vr = compRec(graph,receiver,freq,nphase,whichRec,RecVal)         # Receiver voltage (nphase x 1) 
        for equind in range(nphase):                      # Vr = Vr in
            VIMat[Rowcount,-2*nphase+equind] = 1        
            b[Rowcount] = Vr[equind]
            Rowcount += 1
    elif whichRec == 'I':
        Ir = compRec(graph,receiver,freq,nphase,whichRec,RecVal)         # Receiver current (nphase x 1) 
        for equind in range(nphase):                      # Ir = Ir in
            VIMat[Rowcount,-nphase+equind] = 1
            b[Rowcount] = Ir[equind]        
            Rowcount += 1
    else:
        print('Error: RecEq must be Z, V, or I.')
     
    # Add transmitter/receiver variables to list
    for vlc in range(2*nphase): 
        VarList[vlc] = transmitter                        # Transmitter and receiver nodes are the beginning and end of the VarList, respectively
        VarList[-(vlc+1)] = receiver
    
    # Iterate through nodes and add relevant equations
    # Each node has nphase voltages and nphase currents plus nphase * number of forward paths out (if >1)
    # Ordering is: main path last, smaller number first not on main path 
    # While arbitrary, it will be consistent as long as we track NodeOrder    

    nodecount = 0                        # Tracks index of nodes visited
    varlistcount = 2*nphase              # Tracks which node corresponds to which node in ordering V, I1, I2, ...
    curnode = transmitter                # Start at transmitter
    parentnode = copy.copy(curnode)      # Track parent node
    NodeOrder[nodecount] = curnode       # Book keeping, first visited node is transmitter
    
    # Traverse network, adding applicable equations when new nodes are encountered
    while curnode != receiver:
        
        # If the node has yet to be encountered
        if curnode not in NodeOrder:

            # Add to tracked list and increment number of encountered nodes
            nodecount += 1
            NodeOrder[nodecount] = curnode
            
            # Identify which equations to add based on node properties
            # Although some code is redundant, we leave it to explicitly address all cases for readability
            if len(forward[curnode]) == 0:

                if (outlet_type[curnode] == 'Series outlet'):

                    # It's a load that ends the current subpath
                    
                    # Compute ABCD for wire
                    ABCD = compABCDwireto(graph,curnode,nphase,freq)          # Compute ABCD from parent to curnode

                    # Calculate load impedance
                    ZL = compLoadImp(graph,curnode,nphase,freq)               # Compute load impedance associated with curnode

                    # Add impedance relation equations
                    for equind in range(nphase):
                        for equind2 in range(nphase):
                            VIMat[Rowcount,varlistcount+nphase+equind2] = -ZL[equind,equind2]   # Vi = ZLi*Ii
                        VIMat[Rowcount,varlistcount+equind] = 1
                        Rowcount += 1

                    # Add ABCD equations from parent to child
                    for abcdind in range(2*nphase):                           # ABCD Eqs
                        VIMat[Rowcount,varlistcount+abcdind] = 1
                        for abcdind2 in range(nphase):
                            VIMat[Rowcount,np.where(VarList==parentnode)[0][0]+abcdind2] = -ABCD[abcdind,abcdind2] 
                        for abcdind2 in range(nphase,2*nphase):
                            VIMat[Rowcount,np.where(VarList==parentnode)[0][0]+WhichChild+abcdind2-nphase] = -ABCD[abcdind,abcdind2] 
                        Rowcount += 1
                   
                    # Update list of variable ordering
                    for vlc in range(2*nphase):                               # Add new variables for loads
                        VarList[varlistcount] = curnode
                        varlistcount += 1

                else:
            
                    # Nothing should have an empty forward unless it ends at a load
                    print('Error: forward path empty for non-(series)load, non-receiver node.')

            elif len(forward[curnode]) == 1:

                if curnode in mainpath:
 
                    # This node does not branch and is on main path, need only ABCD

                    # Obtain wire ABCD
                    ABCD = compABCDwireto(graph,curnode,nphase,freq)
                   
                    # Add parent to child ABCD equations
                    for abcdind in range(2*nphase):                             # ABCD Eqs
                        VIMat[Rowcount,varlistcount+abcdind] = 1
                        for abcdind2 in range(nphase):
                            VIMat[Rowcount,np.where(VarList==parentnode)[0][0]+abcdind2] = -ABCD[abcdind,abcdind2] 
                        for abcdind2 in range(nphase,2*nphase):
                            VIMat[Rowcount,np.where(VarList==parentnode)[0][0]+WhichChild+abcdind2-nphase] = -ABCD[abcdind,abcdind2] 
                        Rowcount += 1

                    # If not parallel outlet apply only Kirchhoff's
                    if (outlet_type[curnode] != 'Parallel outlet'):
                        for sIind in range(nphase):
                            VIMat[Rowcount,varlistcount+nphase+sIind] = 1                  # Sum Iin = sum Iout
                            VIMat[Rowcount,varlistcount+2*nphase+sIind] = -1               # For one in one out
                            if (outlet_type[forward[curnode][0]] == 'Series outlet') or (outlet_type[forward[curnode][0]] == 'Parallel outlet'):
                                VIMat[Rowcount,varlistcount+sIind] += -2*np.pi*1j*freq*Cp  # Add parasitic capacitance if outlet (impossible here)  
                            Rowcount += 1
                    else: 
                        
                        # We're at a parallel outlet

                        # Acquire ABCD parameters corresponding to parallel outlet load
                        ABCDpol = calculate_ABCD_parallel_outlet_load(nphase,1,freq,graph,curnode) 

                        # Implement Vin = Vout and Iin = sum Iout via adjusted Vin and Iin
                        for sIind in range(nphase):
                            VIMat[Rowcount,varlistcount+nphase+sIind] = 1                  # Sum Iin = sum Iout
                            for VindI in range(nphase):
                                VIMat[Rowcount,varlistcount+VindI] = ABCDpol[nphase+sIind,VindI] # Iin2 = CVin+Iin
                            
                            VIMat[Rowcount,varlistcount+2*nphase+sIind] = -1               # For one in one out
                            if (outlet_type[forward[curnode][0]] == 'Series outlet') or (outlet_type[forward[curnode][0]] == 'Parallel outlet'):
                                VIMat[Rowcount,varlistcount+sIind] += -2*np.pi*1j*freq*Cp  # Add parasitic capacitance if outlet (impossible here)  
                            Rowcount += 1

                    # Update variable list
                    for vlc in range(3*nphase):                                            # Add new variables 
                        VarList[varlistcount] = curnode
                        varlistcount += 1

                else:

                    # This node is on a branch but only has one forward

                    # Compute ABCD for wire
                    ABCD = compABCDwireto(graph,curnode,nphase,freq)

                    # Add ABCD equations from parent node to child node
                    for abcdind in range(2*nphase):                             # ABCD Eqs
                        VIMat[Rowcount,varlistcount+abcdind] = 1
                        for abcdind2 in range(nphase):
                            VIMat[Rowcount,np.where(VarList==parentnode)[0][0]+abcdind2] = -ABCD[abcdind,abcdind2] 
                        for abcdind2 in range(nphase,2*nphase):
                            VIMat[Rowcount,np.where(VarList==parentnode)[0][0]+WhichChild+abcdind2-nphase] = -ABCD[abcdind,abcdind2] 
                        Rowcount += 1

                    # If not a parallel outlet, use Kirchhoff's
                    if (outlet_type[curnode] != 'Parallel outlet'):
                        for sIind in range(nphase):
                            VIMat[Rowcount,varlistcount+nphase+sIind] = 1                  # Sum Iin = sum Iout
                            VIMat[Rowcount,varlistcount+2*nphase+sIind] = -1               # For one in one out
                            if (outlet_type[forward[curnode][0]] == 'Series outlet') or (outlet_type[forward[curnode][0]] == 'Parallel outlet'):
                                VIMat[Rowcount,varlistcount+sIind] += -2*np.pi*1j*freq*Cp  # Add parasitic capacitance if outlet   
                            Rowcount += 1
                    else: 
                        
                        # We're at a parallel outlet

                        # Acquire ABCD parameters corresponding to parallel outlet load
                        ABCDpol = calculate_ABCD_parallel_outlet_load(nphase,1,freq,graph,curnode) 

                        # Implement Vin = Vout and Iin = sum Iout via adjusted Vin and Iin
                        for sIind in range(nphase):
                            VIMat[Rowcount,varlistcount+nphase+sIind] = 1                  # Sum Iin = sum Iout
                            for VindI in range(nphase):
                                VIMat[Rowcount,varlistcount+VindI] = ABCDpol[nphase+sIind,VindI] # Iin2 = CVin+Iin
                            
                            VIMat[Rowcount,varlistcount+2*nphase+sIind] = -1               # For one in one out
                            if (outlet_type[forward[curnode][0]] == 'Series outlet') or (outlet_type[forward[curnode][0]] == 'Parallel outlet'):
                                VIMat[Rowcount,varlistcount+sIind] += -2*np.pi*1j*freq*Cp  # Add parasitic capacitance if outlet (impossible here)  
                            Rowcount += 1
                    
                    # Update variable list
                    for vlc in range(3*nphase):                                            # Add new variables 
                        VarList[varlistcount] = curnode
                        varlistcount += 1

            elif curnode in mainpath:
      
                # The current node splits into multiple others and is on main path

                # Calculate transmission parameters from parent to child
                ABCD = compABCDwireto(graph,curnode,nphase,freq)

                # Add ABCD equations to system
                for abcdind in range(2*nphase):                             # ABCD Eqs
                        VIMat[Rowcount,varlistcount+abcdind] = 1
                        for abcdind2 in range(nphase):
                            VIMat[Rowcount,np.where(VarList==parentnode)[0][0]+abcdind2] = -ABCD[abcdind,abcdind2] 
                        for abcdind2 in range(nphase,2*nphase):
                            VIMat[Rowcount,np.where(VarList==parentnode)[0][0]+WhichChild+abcdind2-nphase] = -ABCD[abcdind,abcdind2] 
                        Rowcount += 1
                
                # If not parallel outlet, apply Kirchhoffs to all children
                if (outlet_type[curnode] != 'Parallel outlet'):
                    for sIind in range(nphase):
                        VIMat[Rowcount,varlistcount+nphase+sIind] = 1                      # Sum Iin = sum Iout
                        for cofc in range(len(forward[curnode])):
                            VIMat[Rowcount,varlistcount+(2+cofc)*nphase+sIind] = -1
                        for otind in range(len(forward[curnode])):
                            if (outlet_type[forward[curnode][otind]] == 'Series outlet') or (outlet_type[forward[curnode][otind]] == 'Parallel outlet'):
                                VIMat[Rowcount,varlistcount+sIind] += -2*np.pi*1j*freq*Cp  # Add parasitic capacitance if outlet 
                                break  
                        Rowcount += 1

                else: 
                    
                    # It is a parallel outlet

                    # Acquire ABCD parameters corresponding to parallel outlet load
                    ABCDpol = calculate_ABCD_parallel_outlet_load(nphase,1,freq,graph,curnode) 

                    # Implement Vin = Vout and Iin = sum Iout via adjusted Vin and Iin
                    for sIind in range(nphase):
                        VIMat[Rowcount,varlistcount+nphase+sIind] = 1                        # Sum Iin = sum Iout
                        
                        for VindI in range(nphase):
                            VIMat[Rowcount,varlistcount+VindI] = ABCDpol[nphase+sIind,VindI] # Iin2 = CVin+Iin

                        for cofc in range(len(forward[curnode])):
                            VIMat[Rowcount,varlistcount+(2+cofc)*nphase+sIind] = -1
                        for otind in range(len(forward[curnode])):
                            if (outlet_type[forward[curnode][otind]] == 'Series outlet') or (outlet_type[forward[curnode][otind]] == 'Parallel outlet'):
                                VIMat[Rowcount,varlistcount+sIind] += -2*np.pi*1j*freq*Cp    # Add parasitic capacitance if outlet 
                                break  
                        Rowcount += 1

                # Update variable tracking list
                for vlc in range((len(forward[curnode])+2)*nphase):                   # V + I in + n I out
                    VarList[varlistcount] = curnode
                    varlistcount += 1

            elif curnode not in mainpath:

                # The current node splits into multiple others and is not on main path  

                # Determine ABCD parameters
                ABCD = compABCDwireto(graph,curnode,nphase,freq)

                # Add to matrix system
                for abcdind in range(2*nphase):                             # ABCD Eqs
                        VIMat[Rowcount,varlistcount+abcdind] = 1
                        for abcdind2 in range(nphase):
                            VIMat[Rowcount,np.where(VarList==parentnode)[0][0]+abcdind2] = -ABCD[abcdind,abcdind2] 
                        for abcdind2 in range(nphase,2*nphase):
                            VIMat[Rowcount,np.where(VarList==parentnode)[0][0]+WhichChild+abcdind2-nphase] = -ABCD[abcdind,abcdind2] 
                        Rowcount += 1

                # Subject children to Kirchhoff's law if not parallel outlet
                if (outlet_type[curnode] != 'Parallel outlet'):
                    for sIind in range(nphase):
                        VIMat[Rowcount,varlistcount+nphase+sIind] = 1                      # Sum Iin = sum Iout
                        for cofc in range(len(forward[curnode])):
                            VIMat[Rowcount,varlistcount+(2+cofc)*nphase+sIind] = -1
                        for otind in range(len(forward[curnode])):
                            if (outlet_type[forward[curnode][otind]] == 'Series outlet') or (outlet_type[forward[curnode][otind]] == 'Parallel outlet'):
                                VIMat[Rowcount,varlistcount+sIind] += -2*np.pi*1j*freq*Cp  # Add parasitic capacitance if outlet 
                                break  
                        Rowcount += 1

                else: 
                    
                    # It is a parallel outlet

                    # Acquire ABCD parameters corresponding to parallel outlet load
                    ABCDpol = calculate_ABCD_parallel_outlet_load(nphase,1,freq,graph,curnode) 

                    # Implement Vin = Vout and Iin = sum Iout via adjusted Vin and Iin
                    for sIind in range(nphase):
                        VIMat[Rowcount,varlistcount+nphase+sIind] = 1                        # Sum Iin = sum Iout
                        
                        for VindI in range(nphase):
                            VIMat[Rowcount,varlistcount+VindI] = ABCDpol[nphase+sIind,VindI] # Iin2 = CVin+Iin

                        for cofc in range(len(forward[curnode])):
                            VIMat[Rowcount,varlistcount+(2+cofc)*nphase+sIind] = -1
                        for otind in range(len(forward[curnode])):
                            if (outlet_type[forward[curnode][otind]] == 'Series outlet') or (outlet_type[forward[curnode][otind]] == 'Parallel outlet'):
                                VIMat[Rowcount,varlistcount+sIind] += -2*np.pi*1j*freq*Cp    # Add parasitic capacitance if outlet 
                                break  
                        Rowcount += 1

                # Track variable order
                for vlc in range((len(forward[curnode])+2)*nphase):                   # V + I in + n I out
                    VarList[varlistcount] = curnode
                    varlistcount += 1

        # Determine next node that will be visited

        tmpforw = np.sort(forward[curnode])              # sort forward nodes
        tmpforw2 = copy.copy(tmpforw)                    # to track WhichChild      
        tmpforw2 = list(set(tmpforw2)-set(mainpath))     # remove nodes on main path
        tmpforw = list(set(tmpforw)-set(NodeOrder))      # remove nodes we've already visited
        tmpforw2.sort()                                  # re-sort since list unsorts
        tmpforw.sort()                                   # re-sort since list unsorts
        parentnode = copy.copy(curnode)                  # save parent node before moving on
        WhichChild = -867530.9                           # this indexes the current variable to each split

        if len(tmpforw) == 1:

            # Only one unvisited node                  

            if curnode == transmitter:

                WhichChild = nphase                      # Assumes one line from transmitter

            elif curnode in mainpath:

                WhichChild = nphase*(len(tmpforw2)+2)    # Take current associated with last branch or main path 

            else:

                WhichChild = nphase*(len(tmpforw2)+1)

            WireID[WIDC,0] = curnode
            curnode = tmpforw[0]
            WireID[WIDC,1] = curnode
            WCI[WIDC,0] = WhichChild
            WIDC += 1

        elif len(tmpforw) > 1:

            # Multiple unvisited nodes

            WireID[WIDC,0] = curnode

            if curnode not in mainpath:

                curnode = tmpforw[0]
                WhichChild = nphase*(tmpforw2.index(curnode)+2)

            elif tmpforw[0] in mainpath:

                curnode = tmpforw[1]                     # This must be smallest variable due to sort
                WhichChild = nphase*(tmpforw2.index(curnode)+2)

            else:

                curnode = tmpforw[0]
                WhichChild = nphase*(tmpforw2.index(curnode)+2) 

            WireID[WIDC,1] = curnode
            WCI[WIDC,0] = WhichChild
            WIDC += 1

        elif len(tmpforw) == 0:

            # All child nodes visited, but receiver not reached
         
            curnode = backward[curnode][0]               # Return to the parent node

    # Add final equations for receiver
    curnode = receiver
    parentnode = backward[curnode][0]
    WhichChild = (len(forward[parentnode])+1)*nphase     # Receiver must be last child in accordance with ruleset
    if parentnode == transmitter:
        WhichChild -= nphase                             # Parent of receiver is transmitter in single wire case, which does not have children, so use transmitter current directly
    
    # Compute ABCD matrices and add equations to complete matrix system
    ABCD = compABCDwireto(graph,curnode,nphase,freq)  
    for abcdind in range(2*nphase):                      # ABCD Eqs
        VIMat[Rowcount,varlistcount+abcdind] = 1
        for abcdind2 in range(nphase):
            VIMat[Rowcount,np.where(VarList==parentnode)[0][0]+abcdind2] = -ABCD[abcdind,abcdind2] 
        for abcdind2 in range(nphase,2*nphase):
            VIMat[Rowcount,np.where(VarList==parentnode)[0][0]+WhichChild+abcdind2-nphase] = -ABCD[abcdind,abcdind2] 
        Rowcount += 1

    # For completionist sake, add the receiver as the final encountered node
    NodeOrder[-1] = receiver                             


    # Solve the sparse inverse problem

    VI = spsolve(csc_matrix(VIMat), b)


    return NodeOrder, VarList, VIMat, b, Rowcount, VI, WireID, WCI
       
