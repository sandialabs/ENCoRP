import sys
import time
import numpy as np
from compVInodes_Nwire import compVInodes_Nwire
from compVInodes_Functions import compABCDwireto, compLoadImp, compRec, compFullABCD
import scipy.io
from AssembleData import AssembleData

############################################################
##### Solve for Spatially Resolved Voltage and Current #####
############################################################

# Extract and structure data from GUI
graph,allfreq,nphase,nwires,wres,inpVI,UseI,whichRec,RecVal,nodexyz = AssembleData()

# For debugging purposes
printon = False

# Compute point-to-point transmission parameters
ABCD = compFullABCD(graph,allfreq,nphase)

# Details for saving data
saveResults = True
saveFile = '../Data/VIOut.mat'

# Define input voltage or current
VIin = inpVI*np.ones((nphase),dtype=np.cdouble)
if printon:
    print('VIin = ',VIin)

# Initialize output structure
nf = len(allfreq)                           # Number of frequencies
nw = len(graph['forward'])-1                # Number of wires
Out = np.empty(shape=(nf,nw),dtype=object)  # Empty nf x nw celled structure

# If Spice, initialize additional outputs
if graph['Spice'] == 1: # Single-phase case restricted
    NodeVI = -np.ones((len(allfreq),4*nphase+3*nphase*(nw-1)),dtype = np.cdouble)
    ListVI = -np.ones((len(allfreq),4*nphase+3*nphase*(nw-1)),dtype = np.cdouble)

# Compute spatially resolved voltages and currents at each frequency
for i in range(len(allfreq)):

    # Pull i-th frequency
    freq = np.array([allfreq[i]])
    graph['Count'] = i

    # Solve for nodal VI at i-th frequency
    NodeOrder, VarList, VIMat, b, Rowcount, VI, WireID, WCI = compVInodes_Nwire(graph,freq,nwires,UseI,VIin,whichRec,RecVal)

    # Add outputs if SPICE on
    if graph['Spice'] == 1:
        NodeVI[i] = VI
        ListVI[i] = VarList

    # Print everything if interested in per-frequency results (see debugging above)
    if printon:
        print('\nNodeOrder: ',NodeOrder)
        print('VarList: ',VarList)
        print('VIMat: ',VIMat)
        print('b: ',b)
        print('Rowcount: ',Rowcount)
        print('VI: ',VI)
        print('WCI: ',WCI)
        print('Sparsity: ',(1-np.count_nonzero(VIMat)/VIMat.shape[0]/VIMat.shape[1])*100,'%')
        print('Freq: ',freq,' All Freqs: ',allfreq)
        print('WireID: ',WireID+1)
        print('Output Size: ',OutputVI.shape) 
    
    # Propagate nodal solutions along wires
    for ind in range(len(graph['forward'])-1):                     # For each wire
        Vind = np.where(VarList == WireID[ind,0])[0][0]            # Find the wire start index
        VI0_V = VI[Vind:Vind+nphase]                               # Store the voltages
        VI0_I = VI[Vind+WCI[ind][0]:Vind+WCI[ind][0]+nphase]       # Store the corresponding currents (index for correct "child")
        VI0 = np.zeros((2*nphase),dtype = np.cdouble)                 
        VI0[:nphase] = VI0_V                                       # Store VI in vec
        VI0[nphase:] = VI0_I                                       
       
        Full_L = graph['wire_length'][WireID[ind,1]]               # Pull wire length

        # Initialize matrix for all points along the wire for all phases for voltage and current
        npoints = int(np.ceil(Full_L/wres))
        OutputVI = np.zeros((npoints+1,2*nphase),dtype = np.cdouble)

        OutputVI[0,:] = VI0                                            # VI at start are unchanged
        if npoints > 1.5:                                              # If the wire length exceeds 1 wres
            graph['wire_length'][WireID[ind,1]] = wres                 # Set wire length to wres
            ABCDtmp = compABCDwireto(graph,WireID[ind,1],nphase,freq)  # Obtain ABCD for this wire for a length of wire res
            VItmp = np.matmul(ABCDtmp,VI0)                             # Compute the second VI
            OutputVI[1,:] = VItmp                                      # Store second VI (wres along wire)
            for L in range(2,npoints):                                 # Until you reach end-1 of wire
                VItmp = np.matmul(ABCDtmp,VItmp)                       # Compute VI
                OutputVI[L,:] = VItmp                                  # Store
            graph['wire_length'][WireID[ind,1]] = Full_L               # Reset wire length to original value
        ABCDtmp = compABCDwireto(graph,WireID[ind,1],nphase,freq)      # Compute full length ABCD for last point
        VItmp = np.matmul(ABCDtmp,VI0)                                 # Compute VI
        OutputVI[npoints] = VItmp                                      # Store final VI

        # Store this matrix for the given frequency and wire
        Out[i,ind] = OutputVI

# Save results as MAT file:
if saveResults:
    if graph['Spice'] == 1:
        scipy.io.savemat(saveFile,mdict={'Freq': allfreq, 'WireID': WireID+1, 'VI': Out, 'WireRes': wres, 'Nodes': nodexyz, 'ABCD': ABCD, 'NodeVI': NodeVI, 'ListVI': ListVI})
    else:
        scipy.io.savemat(saveFile,mdict={'Freq': allfreq, 'WireID': WireID+1, 'VI': Out, 'WireRes': wres, 'Nodes': nodexyz, 'ABCD': ABCD})
