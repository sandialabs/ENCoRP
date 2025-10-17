import numpy as np
from numpy.linalg import inv, cond

####################################################
#####            Full Network ABCD             #####
####################################################

def compFullABCD(graph,freq,nphase):

    # Set up to compute network ABCD
    main_path = graph['path']
    starting_node = graph['transmitter']

    # Solve path starting from transmitter
    try:
        N_f = len(freq)
    except:
        freq = [freq]
        N_f = 1

    # Compute transmission parameters for each frequency
    ABCD = np.empty((N_f,2*nphase,2*nphase),dtype=np.cdouble)
    for freqind in range(N_f):
        graph['Count'] = freqind
        ABCD0 = solve_path(starting_node,graph,[freq[freqind]],main_path,nphase)
        ABCD[freqind,:,:] = np.squeeze(ABCD0)

    # Reshape
    ABCD = ABCD.transpose(1,2,0)

    return ABCD
    
    
def solve_path(current_node,graph,freq,main_path,nphase):
    
    # Pull forward list details for current node
    forward_nodes = graph['forward'][current_node]
    num_forward = len(forward_nodes)

    # Compute ABCD matrix with wire for current node
    ABCD_loading = calculate_ABCD_load_and_wire(current_node,graph,freq,nphase)
    
    # If no forward nodes
    if num_forward == 0:   

        # Terminating, return ABCD
        ABCD = ABCD_loading

    # If one forward node
    elif num_forward == 1:    

        # Solve for child node and combine
        ABCD_forward = solve_path(forward_nodes[0],graph,freq,main_path,nphase)
        ABCD = np.matmul(ABCD_forward, ABCD_loading)

    # Multiple forward nodes
    else:           

        # Reorder forward nodes so element on main path is first
        forward_nodes = reorder_forward_nodes(forward_nodes,main_path)

        # For each forward node, calculate ABCD matrix
        ABCD_forward_nodes = []
        for forward_node in forward_nodes:
            ABCD_f = solve_path(forward_node,graph,freq,main_path,nphase)
            ABCD_forward_nodes.append(ABCD_f)

        # Combine ABCD matrices
        ABCD_forward = combine_parallel_ABCD_matrices(ABCD_forward_nodes,nphase)
        ABCD = np.matmul(ABCD_forward, ABCD_loading)

    return ABCD


def reorder_forward_nodes(forward_nodes,main_path):
    
    # Determines which forward nodes are on main path
    intersection = set(main_path).intersection(forward_nodes) 

    # If any on main path, reorder
    if len(intersection) > 0:   

        # Pull index of main path node
        intersection_element = intersection.pop()

        # Move main path node to beginning
        forward_nodes.insert(0, forward_nodes.pop(forward_nodes.index(intersection_element)))

    return forward_nodes


def combine_parallel_ABCD_matrices(Phi_matrix_list,n):

    # Combine all elements in parallel using sum inv(B)*A into factor
    factor = 0            

    # Start at 1 to skip main path
    for i in range(1, len(Phi_matrix_list)):

        # Pull i-th ABCD 4D matrix
        Phi_matrix = Phi_matrix_list[i]

        # Add to factor inv(B)*A
        factor += mult(inv(Phi_matrix[:,:,0:n,n:2*n]), Phi_matrix[:,:,0:n,0:n])

    # Pull main path ABCD or arbitrary branch
    Phi_1 = Phi_matrix_list[0]

    # Combine
    Phi = Phi_1.copy()
    Phi[:,:,0:n,0:n]   = Phi_1[:,:,0:n,0:n]   + mult(Phi_1[:,:,0:n,n:2*n], factor)
    Phi[:,:,n:2*n,0:n] = Phi_1[:,:,n:2*n,0:n] + mult(Phi_1[:,:,n:2*n,n:2*n], factor)

    return Phi


def calculate_ABCD_load_and_wire(node_number,graph,freq,n):

    # Pull nodes information
    load_type = graph['load_type'][node_number]
    outlet_type = graph['outlet_type'][node_number]
    load_parameters = graph['load_parameters'][node_number]
    wire_length = graph['wire_length'][node_number]
    count = graph['Count']

    # Calculate combined ABCD matrix of load and backward wire
    N_f = len(freq)
        
    # Loop over frequencies and store
    Phi = np.empty((N_f, 1, 2*n, 2*n), dtype=np.cdouble)
    for freqind in range(N_f):

        # If at transmitter
        if wire_length == None: 
    
            # Set as identity matrix
            Phi_wire = np.zeros((1, 1, 2*n, 2*n), dtype=np.cdouble)
            for i in range(2*n):
                Phi_wire[:,:,i,i] = 1
    
        # Not at transmitter, proceed accordingly
        else:
    
            # Compute ABCD of wire
            Phi_wire = np.empty((1, 1, 2*n, 2*n), dtype=np.cdouble)
            Phi_wire[:,:,:,:] = compABCDwireto(graph,node_number,n,freq)
    
        # If no load, do not need to add load ABCD
        if load_parameters is None:                                 
            Phi[freqind,:,:,:] = Phi_wire
    
        # Series outlet or receiver
        elif outlet_type == 'Series outlet' or outlet_type == 'Receiver':
    
            # Obtain impedance matrix 
            Y_L = load_admittance_matrix(load_type, load_parameters, freq[freqind], n, 1, count)
            Z_L = inv(Y_L)
    
            # Combine with wire
            Phi[freqind,:,:,:] = Phi_wire
            Phi[freqind,:,0:n,0:n]   = Phi_wire[:,:,0:n,0:n]   - mult( Z_L, Phi_wire[:,:,n:2*n,0:n])
            Phi[freqind,:,0:n,n:2*n] = Phi_wire[:,:,0:n,n:2*n] - mult( Z_L, Phi_wire[:,:,n:2*n,n:2*n])
    
        # Parallel outlet
        elif outlet_type == 'Parallel outlet':
    
            # Obtain admittance matrix
            Y_L = load_admittance_matrix(load_type, load_parameters, freq[freqind], n, 1, count)
    
            # Combine with wire
            Phi[freqind,:,:,:] = Phi_wire
            Phi[freqind,:,n:2*n,0:n]   = - mult(Y_L, Phi_wire[:,:,0:n,0:n])   + Phi_wire[:,:,n:2*n,0:n]
            Phi[freqind,:,n:2*n,n:2*n] = - mult(Y_L, Phi_wire[:,:,0:n,n:2*n]) + Phi_wire[:,:,n:2*n,n:2*n]

    return Phi

####################################################
#####             WIRE FUNCTIONS               #####
####################################################

def compABCDwireto(graph,child,nphase,freq):

    # Load relevant parameters     
    wirelen = graph['wire_length'][child]
    wirerad = graph['wire_radius'][child]
    wiredis = graph['wire_distance'][child]
    wiretyp = graph['wire_type'][child]
    Spice = graph['Spice']
    sigmaNew = graph['Sigma'][child]
    epsNew = graph['Epsilon'][child]
 
    # Calculate the transmission matrix for the wire
    phi = calculate_wire_matrix(wirelen, wirerad, wiredis, wiretyp, freq, nphase, len(freq), Spice, sigmaNew, epsNew)

    return phi[0] 


def calc_wire_property_matrices(l, r_w, d0, wire_type, f, n, N_f, Spice, sigmaNew, epsNew):

    # If Spice comparison, return fixed RLCG
    if Spice == 1:
        if n != 1:
            print('Warning: Spice comparison only valid for single-phase networks. Defaulting to computed R, L, C, G.')
        else: # n = 1
            R = np.zeros((N_f,1,1))
            L = np.zeros((N_f,1,1))
            C = np.zeros((N_f,1,1))
            G = np.zeros((N_f,1,1))
            for indf in range(N_f):
                R[indf,0,0] = 1e-2  # These values are near the computed
                L[indf,0,0] = 5e-7  # for conventional values at low-mid
                C[indf,0,0] = 8e-11 # frequencies.
                G[indf,0,0] = 0
            return R, L, C, G

    # Define wire distances based on geometry
    dij = np.zeros((n+1,n+1))
    if wire_type == 'line':        # 2 wires, no other choice
        dij[0,1] = d0
    elif wire_type == 'triangle':  # 3 wires equally spaced in triangle
        dij[0,1] = d0
        dij[0,2] = d0
        dij[1,2] = d0
    elif wire_type == 'square':    # 4 wires in a square
        dij[0,1] = d0
        dij[0,2] = np.sqrt(2)*d0
        dij[0,3] = d0
        dij[1,2] = d0
        dij[1,3] = np.sqrt(2)*d0
        dij[2,3] = d0

    # Properties for wire materials (Default: copper/PVC) (Length n+1)
    mu = 4 * np.pi * 1e-7     # Vacuum permeability
    sigma = sigmaNew          # Default: 5.85e7, Conductor conductivity
    epsilon = epsNew          # Default: 3.6*8.854e-12, Insulator permittivity

    # Skin depth
    delta = np.zeros((N_f,n+1))
    for sdi in range(n+1):
        if N_f > 1:
            delta[:,sdi] = 1/np.sqrt(np.pi*mu*sigma[sdi]*f)
        else:
            delta[:,sdi] = 1/np.sqrt(np.pi*mu*sigma[sdi]*f[0])

    # Resistance matrix 
    R = np.zeros((N_f, n, n))

    # Compute r0 with respect to reference conductor (i.e. ground/neutral)
    for i in range(N_f):     
        if r_w[0] <= 2*delta[i,0]:   
            res0 = 1/(sigma[0]*np.pi*r_w[0]**2)
        else:
            res0 = 1/(2*r_w[0])*np.sqrt(mu*f[i]/(np.pi*sigma[0])) 

        # Size of matrix is a function of nphase
        for ind1 in range(n):
            for ind2 in range(n):
                if ind1 == ind2:
                    # Calculate ri with respect to non-reference conductor i
                    if r_w[ind1+1] <= 2*delta[i,ind1+1]:   
                        resi = 1/(sigma[ind1+1]*np.pi*r_w[ind1+1]**2)
                    else:
                        resi = 1/(2*r_w[ind1+1])*np.sqrt(mu*f[i]/(np.pi*sigma[ind1+1])) 
                    R[i,ind1,ind2] = res0 + resi
                else:
                    R[i,ind1,ind2] = res0

    # PUL inductance matrix 
    L = np.zeros((N_f,n,n))
    for ind1 in range(n):
        for ind2 in range(ind1+1):
            if ind1 == ind2:
                L[:,ind2,ind1] = mu/(2*np.pi) * np.log(dij[0,ind2+1]*dij[0,ind1+1]/r_w[0]/r_w[ind1+1])
            else:
                L[:,ind2,ind1] = mu/(2*np.pi) * np.log(dij[0,ind2+1]*dij[0,ind1+1]/r_w[0]/dij[ind2+1,ind1+1])
    for indf in range(N_f):
        L[indf,:,:] += np.triu(L[indf,:,:],k=1).T

    # PUL capacitance matrix
    C = np.zeros((N_f,n,n))
    
    # Compute C
    for indc in range(N_f):
        C[indc,:,:] = mu * epsilon * inv(L[indf,:,:])

    # PUL conductance matrix 
    G = np.zeros((N_f,n,n))

    return R, L, C, G


def calculate_wire_matrix(l, r_w, d0, wire_type, f, n, N_f, Spice, sigmaNew, epsNew):

    # Compute R, L, C, G
    R, L, C, G = calc_wire_property_matrices(l, r_w, d0, wire_type, f, n, N_f, Spice, sigmaNew, epsNew)
    
    # Calculate Z, Y, and YZ
    Z = R + 2*np.pi*1j * np.expand_dims(np.expand_dims(f,1),2) * L
    Y = G + 2*np.pi*1j * np.expand_dims(np.expand_dims(f,1),2) * C
    YZ = np.matmul(Y,Z)

    # Compute gamma and T
    gam2, T = np.linalg.eig(YZ)
    gam = np.zeros((N_f,n,n), dtype=np.cdouble)
    for i in range(n):
        gam[:,i,i] = np.sqrt(gam2[:,i])

    # Calculate Zc and some inverses for wire definition
    YZ_sqrt = mult(T, gam, inv(T))
    Z_c = np.matmul(inv(Y), YZ_sqrt)
    Y_c = inv(Z_c)
    Y_inv = inv(Y)
    T_inv = inv(T)

    # Compute cosh(gamma * L) & sinh(gamma * L)
    cosh_gamL = np.zeros((N_f,n,n), dtype=np.cdouble)
    sinh_gamL = np.zeros((N_f,n,n), dtype=np.cdouble)
    for i in range(n):
        cosh_gamL[:,i,i] = np.cosh(gam[:,i,i] * l)
        sinh_gamL[:,i,i] = np.sinh(gam[:,i,i] * l)

    # Calculate transmission parameter blocks
    phi11 = mult(Y_inv, T, cosh_gamL, T_inv, Y)         
    phi12 = - mult(Z_c, T, sinh_gamL, T_inv)
    phi21 = - mult(T, sinh_gamL, T_inv, Y_c)
    phi22 = mult(T, cosh_gamL, T_inv)

    # Combine into transmission matrix
    Phi = np.zeros((N_f, 2*n, 2*n), dtype=np.cdouble)
    Phi[:,0:n,0:n] = phi11
    Phi[:,0:n,n:2*n] = phi12
    Phi[:,n:2*n,0:n] = phi21
    Phi[:,n:2*n,n:2*n] = phi22

    return Phi


####################################################
#####             LOAD FUNCTIONS               #####
####################################################


def calculate_ABCD_parallel_outlet_load(n,N_f,f,graph,curnode):

    # Extract load type and parameters
    load_type = graph['load_type'][curnode]
    load_parameters = graph['load_parameters'][curnode]
    count = graph['Count']

    # Set as identity matrix since wire already incorporated
    Phi = np.identity(2*n, dtype=np.cdouble)              
    
    # Compute admittance matrix
    Y_L = load_admittance_matrix(load_type, load_parameters, f, n, N_f, count)

    # Set up transmission parameter matrix
    Phi[n:2*n,0:n] = -Y_L

    return Phi


def compLoadImp(graph,child,nphase,freq):
   
    # Pull load information
    load_type = graph['load_type'][child]
    load_parameters = graph['load_parameters'][child]
    count = graph['Count']

    # Compute load admittance matrix
    Y_L = load_admittance_matrix(load_type, load_parameters, freq, nphase, len(freq), count)
    Z_L = inv(Y_L)

    return Z_L  


def calculate_load_impedance(load_type, load_parameters, lineA, lineB, f, n, N_f):
    
    # Define omega and s and capacitor
    omega = 2*np.pi*f
    s = omega*1j
    
    # Determine if line to line or line to ground impedance
    if ((lineA == 0) or (lineB == 0)) and (n != 1):
        lineType = 'LineToGround'
    else:
        lineType = 'LineToLine'

    # Calculate load impedance in accordance with selection

    # NONE (grounded/zero impedance)
    if load_type is None:
        capacitor = load_parameters[0]
        if lineType == 'LineToLine':
            Z = 1e-16 * np.ones((N_f), dtype=np.cdouble)
        elif lineType == 'LineToGround':
            Z = 1e-16 * np.ones((N_f), dtype=np.cdouble) 

    # OFF (open-circuit/"infinite" impedance)
    elif load_type == 'Off':
        capacitor = load_parameters[0]
        if lineType == 'LineToLine':
            Z = 1e6 * np.ones((N_f), dtype=np.cdouble)
        elif lineType == 'LineToGround': 
            Z = 1e6 * np.ones((N_f), dtype=np.cdouble)

    # CONSTANT    
    elif load_type == 'Constant':
        R, capacitor = load_parameters
        if lineType == 'LineToLine':
            Z = R * np.ones((N_f), dtype=np.cdouble)
        elif lineType == 'LineToGround':
            Z = 1 / (s * capacitor)

    # MOTOR
    elif load_type == 'Motor':
        L, R, C, Cg, r, capacitor = load_parameters
        if lineType == 'LineToLine':
            Z = 1/3 * (L*s + r) / ( (L*s + r)*(C+Cg/2)*s + (L*s+r)/R + 1 )
        elif lineType == 'LineToGround':
            Z = 1/3 * 1 / (s*Cg + 1/(s*Cg + R*(s*L+r)/(s*C*(s*L+r)*R + s*L+r + R) ) )

    # DOUBLE RLC
    elif load_type == 'Double RLC':
        R_ser, omega_ser, zeta_ser, R_par, omega_par, zeta_par, dev13, dev23, capacitor = load_parameters
        if lineType == 'LineToLine':
            omega_par_ratio = omega / omega_par
            omega_ser_ratio = omega / omega_ser
            Z_par = 1j * omega_par_ratio * 2 * R_par * zeta_par / (1 + 1j * omega_par_ratio * 2 * zeta_par - omega_par_ratio**2)
            Z_ser = R_ser + 1j * 2 * R_ser * zeta_ser * (omega_ser_ratio - 1/omega_ser_ratio)
            Z = Z_par + Z_ser
            if ((lineA == 1) and (lineB == 2)) or ((lineA == 2) and (lineB == 1)):
                pass
            elif ((lineA == 1) and (lineB == 3)) or ((lineA == 3) and (lineB == 1)):
                Z = Z * (1 + dev13)
            elif ((lineA == 2) and (lineB == 3)) or ((lineA == 3) and (lineB == 2)):
                Z = Z * (1 + dev23)
        elif lineType == 'LineToGround': 
            Z = 1 / (s * capacitor)

    # SERIES RLC 
    elif load_type == 'Series RLC':
        R, L, C, capacitor = load_parameters
        if lineType == 'LineToLine':
            if R < 1e-15 and L < 1e-15 and C < 1e-15:
                Z = 1e-16 * np.ones((N_f), dtype=np.cdouble)
            elif R > 1e-15 and L < 1e-15 and C < 1e-15:
                Z = R
            elif R < 1e-15 and L > 1e-15 and C < 1e-15:
                Z = s*L
            elif R < 1e-15 and L < 1e-15 and C > 1e-15:
                Z = 1/(s*C)
            elif R > 1e-15 and L > 1e-15 and C < 1e-15:
                Z = R + s*L
            elif R > 1e-15 and L < 1e-15 and C > 1e-15:
                Z = R + 1/(s*C)
            elif R < 1e-15 and L > 1e-15 and C > 1e-15:
                Z = s*L + 1/(s*C)
            else:
                Z = R + s*L + 1/(s*C)
        elif lineType == 'LineToGround':
            Z = 1 / (s * capacitor)
            
    # PARALLEL RLC
    elif load_type == 'Parallel RLC':
        R, L, C, capacitor = load_parameters
        if lineType == 'LineToLine':
            if R < 1e-15 and L < 1e-15 and C < 1e-15:
                Z = 1e6 * np.ones((N_f), dtype=np.cdouble)
            elif R > 1e-15 and L < 1e-15 and C < 1e-15:
                Z = 1 / (1/R)
            elif R < 1e-15 and L > 1e-15 and C < 1e-15:
                Z = 1 / (1/(s*L))
            elif R < 1e-15 and L < 1e-15 and C > 1e-15:
                Z = 1 / (s*C)
            elif R > 1e-15 and L > 1e-15 and C < 1e-15:
                Z = 1 / (1/R + 1/(s*L))
            elif R > 1e-15 and L < 1e-15 and C > 1e-15:
                Z = 1 / (1/R + s*C)
            elif R < 1e-15 and L > 1e-15 and C > 1e-15:
                Z = 1 / (1/(s*L) + s*C)
            else:
                Z = 1 / (1/R + 1/(s*L) + s*C)
        elif lineType == 'LineToGround':
            Z = 1 / (s * capacitor)

    return Z


def load_admittance_matrix(load_type, load_parameters, f, n, N_f, count):

    # Ensure N_f = 1
    if N_f != 1:
        print('Error: expected Nf in LAM to be 1.')
        exit()

    # For built-in loads use the above "library"
    if load_type != 'Custom':

        # Compute line to line and line to ground impedances
        Zij = np.zeros((n+1,n+1), dtype=np.cdouble)
        for ind1 in range(n):
            for ind2 in range(ind1+1,n+1):
                Zij[ind1,ind2] = calculate_load_impedance(load_type, load_parameters, ind1, ind2, f, n, N_f)
      
        # Calculate admittance matrix
        Y_L = np.zeros((n,n), dtype=np.cdouble)
        for ind1 in range(n):
            for ind2 in range(ind1+1):
                if ind1 == ind2:
                    for ind3 in range(ind1+1):
                        Y_L[ind1,ind1] += 1/Zij[ind3,ind1+1]
                    for ind3 in range(ind1+2,n+1):
                        Y_L[ind1,ind1] += 1/Zij[ind1+1,ind3]
                else:
                    Y_L[ind2,ind1] = -1/Zij[ind2+1,ind1+1]

        Y_L += np.triu(Y_L,k=1).T

    else:
        
        # Custom load, use user input to define Z matrix
        userinput = load_parameters

        # Pull data corresponding to i-th freq
        curdata = userinput[count]
        
        # Check if full matrix or L2L/L2G components
        if len(curdata) == sum(range(1,n+1)):
            # L2L and L2G components, set and construct admittance matrix
            # Compute line to line and line to ground impedances
            Zij = np.zeros((n+1,n+1), dtype=np.cdouble)
            tmpc = 0
            for ind1 in range(n):
                for ind2 in range(ind1+1,n+1):
                    Zij[ind1,ind2] = curdata[tmpc]
                    tmpc += 1
      
            # Calculate admittance matrix
            Y_L = np.zeros((n,n), dtype=np.cdouble)
            for ind1 in range(n):
                for ind2 in range(ind1+1):
                    if ind1 == ind2:
                        for ind3 in range(ind1+1):
                            Y_L[ind1,ind1] += 1/Zij[ind3,ind1+1]
                        for ind3 in range(ind1+2,n+1):
                            Y_L[ind1,ind1] += 1/Zij[ind1+1,ind3]
                    else:
                        Y_L[ind2,ind1] = -1/Zij[ind2+1,ind1+1]
           
            Y_L += np.triu(Y_L,k=1).T
           
        else:

            # Full Z matrix, reshape and invert (row ~ "reading order")
            Zmat = np.zeros((n,n), dtype = np.cdouble)
            tmpc = 0
            for zind1 in range(n):
                for zind2 in range(n):
                    Zmat[zind1,zind2] = curdata[tmpc]
                    tmpc += 1
            if cond(Zmat) > 1e16:
                print('Error: user supplied load impedance matrix is ill-conditioned.')
                exit()
 
            Y_L = inv(Zmat)

    ## Check custom load resulting matrix for ill-conditioning
    #if cond(Y_L) > 1e16 and np.min(abs(Y_L)) > 1e16:
    #    print('Error with CLM:',cond(Y_L),Y_L,inv(Y_L))

    return Y_L


###################################################
#####             RECEIVER  MODEL             #####
###################################################


def compRec(graph, receiver, freq, nphase, whichRec, RecVal):
  
    # Call receiver model to obtain impedance matrix or rec voltages or current
    if whichRec == 'Z':
        Z_R = compLoadImp(graph,receiver,nphase,freq)
        return Z_R
    elif whichRec == 'V' or whichRec == 'I':
        return RecVal


###################################################
#####         SUPPLEMENTARY FUNCTIONS         #####
###################################################


def mult(*arg):
    prod = arg[0]
    for i in range(1,len(arg)):
        prod = np.matmul(prod,arg[i])
    return prod

