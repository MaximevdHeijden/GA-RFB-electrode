import openpnm as op
import numpy as np
from inputDict_GA_V import input_dict as cnst

def validate_face_labels(net):
    r"""
    This function validates if the network has the correct face labels.    
    Usage example: validate_face_labels(net=net_c)
    """
# Update face labels if the network is extracted using the SNOW algorithm:
    if 'pore.left' and 'pore.front' and 'pore.top' in net:
        print('Network contains the correct labels, please make sure that the'
              ' labels are consistent with:\n'
                  '\'left\' and \'right\' - x-dimension\n'
                  '\'front\' and \'back\' - y-dimension\n'
                  '\'top\' and \'bottom\' - z-dimension\n')
    else: 
        raise KeyError('Network does not contain the correct boundary pore labels.\n'
                       'Please assign the following network labels:\n'
                       '\'left\' and \'right\' - x-dimension\n'
                       '\'front\' and \'back\' - y-dimension\n'
                       '\'top\' and \'bottom\' - z-dimension\n')


def delete_double_boundary_pores(net):
    r'''
    This function deletes boundary pores that are both on the boundary for the
    advection-diffusion algirhtm and for ohmic conductiction. 
    # NOTE THIS IS DIFFERENT THAT IN THE CASE FOR EXTRACTED NETWORKS! 
    The SNOW algorithm adds the labels 'front' and 'back' to the x-dimension 
                                        'left' and 'right' to the y-dimension                             
                                       'top' and 'bottom' to the z-dimension
    '''
    W_dim = cnst['width_dimension']
    L_dim = cnst['length_dimension']
    
    if W_dim == 0 and L_dim == 1 or  W_dim == 1 and L_dim == 0:
        ps1 = net.pores(['front', 'left'], mode = 'and')
        ps2 = net.pores(['front', 'right'], mode = 'and')
        ps3 = net.pores(['back', 'right'], mode = 'and')
        ps4 = net.pores(['back', 'left'], mode = 'and')        
        op.topotools.trim(network=net, pores=[ps1, ps2, ps3, ps4])
        
    elif W_dim == 0 and L_dim == 2 or L_dim == 0 and W_dim == 2:
        ps1 = net.pores(['front', 'top'], mode = 'and')
        ps2 = net.pores(['front', 'bottom'], mode = 'and')
        ps3 = net.pores(['back', 'top'], mode = 'and')
        ps4 = net.pores(['back', 'bottom'], mode = 'and')       
        op.topotools.trim(network=net, pores=[ps1, ps2, ps3, ps4])
    
    elif W_dim == 1 and L_dim == 2 or W_dim == 2 and L_dim == 1 :
        ps1 = net.pores(['top', 'left'], mode = 'and')
        ps2 = net.pores(['bottom', 'right'], mode = 'and')
        ps3 = net.pores(['top', 'right'], mode = 'and')
        ps4 = net.pores(['bottom', 'left'], mode = 'and')               
        op.topotools.trim(network=net, pores=[ps1, ps2, ps3, ps4])
    
    
def assign_boundary_poresFF(net, W_dim, L_dim):
    r"""
    This function assigns the right labels to the boundary pores for the given network.
    
    The cubic network adds the labels 'left' and 'right' to the y-dimension
                                       'front' and 'back' to the x-dimension
                                       'top' and 'bottom' to the z-dimension
    Because we invert the network in the thickness, we can assign the same boundary
    pores for the membrane in both the anode and the cathode in this step.
    ."""
    if W_dim == 0:
        net['pore.membrane'] = net['pore.front']
        net['pore.current_collector'] = net['pore.back']
    elif W_dim == 1:
        net['pore.membrane'] = net['pore.left']
        net['pore.current_collector'] = net['pore.right']
    elif W_dim == 2:
        net['pore.membrane'] = net['pore.top']
        net['pore.current_collector'] =  net['pore.bottom']
    
    if L_dim == 0:
        net['pore.flow_inlet'] = net['pore.front']
        net['pore.flow_outlet'] = net['pore.back']
    elif L_dim == 1:
        net['pore.flow_inlet'] = net['pore.left']
        net['pore.flow_outlet'] = net['pore.right']
    elif L_dim == 2:
        net['pore.flow_inlet'] = net['pore.bottom']
        net['pore.flow_outlet'] = net['pore.top']
        
        
def assign_boundary_poresIDFF(net, W_dim, L_dim, H_dim, H, H_min, L_min, W_min, L, W):
    r"""
    This function assigns the right labels to the boundary pores for the given 
    network with an interdigitated flow field
    
    The SNOW algorithm adds the labels 'left' and 'right' to the x-dimension
                                       'front' and 'back' to the y-dimension
                                       'top' and 'bottom' to the z-dimension
    Because we invert the network in the thickness, we can assign the same boundary
    pores for the membrane in both the anode and the cathode in this step. 
    
    The inlet and outlet is assigned with boundary conditions or via flow field (FF) pores:

    The network will be sliced in 3 parts via the height along the length and 
    thickness of the electrode:     Inlet channel - rib - Outlet channel
    This will be done by assigning masks based on the height (pore y-coordinate)
    domain and the by selecting all pores on the left boundaries.
        Inlet channel:  H_min < Pore y-coordinate <= Quarter_heigth + H_min
        rib:            H_min + Quarter_heigth < Pore y-coordinate <= H_min + 3 * Quarter_heigth
        Outlet channel: H_min + 3 * Quarter_heigth < Pore y-coordinate <= H_min + 4 * Quarter_heigth
    The boundary pores of the inlet channel will be assigned the label: 'pore.flow_inlet'.
    The boundary pores of the outlet channel will be assigned the label: 'pore.flow_outlet'
        
    The boundary pores of the rib will be assigned the label: 'pore.current_collector'
    ."""
    
    total_heigth = H - H_min
    Quarter_heigth = total_heigth/4
    
    Heigth_0 = H_min
    Heigth_1 = H_min + Quarter_heigth
    Heigth_2 = H_min + 3 * Quarter_heigth
    Heigth_3 = H_min + 4* Quarter_heigth

    # Masking: note "xxx_mask_boolean" and "xxx_mask" contain ALL pores at this height -> we need to
    # select only those present at the left surface of the electrode (this excludes the other boundary pores like bottom and top. 
    # The pores at the left AND at the specefic heights are then selected in 'xxx_mask_left':
    Inlet_mask_boolean = np.logical_and(net['pore.coords'][net.pores(), H_dim] >= Heigth_0,
                                net['pore.coords'][net.pores(), H_dim] <= Heigth_1)
    Inlet_mask = np.where(Inlet_mask_boolean)[0]
    Inlet_heigth_pores = net.pores('all')[Inlet_mask]
    net.set_label(label='height_1', pores = Inlet_heigth_pores)
    Inlet_mask_left = net.pores(['height_1', 'bottom'], mode='and')
    Inlet_pores = net.pores('all')[Inlet_mask_left]
    net.set_label(label = 'flow_inlet', pores=Inlet_pores)

    Rib_mask_boolean = np.logical_and(net['pore.coords'][net.pores(), H_dim] > Heigth_1,
                                net['pore.coords'][net.pores(), H_dim] <= Heigth_2)
    Rib_mask = np.where(Rib_mask_boolean)[0]
    Rib_heigth_pores = net.pores('all')[Rib_mask]
    net.set_label(label='height_2', pores = Rib_heigth_pores)
    Rib_mask_left = net.pores(['height_2', 'bottom'], mode='and')
    Rib_pores = net.pores('all')[Rib_mask_left]
    net.set_label(label = 'pore.current_collector', pores=Rib_pores)

    Outlet_mask_boolean = np.logical_and(net['pore.coords'][net.pores(), H_dim] > Heigth_2,
                                net['pore.coords'][net.pores(), H_dim] <= Heigth_3)
    Outlet_mask = np.where(Outlet_mask_boolean)[0]
    Outlet_heigth_pores = net.pores('all')[Outlet_mask]
    net.set_label(label='height_3', pores = Outlet_heigth_pores)
    Outlet_mask_left = net.pores(['height_3', 'bottom'], mode='and')
    Outlet_pores = net.pores('all')[Outlet_mask_left]
    net.set_label(label = 'pore.flow_outlet', pores=Outlet_pores)
               
    net['pore.membrane'] = net['pore.top']


def add_throat_surface_area_to_pores(net):
    r"""
    This function updatse the internal surface area of every pore with half the area of 
    the connecting throats to account for surface area of the throats.
    
    Usage example: add_throat_surface_area_to_pores(net=net_c)
    """
    net['throat.surface_area'] = op.models.geometry.throat_surface_area.cylinder(net, throat_diameter='throat.diameter', throat_length='throat.length')

    for pore in net.Ps:
        connected_throats = net.find_neighbor_throats(pores=pore)
        net['pore.surface_area'][pore] += 1 / 2 * np.sum(net['throat.surface_area'][connected_throats])
        
        
def bv_rate_constant_oc_c(c, eta, Ai_c, rp_c):
    r"""
    Calculates A2 in linear kinetics for OhmicConduction algorithm in the cathode.
    """
    c1 = cnst['j0_c'] * Ai_c * c / cnst['conc_ref_c']
    c2 = cnst['j0_c'] / (cnst['F'] * cnst['conc_ref_c'] * cnst['D_c'] / rp_c)
    c3 = cnst['j0_c'] / (cnst['F'] * cnst['conc_ref_c'] * cnst['D_a'] / rp_c)
    arg1 = cnst['alpha_a_c'] * cnst['val_c'] * cnst['F'] / (cnst['R_g'] * cnst['T']) * eta
    arg2 = -cnst['alpha_c_c'] * cnst['val_c'] * cnst['F'] / (cnst['R_g'] * cnst['T']) * eta
    return (c1 * (np.exp(arg1) - np.exp(arg2)))/(1 + c2 * np.exp(arg1) + c3 * np.exp(arg2))


def bv_rate_constant_oc_a(c, eta, Ai_a, rp_a):
    r"""
    Calculates A2 in linear kinetics for OhmicConduction algorithm in the anode.
    """
    c1 = cnst['j0_a'] * Ai_a * c / cnst['conc_ref_a']
    c2 = cnst['j0_a'] / (cnst['F'] * cnst['conc_ref_a'] * cnst['D_c'] / rp_a)
    c3 = cnst['j0_a'] / (cnst['F'] * cnst['conc_ref_a'] * cnst['D_a'] / rp_a)
    arg1 = -cnst['alpha_a_a'] * cnst['val_a'] * cnst['F'] / (cnst['R_g'] * cnst['T']) * eta
    arg2 = cnst['alpha_c_a'] * cnst['val_a'] * cnst['F'] / (cnst['R_g'] * cnst['T']) * eta
    return (c1 * (np.exp(arg1) - np.exp(arg2)))/(1 + c2 * np.exp(arg1) + c3 * np.exp(arg2))


def bv_rate_derivative_oc_c(conc_c, eta_c, Ai_c, rp_c):
    c1 = cnst['j0_c'] * Ai_c * conc_c / (cnst['conc_ref_c'])
    c2 = cnst['j0_c'] / (cnst['F'] * cnst['conc_ref_c'] * cnst['D_c'] / rp_c)
    c3 = cnst['j0_c'] / (cnst['F'] * cnst['conc_ref_c'] * cnst['D_a'] / rp_c)
    
    m1 = cnst['alpha_a_c'] * cnst['val_c'] * cnst['F'] / (cnst['R_g'] * cnst['T'])
    m2 = cnst['alpha_c_c'] * cnst['val_c'] * cnst['F'] / (cnst['R_g'] * cnst['T'])
    
    nom1 = -c3*m2*np.exp(m1*eta_c-m2*eta_c)
    nom2 = -c2*m2*np.exp(m1*eta_c-m2*eta_c)
    nom3 = -c2*m1*np.exp(m1*eta_c-m2*eta_c)
    nom4 = -c3*m1*np.exp(m1*eta_c-m2*eta_c)
    nom5 = -m2*np.exp(-m2*eta_c)
    nom6 = -m1*np.exp(m1*eta_c)
    
    den1 = c2*np.exp(m1*eta_c)
    den2 = c3*np.exp(-m2*eta_c)
    
    nom = c1*(nom1+nom2+nom3+nom4+nom5+nom6)
    den = (1+den1+den2)**2
    
    return nom/den
    

def bv_rate_derivative_oc_a(conc_a, eta_a, Ai_a, rp_a):
    c1 = cnst['j0_a'] * Ai_a * conc_a / cnst['conc_ref_a']
    c2 = cnst['j0_a'] / (cnst['F'] * cnst['conc_ref_a'] * cnst['D_c'] / rp_a)
    c3 = cnst['j0_a'] / (cnst['F'] * cnst['conc_ref_a'] * cnst['D_a'] / rp_a)
    
    m1 = cnst['alpha_a_a'] * cnst['val_a'] * cnst['F'] / (cnst['R_g'] * cnst['T'])
    m2 = cnst['alpha_c_a'] * cnst['val_a'] * cnst['F'] / (cnst['R_g'] * cnst['T'])
    
    nom1 = c2*m1*np.exp(m2*eta_a-m1*eta_a)
    nom2 = c3*m1*np.exp(m2*eta_a-m1*eta_a)
    nom3 = c3*m2*np.exp(m2*eta_a-m1*eta_a)
    nom4 = c2*m2*np.exp(m2*eta_a-m1*eta_a)
    nom5 = m1*np.exp(-m1*eta_a)
    nom6 = m2*np.exp(m2*eta_a)
    
    den1 = c2*np.exp(-m1*eta_a)
    den2 = c3*np.exp(m2*eta_a)
    
    nom = c1*(nom1+nom2+nom3+nom4+nom5+nom6)
    den = (1+den1+den2)**2 

    return nom/den 


def bv_rate_constant_ad_c(eta, Ai_c, rp_c):
    r"""
    Calculate A1 in linear kinetics for AdvectionDiffusion algorithm in the cathode.
    """
    # c = 1.0 is a workaround so that A1 is rate "constant" not the actual rate
    return bv_rate_constant_oc_c(c=1, eta=eta, Ai_c=Ai_c, rp_c=rp_c) / (cnst['F'] * cnst['val_c'])


def bv_rate_constant_ad_a(eta, Ai_a, rp_a):
    r"""
    Calculates A1 in linear kinetics for AdvectionDiffusion algorithm in the anode.
    """
    # c = 1.0 is a workaround so that A1 is rate "constant" not the actual rate
    return bv_rate_constant_oc_a(c=1, eta=eta, Ai_a=Ai_a, rp_a=rp_a) / (cnst['F'] * cnst['val_a'])


def rel_error(current_new, current):
    r"""
    Calculates the relative error of the current estimation with the previous estimation."""
    # Check if it's safe to calculate relative error (division by 0)
    if current != 0.0:
        rel_err = abs((current_new - current) / current)
        
        # Solution has converged, but relative error --> infinity (division by 0)
    elif current_new == current == 0.0:
        rel_err = 0.0
        
        # Solution has not converged, set relative error to an arbitrary high value
    else:
        rel_err = 1e3
    return rel_err


def find_eta_act(eta, cell, i_actual, Ai):
    r'''
    find_eta_act is passed onto the least squares optimization to compute the activation 
    overpotential required for a given current. It is used to separate the contribution
    of concentration and activation overpotential on the cell performance. 
    
    Parameters
    ----------
    eta : Guess for the activation overpotential [V].
    cell : cell type (anode or cathode).
    i_actual : current that the pore should generate.
    Ai : Internal surface area of the considered pore.
    p : Internal pore radius of the considered pore.
            
    Returns
    -------
    err : relative error between the calculated current and i_actual
    for a guessed overpotential 
    '''
    if cell == 'anode':
        # Compute the current that can be generated without concentration effects
        c1 = cnst['j0_a'] * Ai * cnst['conc_in_a'] / cnst['conc_ref_a']
        arg1 = -cnst['alpha_a_a'] * cnst['val_a'] * cnst['F'] / (cnst['R_g'] * cnst['T']) * eta
        arg2 = cnst['alpha_c_a'] * cnst['val_a'] * cnst['F'] / (cnst['R_g'] * cnst['T']) * eta        
        
        ideal_current = c1 * (np.exp(arg1) - np.exp(arg2))
        
        # Calculate the error between the ideal current and the actual current
        err = (ideal_current - i_actual) / i_actual          

    elif cell == 'cathode':
        # Compute the current that can be generated without concentration effects
        c1 = cnst['j0_c'] * Ai * cnst['conc_in_c'] / cnst['conc_ref_c']
        arg1 = cnst['alpha_a_c'] * cnst['val_c'] * cnst['F'] / (cnst['R_g'] * cnst['T']) * eta
        arg2 = -cnst['alpha_c_c'] * cnst['val_c'] * cnst['F'] / (cnst['R_g'] * cnst['T']) * eta
        
        ideal_current = c1 * (np.exp(arg1) - np.exp(arg2))
        
        err = (ideal_current - i_actual) / i_actual
                             
    return err


def ideal_I_c(eta, Ai):
    
     # Compute the current that can be generated without concentration effects
     c1 = cnst['j0_a'] * Ai * cnst['conc_in_a'] / cnst['conc_ref_a']
     arg1 = -cnst['alpha_a_a'] * cnst['val_a'] * cnst['F'] / (cnst['R_g'] * cnst['T']) * eta
     arg2 = cnst['alpha_c_a'] * cnst['val_a'] * cnst['F'] / (cnst['R_g'] * cnst['T']) * eta        
        
     ideal_current = c1 * (np.exp(arg1) - np.exp(arg2))
     
     return ideal_current
        
 
def ideal_I_a(eta, Ai):
    
     # Compute the current that can be generated without concentration effects
     c1 = cnst['j0_c'] * Ai * cnst['conc_in_c'] / cnst['conc_ref_c']
     arg1 = cnst['alpha_a_c'] * cnst['val_c'] * cnst['F'] / (cnst['R_g'] * cnst['T']) * eta
     arg2 = -cnst['alpha_c_c'] * cnst['val_c'] * cnst['F'] / (cnst['R_g'] * cnst['T']) * eta
        
     ideal_current = c1 * (np.exp(arg1) - np.exp(arg2))
     
     return ideal_current