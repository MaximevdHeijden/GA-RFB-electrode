import pdb                                                  
import openpnm as op
import numpy as np
import time
import openpyxl 
import random
from datetime import datetime
import os
import scipy.optimize as optimize
import math
from numpy import pi as _pi
op.Workspace().clear()

def algorithm(net_c, net_a, eta_c, eta_a, generation, Extracted, FF, ReloadData, reload, inputDict, cf, i):
    '''Algorithm() computes the coupled mass and charge transport within both electrodes of a redox flow battery.'''
    
    # Load in the input dictionary
    param = inputDict.input_dict     
    # Retrieve the geometry of the networks
    geom_c = net_c.project['geo_01']
    geom_a = net_a.project['geo_01']
        
    # Define and calculate macro-network properties:
    # Electrode properties. Pores in the cathode and anode:
    Ps_c = net_c.pores()  
    Ps_a = net_a.pores()  
    # Define electrode dimensions [x=0, y=1, z=2]:
    H_dim = param['height_dimension']                   # Sides of the electrode
    L_dim = param['length_dimension']                   # Flow field direction
    W_dim = param['width_dimension']                    # Thickness: current collector -> membrane
    # Calculate the length, height and width of the electrodes:
    L_c = np.amax(net_c['pore.coords'][Ps_c, L_dim])  
    H_c = np.amax(net_c['pore.coords'][Ps_c, H_dim])
    H_c_min = np.amin(net_c['pore.coords'][net_c.pores(), H_dim])
    L_c_min = np.amin(net_c['pore.coords'][net_c.pores(), L_dim])
    W_c_min = np.amin(net_c['pore.coords'][net_c.pores(), W_dim])
    W_c = np.amax(net_c['pore.coords'][Ps_c, W_dim])  
    L_a = np.amax(net_a['pore.coords'][Ps_a, L_dim])
    H_a = np.amax(net_a['pore.coords'][Ps_a, H_dim])
    H_a_min = np.amin(net_a['pore.coords'][net_a.pores(), H_dim])
    W_a = np.amax(net_a['pore.coords'][Ps_a, W_dim])    
    L_a_min = np.amin(net_c['pore.coords'][net_c.pores(), L_dim])
    W_a_min = np.amin(net_c['pore.coords'][net_c.pores(), W_dim])
    
    # Assign the labels 'membrane', 'current collector' 'flow inlet' and 'flow outlet' to the correct boundary pores:
    if FF == 0:
        cf.assign_boundary_poresFF(net_c, W_dim, L_dim)
        cf.assign_boundary_poresFF(net_a, W_dim, L_dim)    
    elif FF == 1:
        cf.assign_boundary_poresIDFF(net_c, W_dim, L_dim, H_dim, H_c, H_c_min, L_c_min, W_c_min, L_c, W_c)
        cf.assign_boundary_poresIDFF(net_a, W_dim, L_dim, H_dim, H_a, H_a_min, L_a_min, W_a_min, L_a, W_a) 
    # Position networks next to each other (left: anode, middle: fictitious membrane, right: cathode)
    # Compute area's of the electrodes:
    A_ext_c = L_c * H_c                                 # Cathode area to compute external current density [m2]
    A_ext_a = L_a * H_a                                 # Anode area [m2]
    mem_area = A_ext_c                                  # Fictitious membrane area [m2]
    if FF == 0:
        A_in_c = W_c * H_c                              # Cathode inlet area [m2]
        A_in_a = W_a * H_a                              # Anode inlet area [m2]
    elif FF == 1:
        A_in_c = W_c * L_c                              # Cathode inlet area [m2]
        A_in_a = W_a * L_a                              # Anode inlet area [m2]
        
    # Invert anodic network in the width to obtain the same pores at the membrane side for both electrodes:
        # NOTE: The boundaries of the network are also reversed in this direction, see cf.assign_boundary_pores().
    inversion_factor = net_a['pore.coords'][:, W_dim] - min(net_a['pore.coords'][:, W_dim])
    net_a['pore.coords'][:, W_dim]  = max(net_a['pore.coords'][:, W_dim]) - inversion_factor 
     
    # Obtain electrochemically active surface area (= geometric internal surface area):
    Ai_c = net_c['pore.surface_area']  
    Ai_a = net_a['pore.surface_area'] 
    # Pore radii:
    rp_c = geom_c['pore.diameter'] / 2  
    rp_a = geom_a['pore.diameter'] / 2  

    # Create phase and physics objects:
    # Anolyte
    anolyte = op.phases.Water(network=net_a, name='anolyte')
    anolyte['pore.electrical_conductivity'] = param['anolyte_conductivity']     
    anolyte['pore.density'] = param['anolyte_density']                          
    anolyte['pore.diffusivity'] = param['D_a']                                  
    anolyte['pore.viscosity'] = param['anolyte_viscosity']                      
    # Catholyte
    catholyte = op.phases.Water(network=net_c, name='catholyte')
    catholyte['pore.electrical_conductivity'] = param['catholyte_conductivity'] 
    catholyte['pore.density'] = param['catholyte_density']                      
    catholyte['pore.diffusivity'] = param['D_c']                                
    catholyte['pore.viscosity'] = param['catholyte_viscosity']                  
    # Physics:
    phys_a = op.physics.Standard(network=net_a, geometry=geom_a, phase=anolyte)     
    phys_c = op.physics.Standard(network=net_c, geometry=geom_c, phase=catholyte)   

    # Regenerate the physics and phase models
    phys_c.regenerate_models()
    phys_a.regenerate_models()
    anolyte.regenerate_models()
    catholyte.regenerate_models()
    
    # Create physics for boundary pores (there should be no resistance to flow or charge transport in the boundary pores):
    if Extracted == 0 or Extracted == 2:    
        boundary_throats_c = net_c.find_neighbor_throats(net_c.pores('surface'))
        boundary_throats_a = net_a.find_neighbor_throats(net_a.pores('surface'))        
        phys_c['throat.ad_dif_conductance'][boundary_throats_c] = 1e4 * phys_c['throat.ad_dif_conductance'].max()
        phys_c['throat.electrical_conductance'][boundary_throats_c] = 1e4 * phys_c['throat.electrical_conductance'].max()        
        phys_a['throat.ad_dif_conductance'][boundary_throats_a] = 1e4 * phys_a['throat.ad_dif_conductance'].max()
        phys_a['throat.electrical_conductance'][boundary_throats_a] = 1e4 * phys_a['throat.electrical_conductance'].max()
    
    # Load in parameters for mass and charge transport:
    # Active species / kinetic parameters:
    conc_in_a = param['conc_in_a']      
    conc_in_c = param['conc_in_c']            
    # Cell potential parameters:
    E_red_a = param['E_red_a']          
    E_red_c = param['E_red_c']          
    E_cell_final = param['E_cell_final']                # Final value of the cell voltage range [V] # (E_cell < 0 --> Charging, E_cell > 0 --> Discharging)    
    # Potential parameters
    E_0 = E_red_c - E_red_a                             # Open circuit voltage [V]
    E_cell_vec = [E_cell_final]                         # Cell voltage range for which the algorithm calculates    
    # Membrane properties:
    R_mem = param['membrane_resistivity']               # Flow-through ionic membrane resistivity [Ohm m2] 
    res_mem = R_mem / mem_area                          # Membrane resistance [Ohm]

    # Run the Stokes Flow algorithms for electrolyte transport:
    sf_c = op.algorithms.StokesFlow(network=net_c, phase=catholyte)
    sf_a = op.algorithms.StokesFlow(network=net_a, phase=anolyte)
    
    # Neumann inlet boundary condition:
    v_in_c = param['catholyte_inlet_velocity']
    v_in_a = param['anolyte_inlet_velocity']     
    Q_in_c = v_in_c * A_in_c                            # total flow rate entering the network [m3/s]
    Q_in_a = v_in_a * A_in_a 
       
    # Find the inlet throats and compute the inlet area:
    inlet_throats_c = net_c.find_neighbor_throats(pores=net_c.pores('flow_inlet'))
    inlet_throats_a = net_a.find_neighbor_throats(pores=net_a.pores('flow_inlet'))    
    total_area_inlet_c = net_c['throat.area'][inlet_throats_c].sum()
    total_area_inlet_a = net_a['throat.area'][inlet_throats_a].sum()

    # Initial guess for the Stokes flow algorithm: 
    # Assign inlet boundary conditions to pores:
    for pore in net_c.pores('flow_inlet'):
        throat = net_c.find_neighbor_throats(pores=pore)[0]             # Boundary pores are connected to one neighboring pore via only one throat
        sf_c.set_rate_BC(values=Q_in_c * net_c['throat.area'][throat] / total_area_inlet_c, pores=pore)     # Neumann inlet boundary condition scaled via the throat connected to the first internal pore
    for pore in net_a.pores('flow_inlet'):
        throat = net_a.find_neighbor_throats(pores=pore)[0] 
        sf_a.set_rate_BC(values=Q_in_a * net_a['throat.area'][throat] / total_area_inlet_a, pores=pore)    
    # Dirichlet outlet boundary condition
    Pout = 0                                            # Pressure outlet boundary condition [Pa]
    sf_c.set_value_BC(values=Pout, pores=net_c.pores('flow_outlet'))    # Dirichlet outlet boundary condition
    sf_a.set_value_BC(values=Pout, pores=net_a.pores('flow_outlet'))    

    # Run electrolyte transport algorithms:
    x, y = sf_c.run()
    # If the algoritgm does not converge, stop the function for this network and give it a bad fitness definition:
    if y == False:
        print(y)
        P_pump = 1000
        P_el = 0
        P_el_ideal = 0
        P_theory = 0
        current_sum = 0
        eta_c = eta_c
        eta_a = eta_a
        return P_pump, P_el, P_el_ideal, P_theory, current_sum, eta_c, eta_a
    x, y = sf_a.run()
    # If the algoritgm does not converge, stop the function for this network and give it a bad fitness definition:
    if y == False:
        print(y)
        P_pump = 1000
        P_el = 0
        P_el_ideal = 0
        P_theory = 0
        current_sum = 0
        eta_c = eta_c
        eta_a = eta_a
        return P_pump, P_el, P_el_ideal, P_theory, current_sum, eta_c, eta_a
    
    # Update phases with calculated pressure field:
    catholyte.update(sf_c.results())
    anolyte.update(sf_a.results())

    # Find the inlet pressure, so that the target flow rate is reached:
    def inlet_pressure_c(P_in, Q_desired, inlet_throats, inlet_pores, net, phase, phys, eta_a, eta_c):
        sf_c.set_value_BC(values=P_in, pores=inlet_pores)           # Dirichlet boundary inlet condition
        x, y = sf_c.run()  
        # If the algoritgm does not converge, stop the function for this network and give it a bad fitness definition:
        if y == False:
            print(y)
            P_pump = 1000
            P_el = 0
            P_el_ideal = 0
            P_theory = 0
            current_sum = 0
            eta_c = eta_c
            eta_a = eta_a
            return P_pump, P_el, P_el_ideal, P_theory, current_sum, eta_c, eta_a
    
        phase.update(sf_c.results())        
        pore_1 = net['throat.conns'][inlet_throats][:, 0]
        pore_2 = net['throat.conns'][inlet_throats][:, 1]        
        delta_P = abs(phase['pore.pressure'][pore_1] - phase['pore.pressure'][pore_2])
        Q_tot = (phys['throat.hydraulic_conductance'][inlet_throats] * delta_P).sum()
        
        return (Q_tot - Q_desired)
    def inlet_pressure_a(P_in, Q_desired, inlet_throats, inlet_pores, net, phase, phys, eta_a, eta_c):
        sf_a.set_value_BC(values=P_in, pores=inlet_pores)           # Dirichlet boundary inlet condition
        x, y = sf_a.run()   
        # If the algoritgm does not converge, stop the function for this network and give it a bad fitness definition:
        if y == False:
            print(y)
            P_pump = 1000
            P_el = 0
            P_el_ideal = 0
            P_theory = 0
            current_sum = 0
            eta_c = eta_c
            eta_a = eta_a
            return P_pump, P_el, P_el_ideal, P_theory, current_sum, eta_c, eta_a
        
        phase.update(sf_a.results())       
        pore_1 = net['throat.conns'][inlet_throats][:, 0]
        pore_2 = net['throat.conns'][inlet_throats][:, 1]        
        delta_P = abs(phase['pore.pressure'][pore_1] - phase['pore.pressure'][pore_2])
        Q_tot = (phys['throat.hydraulic_conductance'][inlet_throats] * delta_P).sum()
        
        return (Q_tot - Q_desired)
    
    # Set initial guesses to the mean pressure found above:
    x0_c = catholyte['pore.pressure'][net_c.pores('flow_inlet')].mean()
    x0_a = anolyte['pore.pressure'][net_a.pores('flow_inlet')].mean()    
    # Delete the initial guess boundary condition:
    sf_c['pore.bc_rate'] = np.nan
    sf_a['pore.bc_rate'] = np.nan
    
    # Find the pressure at the inlet at which the total flow rate matches the desired flow rate:
    optimize.fsolve(func=inlet_pressure_c, x0=x0_c, args=(Q_in_c, inlet_throats_c, net_c.pores('flow_inlet'), net_c, catholyte, phys_c, eta_a, eta_c))
    optimize.fsolve(func=inlet_pressure_a, x0=x0_a, args=(Q_in_a, inlet_throats_a, net_a.pores('flow_inlet'), net_a, anolyte, phys_a, eta_a, eta_c))    
    # The pressure drop is the difference between the pores connected to the boundary pores at the inlet and outlet:
    pores_in_c = net_c.find_neighbor_pores(pores=net_c.pores('flow_inlet'))
    pores_out_c = net_c.find_neighbor_pores(pores=net_c.pores('flow_outlet'))   
    dP_c = catholyte['pore.pressure'][pores_in_c].mean() - catholyte['pore.pressure'][pores_out_c].mean() # Pressure drop cathode    
    pores_in_a = net_a.find_neighbor_pores(pores=net_a.pores('flow_inlet'))
    pores_out_a = net_a.find_neighbor_pores(pores=net_a.pores('flow_outlet'))    
    dP_a = anolyte['pore.pressure'][pores_in_a].mean() - anolyte['pore.pressure'][pores_out_a].mean() # Pressure drop cathode
    
    # Set up the advection-diffusion algorithms for species mass transport:
    ad_a = op.algorithms.AdvectionDiffusion(network=net_a, phase=anolyte)
    ad_c = op.algorithms.AdvectionDiffusion(network=net_c, phase=catholyte)
    
    # Boundary conditions:
        # NOTE: Boundary condition for charge transport at the membrane is a function of the current and can therefore be found in the iteration loop.
        # NOTE: OpenPNM automatically sets all other boundaries to no flux, i.e. dV/dx = 0.
    # Compute ad-dif conductance:
    mod = op.models.physics.ad_dif_conductance.ad_dif
    phys_a.add_model(propname='throat.ad_dif_conductance', model=mod, s_scheme='powerlaw')
    phys_c.add_model(propname='throat.ad_dif_conductance', model=mod, s_scheme='powerlaw')    
    ad_a.set_value_BC(values=conc_in_a, pores=net_a.pores('flow_inlet'))
    ad_c.set_value_BC(values=conc_in_c, pores=net_c.pores('flow_inlet'))    
    ad_a.set_outflow_BC(pores=net_a.pores('flow_outlet'))
    ad_c.set_outflow_BC(pores=net_c.pores('flow_outlet'))    
    
    # Source term (it should always be a 'string' loaded into the model parameters, no 0):
        # NOTE: regen_mode = "deferred" bypasses OpenPNM warning that concentration is yet undefined.
        # NOTE: Function can be linearized later on to improve stability of the algorithm.
        # NOTE: 'pore.current_A2 is also calculated in the iterative loop (== amount of current generated in every pore)
    source_term = op.models.physics.generic_source_term.linear    
    phys_a['pore.rxn_A2'] = 0.0
    phys_c['pore.rxn_A2'] = 0.0    
    anolyte.add_model(propname='pore.butler_volmer', model=source_term, A1='pore.rxn_A1', A2="pore.rxn_A2", X='pore.concentration', regen_mode="deferred")    
    catholyte.add_model(propname='pore.butler_volmer', model=source_term, A1='pore.rxn_A1', A2="pore.rxn_A2",X='pore.concentration', regen_mode="deferred")    
    ad_a.set_source(propname='pore.butler_volmer', pores=net_a.pores('internal'))
    ad_c.set_source(propname='pore.butler_volmer', pores=net_c.pores('internal'))   
    # Set up the ohmic conduction algorithms for ionic charge transport:
    oc_a = op.algorithms.OhmicConduction(network=net_a, phase=anolyte)
    oc_c = op.algorithms.OhmicConduction(network=net_c, phase=catholyte)
    # Source term:
    phys_a['pore.current_catholyte_A1'] = 0.0
    phys_c['pore.current_anolyte_A1'] = 0.0    
    anolyte.add_model(propname='pore.proton_anolyte', model=source_term, A1='pore.current_anolyte_A1', A2='pore.current_anolyte_A2', X='pore.voltage', regen_mode="deferred")    
    catholyte.add_model(propname='pore.proton_catholyte', model=source_term, A1='pore.current_catholyte_A1', A2='pore.current_catholyte_A2', X='pore.voltage', regen_mode="deferred")    
    oc_a.set_source(propname='pore.proton_anolyte', pores=net_a.pores('internal'))
    oc_c.set_source(propname='pore.proton_catholyte', pores=net_c.pores('internal'))

    # Compute the coupled charge and mass transfer equations:
    # Initialize overpotential and concentration guesses:   
    if ReloadData == 0:
        if generation == 0  and i == 0:
            eta0_a = 0
            eta0_c = 0
        else:
            eta0_a = eta_a                              # Initial guess for overpotential
            eta0_c = eta_c  
    elif ReloadData == 1:
        if generation == reload:
            eta0_a = 0
            eta0_c = 0
        else:
            eta0_a = eta_a                              # Initial guess for overpotential
            eta0_c = eta_c       
    
    conc0_a = conc_in_a                                 # Initialize concentration in pores to inlet concentration
    conc0_c = conc_in_a       
    rel_tol = param['rel_tol']      
    abs_tol = param['abs_tol']      
    max_iter = param['max_iter']    
    omega = param['omega']         
    
    # Loop over the different cell voltages to obtain the polarization curve:
    for E_cell in E_cell_vec:
        # Update guesses for the overpotential and concentration based on the previous result
        eta_c = eta0_c  
        eta_a = eta0_a        
        conc_c = conc0_c 
        conc_a = conc0_a        
        # The current estimation is reset, to reset the relative error
        current_ad_a = 0.0 
        current_oc_a = 0.0
        current_ad_c = 0.0
        current_oc_c = 0.0    
        # The solid cathode operates at a potential of Vcell, the anode at 0 V:
        V_c = E_cell - E_0 - eta_c                      # Voltage in the liquid phase [V]
        V_a = 0 - eta_a        
        anolyte['pore.voltage'] = V_a
        catholyte['pore.voltage'] = V_c
        
        # Iterative loop to calculate coupled charge and mass transfer equations:
        for itr in range(max_iter):                            
            # Cathodic half cell:
            # Update ad reaction term with new values for eta_c:
            phys_c['pore.rxn_A1'] = cf.bv_rate_constant_ad_c(eta_c, Ai_c, rp_c)    
            # Solve the advection-diffusion-reaction equation (concentration field):
            x, y = ad_c.run()   
            # If the algoritgm does not converge, stop the function for this network and give it a bad fitness definition:
            if y == False:
                print(y)
                P_pump = 1000
                P_el = 0
                P_el_ideal = 0
                P_theory = 0
                current_sum = 0
                eta_c = eta_c
                eta_a = eta_a
                return P_pump, P_el, P_el_ideal, P_theory, current_sum, eta_c, eta_a
            # Update the value for the concentration:
            conc_new_c = ad_c['pore.concentration']
            conc_c = conc_new_c * omega + conc_c * (1 - omega)
            catholyte['pore.concentration'] = conc_c
            
            # Compute current according to the reaction of species:
                # NOTE: we only use the internal pores for the relative error calculation, as both the ohmic conduction and advection-diffusion algorithm use different boundary pores.
            current_estimation_ad_c = cf.bv_rate_constant_oc_c(conc_c, eta_c, Ai_c, rp_c)[net_c.pores('internal')].sum()           
            # Compute oc reaction term with new values for conc_c:
            drdv = cf.bv_rate_derivative_oc_c(conc_c, eta_c, Ai_c, rp_c)
            r = cf.bv_rate_constant_oc_c(conc_c, eta_c, Ai_c, rp_c)                
            phys_c['pore.current_catholyte_A1'] = drdv
            phys_c['pore.current_catholyte_A2'] = r - drdv * V_c
            
            # Couple the potential field in both half cells via a continuum averaged membrane approach:
                # NOTE: This only works when two identical networks are used (networks with identical membrane pores)
            V_membrane_pores_a = anolyte['pore.voltage'][net_a.pores('membrane') ]
            V_membrane_pores_c = V_membrane_pores_a + res_mem * current_estimation_ad_c        
            # Update the membrane boundary condition for the cathode:
            oc_c.set_value_BC(values=V_membrane_pores_c, pores=net_c.pores('membrane'))
            
            # Solve the ohmic conduction equation (potential field) for the liquid electrolyte in the porous cathode
            x, y = oc_c.run()
            # If the algoritgm does not converge, stop the function for this network and give it a bad fitness definition:
            if y == False:
                print(y)
                P_pump = 1000
                P_el = 0
                P_el_ideal = 0
                P_theory = 0
                current_sum = 0
                eta_c = eta_c
                eta_a = eta_a
                return P_pump, P_el, P_el_ideal, P_theory, current_sum, eta_c, eta_a
            # Update the value for the liquid potential:
            V_new_c = oc_c['pore.voltage']
            V_c = V_new_c * omega + V_c * (1 - omega)
            catholyte['pore.voltage'] = V_c            
            # Compute new cathodic overpotential:
            eta_c = E_cell - V_c - E_0            
            # Compute current according to the transport of species:
            current_estimation_oc_c = cf.bv_rate_constant_oc_c(conc_c, eta_c, Ai_c, rp_c)[net_c.pores('internal')].sum()
            
            # Anodic half cell
            # Update ad reaction term with new values for eta_a:
            phys_a['pore.rxn_A1'] = cf.bv_rate_constant_ad_a(eta_a, Ai_a, rp_a)    
            # Solve the advection-diffusion-reaction equation (concentration field)
            x, y = ad_a.run()   
            # If the algoritgm does not converge, stop the function for this network and give it a bad fitness definition:
            if y == False:
                print(y)
                P_pump = 1000
                P_el = 0
                P_el_ideal = 0
                P_theory = 0
                current_sum = 0
                eta_c = eta_c
                eta_a = eta_a
                return P_pump, P_el, P_el_ideal, P_theory, current_sum, eta_c, eta_a
            # Update the value for the concentration
            conc_new_a = ad_a['pore.concentration']
            conc_a = conc_new_a * omega + conc_a * (1 - omega)
            anolyte['pore.concentration'] = conc_a
            
            # Compute current according to the reaction of species
            current_estimation_ad_a = cf.bv_rate_constant_oc_a(conc_a, eta_a, Ai_a, rp_a)[net_a.pores('internal')].sum()            
            # Compute oc (linearized) reaction term with new values for conc_c
                # NOTE: the current generation is a source instead of a sink term, so minus bv_rate_constant_oc_anode
            drdv = cf.bv_rate_derivative_oc_c(conc_a, eta_a, Ai_a, rp_a)
            r = -cf.bv_rate_constant_oc_a(conc_a, eta_a, Ai_a, rp_a)                
            phys_a['pore.current_anolyte_A1'] = drdv
            phys_a['pore.current_anolyte_A2'] = r - drdv * V_a
                            
            # Couple the potential field in both half cells via a continuum averaged membrane approach
            V_membrane_pores_c = catholyte['pore.voltage'][net_c.pores('membrane')]
            V_membrane_pores_a = V_membrane_pores_c - res_mem * current_estimation_ad_a    
            # Update the membrane boundary condition for the anode
            oc_a.set_value_BC(values=V_membrane_pores_a, pores=net_a.pores('membrane'))
    
            # Solve the ohmic conduction equation (potential field) for the liquid electrolyte in the porous cathode
            x, y = oc_a.run()     
            # If the algoritgm does not converge, stop the function for this network and give it a bad fitness definition:
            if y == False:
                print(y)
                P_pump = 1000
                P_el = 0
                P_el_ideal = 0
                P_theory = 0
                current_sum = 0
                eta_c = eta_c
                eta_a = eta_a
                return P_pump, P_el, P_el_ideal, P_theory, current_sum, eta_c, eta_a
            # Update voltage
            V_new_a = oc_a['pore.voltage']
            V_a = V_new_a * omega + V_a * (1 - omega)
            anolyte['pore.voltage'] = V_a           
            # Compute new anodic overpotentials
            eta_a = E_0 - V_a                
            # Compute current according to the transport of species
            current_estimation_oc_a = cf.bv_rate_constant_oc_a(conc_a, eta_a, Ai_a, rp_a)[net_a.pores('internal')].sum()
            
            # Convergence and updating
            # Calculate the relative error with the previous solution
            rel_err_current_ad_a = cf.rel_error(current_estimation_ad_a, current_ad_a)
            rel_err_current_oc_a = cf.rel_error(current_estimation_oc_a, current_oc_a)
            rel_err_current_ad_c = cf.rel_error(current_estimation_ad_c, current_ad_c)
            rel_err_current_oc_c = cf.rel_error(current_estimation_oc_c, current_oc_c)                
            # Store the previous solution
            current_ad_a = current_estimation_ad_a
            current_oc_a = current_estimation_oc_a
            current_ad_c = current_estimation_ad_c
            current_oc_c = current_estimation_oc_c    
            # Calculate the absolute error between the current density found in the anodic and cathodic compartment in [A cm-2] (1e4 = conversion m2 to cm2)
            abs_err_current_cat_an = abs((current_oc_a / A_ext_a / 1e4 - current_oc_c / A_ext_c / 1e4))    
            # Check if found progression of the solution is within tolerances
            convergence_ad_a = rel_err_current_ad_a < rel_tol
            convergence_oc_a = rel_err_current_oc_a < rel_tol
            convergence_ad_c = rel_err_current_ad_c < rel_tol
            convergence_oc_c = rel_err_current_oc_c < rel_tol
            convergence_a_c = abs_err_current_cat_an < abs_tol
                     
            # Convergence criterion:
            if (convergence_ad_a and convergence_oc_a and convergence_ad_c and convergence_oc_c and convergence_a_c) or itr == max_iter-1:
                #print(f"CONVERGED at {round(E_cell,2)} [V] with {current_estimation_ad_c/A_ext_c/1e4:.5f} [A/cm2]!")                            
                # Post-treatment:
                # Current in every pore:
                anolyte['pore.current'] = cf.bv_rate_constant_oc_a(conc_a, eta_a, Ai_a, rp_a)
                catholyte['pore.current'] = cf.bv_rate_constant_oc_c(conc_c, eta_c, Ai_c, rp_c)               
                # Ohmic overpotential contribution:
                catholyte['pore.ohmic_overpotential'] = V_c - V_membrane_pores_c.mean()
                anolyte['pore.ohmic_overpotential'] = V_a  - V_membrane_pores_a.mean()
                eta_ohm_mem = V_membrane_pores_c.mean()-V_membrane_pores_a.mean()                
                # Activation overpotential contribution:
                eta_act_c = np.zeros(net_c.Np)
                eta_act_a = np.zeros(net_a.Np)                    
                for pore in range(net_c.Np):
                    eta_act_c[pore] = optimize.fsolve(func = cf.find_eta_act, x0=eta_c[pore], args=('cathode', catholyte['pore.current'][pore], Ai_c[pore])) 
                for pore in range(net_a.Np):
                    eta_act_a[pore] = optimize.fsolve(func = cf.find_eta_act, x0=eta_a[pore], args=('anode', anolyte['pore.current'][pore], Ai_a[pore]))                 
                catholyte['pore.activation_overpotential'] = eta_act_c
                anolyte['pore.activation_overpotential'] = eta_act_a                   
                # Concentration overpotential contribution
                catholyte['pore.concentration_overpotential'] = eta_c-catholyte['pore.activation_overpotential'] 
                anolyte['pore.concentration_overpotential'] = eta_a-anolyte['pore.activation_overpotential']
                                       
                # Compute total current:
                    # NOTE: Boundary pores are excluded, since these are artificial pores added by the SNOW algorithm, that are not part of the real network
                current_sum_a = anolyte['pore.current'][net_a.pores('internal')].sum()
                current_sum_c = catholyte['pore.current'][net_c.pores('internal')].sum()
                current_sum = (abs(current_sum_a)+abs(current_sum_c))/2               
                # Calculate maximum achieveable (ideal) current 
                    # NOTE: In the symmetrical flow cell activation potential defines maximum current
                catholyte['pore.id_current'] = cf.ideal_I_c(eta_c, Ai_c)
                anolyte['pore.id_current'] = cf.ideal_I_a(eta_a, Ai_a)                              
                ideal_I_c_sum = catholyte['pore.id_current'][net_c.pores('internal')].sum()
                ideal_I_a_sum = anolyte['pore.id_current'][net_a.pores('internal')].sum()                
                ideal_I = (abs(ideal_I_c_sum) + abs(ideal_I_a_sum))/2
                                
                # Update initial guesses for next cell voltage:
                eta_0a = eta_a
                eta_0c = eta_c
                conc_0a = conc_a
                conc_0c = conc_c
                break
    
    # Calculate parameters necessary for the fitness definition:
    P_pump_ideal = dP_c*Q_in_c
    P_el = abs(E_cell * current_sum)
    P_el_ideal = E_cell * ideal_I
    P_theory = 1.26 * current_sum

    eta_c = eta_c.mean()
    eta_a = eta_a.mean()
    
    return P_pump_ideal, P_el, P_el_ideal, P_theory, current_sum, eta_c, eta_a
    

def determine_fitness(P_pump_ideal, P_el, P_el_ideal, P_theory, eta_pump = 0.9):
    '''Fitness(net) evaluates the fitness of a given network by comparing the obtained power output with the required pumping power '''
    P_pump = P_pump_ideal/eta_pump                          # Real required pumping power [W]
    fitness = (P_theory - P_el)/(P_theory + P_pump)         # Fitness of the network [-]
    
    return fitness


def initialize_population(population_size, network_shape, spacing, connectivity, network_shape_2, num_points, ws, throat_condition, electrode_name, mutation_chance, mutation_range, ref_pore_diameter, Extracted, folder):
    '''Initalize_population initalizes a set of random generately pore size distributions for the supplied network_template.
    A random pore value for each pore is chosen between 0 < value < maximum allowed pore size'''    
    # Create empty dictionaries to store networks and geometries
    net = {}
    geo = {}
    proj = {}
           
    if Extracted == 2:
        # Create a Voronoi network:
        net_init=op.network.Voronoi(shape=network_shape_2, num_points=num_points) 
        op.utils.Workspace().save_project(project=net_init.project, filename='./Genetic_algorithm' + '/' + folder + '/Voronoi_initial')
        op.utils.Workspace().close_project(net_init.project)
    
    # Load and assign respectively every imported network to a project 
    for i in range(population_size):
        if Extracted == 1:
            proj [i] = op.utils.Workspace().load_project(filename=electrode_name)
            # Assign the imported network to a project
            net[i] = proj[i].network
            # Assign geometry to the imported network
            geo[i] = proj[i]['geo_01']    
        if Extracted == 0:
            # Create a new project
            proj[i] = ws.new_project(name='prj_' +str(i))
            # Load the template
            net[i] = op.network.Cubic(shape=network_shape, spacing=spacing, connectivity=connectivity, project=proj[i])  
        if Extracted == 2:
            # Load the Voronoi network:
            proj[i] = op.utils.Workspace().load_project(filename='./Genetic_algorithm' + '/' + folder + '/Voronoi_initial')
            net[i] = proj[i].network             
            
        # Assign every key to a newly created project    
        if Extracted == 1:
            # Create boundary pores
            net[i]['pore.surface']=[False]*len(net[i]['pore.coords'])
            net[i]['pore.surface'][np.where(net[i]['pore.internal']==False)[0]]=True
            net[i]['throat.surface']=[False]*len(net[i]['throat.all'])
            net[i]['throat.surface'][np.where(net[i]['throat.internal']==False)[0]]=True 
        else:
            # Create boundary pores
            net[i]['pore.internal'][net[i].pores('surface')] = False
            net[i]['throat.internal'][net[i].throats('surface')] = False
        if Extracted == 0:
            # Delete all pores that are on the corners, so that you have no shared inflow and ionic conductivity boundary pores. This assumes that ionic transport takes place from left -> right. And hydraulic transport from front -> back        
            op.topotools.trim(network=net[i], throats=net[i].throats('surface'))
            health = net[i].check_network_health()
            op.topotools.trim(network=net[i], pores=health['trim_pores'])
            health = net[i].check_network_health()
            # Add geometry for internal pores and boundary pores
            geo[i] = op.geometry.GenericGeometry(network=net[i], pores=net[i].Ps, throats=net[i].Ts, project=proj[i])   
            geo[i].add_model(propname='pore.max_size', model=op.models.geometry.pore_size.largest_sphere,iters=10)
            geo[i].add_model(propname='pore.seed', model=op.models.misc.random, element='pore', num_range=[0.2, 0.7], seed=None) 
            geo[i].add_model(propname='pore.diameter', model=op.models.misc.product, prop1='pore.max_size', prop2='pore.seed')
        if Extracted == 2:
            # Add geometry for internal pores and boundary pores
            geo[i] = op.geometry.GenericGeometry(network=net[i], pores=net[i].Ps, throats=net[i].Ts, project=proj[i])  
            # geo[i]['pore.start_size'] = ref_pore_diameter
            geo[i].add_model(propname='pore.seed', model=op.models.misc.random, element='pore', num_range=[0.2, 0.7], seed=None) 
            geo[i].add_model(propname='pore.max_size', model=op.models.geometry.pore_size.largest_sphere,iters=10)
            geo[i].add_model(propname='pore.diameter', model=op.models.misc.product, prop1='pore.max_size', prop2='pore.seed')                               
    
    # Add the remaining geometry models:
    geo = add_geometry_models(net=net, geo=geo, throat_condition=throat_condition, spacing=spacing, Extracted=Extracted)
    
    if Extracted == 1:
        # Define the surface area of the extracted network with the updated geometry and calculate the surface area scaling factor:        
        # Perform mutation on every parent network        
        net, geo = perform_mutation(net=net, geo=geo, mutation_chance=mutation_chance, mutation_range=mutation_range, Extracted=Extracted)
            
    return net, geo, proj
  
    
def select_mating_pool(fitness_of_generation, num_parents, net, geo, proj):
    '''select_mating_pool selects the parents that are most fit for reproduction for the next generation of networks             
    Returns: parent_networks  Dictionairy with the networks used as parents for the next generation.'''    
    # Find the indices of the maximum values for the fitness of the current generation
    ordered_network_list = fitness_of_generation.argsort()          # Orders the indices of the fitness list from small to large.
    parent_indices = ordered_network_list[-num_parents:][::-1]      # Take 'num_parents' amount of networks that have the best fitness value and order the list from best -> worst.    
    parent_networks = {}
    parent_geometries = {}
    parent_projects = {}
    for i, parent in enumerate(parent_indices):
        parent_networks[i] = net[parent] 
        parent_geometries[i] = geo[parent]
        parent_projects[i] = proj[parent]
    
    return parent_networks, parent_geometries, parent_projects


def perform_crossover(parent_networks, offspring_size, generation, Extracted, folder):
    '''perform_crossover performs the crossover of parents to form a new pore size distribution in the 'child'. For each pore, the pore size is randomly chosen as one of the pore sizes of the parent network. 
    Returns:
    new_net : Dictionary of child networks that are produced
    new_geo : Dictionary of geometries that belong to the child networks
    The new networks only contain the new value for the pore diameter. Other geometry values are filled in after possible mutations in perform_mutation().'''
    # Create empty dictionaries for the child networks and geometries.
    new_net = {}
    new_geo = {}
    new_proj = {}
    
    # Select two parents that will mate 
    for i in range(offspring_size):     
        # Load individuals:
        new_proj[i] = op.utils.Workspace().load_project(filename='./Genetic_algorithm' + '/' + folder + '/j0_1_generation_' + str(generation) + '_individual_' + str(i))
        new_net[i] = new_proj[i].network
        new_geo[i] = new_proj[i]['geo_01']
        # Select two random parents from the parents list
        parent_1_int = np.random.randint(low=0,high=len(parent_networks))
        parent_2_int = parent_1_int        
        while parent_2_int == parent_1_int:
            parent_2_int = np.random.randint(low=0,high=len(parent_networks))        
        parent_1 = parent_networks[parent_1_int]
        parent_2 = parent_networks[parent_2_int]                
        # Select the locus for single-point crossover
        crossover_point = np.random.randint(low=1,high=len(new_net[i].pores('internal')))
        
        # Construct new pore diameter array
        if Extracted == 0 or Extracted == 1:
            pore_diameter = np.zeros(len(new_net[i].pores('internal'))) 
            # Inheriting of pore diameters:
            pore_diameter[:crossover_point] = parent_1['pore.diameter'][parent_1.pores('internal')][:crossover_point]
            pore_diameter[crossover_point:] = parent_2['pore.diameter'][parent_2.pores('internal')][crossover_point:] 
            if Extracted == 0:
                new_geo[i]['pore.max_size'] = parent_1['pore.max_size']
            pore_diameter_total = np.zeros(new_net[i].Np)
            pore_diameter_total[new_net[i].pores('surface')] = parent_networks[0]['pore.diameter'][parent_networks[0].pores('surface')]
            pore_diameter_total[new_net[i].pores('internal')] = pore_diameter
            new_geo[i]['pore.diameter'] = pore_diameter_total
        elif Extracted == 2:
            pore_diameter = np.zeros(len(new_net[i].pores('all'))) 
            # Inheriting of pore diameters:
            pore_diameter[:crossover_point] = parent_1['pore.diameter'][parent_1.pores('all')][:crossover_point]
            pore_diameter[crossover_point:] = parent_2['pore.diameter'][parent_2.pores('all')][crossover_point:]         
            pore_diameter_total = np.zeros(new_net[i].Np)        
            pore_diameter_total[new_net[i].pores('all')] = pore_diameter        
            new_geo[i]['pore.diameter'] = pore_diameter_total
       
    return new_net, new_geo    


def perform_crossover2(parent_networks, offspring_size, spacing, generation, Extracted, folder):
    '''perform_crossover performs the crossover of parents to form a new pore size distribution in the 'child'. For each pore, the pore size is randomly chosen as one of the pore sizes of the parent network. 
    Returns:
    new_net : Dictionary of child networks that are produced
    new_geo : Dictionary of geometries that belong to the child networks
    The new networks only contain the new value for the pore diameter. Other geometry values are filled in after possible mutations in perform_mutation().'''  
    # Create empty dictionaries for the child networks and geometries.
    new_net = {}
    new_geo = {}
    new_proj={}
    parent_copies = {}

    # Select two parents that will mate
    for i in range(offspring_size):
        op.Workspace().clear()
        ws = op.Workspace()
        for z in range(len(parent_networks)):
            op.utils.Workspace().load_project(filename='./Genetic_algorithm' + '/' + folder + '/j0_1_generation_' + str(generation) + '_parent_' + str(z), overwrite=True)
        keys = list(ws.keys())
        for z in range(len(parent_networks)):
            new_proj[z]=ws[keys[z]]
            parent_copies[z] = new_proj[z].network
            new_geo[z]=new_proj[z]['geo_01']
    
        # Select two random parents from the parents list
        parent_1_int = np.random.randint(low=0,high=len(parent_copies))
        parent_2_int = parent_1_int        
        while parent_2_int == parent_1_int:
            parent_2_int = np.random.randint(low=0,high=len(parent_copies))        
        parent_1 = parent_copies[parent_1_int]
        parent_geo_1 = new_geo[parent_1_int]
        parent_proj_1 = new_proj[parent_1_int]
        parent_2 = parent_copies[parent_2_int]
        parent_geo_2 = new_geo[parent_2_int]
        parent_proj_2 = new_proj[parent_2_int]
        
        # Select the locus for single-point crossover
        plane = 1
        
        crossover_coordinate = np.mean(parent_1['pore.coords'][:,plane])  
        parent_1['pore.to_be_removed']=len(parent_1['pore.coords'])*[False]
        parent_1['pore.to_be_removed'][np.where(parent_1['pore.coords'][:,plane] >= crossover_coordinate)] = True
        parent_2['pore.to_be_removed']=len(parent_2['pore.coords'])*[False]
        parent_2['pore.to_be_removed'][np.where(parent_2['pore.coords'][:,plane] < crossover_coordinate)] = True        
        parent_1['pore.stitch_1']=len(parent_1['pore.coords'])*[False]
        parent_2['pore.stitch_2']=len(parent_2['pore.coords'])*[False]       
        bond_count = np.zeros((len(np.where(parent_1['pore.to_be_removed']==False)[0])))
        bond_pore =np.zeros((len(np.where(parent_1['pore.to_be_removed']==False)[0]),3))
        
        for k in range(len(np.where(parent_1['pore.to_be_removed'] == False)[0])):
            neighbors=parent_1.find_neighbor_pores(pores=(np.where(parent_1['pore.to_be_removed'] == False)[0][k]))
            if any(parent_1['pore.to_be_removed'][neighbors]==True):
                parent_1['pore.stitch_1'][np.where(parent_1['pore.to_be_removed']==False)[0][k]]=  True
                bond_count[k] =sum(parent_1['pore.to_be_removed'][neighbors])
                bond_pore[k,:] = parent_1['pore.coords'][np.where(parent_1['pore.to_be_removed'] == False)[0][k]] 
            else: 
                continue
            del neighbors        
        for k in range(len(np.where(parent_2['pore.to_be_removed'] == False)[0])):
            neighbors=parent_2.find_neighbor_pores(pores=(np.where(parent_2['pore.to_be_removed'] == False)[0][k]))
            if any(parent_2['pore.to_be_removed'][neighbors]==True):
                parent_2['pore.stitch_2'][np.where(parent_2['pore.to_be_removed']==False)[0][k]]=  True
            else: 
                continue
            del neighbors        
        op.topotools.trim(network=parent_1, pores=np.where(parent_1['pore.to_be_removed'] == True)[0])
        op.topotools.trim(network=parent_2, pores=np.where(parent_2['pore.to_be_removed'] == True)[0])               
        op.topotools.stitch(network=parent_1, donor=parent_2, P_network=parent_1.pores('stitch_1'), P_donor=parent_2.pores('stitch_2'), method='nearest')            
        bond_count = [i for i in bond_count if i != 0]
        bond_count = ([int(bond_count) for bond_count in bond_count])
        bond_pore = bond_pore[~np.all(bond_pore == 0, axis=1)]
        
        for k in range(len(bond_pore)):
            if bond_count[k] > 1:
                near_pores = parent_1.find_nearby_pores(pores=np.where((parent_1['pore.coords'][:, None] == bond_pore[k]).all(-1).any(-1) == True)[0][0], r=parent_1['pore.diameter'].max()*2)[0]
                near_pores_ind = np.where(parent_1['pore.geo_02'][near_pores]==True)[0]
                near_pores = near_pores[near_pores_ind]
                pore_dist = ((bond_pore[k,0]-parent_1['pore.coords'][near_pores][:,0])**2+(bond_pore[k,1]-parent_1['pore.coords'][near_pores][:,1])**2+ (bond_pore[k,2]-parent_1['pore.coords'][near_pores][:,2])**2)**(1/2.0)
                pore_dist_key = np.argsort(pore_dist)
                for bonds in range(len(pore_dist_key)-1):
                    if (bonds+1) < (bond_count[k]):
                        op.topotools.connect_pores(network=parent_1, pores1=np.where((parent_1['pore.coords'][:, None] == bond_pore[k]).all(-1).any(-1) == True), pores2=near_pores[pore_dist_key[bonds+1]], labels=['connected'])        
        geo_info = {}
        geo_keys = list(parent_geo_2.keys())
        geo_plen = len(parent_geo_2['pore.all'])
        geo_tlen = len(parent_geo_2['throat.all'])       
        for k in range(len(geo_keys)):
            geo_info[k] = parent_geo_2[geo_keys[k]]
        
        # Add elements from 2nd to the 1st geometry
        parent_geo_2._set_locations(element='pores', indices=np.where(parent_1['pore.geo_01'] == False)[0], mode='drop')
        parent_geo_1._set_locations(element='pores', indices=np.where(parent_1['pore.geo_01'] == False)[0], mode='add')        
        
        if Extracted == 0 or Extracted == 2:
            for k in range(0,7):
                parent_geo_1[geo_keys[k]][-geo_plen:] = geo_info[k]            
            parent_geo_2._set_locations(element='throats', indices=np.where(parent_1['throat.geo_01'] == False)[0], mode='drop')
            parent_geo_1._set_locations(element='throats', indices=np.where(parent_1['throat.geo_01'] == False)[0], mode='add')        
            for k in range(8,19):
                parent_geo_1[geo_keys[k]][-geo_tlen:] = geo_info[k]
        elif Extracted == 1:
            for k in range(0,10):
                parent_geo_1[geo_keys[k]][-geo_plen:] = geo_info[k]            
            parent_geo_2._set_locations(element='throats', indices=np.where(parent_1['throat.geo_01'] == False)[0], mode='drop')
            parent_geo_1._set_locations(element='throats', indices=np.where(parent_1['throat.geo_01'] == False)[0], mode='add')        
            for k in range(11,29):
                parent_geo_1[geo_keys[k]][-geo_tlen:] = geo_info[k]  
        
        if Extracted == 0 or Extracted == 2:
               parent_geo_1.regenerate_models(['pore.volume', 'throat.max_size', 'throat.length', 'throat.endpoints', 'throat.cross_sectional_area',
                                           'throat.conduit_lengths', 'throat.area', 'throat.volume','throat.surface_area', 'pore.area', 'pore.surface_area'])
        else:
               parent_geo_1.regenerate_models(['pore.volume', 'pore.area', 'throat.length', 'throat.volume', 'throat.endpoints', 'throat.centroid', 'throat.conduit_lengths', 'throat.cross_sectional_area',
                                           'throat.area', 'throat.surface_area', 'pore.surface_area',  'throat.max_size', 'pore.equivalent_diameter', 'pore.centroid', 'pore.label', 
                                           'pore.region_volume', 'throat.total_length', 'throat.equivalent_diameter', 'throat.inscribed_diameter', 'throat.direct_length', 'throat.perimeter',
                                           'pore.extended_diameter', 'pore.inscribed_diameter'])      
               parent_geo_1['pore.centroid'] = parent_geo_1['pore.coords'] 
               
        new_net[i] = parent_1
        new_geo[i] = parent_geo_1        
        parent_proj_1.clear(objtype=['phase', 'algorithm', 'physics'])
        parent_1['throat.internal'][np.where(parent_1['throat.stitched']==True)[0]] = True
        del(parent_proj_1[2], parent_1['throat.geo_02'], parent_1['pore.geo_02'], parent_1['pore.stitch_1'], parent_1['pore.stitch_2'], parent_1['throat.stitched'], parent_1['pore.to_be_removed'])
        
        health = parent_1.check_network_health()
        health_keys=list(health.keys())
        if bool(health[health_keys[2]]) == True:
            op.topotools.trim(network=parent_1, pores=health['trim_pores'])
        #if bool(health[health_keys[3]]) == True:
        #    op.topotools.trim(network=parent_1, throats=health['duplicate_throats'])        
        op.utils.Workspace().save_project(project=parent_1.project, filename='./Genetic_algorithm' + '/' + folder + '/j0_1_generation_' + str(generation) + '_parent1_' + str(i))       
    op.Workspace().clear()
    
    for i in range(offspring_size):
        op.utils.Workspace().load_project(filename='./Genetic_algorithm' + '/' + folder + '/j0_1_generation_' + str(generation) + '_parent1_' + str(i), overwrite=False)
    keys = list(ws.keys())
    for i in range(offspring_size):
            new_proj[i]=ws[keys[i]]
            new_net[i] = new_proj[i].network
            new_geo[i]=new_proj[i]['geo_01']
            os.remove('./Genetic_algorithm' + '/' + folder + "/j0_1_generation_" + str(generation) + "_parent1_" + str(i) + ".pnm")     
    
    return new_net, new_geo   


def perform_mutation(net, geo, mutation_chance, mutation_range, Extracted):
    '''perform_mutation performs the mutation operation on the new generation of networks.
    For every pore, a chance of size 'mutation_chance' is allowed to deviate the pore size with a random value that lies between 1-mutation_range and 1+mutation_range.
    Returns: net:  The dictionairy containing all the networks. & geo:  Dictionary of geometries that belong to the networks.'''
    for i in range(len(net)):        
        # Empty copied array for pore size
        pore_diameter = geo[i]['pore.diameter'].copy()
        geo[i].pop('pore.diameter')      
        for pore in net[i].pores('internal'):
            # Generate a random number and let a mutation occur when it is smaller than a set mutation chance.
            if np.random.random() <= mutation_chance:
                # Draw a random number between the mutation limits of +- mutation_fraction % above or below the current pore diameter
                mutation_factor = np.random.uniform((1-mutation_range), (1+mutation_range))                
                # Perform mutation on pore
                pore_diameter[pore] = mutation_factor * pore_diameter[pore]                
                geo[i]['pore.diameter'] = pore_diameter
                # Check if the new pore diameter is allowed (within bounds)
                correction_factor = 0.999   

                def overlap_check(net, geo):
                # Define overlap as False
                    overlap = True
                    while overlap == True:   
                        overlap = False  
                        # Find the neighboring pores:
                        neighbors_overlap = net.find_nearby_pores(pores=pore, r=geo['pore.diameter'][pore], flatten=True)
                        # Find if the new pore overlaps with existing pores and otherwise delete the location of the new pore:
                        for i in neighbors_overlap:
                          # Find the coordinates of the new pores and the neighboring pores:
                          coords_new = net['pore.coords'][pore]
                          coords = net['pore.coords'][i]
                          # Calculate the center-to-center distance:
                          center_to_center = np.sqrt((coords_new[0]-coords[0])**2 + (coords_new[1]-coords[1])**2 + (coords_new[2]-coords[2])**2)
                          # Check if the new pore overlaps a neighboring pore:
                          if center_to_center < (geo['pore.diameter'][pore]/2 + geo['pore.diameter'][i]/2):
                              pore_diameter[pore] = correction_factor*net['pore.diameter'][pore] 
                              overlap = True
                              break
                          else:
                              overlap = False
                    return net, geo 
                    
                # Check if the new pore overlaps a neighboring pore:
                if mutation_factor > 1:
                    net[i], geo[i] = overlap_check(net=net[i], geo=geo[i])

    return net, geo
                

def split_merge(net, geo, spacing, throat_condition, M_S_chance, M_S_condition, Extracted, ref_network, inputDict): 
    """This function creates two layers of boundaries at the side of the current collector and the membrane side respectively. 
    At the current collector pores with smaller diameter than the predetermined value should merge together, 
    while at the membrane pores with larger diameter than the predetermined value should split."""        
    # Load in the input dictionary
    param = inputDict.input_dict  
    # Electrode properties. Pores in the cathode and anode:
    Ps = net[0].pores()  
    # Define electrode dimensions [x=0, y=1, z=2]:
    if Extracted == 0 or Extracted == 2:
        H_dim = param['height_dimension']                    # Sides of the electrode
        L_dim = param['length_dimension']                    # Flow field direction
        W_dim = param['width_dimension']                     # Thickness: current collector -> membrane
        # Calculate the length, height and width of the electrodes:
        L_max = np.amax(ref_network['pore.coords'][np.where(ref_network['pore.internal'])][:, L_dim]) 
        H_max = np.amax(ref_network['pore.coords'][np.where(ref_network['pore.internal'])][:, H_dim]) 
        W_max = np.amax(ref_network['pore.coords'][np.where(ref_network['pore.internal'])][:, W_dim]) 
        L_min = np.amin(ref_network['pore.coords'][np.where(ref_network['pore.internal'])][:, L_dim]) 
        H_min = np.amin(ref_network['pore.coords'][np.where(ref_network['pore.internal'])][:, H_dim])
        W_min = np.amin(ref_network['pore.coords'][np.where(ref_network['pore.internal'])][:, W_dim]) 
    elif Extracted == 1:
        H_dim = param['height_dimension']                    # Sides of the electrode
        L_dim = param['length_dimension']                    # Flow field direction
        W_dim = param['width_dimension']                     # Thickness: current collector -> membrane
        # Calculate the length, height and width of the electrodes:
        H_max = np.amax(ref_network['pore.coords'][np.where(ref_network['pore.internal'])][:, L_dim]) 
        W_max = np.amax(ref_network['pore.coords'][np.where(ref_network['pore.internal'])][:, H_dim]) 
        L_max = np.amax(ref_network['pore.coords'][np.where(ref_network['pore.internal'])][:, W_dim]) 
        H_min = np.amin(ref_network['pore.coords'][np.where(ref_network['pore.internal'])][:, L_dim]) 
        W_min = np.amin(ref_network['pore.coords'][np.where(ref_network['pore.internal'])][:, H_dim])
        L_min = np.amin(ref_network['pore.coords'][np.where(ref_network['pore.internal'])][:, W_dim])         
    
    for num in range(len(net)):  
        # Determine if a pore transforms, and if yes, whether merges or splits     
        # Create an array of pore coordinates and zeros for further processing
        merge_pores=np.zeros((len(net[num]['pore.internal']),3))
        # Create an array of pore coordinates and zeros for further processing
        split_pores=np.zeros((len(net[num]['pore.internal']),3))           
        # Loop through all pores, if conditions are met, pore coordinates will be stored        
        for k in range(len(np.where(net[num]['pore.internal']==True)[0])):
         # Transformation chance is a given 0.5% 
           if random.random() > M_S_chance:
             if random.random() > M_S_condition:
               # Store coordinates of to be merged pores           
               merge_pores[k]=net[num]['pore.coords'][np.where(net[num]['pore.internal']==True)[0][k]]
             else:   
               # Store coordinates of to be splitted pores           
               split_pores[k]=net[num]['pore.coords'][np.where(net[num]['pore.internal']==True)[0][k]]
           else:
             # If criteria are not met store a 0
             merge_pores[k]=0
             split_pores[k]=0   
        # Erase all zero values from the array     
        merge_pores = merge_pores[~np.all(merge_pores == 0, axis=1)]       
        # 3D array that stores pore coordinates to be splitted, in case no splitting is necessary coordinates are 0         
        split_pores = split_pores[~np.all(split_pores == 0, axis=1)]
 
        #  Merging pores 
        # If length of merge_pores is zero, no merging will take place
        if len(merge_pores) == 0:
            pass       
        else:
          # Loop through all the mergeable pores 
          for k in range(len(merge_pores)):            
              # Find all the neighboring pores next to the selected one
              neighbors=net[num].find_neighbor_pores(pores=np.where((net[num]['pore.coords'][:, None] == merge_pores[k]).all(-1).any(-1) == True))                  
              for z in range(len(neighbors)):
                      # Find a pore with the smallest diameter among the neighbors 
                      to_be_merged_2=neighbors[np.argsort(net[num]['pore.diameter'][neighbors])[z]]
                      # Check whether this pore on the surface is or not
                      if net[num]['pore.surface'][to_be_merged_2]==True:
                          pass
                      # Break the for loop if requirement fulfilled
                      else:  
                          break
              # Create an array for the to-be-merged pore
              to_be_merged=np.append(np.where((net[num]['pore.coords'][:, None] == merge_pores[k]).all(-1).any(-1) == True), to_be_merged_2)
              if len(to_be_merged) < 2:
                  break
              else:
                  pass
              # Calculate the sum of volume for the two to-be-merged pores
              to_be_merged_vol=(net[num]['pore.diameter'][to_be_merged[0]])**3*4*np.pi/3+(net[num]['pore.diameter'][to_be_merged[1]])**3*4*np.pi/3
              # Merge the pores and store them under a label
              op.topotools.merge_pores(network=net[num], pores=[to_be_merged], geo=geo[num], labels=['merged']) 
        
              # Add the merged pore to the geometry
              geo[num]._set_locations(element='pores', indices=len(net[num]['pore.coords'])-1, mode='add')  
              # Calculate new diameter
              geo[num]['pore.diameter'][-1]=np.cbrt(to_be_merged_vol/4*3/np.pi)
              
              # Function to prevent overlap of the new pore with any other pore:
              def overlap_check_M(net, geo, to_be_merged_vol, H_min, H_max, W_max, W_min, L_max, L_min):
                  # Define overlap as False
                  overlap = True
                  while overlap == True:   
                      overlap = False  
                      # Find the neighboring pores:
                      neighbors_overlap = net.find_nearby_pores(pores=len(net['pore.coords'])-1, r=geo['pore.diameter'].max(), flatten=True)
                      # Find if the new pore overlaps with existing pores and otherwise delete the location of the new pore:
                      for i in neighbors_overlap:
                        # Find the coordinates of the new pores and the neighboring pores:
                        coords_new = net['pore.coords'][-1]
                        coords = net['pore.coords'][i]
                        # Calculate the center-to-center distance:
                        center_to_center = np.sqrt((coords_new[0]-coords[0])**2 + (coords_new[1]-coords[1])**2 + (coords_new[2]-coords[2])**2)
                        # Check if the new pore overlaps a neighboring pore:
                        if center_to_center < (geo['pore.diameter'][-1]/2 + geo['pore.diameter'][i]/2):
                            # Remove the location of the new pore when it overlaps another pore:
                            x = coords_new[0]+[-1,1][random.randrange(2)]*random.uniform(0, 1)*geo['pore.diameter'].mean()/2
                            y = coords_new[1]+[-1,1][random.randrange(2)]*random.uniform(0, 1)*geo['pore.diameter'].mean()/2
                            z = coords_new[2]+[-1,1][random.randrange(2)]*random.uniform(0, 1)*geo['pore.diameter'].mean()/2
                            # Check if the new coordinates are within the network dimensions:
                            if not H_min < x:
                                x = H_min + geo['pore.diameter'].mean()/2
                            if not H_max > x:
                                x = H_max - geo['pore.diameter'].mean()/2
                            if not L_min < y:
                                y = L_min + geo['pore.diameter'].mean()/2
                            if not L_max > y:
                                y = L_max - geo['pore.diameter'].mean()/2
                            if not W_min < z:
                                z = W_min + geo['pore.diameter'].mean()/2
                            if not W_max > z:
                                z = W_max - geo['pore.diameter'].mean()/2                               
                            # Update the pore coordinates:
                            net['pore.coords'][-1] = [x, y, z]
                            overlap = True
                            break
                        else:
                            overlap = False
                  return net, geo    
             
              # Check if the new pore overlaps a neighboring pore:
              net[num], geo[num] = overlap_check_M(net=net[num], geo=geo[num], to_be_merged_vol=to_be_merged_vol, H_min=H_min, H_max=H_max, W_max=W_max, W_min=W_min, L_max=L_max, L_min=L_min)
              # Add new throats to the geometry
              geo[num]._set_locations(element='throats', indices=np.where(net[num]['throat.geo_01']==False) , mode='add')           
              # Find the pores to which the new merged pore is connected to
              neighbors_connected = net[num].find_neighbor_pores(pores=len(net[num]['pore.coords'])-1)
              neighbors_connected_2 = neighbors_connected.astype(np.float32)
              # Define the thoat diameter for the new throat connections based on the smallest connected pore diameter
              for p in range(len(neighbors_connected_2)):
                  # print(p)
                  neighbor_diameter = net[num]['pore.diameter'][neighbors_connected][p]
                  if neighbor_diameter > net[num]['pore.diameter'][-1]:
                      throat_diameter_check = net[num]['pore.diameter'][-1] * throat_condition
                  else:
                    throat_diameter_check = neighbor_diameter * throat_condition
                  connecting_throat_1 = net[num].find_connecting_throat(len(net[num]['pore.coords'])-1, neighbors_connected[p])
                  connecting_throat = connecting_throat_1[0]
                  geo[num]['throat.diameter'][connecting_throat] =  throat_diameter_check
              
              # Delete existing to_be_merged arrays
              del to_be_merged, to_be_merged_vol                          
              #  Regenerate all the models 
              if Extracted == 0 or Extracted == 2:
                   geo[num].regenerate_models(['pore.volume', 'throat.max_size', 'throat.length', 'throat.endpoints', 'throat.cross_sectional_area',
                                                   'throat.conduit_lengths', 'throat.area', 'throat.volume','throat.surface_area', 'pore.area', 'pore.surface_area'])
              else:
                   geo[num].regenerate_models(['pore.volume', 'pore.area', 'throat.length', 'throat.volume', 'throat.endpoints', 'throat.centroid', 'throat.conduit_lengths', 'throat.cross_sectional_area',
                                                   'throat.area', 'throat.surface_area', 'pore.surface_area',  'throat.max_size', 'pore.equivalent_diameter', 'pore.centroid', 'pore.label', 
                                                   'pore.region_volume', 'throat.total_length', 'throat.equivalent_diameter', 'throat.inscribed_diameter', 'throat.direct_length', 'throat.perimeter',
                                                   'pore.extended_diameter', 'pore.inscribed_diameter'])      
                   geo[num]['pore.centroid'] = geo[num]['pore.coords']
              net[num]['pore.internal'][net[num].pores('merged')] = True
              net[num]['throat.internal'][net[num].throats('merged')] = True            
          # Delete variables:
          del neighbors, to_be_merged_2, neighbors_connected, neighbors_connected_2, neighbor_diameter, throat_diameter_check, connecting_throat_1, connecting_throat
   
        # Splitting pores          
        # If length of split_pores is zero, no merging will take place 
        if len(split_pores) == 0:
            pass
        else: 
            for k in range(len(split_pores)):
               if not np.where((net[num]['pore.coords'][:, None] == split_pores[k]).all(-1).any(-1) == True):
                  continue
               else:
                  pass
                  # Number of new pores created
                  num_pores=2
                  # Store new pore coordinates within new_pores array
                  new_pores=np.zeros(shape=(num_pores,3))                       
                  # Define the coordinates of the new pores:
                  for i in range(num_pores):
                       new_pores[i]=[split_pores[k][0]+[-1,1][random.randrange(2)]*random.uniform(0, 1)*geo[num]['pore.diameter'].mean()/2,
                       split_pores[k][1]+[-1,1][random.randrange(2)]*random.uniform(0, 1)*geo[num]['pore.diameter'].mean()/2,
                       split_pores[k][2]+[-1,1][random.randrange(2)]*random.uniform(0, 1)*geo[num]['pore.diameter'].mean()/2]
                       # Check if the new coordinates are within the network dimensions:
                       if not H_min < new_pores[i][0]:
                           new_pores[i][0] = H_min + geo[num]['pore.diameter'].mean()/2
                       if not H_max > new_pores[i][0]:
                           new_pores[i][0] = H_max - geo[num]['pore.diameter'].mean()/2
                       if not L_min < new_pores[i][1]:
                           new_pores[i][1] = L_min + geo[num]['pore.diameter'].mean()/2
                       if not L_max > new_pores[i][1]:
                           new_pores[i][1] = L_max - geo[num]['pore.diameter'].mean()/2
                       if not W_min < new_pores[i][2]:
                           new_pores[i][2] = W_min + geo[num]['pore.diameter'].mean()/2
                       if not W_max > new_pores[i][2]:
                           new_pores[i][2] = W_max - geo[num]['pore.diameter'].mean()/2 

                  # Extend network with previously defined pores    
                  op.topotools.extend(network=net[num], pore_coords=new_pores, labels=['splitted'])     
                  # Volume of splitted pore
                  split_volume=(geo[num]['pore.diameter'][np.where((net[num]['pore.coords'][:, None] == split_pores[k]).all(-1).any(-1) == True)])**3*4/3*np.pi
                  # Move on if split_volume empty
                  if len(split_volume) == 0:
                        op.topotools.trim(net[num], pores=np.argwhere(np.isnan(geo[num]['pore.diameter'])))
                        continue
                  # Find the neighbors of the to be splitted pore:  
                  neighbors_old = net[num].find_neighbor_pores(np.where((net[num]['pore.coords'][:, None] == split_pores[k]).all(-1).any(-1) == True)) 
                  # Add the new pores to the geometry:
                  geo[num]._set_locations(element='pores', indices=len(net[num]['pore.all'])-2, mode='add')  
                  geo[num]._set_locations(element='pores', indices=len(net[num]['pore.all'])-1, mode='add') 
                  # Define the diameter of the new pored:
                  geo[num]['pore.diameter'][-2]=np.cbrt(split_volume/num_pores*3/4/np.pi)
                  geo[num]['pore.diameter'][-1]=np.cbrt(split_volume/num_pores*3/4/np.pi)
                  
                  # Function to prevent overlap of the new pores with any other pore:
                  def overlap_check_S(j, net, geo, split_volume, H_min, H_max, W_max, W_min, L_max, L_min, new_pores):                   
                     # Condition that checks for overlap:
                     overlap = True
                     # Add the splitted pore to the geometry
                     while overlap == True:
                            overlap = False
                            # Find the neighboring pores:
                            neighbors_overlap = net.find_nearby_pores(pores=len(net['pore.coords'])-num_pores+j, r=geo['pore.diameter'].max(), flatten=True)
                            # Find if the new pore overlaps with existing pores and otherwise delete the location of the new pore:
                            for i in neighbors_overlap:
                                # Find the coordinates of the new pores and the neighboring pores:
                                coords_new = net['pore.coords'][-num_pores+j]
                                coords = net['pore.coords'][i]
                                # Calculate the center-to-center distance:
                                center_to_center = np.sqrt((coords_new[0]-coords[0])**2 + (coords_new[1]-coords[1])**2 + (coords_new[2]-coords[2])**2)
                                # Check if the new pore overlaps a neighboring pore:
                                if center_to_center < (geo['pore.diameter'][-num_pores+j]/2 + geo['pore.diameter'][i]/2):
                                   # Remove the location of the new pore when it overlaps another pore:                              
                                   x = coords_new[0]+[-1,1][random.randrange(2)]*random.uniform(0, 1)*geo['pore.diameter'].mean()/2
                                   y = coords_new[1]+[-1,1][random.randrange(2)]*random.uniform(0, 1)*geo['pore.diameter'].mean()/2
                                   z = coords_new[2]+[-1,1][random.randrange(2)]*random.uniform(0, 1)*geo['pore.diameter'].mean()/2
                                   # Check if the new coordinates are within the network dimensions:
                                   if not H_min < x:
                                        x = H_min + geo['pore.diameter'].mean()/2
                                   if not H_max > x:
                                        x = H_max - geo['pore.diameter'].mean()/2
                                   if not L_min < y:
                                        y = L_min + geo['pore.diameter'].mean()/2
                                   if not L_max > y:
                                        y = L_max - geo['pore.diameter'].mean()/2
                                   if not W_min < z:
                                        z = W_min + geo['pore.diameter'].mean()/2
                                   if not W_max > z:
                                        z = W_max - geo['pore.diameter'].mean()/2     
                                   # Update the pore coordinates:
                                   net['pore.coords'][-num_pores+j] = [x, y, z]
                                   new_pores[j] = [x, y, z]
                                   # Set the overlap condition back to True:
                                   overlap = True
                                   break
                                else: 
                                   overlap = False                                                            
                     return net, geo, new_pores

                  # Check if the new pore overlaps a neighboring pore:
                  for j in range(num_pores):
                      if i == num_pores-1:      
                         break
                      net[num], geo[num], new_pores = overlap_check_S(j=j, net=net[num], geo=geo[num], split_volume=split_volume, H_min=H_min, H_max=H_max, W_max=W_max, W_min=W_min, L_max=L_max, L_min=L_min, new_pores=new_pores)

                  # Connect the two new pores:
                  op.topotools.connect_pores(network=net[num], pores1=len(net[num]['pore.coords'])-2, pores2=len(net[num]['pore.coords'])-1, labels=['splitted'])
                  geo[num]._set_locations(element='throats', indices=len(net[num]['throat.all'])-1, mode='add')                 
                  # Define the new throat connections to the new pores:
                  neighbors_connected = neighbors_old.astype(np.float32)
                  for p in range(len(neighbors_connected)):
                      # Find the pore coordinates:
                      coords_1 = net[num]['pore.coords'][-2]
                      coords_2 = net[num]['pore.coords'][-1]
                      coords_neighbor = net[num]['pore.coords'][neighbors_old][p]
                      # Calculate the center-to-center distance:
                      center_to_center_1 = np.sqrt((coords_1[0]-coords_neighbor[0])**2 + (coords_1[1]-coords_neighbor[1])**2 + (coords_1[2]-coords_neighbor[2])**2)
                      center_to_center_2 = np.sqrt((coords_2[0]-coords_neighbor[0])**2 + (coords_2[1]-coords_neighbor[1])**2 + (coords_2[2]-coords_neighbor[2])**2)
                      # Find the new pore to which the old neighbor should be connected to:
                      if center_to_center_1 > center_to_center_2:
                          op.topotools.connect_pores(network=net[num], pores1=len(net[num]['pore.coords'])-2, pores2=neighbors_old[p], labels=['splitted'])
                      else:
                          op.topotools.connect_pores(network=net[num], pores1=len(net[num]['pore.coords'])-1, pores2=neighbors_old[p], labels=['splitted'])
                      # Add the throats:  
                      geo[num]._set_locations(element='throats', indices=len(net[num]['throat.all'])-1, mode='add')
                      
                  # Find the diameters:
                  diameter_1 = net[num]['pore.diameter'][-2]
                  diameter_2 = net[num]['pore.diameter'][-1]
                  # Calculate the throat diameter of the throat connecting the two new pores:
                  if diameter_1 > diameter_2:
                      throat_diameter1_check = diameter_2 * throat_condition
                  else:
                      throat_diameter1_check = diameter_1 * throat_condition
                  connecting_throat1_1 = net[num].find_connecting_throat(len(net[num]['pore.coords'])-2, len(net[num]['pore.coords'])-1)
                  connecting_throat1 = connecting_throat1_1[0]
                  geo[num]['throat.diameter'][connecting_throat1] =  throat_diameter1_check                 
                  # Calculate the throat diameters of the other throats connecting one of the new pores to an old neighbor of the splitted pore:
                  for p in range(len(neighbors_connected)):                     
                       neighbor_diameter = net[num]['pore.diameter'][neighbors_old][p]
                       connecting_throat1_2 = net[num].find_connecting_throat(len(net[num]['pore.coords'])-1, neighbors_old[p])
                       connecting_throat2_2 = net[num].find_connecting_throat(len(net[num]['pore.coords'])-2, neighbors_old[p])
                       # Loop to find the correct pore to which the new throat is connected to:
                       if connecting_throat1_2 == [None]:
                            if neighbor_diameter > net[num]['pore.diameter'][-2]:
                                throat_diameter_check = net[num]['pore.diameter'][-2] * throat_condition
                            else:
                                throat_diameter_check = neighbor_diameter * throat_condition
                            connecting_throat_1 = net[num].find_connecting_throat(len(net[num]['pore.coords'])-2, neighbors_old[p])
                            connecting_throat = connecting_throat_1[0]
                            geo[num]['throat.diameter'][connecting_throat] =  throat_diameter_check
                       elif connecting_throat2_2 == [None]:
                            if neighbor_diameter > net[num]['pore.diameter'][-1]:
                                throat_diameter_check = net[num]['pore.diameter'][-1] * throat_condition
                            else:
                                throat_diameter_check = neighbor_diameter * throat_condition
                            connecting_throat_1 = net[num].find_connecting_throat(len(net[num]['pore.coords'])-1, neighbors_old[p])
                            connecting_throat = connecting_throat_1[0]
                            geo[num]['throat.diameter'][connecting_throat] =  throat_diameter_check

                  # Trim splitted pore:
                  op.topotools.trim(network=net[num], pores=[np.where((net[num]['pore.coords'][:, None] == split_pores[k]).all(-1).any(-1) == True)]) 
                  net[num]['pore.internal'][net[num].pores('splitted')] = True
                  net[num]['throat.internal'][net[num].throats('splitted')] = True                 
                  # Regenerate models:
                  if Extracted == 0 or Extracted == 2:
                       geo[num].regenerate_models(['pore.volume', 'pore.max_size', 'throat.max_size', 'throat.length', 'throat.endpoints', 'throat.cross_sectional_area',
                                                   'throat.conduit_lengths', 'throat.area', 'throat.volume','throat.surface_area', 'pore.area', 'pore.surface_area'])
                  else:
                       geo[num].regenerate_models(['pore.volume', 'pore.area', 'throat.length', 'throat.volume', 'throat.endpoints', 'throat.centroid', 'throat.conduit_lengths', 'throat.cross_sectional_area',
                                                   'throat.area', 'throat.surface_area', 'pore.surface_area',  'throat.max_size', 'pore.equivalent_diameter', 'pore.centroid', 'pore.label', 
                                                   'pore.region_volume', 'throat.total_length', 'throat.equivalent_diameter', 'throat.inscribed_diameter', 'throat.direct_length', 'throat.perimeter',
                                                   'pore.extended_diameter', 'pore.inscribed_diameter'])      
                       geo[num]['pore.centroid'] = geo[num]['pore.coords']                         
                  # Delete variables:
                  del split_volume, neighbors_old, new_pores, neighbors_connected, coords_1, coords_2, coords_neighbor, center_to_center_1, center_to_center_2, diameter_1, diameter_2, throat_diameter1_check, connecting_throat1_1, connecting_throat1, neighbor_diameter, connecting_throat1_2 , connecting_throat2_2, throat_diameter_check, connecting_throat_1, connecting_throat    

        net[num]['throat.conns']=net[num]['throat.conns'].astype(int)        
    # Delete variables:
    del merge_pores, split_pores
    
    return net, geo


def check_volume(net, geo, ref_network, Extracted, spacing):    
    for i in range(len(net)): 
        # print('Checking overall pore volume of network:', i+1)
        norm_factor = 1.005        
        # Check if the overall pore volume must be scaled up or down:
        while (net[i]['pore.volume'][net[i].pores('internal')].sum() + net[i]['throat.volume'][net[i].throats('internal')].sum()) > 1.01 * (ref_network['pore.volume'][ref_network.pores('internal')].sum() + ref_network['throat.volume'][ref_network.throats('internal')].sum()) or (net[i]['pore.volume'][net[i].pores('internal')].sum() + net[i]['throat.volume'][net[i].throats('internal')].sum()) < 0.99 * (ref_network['pore.volume'][ref_network.pores('internal')].sum() + ref_network['throat.volume'][ref_network.throats('internal')].sum()):    
            if (net[i]['pore.volume'][net[i].pores('internal')].sum() + net[i]['throat.volume'][net[i].throats('internal')].sum()) > (ref_network['pore.volume'][ref_network.pores('internal')].sum() + ref_network['throat.volume'][ref_network.throats('internal')].sum()):
                # Decrease the pore diameters with a constant multiplication until the overall pore volume matches up to reasonable accuracy.
                while (net[i]['pore.volume'][net[i].pores('internal')].sum() + net[i]['throat.volume'][net[i].throats('internal')].sum()) > 1.01 * (ref_network['pore.volume'][ref_network.pores('internal')].sum() + ref_network['throat.volume'][ref_network.throats('internal')].sum()):
                    # Decrease pore diameter with constant multiplication:
                    geo[i]['pore.diameter'][net[i].pores('internal')] /= norm_factor                
                    if Extracted == 0:
                        geo[i]['pore.diameter'][net[i].pores('surface')] = 1/50 * spacing
                    # Recalculate volume
                    geo[i].regenerate_models(['throat.max_size', 'throat.diameter', 'pore.volume', 'throat.endpoints', 'throat.volume', 'throat.length'])                
            if (net[i]['pore.volume'][net[i].pores('internal')].sum() + net[i]['throat.volume'][net[i].throats('internal')].sum()) < (ref_network['pore.volume'][ref_network.pores('internal')].sum() + ref_network['throat.volume'][ref_network.throats('internal')].sum()):
                # Increase the pore diameters with a constant multiplication until the overall pore volume matches up to reasonable accuracy.
                while (net[i]['pore.volume'][net[i].pores('internal')].sum() + net[i]['throat.volume'][net[i].throats('internal')].sum()) < 0.99 * (ref_network['pore.volume'][ref_network.pores('internal')].sum() + ref_network['throat.volume'][ref_network.throats('internal')].sum()):
                    # Increase pore diameter with constant multiplication:
                    geo[i]['pore.diameter'][net[i].pores('internal')] *= norm_factor
                    if Extracted == 0:
                        geo[i]['pore.diameter'][net[i].pores('surface')] = 1/50 * spacing           
                    # Recalculate volume
                    geo[i].regenerate_models(['throat.max_size', 'throat.diameter', 'pore.volume', 'throat.endpoints', 'throat.volume', 'throat.length'])                       
        
        geo[i].regenerate_models(['throat.cross_sectional_area',
                                              'throat.conduit_lengths', 'throat.area', 'throat.surface_area', 'pore.area', 'pore.surface_area'])
        
    return net, geo    


def add_geometry_models(net, geo, throat_condition, spacing, Extracted):
    for i in range(len(net)):    
        geo[i].add_model(propname='pore.area', model=op.models.geometry.pore_area.sphere, pore_diameter='pore.diameter')
        geo[i].add_model(propname='pore.volume', model=op.models.geometry.pore_volume.sphere, pore_diameter='pore.diameter')    
        geo[i].add_model(propname='throat.max_size', model=op.models.misc.from_neighbor_pores, mode='min', prop='pore.diameter')
        geo[i].add_model(propname='throat.diameter', model=op.models.misc.scaled, factor=throat_condition, prop='throat.max_size')
        geo[i].add_model(propname='throat.endpoints', model=op.models.geometry.throat_endpoints.spherical_pores, pore_diameter='pore.diameter', throat_diameter='throat.diameter')  
        geo[i].add_model(propname='throat.length', model=op.models.geometry.throat_length.piecewise, throat_endpoints='throat.endpoints')
        geo[i].add_model(propname='throat.surface_area', model=op.models.geometry.throat_surface_area.cylinder, throat_diameter='throat.diameter', throat_length='throat.length')
        geo[i].add_model(propname='throat.volume', model=op.models.geometry.throat_volume.cylinder, throat_diameter='throat.diameter', throat_length='throat.length')
        geo[i].add_model(propname='throat.cross_sectional_area', model=op.models.geometry.throat_area.cylinder, throat_diameter='throat.diameter')
        geo[i].add_model(propname='throat.conduit_lengths', model=op.models.geometry.throat_length.conduit_lengths, throat_endpoints='throat.endpoints', throat_length='throat.length')
        geo[i].add_model(propname='throat.area', model=op.models.geometry.throat_area.cylinder, throat_diameter='throat.diameter')       
        geo[i].add_model(propname='pore.surface_area', model=op.models.geometry.pore_surface_area.sphere, pore_diameter='pore.diameter',throat_cross_sectional_area='throat.cross_sectional_area', throat_surface_area='throat.surface_area')
         
        if Extracted == 0:
            geo[i]['pore.diameter'][net[i].pores('surface')] = 1/50 * spacing
        if Extracted == 1:
            geo[i].add_model(propname='pore.region_volume', model=op.models.geometry.pore_volume.sphere, pore_diameter='pore.diameter')
            geo[i].add_model(propname='throat.inscribed_diameter', model=op.models.misc.scaled, factor=throat_condition, prop='throat.max_size')
            geo[i].add_model(propname='throat.equivalent_diameter', model=op.models.misc.scaled, factor=throat_condition, prop='throat.max_size')
            geo[i].add_model(propname='throat.total_length', model=op.models.geometry.throat_length.piecewise, throat_endpoints='throat.endpoints')
            geo[i].add_model(propname='throat.direct_length', model=op.models.geometry.throat_length.piecewise, throat_endpoints='throat.endpoints')
            geo[i].add_model(propname='throat.perimeter', model=op.models.geometry.throat_perimeter.cylinder, throat_diameter='throat.diameter')  
            geo[i].add_model(propname='pore.equivalent_diameter', model=op.models.geometry.pore_size.equivalent_diameter, pore_volume='pore.volume', pore_shape='sphere')
            geo[i].add_model(propname='pore.extended_diameter', model=op.models.geometry.pore_size.equivalent_diameter, pore_volume='pore.volume', pore_shape='sphere')
            geo[i].add_model(propname='pore.inscribed_diameter', model=op.models.geometry.pore_size.equivalent_diameter, pore_volume='pore.volume', pore_shape='sphere') 
            geo[i].add_model(propname='pore.label', model=op.models.geometry.pore_label) 
            f_throat_centroids = throat_centroid
            geo[i].add_model(propname = 'throat.centroid', model = f_throat_centroids)    
        
    return geo


def add_geometry_models_ref(net, geo, throat_condition, spacing, Extracted):      
    geo.add_model(propname='pore.area', model=op.models.geometry.pore_area.sphere, pore_diameter='pore.diameter')
    geo.add_model(propname='pore.volume', model=op.models.geometry.pore_volume.sphere, pore_diameter='pore.diameter')    
    geo.add_model(propname='throat.max_size', model=op.models.misc.from_neighbor_pores, mode='min', prop='pore.diameter')
    if Extracted == 1 or Extracted == 2:
        geo.add_model(propname='throat.diameter', model=op.models.misc.scaled, factor=throat_condition, prop='throat.max_size')
    geo.add_model(propname='throat.endpoints', model=op.models.geometry.throat_endpoints.spherical_pores, pore_diameter='pore.diameter', throat_diameter='throat.diameter')  
    geo.add_model(propname='throat.length', model=op.models.geometry.throat_length.piecewise, throat_endpoints='throat.endpoints')
    geo.add_model(propname='throat.surface_area', model=op.models.geometry.throat_surface_area.cylinder, throat_diameter='throat.diameter', throat_length='throat.length')
    geo.add_model(propname='throat.volume', model=op.models.geometry.throat_volume.cylinder, throat_diameter='throat.diameter', throat_length='throat.length')
    geo.add_model(propname='throat.cross_sectional_area', model=op.models.geometry.throat_area.cylinder, throat_diameter='throat.diameter')
    geo.add_model(propname='throat.conduit_lengths', model=op.models.geometry.throat_length.conduit_lengths, throat_endpoints='throat.endpoints', throat_length='throat.length')
    geo.add_model(propname='throat.area', model=op.models.geometry.throat_area.cylinder, throat_diameter='throat.diameter')       
    geo.add_model(propname='pore.surface_area', model=op.models.geometry.pore_surface_area.sphere, pore_diameter='pore.diameter',throat_cross_sectional_area='throat.cross_sectional_area', throat_surface_area='throat.surface_area')
     
    if Extracted == 0:
        geo['pore.diameter'][net.pores('surface')] = 1/50 * spacing
    if Extracted == 1:
        geo.add_model(propname='pore.region_volume', model=op.models.geometry.pore_volume.sphere, pore_diameter='pore.diameter')
        geo.add_model(propname='throat.inscribed_diameter', model=op.models.misc.scaled, factor=throat_condition, prop='throat.max_size')
        geo.add_model(propname='throat.equivalent_diameter', model=op.models.misc.scaled, factor=throat_condition, prop='throat.max_size')
        geo.add_model(propname='throat.total_length', model=op.models.geometry.throat_length.piecewise, throat_endpoints='throat.endpoints')
        geo.add_model(propname='throat.direct_length', model=op.models.geometry.throat_length.piecewise, throat_endpoints='throat.endpoints')
        geo.add_model(propname='throat.perimeter', model=op.models.geometry.throat_perimeter.cylinder, throat_diameter='throat.diameter')  
        geo.add_model(propname='pore.equivalent_diameter', model=op.models.geometry.pore_size.equivalent_diameter, pore_volume='pore.volume', pore_shape='sphere')
        geo.add_model(propname='pore.extended_diameter', model=op.models.geometry.pore_size.equivalent_diameter, pore_volume='pore.volume', pore_shape='sphere')
        geo.add_model(propname='pore.inscribed_diameter', model=op.models.geometry.pore_size.equivalent_diameter, pore_volume='pore.volume', pore_shape='sphere') 
        geo.add_model(propname='pore.label', model=op.models.geometry.pore_label) 
        f_throat_centroids = throat_centroid
        geo.add_model(propname = 'throat.centroid', model = f_throat_centroids)    
        
    return geo
               
        
def nw_health_diag(net, geo):    
    """ This function checks the network health based on the in-built OpenPNM function 'check_network_health(). 
    After finding faulty network segments (isolated pores, duplicated throats, etc.), trims the corresponding throats and pores, 
    so the algorithm can handle the network and geometry afterwards. """ 
      # Loop through the range of networks
    for i in range(len(net)):
        # Create network health dictionary
        health = net[i].check_network_health()
        # Identify dictionary keys
        health_keys=list(health.keys())
        # Check if 'trim_pores' key contains any information
        if bool(health[health_keys[2]]) == True:
           # Trim the corresponding pores, if yes
           op.topotools.trim(network=net[i], pores=health['trim_pores'])
        # Check if 'duplicate_throats' key contains any information
        if bool(health[health_keys[3]]) == True:
           # Trim the corresponding throats, if yes
           op.topotools.trim(network=net[i], throats=health['duplicate_throats'])
      
    return net, geo  


def throat_centroid(target):
    network = target.project.network
    Ts = network.throats(target.name)
    conns = network['throat.conns']
    coords = network['pore.coords']
    
    return np.mean(coords[conns], axis=1)[Ts]