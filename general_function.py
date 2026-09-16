
#%%
# Jupytext: {"Function for tolerence distribution": {"extension": ".py", "format_name": "light", "format_version": "0.4"}}
import numpy as np

def Gaussian_tolerence(ID, tolerance,std,num_samples):
    nominal = ID
    std_dev = tolerance / std

    return np.random.normal(nominal, std_dev, num_samples)

#%%
# Jupytext: {"Gerneral functions for fluid flow friction": {"extension": ".py", "format_name": "light", "format_version": "0.4"}}
def reynolds_number(Q, D, rho, mu):
    """Calculate Reynold's number for a tube of diameter D with flow rate Q
    Flow rate, Q [m^3/s]
    Diameter, D [m]
    Density, rho []
    Viscosity, mu []"""
    A = np.pi * (D/2)**2
    V = Q / A
    return (rho * V * D) / mu


def attenuate_length(Re, D):
    """Calculate length to establish flow pattern give a Reynold number.
    Tube ID, D [m]"""
    if Re < 2000:
        return 0.13*Re*D  # Laminar flow
    else:
        return 0.0573*D  # Approximate for turbulent flow


def flow_rate_laminar(P_inlet, P_outlet, D, length, mu):
     """Hagen-Poiseuille for Laminar Flow """
     return (np.pi * (P_inlet - P_outlet) * D**4) / (128 * mu * length) # not needed, as the turbulence funtion already cover this laminar flow range with friction factor


def friction_factor(Re):
    """Darcy-Weisbach friction factor assuming a smooth-wall tube
     (Re < 65D/k) k of silcone tube 0.0015mm Re<1.7e5"""
    if Re <= 2000: # Laminar flow
        return 64/Re

    #TODO: Critical regime 2000 < Re < 4000. Perhaps raise a warning ?
    elif 2000 < Re < 100000:
        return 0.3164 / (Re ** 0.25)  # Approximate for turbulent flow
    else: # Is this condition Re>1e5 ever met?
        return 0.0032+0.221/ (Re ** 0.237)  # Approximate for turbulent flow


def effective_f_friction(num,f,length,inner_diameter):
    Delta = length/(num+1)
    x = Delta/inner_diameter
    if x < 10:
        return f*(100+f/(x+1))
    else:
        return f

def InletDrag_coefficient(Re, ID_small, ID_big):
    """ Drag coefficient for tube ID big to small (e.g vessel to tube or tube to barb) """
    if Re < 2000:
        return 1.08-0.4*((ID_small/2)**2/(ID_big/2)**2)  # Laminar flow
    else:
        #return 0.55*(0.617719 + (-0.010875) * np.log(Re))-0.4*((ID_small/2)**2/(ID_big/2)**2)  # Approximate for turbulent flow
        #return (1-(ID_small/2)/(ID_big/2))**2  # Approximate for turbulent flow
        term = (0.5 / 3.7)**1.11 + 6.9 / np.maximum(Re, 1e-12)
        return 1.0 / (-1.8 * np.log10(term))**2   
    
def get_zeta_bend(reynolds, r_over_di):
    """
    Estimates the drag coefficient (zeta_b) for 90-degree bends/elbows.
    
    Parameters:
    -----------
    reynolds : float
        Reynolds number (Re_i).
    r_over_di : float
        The radius-to-diameter ratio (r/di). 
        Available: 2.26, 3.04, 6.53, 11.71.
    """
    # Empirical log-log data points extracted from the curves
    # Re: [10^1, 10^2, 10^3, 10^4, 10^5]
    data = {
        2.26:  [200, 10, 2, 0.5, 0.2],
        3.04:  [150, 12, 2.5, 0.8, 0.3],
        6.53:  [80, 8, 1.8, 0.4, 0.15],
        11.71: [70, 7, 1.5, 0.3, 0.1]
    }
    
    if r_over_di not in data:
        raise ValueError("r/di must be one of: 2.26, 3.04, 6.53, 11.71")
        
    # Use log-log linear interpolation
    log_re = np.log10(reynolds)
    log_re_points = np.log10([10, 100, 1000, 10000, 100000])
    log_zeta_points = np.log10(data[r_over_di])
    
    return 10**np.interp(log_re, log_re_points, log_zeta_points)

# --- Example Usage ---
re_val = 5000
ratio = 3.04
zeta = get_zeta_bend(re_val, ratio)

print(f"For Re={re_val} and r/di={ratio}, zeta_b ≈ {zeta:.3f}")


    
def get_zeta_in(small_ID, large_ID, regime='turbulent', Re=None):
    """
    Calculates the drag coefficient (zeta_in) for an abrupt reduction in tube diameter.
    
    Parameters:
    -----------
    area_ratio : float
        The ratio of the cross-sectional areas (A2/A1). Must be between 0.0 and 1.0.
    regime : str
        Flow regime: 'laminar' or 'turbulent'. Defaults to 'turbulent'.
    Re : float or np.inf, optional
        Reynolds number. Required if regime is 'turbulent'.
        Handles exact lines (2000, 10000, inf) or interpolates between them.
        
    Returns:
    --------
    float : Estimated drag coefficient (zeta_in)
    """

    area_ratio = small_ID/large_ID
    # 1. Input Validation
    if not (0.0 <= area_ratio <= 1.0):
        raise ValueError("area_ratio (A2/A1) must be between 0.0 and 1.0")
        
    regime = regime.lower()
    if regime not in ['laminar', 'turbulent']:
        raise ValueError("regime must be either 'laminar' or 'turbulent'")

    # 2. Laminar Regime Function
    if regime == 'laminar':
        return 1.10 - 0.43 * area_ratio

    # 3. Turbulent Regime Functions
    if Re is None:
        raise ValueError("Reynolds number (Re) must be specified for turbulent regime.")
        
    # Exact curve equations derived from the plot
    zeta_2000 = 0.56 - 0.41 * area_ratio
    zeta_10k  = 0.50 - 0.40 * area_ratio
    zeta_inf  = 0.42 - 0.42 * area_ratio

    # 4. Handle Exact Reynolds Numbers or Interpolation
    if Re <= 2000:
        return zeta_2000
    elif Re == 10000:
        return zeta_10k
    elif Re == np.inf or Re >= 1e7: # Approximating infinity for practical numbers
        return zeta_inf
    
    # Logarithmic interpolation between 2000 and 10,000
    elif 2000 < Re < 10000:
        log_Re = np.log10(Re)
        log_2k = np.log10(2000)
        log_10k = np.log10(10000)
        # Linear weight based on log scale
        weight = (log_Re - log_2k) / (log_10k - log_2k)
        return zeta_2000 + weight * (zeta_10k - zeta_2000)
        
    # Logarithmic interpolation between 10,000 and "Infinity" (Approximated as Re = 1,000,000)
    else: 
        max_practical_re = 1000000
        log_Re = np.log10(min(Re, max_practical_re))
        log_10k = np.log10(10000)
        log_max = np.log10(max_practical_re)
        
        weight = (log_Re - log_10k) / (log_max - log_10k)
        return zeta_10k + weight * (zeta_inf - zeta_10k)                   #!!!!! Compare


# Drag coefficient for tube ID small to big (tube to vessel/ barb to tube)
def OutletDrag_coefficient(ID_small, ID_big):
    return (1-(ID_small**2)/(ID_big**2))**2  # Approximate for turbulent flow
    #return ((ID_big/2)/(ID_small/2)-1)**2  # Approximate for turbulent flow    #!!!!! Compare

# Drag coefficient for the orifice
def OrificeDrag_coefficient(inner_diameter, ID_orifice):
    return 1-1.05*((ID_orifice/inner_diameter)**2 )     #!!!!! Verify

# Drag coefficient for the orifice
def Orifice_thin_Drag_coefficient(inner_diameter, ID_orifice):
    return ((1/0.611)-(ID_orifice**2/inner_diameter**2))**2

def Orifice_thick_Drag_coefficient(inner_diameter, ID_orifice):
    return (1/(0.6+0.4*(ID_orifice**2/inner_diameter**2))-1)**2+(1-(ID_orifice**2/inner_diameter**2))**2


def barb_dp (Q,mu, rho, ID_barb_in,ID_barb_out, ID_in, ID_out, length_barb):
    A_barb = np.pi * (((ID_barb_in+ID_barb_out)/2)/2)**2        
    Re = reynolds_number(Q, ID_in, rho, mu)
     
    # step 2: Driven the liquid from 6mm PA tubing ID2 to the barb ID_1
    V_barb = Q / A_barb
    Re_barb = reynolds_number(Q,(ID_barb_in+ID_barb_out)/2,rho,mu)
    f_barb = friction_factor(Re_barb)
             
    
    # Fittings/barbed
    f_inlet_barb = get_zeta_in(ID_barb_in, ID_in, regime='turbulent', Re = (Re + Re_barb)/2) #InletDrag_coefficient(Re,ID_barb_in,ID_in)   #liquid goes in from 6mm PA tubing to barb
            # step 3: Driven the liquid from the barb ID_2 to the 6.35 multilater tubing
    f_outlet_barb = OutletDrag_coefficient(ID_barb_out, ID_out)   #liquid goes out from barb to the 6.35 multilayer tubing
    
    dp_barb = 1* ((f_inlet_barb + f_outlet_barb  ) * (rho *  V_barb**2 / 2)        
                                 +(f_barb ) * ((length_barb) / ((ID_barb_in+ID_barb_out)/2)) * (rho * V_barb**2 / 2))

    

    return dp_barb


def Lbarb_dp (Q,mu, rho, ID_barb_in,ID_barb_out, ID_in, ID_out, length_barb):
    A_barb = np.pi * (((ID_barb_in+ID_barb_out)/2)/2)**2        
    Re = reynolds_number(Q, ID_in, rho, mu)
     
    # step 2: Driven the liquid from 6mm PA tubing ID2 to the barb ID_1
    V_barb = Q / A_barb
    Re_barb = reynolds_number(Q,(ID_barb_in+ID_barb_out)/2,rho,mu)
    f_barb = friction_factor(Re_barb)
             
    
    # Fittings/barbed
    f_inlet_barb = get_zeta_in(ID_barb_in, ID_in, regime='turbulent', Re = (Re + Re_barb)/2) #InletDrag_coefficient(Re,ID_barb_in,ID_in)   #liquid goes in from 6mm PA tubing to barb
            # step 3: Driven the liquid from the barb ID_2 to the 6.35 multilater tubing
    f_outlet_barb = OutletDrag_coefficient(ID_barb_out, ID_out)   #liquid goes out from barb to the 6.35 multilayer tubing
    
    dp_barb = 1* ((f_inlet_barb + f_outlet_barb + 1.15  ) * (rho *  V_barb**2 / 2)        
                                 +(f_barb ) * ((length_barb) / ((ID_barb_in+ID_barb_out)/2)) * (rho * V_barb**2 / 2))

    

    return dp_barb


"""
#dp_Tconver = 1*1.56* (rho * (2*V)**2 / 2) # suppose there is one T connector combiand MTBE and acetone , Ref VDI p.1070   for later
def Tbarb_dp (Q,mu, rho, ID_barb_in,ID_barb_out, ID_in, ID_out, length_barb):
    A_barb = np.pi * (((ID_barb_in+ID_barb_out)/2)/2)**2        
    
   
    Re = reynolds_number(Q, ID_in, rho, mu)
     
    # step 2: Driven the liquid from 6mm PA tubing ID2 to the barb ID_1
    V_barb = Q / A_barb
    Re_barb = reynolds_number(Q,(ID_barb_in+ID_barb_out)/2,rho,mu)
    f_barb = friction_factor(Re_barb)
             
    
    # Fittings/barbed
    f_inlet_Tbarb = get_zeta_in(ID_barb_in, ID_in, regime='turbulent', Re = (Re + Re_barb)/2) #InletDrag_coefficient(Re,ID_barb_in,ID_in)   #liquid goes in from 6mm PA tubing to barb
            # step 3: Driven the liquid from the barb ID_2 to the 6.35 multilater tubing
    f_outlet_Tbarb = OutletDrag_coefficient(ID_barb_out, ID_out)   #liquid goes out from barb to the 6.35 multilayer tubing
    
    dp_Tbarb =  ((f_inlet_Tbarb*2 + f_outlet_Tbarb + 1.56 ) * (rho *  V_barb**2 / 2)        
                                 +(f_barb ) * ((length_barb) / ((ID_barb_in+ID_barb_out)/2)) * (rho * V_barb**2 / 2))

    return dp_Tbarb

"""
def Tbarb_dp_two_inlets(
    Q1,
    mu1,
    rho1,
    Q2,
    mu2,
    rho2,
    ID_barb_in,
    ID_barb_out,
    ID_in,
    ID_out,
    length_barb,
):
    """Estimate one total T-barb pressure drop for two inlet streams."""
    if Q1 < 0 or Q2 < 0:
        raise ValueError("Inlet flow rates must be non-negative")
    Q_out = Q1 + Q2
    if Q_out <= 0:
        raise ValueError("At least one inlet flow rate must be positive")
    rho_out = (Q1 * rho1 + Q2 * rho2) / Q_out
    mu_out = (Q1 * mu1 + Q2 * mu2) / Q_out


    def inlet_loss(Q, rho, mu):
        reynolds_in = reynolds_number(Q, ID_barb_in, rho, mu)
        
        inlet_coefficient = get_zeta_in(
            ID_barb_in,
            ID_in,
            regime='turbulent',
            Re=reynolds_in,
        )
        dp_inlet = (inlet_coefficient ) * (rho * (Q / (np.pi * (ID_barb_in / 2) ** 2)) ** 2 / 2) 

        return dp_inlet

    area_barb = np.pi * (((ID_barb_in + ID_barb_out) / 2) / 2) ** 2
    diameter_barb = (ID_barb_in + ID_barb_out) / 2
    velocity_barb = Q_out / area_barb
    
    reynolds_barb = reynolds_number(Q_out, diameter_barb, rho_out, mu_out)
    friction_barb = friction_factor(reynolds_barb)

    dp_friction = friction_barb * (length_barb / diameter_barb) * (rho_out * velocity_barb**2 / 2)
    outlet_loss = (OutletDrag_coefficient(ID_barb_out, ID_out)+1.56) * (rho_out * (Q_out / (np.pi * (ID_barb_out / 2) ** 2)) ** 2 / 2)
    
    dp_Tbarb = (inlet_loss(Q1, rho1, mu1) + inlet_loss(Q2, rho2, mu2) + outlet_loss ) + dp_friction


    return {
        'outlet_flow_rate': Q_out,
        'inlet_1_pressure_drop': inlet_loss(Q1, rho1, mu1),
        'inlet_2_pressure_drop': inlet_loss(Q2, rho2, mu2),
        'outlet_density': rho_out,
        'outlet_viscosity': mu_out,
        'total_pressure_drop': dp_Tbarb
    }


def mixer_feed_pressure_drop(
    Q,
    rho,
    mu,
    D_tube,
    D_connection,
    L_connection,
    D_mixer,
    L_mixer,
    K_LT,
):
    """Calculate losses from the inlet tube through and out of the mixer."""
    
    A_connection = np.pi * D_connection**2 / 4
    velocity_connection = Q / A_connection
    reynolds_connection = rho * velocity_connection * D_connection / mu
    friction_connection = friction_factor(reynolds_connection)

    
    dp_tube_connection_mixer= barb_dp (Q,mu, rho, D_connection,D_connection,  D_tube, D_mixer, L_connection)
    dp_mixer_connection_tube= barb_dp (Q,mu, rho, D_connection,D_connection,  D_mixer, D_tube, L_connection)


    area_mixer = np.pi * D_mixer**2 / 4
    velocity_mixer = Q / area_mixer
    reynolds_mixer = rho * velocity_mixer * D_mixer / mu

    fd_mixer = 64 / reynolds_mixer if reynolds_mixer <= 2000 else 0.3164 / reynolds_mixer**0.25

    K_mixer = fluids.mixing.K_motionless_mixer(
        K=K_LT,
        L=L_mixer,
        D=D_mixer,
        fd=fd_mixer,
    )       #K_mixer =K L/T *f_D*L/​D​

    dP_mixer = 16 * rho * velocity_mixer**2 / 2 

    """
    pressure_drop_mixer_branch = (
        contraction_loss
        + connection_friction
        + expansion_loss
        + outlet_contraction_loss
        + outlet_connection_friction
        + outlet_expansion_loss
        + pressure_drop_mixer
    )
     """
    pressure_drop_mixer_branch = (
            dp_tube_connection_mixer
            + dp_mixer_connection_tube
            + dP_mixer
        )

    return {
        'mixer_inlet_pressure_drop': dp_tube_connection_mixer,
        'mixer_inlet_pressure_drop': dp_mixer_connection_tube,
        'mixer_pressure_drop': dP_mixer,
        'pressure_drop_mixer_branch': pressure_drop_mixer_branch    
    } 


def syringe_pressure (T, P_in, P_out, fluid, D_syringe1,pump_speed, V_syringe1, D_luerlock, L_luerlock, D_tube, L_tube):
    if isinstance(fluid, str):
            #if fluid.lower() == 'mixer':
               # chem = mixer
           # else:
                chem = Chemical(fluid, T=T, P=(P_in+P_out)/2)
    else: # Assumes it's a pre-configured Mixture object
            chem = fluid
    rho = chem.rho
    mu = chem.mu


    # All geometry inputs use meters; areas are calculated in square meters.
    A_syringe = np.pi * (D_syringe1/2)**2
    A_luerlock = np.pi * (D_luerlock/2)**2
    A_tube = np.pi * (D_tube/2)**2

    # Plunger speed times area gives the volumetric flow rate from this syringe.
    Q_syringe = pump_speed * A_syringe 
    

    # V_syringe1 is the liquid volume being pushed, so this gives the travel time.
    Time_syringe = V_syringe1 / Q_syringe
    

    # Reynolds number identifies the flow regime used to select each friction factor.
    Re_syringe = reynolds_number(Q_syringe, D_syringe1, rho, mu)
    Re_luerlock = reynolds_number(Q_syringe, D_luerlock, rho, mu)
    Re_tube = reynolds_number(Q_syringe, D_tube, rho, mu)

    f_syringe = friction_factor(Re_syringe)
    f_luerlock = friction_factor(Re_luerlock)
    f_tube = friction_factor(Re_tube)

    # Darcy-Weisbach: dP = f*(L/D)*(rho*v^2/2). Here L = V/A is liquid-column length.
    dp_syringe_liquid = f_syringe * (V_syringe1 /A_syringe/ D_syringe1) * (rho * (Q_syringe/A_syringe)**2 / 2)
    dp_luerlock = f_luerlock * (L_luerlock/ D_luerlock) * (rho * (Q_syringe/A_luerlock)**2 / 2)
    dp_tube = f_tube * (L_tube/ D_tube) * (rho * (Q_syringe/A_tube)**2 / 2)

    k_syringe_luer = InletDrag_coefficient(Re_luerlock, D_luerlock, D_syringe1) #far big surface to small surface
    k_luer_tube = OutletDrag_coefficient(D_luerlock, D_tube)
    # Local contraction/expansion losses use K*(rho*v^2/2), with velocity in the luer lock.
    dp_connection = (k_syringe_luer + k_luer_tube) *  (rho * (Q_syringe/A_luerlock)**2 / 2)


    # Sum the steady-state losses from the syringe outlet to the first T-junction.
    dp_syringe_branch = dp_syringe_liquid +  dp_luerlock + dp_tube + dp_connection
    F_syringe = dp_syringe_branch * A_syringe

    output_arrays = {
            'surface': A_syringe,
            'density': rho,
            'viscosity': mu,
            'flow_rate': Q_syringe,
            'pressure_drop_syringe_branch': dp_syringe_branch,
            'force_syringe_branch': F_syringe,
            'time_finish': Time_syringe,
            'pressure_drop_syringe_liquid': dp_syringe_liquid,
            'pressure_drop_luerlock': dp_luerlock,
            'pressure_drop_tube': dp_tube,
            'pressure_drop_connection': dp_connection
        }

    return output_arrays



#%%
# Jupytext: {"Gernral function for Modules": {"extension": ".py", "format_name": "light", "format_version": "0.4"}}

from thermo import Chemical, Mixture
import numpy as np
import pandas as pd
import scipy.optimize as opt

# --- Centralized Fluid Definitions ---
#

def turbulent_flow_3(P_in, P_out, fluid, T, ID, length, ID2, length2,ID3, length3, ID_vessel, ID_orifice, ID_barb_1,ID_barb_2, num_barb, length_barb):  # if 2 barbs are different...Volume_in_vessel ??? if needed? and tube id change as well.....
    """ Calculate the flow rate under a pressure difference of P_in to P_out
    Consider the total length of tubing, the total number of barb fittings and the orifice"""
    chem = Chemical(fluid, T=T, P=(P_in+P_out)/2)
    rho = chem.rho
    mu = chem.mu

    # Tubing precalculations
    A = np.pi * (ID/2)**2
    A2 = np.pi * (ID2/2)**2   # for PA tubing in the vessel
    A3 = np.pi * (ID3/2)**2   # for PA tubing in the vessel
    ####Height_vessel = Volume_in
    # Barbed fittings precalculation
    A_barb = np.pi * (((ID_barb_1+ID_barb_2)/2)/2)**2
    A_orifice = np.pi * (ID_orifice/2)**2


    def equation(Q):
         # Tubing 1--- 6.35mm ID Multilayer PA
        V = Q / A
        Re = reynolds_number(Q, ID, rho, mu)
        f = friction_factor(Re)
        #f_inlet_vessel = InletDrag_coefficient(Re,ID2,ID_vessel)+1   #### verify why +1 #liquid goes in from vessel to 6mm PA tubing
        dp_tube = f * (length / ID) * (rho * V**2 / 2)

         # Tubing 2--- 6mm ID PA long in
        V2 = Q / A2
        Re2 = reynolds_number(Q, ID2, rho, mu)
        f2 = friction_factor(Re2)
        dp_tube2 = f2 * (length2 / ID2) * (rho * V2**2 / 2)

         # Tubing 3--- 6mm ID PA short in
        V3 = Q / A3
        Re3 = reynolds_number(Q, ID3, rho, mu)
        f3 = friction_factor(Re3)
        dp_tube3 = f3 * (length3 / ID3) * (rho * V3**2 / 2)

        # step 1: Driven the liquid from vessel into the 6mm PA tubing ID2
        f_inlet_vessel = InletDrag_coefficient(Re2,ID2,ID_vessel)   #### verify why +1

        # step 2: Driven the liquid from 6mm PA tubing ID2 to the barb ID_1
        
        V_barb = Q / A_barb
        Re_barb = reynolds_number(Q,(ID_barb_1+ID_barb_2)/2,rho,mu)
        f_barb = friction_factor(Re_barb)
        

        #dp_barb1 = 1* ((f_inlet_barb + f_outlet_barb ) * (rho * V**2 / 2)+1.5 * (rho * V**2 / 2))

        # Fittings/barbed
        f_inlet_barb = get_zeta_in(ID_barb_1, ID2, regime='turbulent', Re = (Re2 + Re_barb)/2) #InletDrag_coefficient(Re2,ID_barb_1,ID2)   #liquid goes in from 6mm PA tubing to barb
        # step 3: Driven the liquid from the barb ID_2 to the 6.35 multilater tubing
        f_outlet_barb = OutletDrag_coefficient(ID_barb_2, ID)   #liquid goes out from barb to the 6.35 multilayer tubing

        dp_barb1 = 1* ((f_inlet_barb + f_outlet_barb + 1.15 ) * (rho *  V_barb**2 / 2)        
                             +(f_barb ) * ((length_barb) / ((ID_barb_1+ID_barb_2)/2)) * (rho * V_barb**2 / 2))


        dp_barb3 = barb_dp(Q,mu, rho,ID_barb_1,ID_barb_2, ID, ID2, length_barb)


        # step 4: Driven the liquid from 6.35mm PAmulti tubing ID to the barb ID_2
        f_inlet_barb2 = get_zeta_in(ID_barb_2, ID, regime='turbulent', Re = (Re + Re_barb)/2) #InletDrag_coefficient(Re,ID_barb_2,ID)
        # step 5: Driven the liquid from the barb ID_1 to the ID2 6mm PA tubing
        f_outlet_barb2 = OutletDrag_coefficient(ID_barb_1, ID2)

        dp_barb2 = 1* ((f_inlet_barb2 + f_outlet_barb2 + 1.15 ) * (rho * V_barb**2 / 2)
                             +(f_barb) * ((length_barb) / ((ID_barb_1+ID_barb_2)/2)) * (rho * V_barb**2 / 2))

        #dp_barb2 = 5* ((f_inlet_barb2 + f_outlet_barb2 ) * (rho * V**2 / 2) +1.5 * (rho * V**2 / 2))


        ##### dp_h = rho * 9.81 * (0.05) need or not???
        #dp_Yconver = 1*0.1* (rho * (2*V)**2 / 2) # suppose there is one Y connector combiand MTBE and acetone , Ref VDI p.1070   for later
        #dp_Tconver = 1*1.56* (rho * (2*V)**2 / 2) # suppose there is one T connector combiand MTBE and acetone , Ref VDI p.1070   for later
        #dp_Lbarb = 1*1.15* (rho * (2*V)**2 / 2)

        # Orifice
        V_orifice = Q / A_orifice
        #f_orifice = OrificeDrag_coefficient(ID,ID_orifice)  # ???
        #dp_orifice = (f_orifice)* (rho * V_orifice**2 / 2)
        #dp_orifice = (rho/2.0) * (Q/(0.62*3.14*(ID_orifice/2)**2))**2 #wht I have this??? verify

        Re_orifice = reynolds_number(Q, ID_orifice, rho, mu)
        f_orifice = friction_factor(Re_orifice)

        f_inlet_orifice = get_zeta_in( ID_orifice, ID2, regime='turbulent', Re = (Re + Re_orifice)/2)
        f_outlet_orifice = OutletDrag_coefficient(ID_orifice, ID)

        dp_orifice = (f_inlet_orifice + f_outlet_orifice + f_orifice) * (rho * V_orifice**2 / 2)

        ####Height_vessel = Volume_in_vessel/(ID_vessel/2)**2/np.pi # Height of the liquid in the vessel ???if need???
        # step 6: Driven the liquid from the ID2 6mm PA tubing to the vessel
        f_outlet_vessel = OutletDrag_coefficient(ID2, ID_vessel)   #liquid goes out from 6.35 multilayer tubing to the vessel
        #dp_in_out = (20+1)* (rho * V**2 / 2)+20* (rho * V**2 / 2) # suppose the inlet from vessel to the tube is very sharp, worst condition, and the outlet vessel is empty, no friction Ref, VDI p.1067 Fig.4
        dp_in_out = (f_outlet_vessel +  f_inlet_vessel)* (rho * V2**2 / 2) # suppose the inlet from vessel to the tube is very sharp, worst condition, and the outlet vessel is empty, no friction Ref, VDI p.1067 Fig.4

        f_bend_tube = 4 * get_zeta_bend(Re, 11.71)
        dp_bend_tube = f_bend_tube * (rho * V **2 / 2)

        # Calculate the final residual (this is what fsolve needs to hit 0)
        return P_in - P_out - (dp_tube + dp_tube2 + dp_tube3 + dp_in_out + dp_barb1 + dp_barb2 + dp_orifice + dp_bend_tube)
        #residual = P_in- P_out - (50 * (rho * V**2 / 2))

        #return P_in- P_out - (50 * (rho * V**2 / 2))
        #return P_in- P_out - (dp_tube + dp_tube2 + dp_tube3 + dp_in_out + dp_barb1 + dp_barb2 + dp_orifice)
          #return P_inlet- P_outlet + dp_h - wash_time* (dp_tube  + dp_in_out + dp_barb+dp_h+dp_Yconver) dp_h???

    Q_guess = 0.001  # Initial guess for flow rate (m³/s)
    #num_valves_guess = 0
    return opt.fsolve(equation, Q_guess)[0]


import fluids
def turbulent_flow_fluids(P_in, P_out, fluid, T, ID, length, ID2, length2,ID3, length3, ID_vessel, ID_orifice, ID_barb_1,ID_barb_2, num_barb, length_barb):  # if 2 barbs are different...Volume_in_vessel ??? if needed? and tube id change as well.....
    """ Calculate the flow rate under a pressure difference of P_in to P_out
    Consider the total length of tubing, the total number of barb fittings and the orifice"""
    chem = Chemical(fluid, T=T, P=(P_in+P_out)/2)
    rho = chem.rho
    mu = chem.mu

    # Tubing precalculations
    A = np.pi * (ID/2)**2
    A2 = np.pi * (ID2/2)**2   # for PA tubing in the vessel
    A3 = np.pi * (ID3/2)**2   # for PA tubing in the vessel
    A_orifice = np.pi * (ID_orifice/2)**2
    ####Height_vessel = Volume_in
    # Barbed fittings precalculation
    A_barb = np.pi * (((ID_barb_1+ID_barb_2)/2)/2)**2


    def calculate_all_physics(Q):
        # Tubing 1--- 6.35mm ID Multilayer PA
        V = Q / A
        Re = fluids.Reynolds(V, ID, rho, mu)
        f = fluids.friction_factor(Re, 0.00045)
        #f_inlet_vessel = InletDrag_coefficient(Re,ID2,ID_vessel)+1   #### verify why +1 #liquid goes in from vessel to 6mm PA tubing
        dp_tube = f * (length / ID) * (rho * V**2 / 2)
        #dp_tube = fluids.one_phase_dP(Q*rho, rho, mu, ID, 0.00045, length)

         # Tubing 2--- 6mm ID PA long in
        V2 = Q / A2
        Re2 = fluids.Reynolds(V2, ID2, rho, mu)
        f2 = fluids.friction_factor(Re2, 0.00045)
        dp_tube2 = f2 * (length2 / ID2) * (rho * V2**2 / 2)
        #dp_tube2 = fluids.one_phase_dP(Q*rho, rho, mu, ID2, 0.00045, length2)

         # Tubing 3--- 6mm ID PA short in
        V3 = Q / A3
        Re3 = fluids.Reynolds(V3, ID3, rho, mu)
        f3 = fluids.friction_factor(Re3, 0.00045)
        dp_tube3 = f3 * (length3 / ID3) * (rho * V3**2 / 2)
        #dp_tube3 = fluids.one_phase_dP(Q*rho, rho, mu, ID3, 0.00045, length3)

        # step 1: Driven the liquid from vessel into the 6mm PA tubing ID2
        #f_inlet_vessel = InletDrag_coefficient(Re2,ID2,ID_vessel)+1   #### verify why +1
        f_inlet_vessel = fluids.entrance_distance(ID2, 0.0001, 0.05,method='Idelchik'  )
        #f_inlet_vessel = fluids.contraction_sharp(ID_vessel, ID2) 


        # step 2: Driven the liquid from 6mm PA tubing ID2 to the barb ID_1
        # Fittings/barbed
        #f_inlet_barb = InletDrag_coefficient(Re2,ID_barb_1,ID2)+1   #liquid goes in from 6mm PA tubing to barb
        f_inlet_barb = fluids.contraction_sharp(ID2, ID_barb_1)
        # step 3: Driven the liquid from the barb ID_2 to the 6.35 multilater tubing
        #f_outlet_barb = OutletDrag_coefficient(ID_barb_2, ID)   #liquid goes out from barb to the 6.35 multilayer tubing
        f_outlet_barb = fluids.diffuser_sharp(ID_barb_2, ID)
        #f_outlet_barb = (1.0 - (ID_barb_2 / ID)**2)**2

      
        V_barb = Q / A_barb
        Re_barb = fluids.Reynolds(V_barb, (ID_barb_1+ID_barb_2)/2, rho, mu)
        f_barb = fluids.friction_factor(Re_barb, 0.00045)
        f_L_bend = fluids.bend_miter(90,ID_barb_1, Re_barb, 0.00045)  # for the bend in the barb, need to verify the Re number and the roughness
        #dp_barb1 = 1* ((f_inlet_barb + f_outlet_barb ) * (rho * V**2 / 2)
                             #+f_barb * ((length_barb) / ((ID_barb_1+ID_barb_2)/2)) * (rho * V**2 / 2))

        #dp_barb1 = 1* ((f_inlet_barb + f_outlet_barb ) * (rho * V**2 / 2) +1.5 * (rho * V**2 / 2))
        dp_barb1 = ( f_inlet_barb + f_outlet_barb + f_L_bend)* (rho * V_barb**2 / 2) + f_barb * ((length_barb) / ((ID_barb_1+ID_barb_2)/2)) * (rho * V_barb**2 / 2)

        # step 4: Driven the liquid from 6.35mm PAmulti tubing ID to the barb ID_2
        #f_inlet_barb2 = InletDrag_coefficient(Re,ID_barb_2,ID)+1
        f_inlet_barb2 = fluids.contraction_sharp(Di1 = ID, Di2 = ID_barb_2, method = 'Rennels')
        # step 5: Driven the liquid from the barb ID_1 to the ID2 6mm PA tubing
        #f_outlet_barb2 = OutletDrag_coefficient(ID_barb_1, ID2)
        f_outlet_barb2 = fluids.diffuser_sharp(ID_barb_1, ID2)
        

        #dp_barb2 = 1* ((f_inlet_barb2 + f_outlet_barb2 ) * (rho * V**2 / 2)
                            # +f_barb * ((length_barb) / ((ID_barb_1+ID_barb_2)/2)) * (rho * V**2 / 2))

       # dp_barb2 = 5* ((f_inlet_barb2 + f_outlet_barb2 ) * (rho * V**2 / 2) +1.5 * (rho * V**2 / 2))
        
        dp_barb2 = ( f_inlet_barb2 + f_outlet_barb2 + f_L_bend)* (rho * V_barb**2 / 2) + f_barb * ((length_barb) / ((ID_barb_1+ID_barb_2)/2)) * (rho * V_barb**2 / 2)


        ##### dp_h = rho * 9.81 * (0.05) need or not???
        #dp_Yconver = 1*0.1* (rho * (2*V)**2 / 2) # suppose there is one Y connector combiand MTBE and acetone , Ref VDI p.1070   for later
        #dp_Tconver = 1*1.56* (rho * (2*V)**2 / 2) # suppose there is one T connector combiand MTBE and acetone , Ref VDI p.1070   for later
        #dp_Lbarb = 1*1.15* (rho * (2*V)**2 / 2)

        # Orifice
        #f_orifice = OrificeDrag_coefficient(ID,ID_orifice)  # ???
        #f_orifice = fluids.C_Reader_Harris_Gallagher(ID, ID_orifice, 0.00045)
        #dp_orifice = (f_orifice)* (rho * V**2 / 2)
        #dp_orifice = fluids.orifice_dP(ID, ID_orifice, Q*rho, rho, mu, 0.00045)
        ####dp_orifice = (rho/2.0) * (Q/(0.62*3.14*(ID_orifice/2)**2))**2 wht I have this??? verify


        V_orifice = Q / A_orifice
        Re_orifice = fluids.Reynolds(V_orifice, ID_orifice, rho, mu)
        f_orifice = fluids.friction_factor(Re_orifice, 0.00045)
            
        f_inlet_orifice = fluids.contraction_sharp(ID2, ID_orifice)       
        f_outlet_orifice= fluids.diffuser_sharp(ID_orifice, ID)   
        
        dp_orifice = ( f_inlet_orifice + f_outlet_orifice)* (rho *  V_orifice**2 / 2) + f_orifice * (length_barb / ID_orifice) * (rho * V_orifice**2 / 2)

        ####Height_vessel = Volume_in_vessel/(ID_vessel/2)**2/np.pi # Height of the liquid in the vessel ???if need???

        # step 6: Driven the liquid from the ID2 6mm PA tubing to the vessel
        #f_outlet_vessel = OutletDrag_coefficient(ID2, ID_vessel)   #liquid goes out from 6.35 multilayer tubing to the vessel
        f_outlet_vessel = fluids.exit_normal()
        #f_outlet_vessel = fluids.diffuser_sharp(ID2, ID_vessel)
        dp_in_out = (f_outlet_vessel +  f_inlet_vessel)* (rho * V2**2 / 2) # suppose the inlet from vessel to the tube is very sharp, worst condition, and the outlet vessel is empty, no friction Ref, VDI p.1067 Fig.4

        f_bend_tube= 4 * fluids.bend_rounded_Miller(Di = ID, angle = 90, bend_diameters=11.71, Re = Re)  # for the bend in the barb, need to verify the Re number and the roughness
        dp_bend_tube = f_bend_tube * (rho * V **2 / 2)


        residual = P_in- P_out - (dp_tube + dp_tube2 + dp_tube3 + dp_in_out + dp_barb1 + dp_barb2 + dp_orifice + dp_bend_tube)
        #residual = P_in- P_out - (50 * (rho * V**2 / 2))

        # Return EVERYTHING as a packaged dictionary
        return {
            "residual": residual,
            "Reynolds_number_tube1": Re,
            "Reynolds_number_tube2": Re2,
            "Reynolds_number_tube3": Re3,
            "f_tube1": f,
            "f_tube2": f2,
            "f_tube3": f3,
            "f_bend_tube": f_bend_tube,
            "f_inlet_vessel": f_inlet_vessel,
            "f_outlet_vessel": f_outlet_vessel,
            "f_inlet_barb": f_inlet_barb,
            "f_outlet_barb": f_outlet_barb,
            "f_barb": f_barb,
            "f_L_bend": f_L_bend,
            "f_inlet_barb2": f_inlet_barb,
            "f_outlet_barb2": f_outlet_barb2,
            "f_orifice": f_orifice,
            "f_inlet_orifice": f_inlet_orifice,
            "f_outlet_orifice": f_outlet_orifice,
            "dp_tube_total": dp_tube + dp_tube2 + dp_tube3,
            "dp_barbs_total": dp_barb1 + dp_barb2,
            "dp_orifice": dp_orifice,
            "dp_bend_tube": dp_bend_tube
        }

       
 # --- STEP 2: Make a tiny wrapper function for the solver ---
    def equation_for_solver(Q):
        return calculate_all_physics(Q)["residual"]

    # --- STEP 3: The Ultimate Bounded Solver ---
    # First, a physics reality check:
    if P_in <= P_out:
        print("⚠️ WARNING: P_in is NOT greater than P_out. Flow is physically impossible!")
        return calculate_all_physics(1e-10) # Return effectively zero flow

    try:
        # root_scalar absolutely guarantees the equation balances to zero.
        # It searches exclusively between 1e-10 m^3/s (0 mL/min) and 0.005 m^3/s (300,000 mL/min)
        result = opt.root_scalar(equation_for_solver, bracket=[1e-10, 0.01], method='brentq')
        
        if result.converged:
            Q_final = result.root
        else:
            print("⚠️ WARNING: Solver stopped early!")
            Q_final = result.root
            
    except ValueError:
        print("❌ SOLVER FAILED: Your available pressure is either too low to overcome static friction, or too high for this bracket!")
        Q_final = 1e-10

    # --- STEP 4: Run the math ONE LAST TIME with the perfectly balanced Q ---
    final_data = calculate_all_physics(Q_final)
    final_data["Flow_Rate_Q"] = Q_final
    
    return final_data



def turbulent_flow_X(P_in, P_out, fluid, T, ID, length, ID2, length2,ID3, length3, ID_vessel, ID_orifice, ID_barb_1,ID_barb_2, num_barb, length_barb):  # if 2 barbs are different...Volume_in_vessel ??? if needed? and tube id change as well.....
    """ Calculate the flow rate under a pressure difference of P_in to P_out
    Consider the total length of tubing, the total number of barb fittings and the orifice"""
    if isinstance(fluid, str):
        #if fluid.lower() == 'mixer':
           # chem = mixer
       # else:
            chem = Chemical(fluid, T=T, P=(P_in+P_out)/2)
    else: # Assumes it's a pre-configured Mixture object
        chem = fluid
    rho = chem.rho
    mu = chem.mu

    # Tubing precalculations
    A = np.pi * (ID/2)**2
    A2 = np.pi * (ID2/2)**2   # for PA tubing in the vessel
    A3 = np.pi * (ID3/2)**2   # for PA tubing in the vessel
    ####Height_vessel = Volume_in
    # Barbed fittings precalculation
    A_barb = np.pi * (((ID_barb_1+ID_barb_2)/2)/2)**2
    A_orifice = np.pi * (ID_orifice/2)**2


    def calculate_all_physics(Q):
        # Tubing 1--- 6.35mm ID Multilayer PA
        V = Q / A
        Re = reynolds_number(Q, ID, rho, mu)
        f = friction_factor(Re)
        #f_inlet_vessel = InletDrag_coefficient(Re,ID2,ID_vessel)+1   #### verify why +1 #liquid goes in from vessel to 6mm PA tubing
        dp_tube = f * (length / ID) * (rho * V**2 / 2)

         # Tubing 2--- 6mm ID PA long in
        V2 = Q / A2
        Re2 = reynolds_number(Q, ID2, rho, mu)
        f2 = friction_factor(Re2)
        dp_tube2 = f2 * (length2 / ID2) * (rho * V2**2 / 2)

         # Tubing 3--- 6mm ID PA short in
        V3 = Q / A3
        Re3 = reynolds_number(Q, ID3, rho, mu)
        f3 = friction_factor(Re3)
        dp_tube3 = f3 * (length3 / ID3) * (rho * V3**2 / 2)

        # step 1: Driven the liquid from vessel into the 6mm PA tubing ID2
        f_inlet_vessel = InletDrag_coefficient(Re2,ID2,ID_vessel)   #### verify why +1

        # step 2: Driven the liquid from 6mm PA tubing ID2 to the barb ID_1
        
        V_barb = Q / A_barb
        Re_barb = reynolds_number(Q,(ID_barb_1+ID_barb_2)/2,rho,mu)
        f_barb = friction_factor(Re_barb)
        

        #dp_barb1 = 1* ((f_inlet_barb + f_outlet_barb ) * (rho * V**2 / 2)+1.5 * (rho * V**2 / 2))

        # Fittings/barbed
        f_inlet_barb = get_zeta_in(ID_barb_1, ID2, regime='turbulent', Re = (Re2 + Re_barb)/2) #InletDrag_coefficient(Re2,ID_barb_1,ID2)   #liquid goes in from 6mm PA tubing to barb
        # step 3: Driven the liquid from the barb ID_2 to the 6.35 multilater tubing
        f_outlet_barb = OutletDrag_coefficient(ID_barb_2, ID)   #liquid goes out from barb to the 6.35 multilayer tubing

        dp_barb1 = 1* ((f_inlet_barb + f_outlet_barb + 1.15 ) * (rho *  V_barb**2 / 2)        
                             +(f_barb ) * ((length_barb) / ((ID_barb_1+ID_barb_2)/2)) * (rho * V_barb**2 / 2))


        dp_barb3 = barb_dp (Q,mu, ID_barb_1,ID_barb_2, ID2, ID, length_barb)
        
        # step 4: Driven the liquid from 6.35mm PAmulti tubing ID to the barb ID_2
        f_inlet_barb2 = get_zeta_in(ID_barb_2, ID, regime='turbulent', Re = (Re + Re_barb)/2) #InletDrag_coefficient(Re,ID_barb_2,ID)
        # step 5: Driven the liquid from the barb ID_1 to the ID2 6mm PA tubing
        f_outlet_barb2 = OutletDrag_coefficient(ID_barb_1, ID2)

        dp_barb2 = ((f_inlet_barb2 + f_outlet_barb2 + 1.15 ) * (rho * V_barb**2 / 2)
                             +(f_barb) * ((length_barb) / ((ID_barb_1+ID_barb_2)/2)) * (rho * V_barb**2 / 2))

        #dp_barb2 = 5* ((f_inlet_barb2 + f_outlet_barb2 ) * (rho * V**2 / 2) +1.5 * (rho * V**2 / 2))


        ##### dp_h = rho * 9.81 * (0.05) need or not???
        #dp_Yconver = 1*0.1* (rho * (2*V)**2 / 2) # suppose there is one Y connector combiand MTBE and acetone , Ref VDI p.1070   for later
        #dp_Tconver = 1*1.56* (rho * (2*V)**2 / 2) # suppose there is one T connector combiand MTBE and acetone , Ref VDI p.1070   for later
        #dp_Lbarb = 1*1.15* (rho * (2*V)**2 / 2)

        # Orifice
        V_orifice = Q / A_orifice
        #f_orifice = OrificeDrag_coefficient(ID,ID_orifice)  # ???
        #dp_orifice = (f_orifice)* (rho * V_orifice**2 / 2)
        #dp_orifice = (rho/2.0) * (Q/(0.62*3.14*(ID_orifice/2)**2))**2 #wht I have this??? verify

        Re_orifice = reynolds_number(Q, ID_orifice, rho, mu)
        f_orifice = friction_factor(Re_orifice)

        f_inlet_orifice = get_zeta_in( ID_orifice, ID2, regime='turbulent', Re = (Re + Re_orifice)/2)
        f_outlet_orifice = OutletDrag_coefficient(ID_orifice, ID)

        dp_orifice = (f_inlet_orifice + f_outlet_orifice + f_orifice) * (rho * V_orifice**2 / 2)

        ####Height_vessel = Volume_in_vessel/(ID_vessel/2)**2/np.pi # Height of the liquid in the vessel ???if need???
        # step 6: Driven the liquid from the ID2 6mm PA tubing to the vessel
        f_outlet_vessel = OutletDrag_coefficient(ID2, ID_vessel)   #liquid goes out from 6.35 multilayer tubing to the vessel
        #dp_in_out = (20+1)* (rho * V**2 / 2)+20* (rho * V**2 / 2) # suppose the inlet from vessel to the tube is very sharp, worst condition, and the outlet vessel is empty, no friction Ref, VDI p.1067 Fig.4
        dp_in_out = (f_outlet_vessel +  f_inlet_vessel)* (rho * V2**2 / 2) # suppose the inlet from vessel to the tube is very sharp, worst condition, and the outlet vessel is empty, no friction Ref, VDI p.1067 Fig.4

        f_bend_tube = 4 * get_zeta_bend(Re, 11.71)
        dp_bend_tube = f_bend_tube * (rho * V **2 / 2)

        # Calculate the final residual (this is what fsolve needs to hit 0)
        residual = P_in - P_out - (dp_tube + dp_tube2 + dp_tube3 + dp_in_out + dp_barb1 + dp_barb2 + dp_orifice + dp_bend_tube)
        #residual = P_in- P_out - (50 * (rho * V**2 / 2))
        # Return EVERYTHING as a packaged dictionary
        return {
            "residual": residual,
            "f_tube1": f,
            "f_tube2": f2,
            "f_tube3": f3,
            "f_bend_tube": f_bend_tube,
            "f_inlet_vessel": f_inlet_vessel,
            "f_outlet_vessel": f_outlet_vessel,
            "f_inlet_barb": f_inlet_barb,
            "f_outlet_barb": f_outlet_barb,
            "f_barb": f_barb,
            "f_L_bend": 1.15,
            "f_inlet_barb2": f_inlet_barb,
            "f_outlet_barb2": f_outlet_barb2,
            "f_orifice": f_orifice,
            "f_inlet_orifice": f_inlet_orifice,
            "f_outlet_orifice": f_outlet_orifice,
            "dp_tube_total": dp_tube + dp_tube2 + dp_tube3,
            "dp_barbs_total": dp_barb1 + dp_barb2,
            "dp_orifice": dp_orifice,
            "dp_bend_tube": dp_bend_tube
        }
    



    # --- STEP 2: Make a tiny wrapper function for the solver ---
    def equation_for_solver(Q):
        return calculate_all_physics(Q)["residual"]

    # --- STEP 3: The Ultimate Bounded Solver ---
    # First, a physics reality check:
    if P_in <= P_out:
        print("⚠️ WARNING: P_in is NOT greater than P_out. Flow is physically impossible!")
        return calculate_all_physics(1e-10) # Return effectively zero flow

    try:
        # root_scalar absolutely guarantees the equation balances to zero.
        # It searches exclusively between 1e-10 m^3/s (0 mL/min) and 0.01 m^3/s (600,000 mL/min)
        result = opt.root_scalar(equation_for_solver, bracket=[1e-10, 0.01], method='brentq')
        
        if result.converged:
            Q_final = result.root
        else:
            print("⚠️ WARNING: Solver stopped early!")
            Q_final = result.root
            
    except ValueError:
        print("❌ SOLVER FAILED: Your available pressure is either too low to overcome static friction, or too high for this bracket!")
        Q_final = 1e-10

    # --- STEP 4: Run the math ONE LAST TIME with the perfectly balanced Q ---
    final_data = calculate_all_physics(Q_final)
    final_data["Flow_Rate_Q"] = Q_final
    
    return final_data


def turbulent_flow_X_numB(P_in, P_out, fluid, T, ID, length, ID2, length2,ID3, length3, ID_vessel, ID_orifice, ID_barb_1,ID_barb_2, num_barb, length_barb):  # if 2 barbs are different...Volume_in_vessel ??? if needed? and tube id change as well.....
    """ Calculate the flow rate under a pressure difference of P_in to P_out
    Consider the total length of tubing, the total number of barb fittings and the orifice"""
    if isinstance(fluid, str):
        #if fluid.lower() == 'mixer':
           # chem = mixer
       # else:
            chem = Chemical(fluid, T=T, P=(P_in+P_out)/2)
    else: # Assumes it's a pre-configured Mixture object
        chem = fluid
    rho = chem.rho
    mu = chem.mu

    # Tubing precalculations
    A = np.pi * (ID/2)**2
    A2 = np.pi * (ID2/2)**2   # for PA tubing in the vessel
    A3 = np.pi * (ID3/2)**2   # for PA tubing in the vessel
    ####Height_vessel = Volume_in
    # Barbed fittings precalculation
    A_barb = np.pi * (((ID_barb_1+ID_barb_2)/2)/2)**2
    A_orifice = np.pi * (ID_orifice/2)**2


    def calculate_all_physics(Q):
        # Tubing 1--- 6.35mm ID Multilayer PA
        V = Q / A
        Re = reynolds_number(Q, ID, rho, mu)
        f = friction_factor(Re)
        #f_inlet_vessel = InletDrag_coefficient(Re,ID2,ID_vessel)+1   #### verify why +1 #liquid goes in from vessel to 6mm PA tubing
        dp_tube = f * (length / ID) * (rho * V**2 / 2)

         # Tubing 2--- 6mm ID PA long in
        V2 = Q / A2
        Re2 = reynolds_number(Q, ID2, rho, mu)
        f2 = friction_factor(Re2)
        dp_tube2 = f2 * (length2 / ID2) * (rho * V2**2 / 2)

         # Tubing 3--- 6mm ID PA short in
        V3 = Q / A3
        Re3 = reynolds_number(Q, ID3, rho, mu)
        f3 = friction_factor(Re3)
        dp_tube3 = f3 * (length3 / ID3) * (rho * V3**2 / 2)

        # step 1: Driven the liquid from vessel into the 6mm PA tubing ID2
        f_inlet_vessel = InletDrag_coefficient(Re2,ID2,ID_vessel)   #### verify why +1

        # step 2: Driven the liquid from 6mm PA tubing ID2 to the barb ID_1
        
        V_barb = Q / A_barb
        Re_barb = reynolds_number(Q,(ID_barb_1+ID_barb_2)/2,rho,mu)
        f_barb = friction_factor(Re_barb)
        

        #dp_barb1 = 1* ((f_inlet_barb + f_outlet_barb ) * (rho * V**2 / 2)+1.5 * (rho * V**2 / 2))

        # Fittings/barbed
        f_inlet_barb = get_zeta_in(ID_barb_1, ID2, regime='turbulent', Re = (Re2 + Re_barb)/2) #InletDrag_coefficient(Re2,ID_barb_1,ID2)   #liquid goes in from 6mm PA tubing to barb
        # step 3: Driven the liquid from the barb ID_2 to the 6.35 multilater tubing
        f_outlet_barb = OutletDrag_coefficient(ID_barb_2, ID)   #liquid goes out from barb to the 6.35 multilayer tubing

        dp_barb1 = num_barb* ((f_inlet_barb + f_outlet_barb + 1.15 ) * (rho *  V_barb**2 / 2)        
                             +(f_barb ) * ((length_barb) / ((ID_barb_1+ID_barb_2)/2)) * (rho * V_barb**2 / 2))


        # step 4: Driven the liquid from 6.35mm PAmulti tubing ID to the barb ID_2
        f_inlet_barb2 = get_zeta_in(ID_barb_2, ID, regime='turbulent', Re = (Re + Re_barb)/2) #InletDrag_coefficient(Re,ID_barb_2,ID)
        # step 5: Driven the liquid from the barb ID_1 to the ID2 6mm PA tubing
        f_outlet_barb2 = OutletDrag_coefficient(ID_barb_1, ID2)

        dp_barb2 =0* ((f_inlet_barb2 + f_outlet_barb2 + 1.15 ) * (rho * V_barb**2 / 2)
                             +(f_barb) * ((length_barb) / ((ID_barb_1+ID_barb_2)/2)) * (rho * V_barb**2 / 2))

        #dp_barb2 = 5* ((f_inlet_barb2 + f_outlet_barb2 ) * (rho * V**2 / 2) +1.5 * (rho * V**2 / 2))


        ##### dp_h = rho * 9.81 * (0.05) need or not???
        #dp_Yconver = 1*0.1* (rho * (2*V)**2 / 2) # suppose there is one Y connector combiand MTBE and acetone , Ref VDI p.1070   for later
        #dp_Tconver = 1*1.56* (rho * (2*V)**2 / 2) # suppose there is one T connector combiand MTBE and acetone , Ref VDI p.1070   for later
        #dp_Lbarb = 1*1.15* (rho * (2*V)**2 / 2)

        # Orifice
        V_orifice = Q / A_orifice
        #f_orifice = OrificeDrag_coefficient(ID,ID_orifice)  # ???
        #dp_orifice = (f_orifice)* (rho * V_orifice**2 / 2)
        #dp_orifice = (rho/2.0) * (Q/(0.62*3.14*(ID_orifice/2)**2))**2 #wht I have this??? verify

        Re_orifice = reynolds_number(Q, ID_orifice, rho, mu)
        f_orifice = friction_factor(Re_orifice)

        f_inlet_orifice = get_zeta_in( ID_orifice, ID2, regime='turbulent', Re = (Re + Re_orifice)/2)
        f_outlet_orifice = OutletDrag_coefficient(ID_orifice, ID)

        dp_orifice = (f_inlet_orifice + f_outlet_orifice + f_orifice) * (rho * V_orifice**2 / 2)

        ####Height_vessel = Volume_in_vessel/(ID_vessel/2)**2/np.pi # Height of the liquid in the vessel ???if need???
        # step 6: Driven the liquid from the ID2 6mm PA tubing to the vessel
        f_outlet_vessel = OutletDrag_coefficient(ID2, ID_vessel)   #liquid goes out from 6.35 multilayer tubing to the vessel
        #dp_in_out = (20+1)* (rho * V**2 / 2)+20* (rho * V**2 / 2) # suppose the inlet from vessel to the tube is very sharp, worst condition, and the outlet vessel is empty, no friction Ref, VDI p.1067 Fig.4
        dp_in_out = (f_outlet_vessel +  f_inlet_vessel)* (rho * V2**2 / 2) # suppose the inlet from vessel to the tube is very sharp, worst condition, and the outlet vessel is empty, no friction Ref, VDI p.1067 Fig.4

        f_bend_tube = 4 * get_zeta_bend(Re, 11.71)
        dp_bend_tube = f_bend_tube * (rho * V **2 / 2)

        # Calculate the final residual (this is what fsolve needs to hit 0)
        residual = P_in - P_out - (dp_tube + dp_tube2 + dp_tube3 + dp_in_out + dp_barb1 + dp_barb2 + dp_orifice + dp_bend_tube)
        #residual = P_in- P_out - (50 * (rho * V**2 / 2))
        # Return EVERYTHING as a packaged dictionary
        return {
            "residual": residual,
            "f_tube1": f,
            "f_tube2": f2,
            "f_tube3": f3,
            "f_bend_tube": f_bend_tube,
            "f_inlet_vessel": f_inlet_vessel,
            "f_outlet_vessel": f_outlet_vessel,
            "f_inlet_barb": f_inlet_barb,
            "f_outlet_barb": f_outlet_barb,
            "f_barb": f_barb,
            "f_L_bend": 1.15,
            "f_inlet_barb2": f_inlet_barb,
            "f_outlet_barb2": f_outlet_barb2,
            "f_orifice": f_orifice,
            "f_inlet_orifice": f_inlet_orifice,
            "f_outlet_orifice": f_outlet_orifice,
            "dp_tube_total": dp_tube + dp_tube2 + dp_tube3,
            "dp_barbs_total": dp_barb1 + dp_barb2,
            "dp_orifice": dp_orifice,
            "dp_bend_tube": dp_bend_tube
        }




    # --- STEP 2: Make a tiny wrapper function for the solver ---
    def equation_for_solver(Q):
        return calculate_all_physics(Q)["residual"]

    # --- STEP 3: The Ultimate Bounded Solver ---
    # First, a physics reality check:
    if P_in <= P_out:
        print("⚠️ WARNING: P_in is NOT greater than P_out. Flow is physically impossible!")
        return calculate_all_physics(1e-10) # Return effectively zero flow

    try:
        # root_scalar absolutely guarantees the equation balances to zero.
        # It searches exclusively between 1e-10 m^3/s (0 mL/min) and 0.01 m^3/s (600,000 mL/min)
        result = opt.root_scalar(equation_for_solver, bracket=[1e-10, 0.01], method='brentq')
        
        if result.converged:
            Q_final = result.root
        else:
            print("⚠️ WARNING: Solver stopped early!")
            Q_final = result.root
            
    except ValueError:
        print("❌ SOLVER FAILED: Your available pressure is either too low to overcome static friction, or too high for this bracket!")
        Q_final = 1e-10

    # --- STEP 4: Run the math ONE LAST TIME with the perfectly balanced Q ---
    final_data = calculate_all_physics(Q_final)
    final_data["Flow_Rate_Q"] = Q_final
    
    return final_data
