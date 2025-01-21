# !! Note!!
# in Matlab code, xh is computed in f() and placed in x(kk+1,:)
# here, xh is computed in g() and placed in y(kk,:)
# consequently, there is a difference between the Matlab and Python implementation in y_{10}
# another differences due to discretization method and sample period choice 

import numpy as np
import pdb #pdb.set_trace()

def DefineParameters():
# Model parameters
    p = {}
      
    p['fc.h_melt']     = 334   # [kJ/kg] latent heat of fusion of water
    p['fc.c_water']    = 4.2   # [kJ/kg.gC] specific heat of water
    p['fc.c_ice']      = 2.2   # [kJ/kg.gC] specific heat of ice
    p['fc.rho_water']  = 1000  # [kg/m3] density of water @ +20 gC
    p['fc.rho_ice']    = 917   # [kg/m3] density of ice @ -4 gC 
 
    # KOUDEOPWEKKING MODELPARAMETERS
    p['ko.pm(1)'] = 50  # [kW] max. nominaal thermisch koelvermogen van compressoren
    p['ko.pm(2)'] = 35  # [gC] condensatietemperatuur waarbij max. nominaal thermisch koelvermogen is gespecificeerd
    p['ko.pm(3)'] = -15 # [gC] verdampingstemperatuur waarbij max. nominaal thermisch koelvermogen is gespecificeerd
    p['ko.pm(4)'] = 0.7 # [(kW th)/(kW el)] fixed compressor rendement
     
    # IJSBANK MODELPARAMETERS
    p['ib.pm(1)'] = 50             #9173 # [kg] mass of water/ice in storage tank (9173 kg after initial filling by Jouke, Sept. '17)
    p['ib.pm(2)'] = 1.23           # [m] effective ice bank radius: true radius minus a few percent for space occupied by glycol tubing, based on jouke's measurements
    p['ib.pm(3)'] = 2              # [m] physical ice bank height : based on jouke's measurements
    p['ib.pm(4)'] = 2*np.pi*p['ib.pm(2)']*p['ib.pm(3)']+2*np.pi*p['ib.pm(2)']**2 # [m2] surface area of storage tank
    p['ib.pm(5)'] = 0.1            #0.2 # [W/m2.gC] specific heat transmission coeff. of storage tank wall
    p['ib.pm(6)'] = 50*6e1         # [W/gC] UA_ib_charge_melted, heat transfer coefficient from ice/water to glycol in case all content is water
    p['ib.pm(7)'] = 50*3e1         # [W/gC] UA_ib_charge_iced, heat transfer coefficient from ice/water to glycol in case all content is ice
    p['ib.pm(8)'] = 15*4e1         # [W/gC] UA_ib_discharge_melted, heat transfer coefficient from ice/water to subcoolerglycol in case all content is water
    p['ib.pm(9)'] = 15*2e1         # [W/gC] UA_ib_discharge_iced, heat transfer coefficient from ice/water to sbucoolerglycol in case all content is ice
    
    # KOUDEOPWEKKING MODELPARAMETERS
    p['brine.pm(1)'] = 3.5     # [kJ/kg.gC] specific heat of used brine
    p['brine.pm(2)'] = 1000    # [kg] amount of brine in storage vessel + piping
    
    # KLIMAATCOMPUTERINSTELLINGEN
    p['KlimaatCompuSettings.persdruk.RegelenOp']                      = 'buitentemperatuur'    #{persdruk setpoint} OR {buitentemperatuur}
    p['KlimaatCompuSettings.persdruk.GewenstePersdruk']               = 32                     #ingestelde condensatietemperatuur wanneer RegelenOp='persdruk setpoint'
    p['KlimaatCompuSettings.persdruk.DifferentieOpBuitentemperatuur'] = 8                      #wanneer RegelenOp='buitentemperatuur' ingestelde condensatietemperatuur wanneeer RegelenOp='persdruk setpoint'
    p['KlimaatCompuSettings.persdruk.MaximaleOpschakeltmeperatuur']   = 36
    p['KlimaatCompuSettings.persdruk.MinimaleAfschakeltemperatuur']   = 22
    p['KlimaatCompuSettings.zuigdruk.regelfunctie']                   = 'vast'                 #{vast} OR {op cel}
    p['KlimaatCompuSettings.zuigdruk.Regelwaarde']                    = -7                     #streefwaarde zuigdruk = regelwaarde + verlies leidingwerk
    p['KlimaatCompuSettings.zuigdruk.VerliesLeidingwerk']             = -3                     #streefwaarde zuigdruk = regelwaarde + verlies leidingwerk
     
    # DOELFUNCTIETUNINGPARAMETERS
    p['pd.AC_max']  = 55        # [kW] max. power uptake from AC grid
    p['pj(1)']      = 10        # [-], L(x,u) = ... + pj(1)*xL_el + ...
    p['pj(2)']      = 1         # [€/h/kW^pj(5)], L(x,u) = ... + pj(2)*xL_elACmax + ...
    p['pj(3)']      = 0.5       # [€/h/kW^pj(6)], L(x,u) = ... + pj(3)*xL_Qib_dsicharge_max
    p['pj(4)']      = 1         # [€/h/kW^pj(7)], L(x,u) = ... + pj(4)*xL_Tbrine
    p['pj(5)']      = 2         # [-], macht in Lel_AC_max, L(x,u) = ... + (abs(Pel_AC)>(pd(1)-1))*(abs(Pel_AC)-(pd(1)-1))^pj(4) + ...
    p['pj(6)']      = 2         # [-], macht in L_Qib_discharge_max
    p['pj(7)']      = 2         # [-], macht in L_Tbrine
    
    
    p['Qcellen_min24Uur'] = 10  #[kW], cooling demand by chambers observed 24 hours ago
          
    return p


def LoadDisturbances(ops):
   
    d              = np.zeros((4,ops['N']))
    
    #buitentemperatuur
    d[0,:]         = 10                                                    # disturbance input, assuming constant Tamb = 10 gC   
    #door PV opgewekte elektriciteit
    d[1,:]         = 0                                                     # disturbance input Pel_PV (kW)
    
    #pdb.set_trace()
    
    d[1,int(12/ops['h']):int(13/ops['h'])]   = 25                        # disturbance input Pel_PV (kW)  
     
    #laagtarief   = 0.01                                                                                                # euro/kWh,   bron: https://www.eneco.nl/energieprijzen/actuele-energieprijzen/ (26-4-2016)
    #hoogtarief   = 0.07                                                                                                # euro/kWh,   bron: https://www.eneco.nl/energieprijzen/actuele-energieprijzen/ (26-4-2016)    
    #e_afname     = laagtarief*np.ones((1,size(ops['t']))) 
    #ih           = find( (mod(treal,7*24)>24) & (mod(treal,7*24)<(6*24)) & (mod(treal,24)>=7) & (mod(treal,24)<23))    # hoogtarief geldt op werkdagen van 7 tot 23 uur
    #e_afname[ih] = hoogtarief                                                                                          # prijs die gebruiker betaalt voor afname van elektrische energie van het net
    #e_levering   = e_afname/2                                                                                          # prijs die gebruiker ontvangt voor teruglevering van elektrische energie aan het net    
    #d[2:4,:]     = [e_afname,e_levering]                                                                               #[prijs el. afname, prijs el. levering] (€/kWh)
    
    d[2:4,:]           = 0.04                                                 #[prijs el. afname, prijs el. levering] (€/kWh)
    d[2,int(10.1/ops['h']):int(11/ops['h'])] = 0.40         
    
    return d


def fd(x,u,d,p,h):
    
#       u   =   control input (dim. 1 * (AantalKoelcellen + 3))
#             u[0] = [0, 1] [-] u_comp, fraction of compressor capacity (scalar)
#             u[1] = [0, 1] [-] u_charge, fraction of icebank charging capacity used
#             u[2] = [0, 1] [-] u_discharge, fraction of icebank discharging capacity used

#       x   =   state [xm, xL]
#             xm[0] = enthalpy content of ice buffer [kJ]
#             xm[1] = brine temperature [gC]
#             xL[0] = costs of AC electricity consumption [€]
#             xL[1] = penalty for AC power uptake above a user-defined set limit Pel_AC_max (typically set at 55 kW) [kW^pj(5).h]
#             xL[2] = penalty on Qib_discharge > 0.1*Pth_comp [kW^pj(6).h]
#             xL[3] = penalty op Tbrine > zuigdrukregelwaarde + 4 gC [kW^pj(7)]

#       y   =   measurement [x, xh]
#             xh[0] = Ticebank, temperatur of water/ice in ice bank temperature [gC]
#             xh[1] = h_icebank, height of liquid level in ice buffer [m]
#             xh[2] = m_ice/ModelParameters.ib.pm(1), ice fraction [-]
#             xh[3] = Pth_comp, thermal compressor capacity [kW]
#             xh[4] = Q_ib_charge, heat flow rate from water/ice to  brine in ice bank [kW]
#             xh[5] = Q_ib_discharge, heat flow rate from liquid refrigerant to  water/ice in ice bank [kW]
#             xh[6] = Pel_AC, electric power uptake from AC power grid [kW]

#       d   =   disturbance input
#             d[0] = ambient temperature [gC]
#             d[1] = elektrisch vermogen opgewekt door PV panelen [kW]
#             d[2] = prijs voor afname van electriciteit uit het net [€/kWh]
#             d[3] = prijs voor levering van electriciteit aan het net [€/kWh]
   
    k1  = f(x,u,d,p)
    k2  = f(x + h/2 * k1,u,d,p)
    k3  = f(x + h/2 * k2,u,d,p)
    k4  = f(x + h * k3,u,d,p)
    x   = x + h/6*(k1 + 2*k2 + 2*k3 + k4)

    dx  = x
    
    return dx

def f(x,u,d,p):
        
    
    Tamb       = d[0]   # [gC] ambient temperature
    Pel_PV     = d[1]   # [kW] elektrisch vermogen opgewekt door PV panelen
    e_afname   = d[2]   # [€/kWh] prijs voor afname van electriciteit uit het net
    e_levering = d[3]   # [€/kWh] prijs voor levering van electriciteit aan het net
    
    
    # KOUDEOPWEKKING: KOELVERMOGEN + ELEKTRISCH VERMOGEN
    [Pth_comp,Pel_comp] = fcn_koudeopwekking(u[0],x[1],Tamb,p)
    
    # IJSBANK: HEAT BALANCE + ICE BANK BRINE OUTLET TEMPERATURE
    if p['KlimaatCompuSettings.persdruk.RegelenOp']=='persdruk setpoint':
        Tc = p['KlimaatCompuSettings.persdruk.GewenstePersdruk']                             # [gC] condensatietemperatuur
    else:
        Tc = max(p['KlimaatCompuSettings.persdruk.MinimaleAfschakeltemperatuur'],  
               min(Tamb+p['KlimaatCompuSettings.persdruk.DifferentieOpBuitentemperatuur'],  
                   p['KlimaatCompuSettings.persdruk.MaximaleOpschakeltmeperatuur']))         # [gC] condensatietemperatuur

    [dHicebankdt, Ticebank, m_ice, h_icebank, Q_ib_charge, Q_ib_discharge, UA_ib_discharge] = fcn_icebank(u[1:3],x[0],Tamb,p, x[1], Tc)

    # BRINE: HEAT BALANCE
    c_brine     = p['brine.pm(1)']              # [kJ/kg.gC] specific heat of used brine
    m_brine     = p['brine.pm(2)']              # [kg] amount of brine in storage vessel + piping
    Qcellen     = p['Qcellen_min24Uur']         # average cooling demand by chambers 24 hours ago
    dTbrinedt   = (3600/(m_brine*c_brine))*(Qcellen+Q_ib_charge-Pth_comp-Q_ib_discharge+0.3*(20-x[1])) # [gC/h] heat balance brine storage tank

    # ELEKTRISCHE VERMOGENS: Pel_AC + Pel_PV = Pel_comp +Pel_others
    Pel_others  = u[1]+u[2]                    # [kW], elektrisch vermogensopname door fans, pompen, licht, wasmachine, hogedrukspuit, etc.
    Pel_AC      = Pel_comp+Pel_others-Pel_PV   # [kW], elektrisch vermogen opgenomen uit (of, indien negatief, geleverd aan) het net
        
    # OBJECTIVE FUNCTIONAL STATE EQUATIONS
    dxL_eldt                = e_afname*Pel_AC*(Pel_AC>0)+e_levering*Pel_AC*(Pel_AC<0)                  # €/h
    
    AC_max                  = p['pd.AC_max']                                                                # [kW] max. power uptake from AC grid
    
    #pdb.set_trace()    
    dxL_elACmaxdt           = (abs(Pel_AC)>(AC_max-1))*(abs(Pel_AC)-(AC_max-1))**p['pj(5)']                  # kW^pj(5)
    dxL_Qib_discharge_maxdt = ((Q_ib_discharge-0.1*Pth_comp)>-1)*(Q_ib_discharge-0.1*Pth_comp+1)**p['pj(6)'] # kW^pj(6)
    dxL_Tbrinedt            = ((x[1]-p['KlimaatCompuSettings.zuigdruk.Regelwaarde']-4)>-1)*((x[1]-p['KlimaatCompuSettings.zuigdruk.Regelwaarde']-4)+1)**p['pj(7)'] # kW^pj(7)
    
   
    ki =  np.array([dHicebankdt,dTbrinedt,dxL_eldt,dxL_elACmaxdt,dxL_Qib_discharge_maxdt,dxL_Tbrinedt])
       

    return ki
    
def g(x,u,d,p):
    
    
    Tamb       = d[0]   # [gC] ambient temperature
    Pel_PV     = d[1]   # [kW] elektrisch vermogen opgewekt door PV panelen

    # KOUDEOPWEKKING: KOELVERMOGEN + ELEKTRISCH VERMOGEN
    [Pth_comp,Pel_comp] = fcn_koudeopwekking(u[0],x[1],Tamb,p)
    
    # IJSBANK: HEAT BALANCE + ICE BANK BRINE OUTLET TEMPERATURE
    if p['KlimaatCompuSettings.persdruk.RegelenOp']=='persdruk setpoint':
        Tc = p['KlimaatCompuSettings.persdruk.GewenstePersdruk']                             # [gC] condensatietemperatuur
    else:
        Tc = max(p['KlimaatCompuSettings.persdruk.MinimaleAfschakeltemperatuur'],  
               min(Tamb+p['KlimaatCompuSettings.persdruk.DifferentieOpBuitentemperatuur'],  
                   p['KlimaatCompuSettings.persdruk.MaximaleOpschakeltmeperatuur']))         # [gC] condensatietemperatuur
       

    [dHicebankdt, Ticebank, m_ice, h_icebank, Q_ib_charge, Q_ib_discharge, UA_ib_discharge] = fcn_icebank(u[1:3],x[0],Tamb,p, x[1], Tc)
    
    # ELEKTRISCHE VERMOGENS: Pel_AC + Pel_PV = Pel_comp +Pel_others
    Pel_others  = u[1]+u[2]                    # [kW], elektrisch vermogensopname door fans, pompen, licht, wasmachine, hogedrukspuit, etc.
    Pel_AC      = Pel_comp+Pel_others-Pel_PV   # [kW], elektrisch vermogen opgenomen uit (of, indien negatief, geleverd aan) het net

    xh = [Ticebank, h_icebank, m_ice/p['ib.pm(1)'], Pth_comp, Q_ib_charge, Q_ib_discharge, Pel_AC]
    
    y  = np.concatenate((x,xh), axis=0)
         
     
    return y
    
        
    return np.concatenate((x,xh))
    
def controller(x,u,d,p,ops):
    
        # CONTROL INPUTS
        u_comp      = 0.24          # [-] fraction of cooling capacity used
        u_charge    = 0             # [-] fraction of icebank charging capacity used (0 = no charging, 1 = max. charging)
        u_discharge = 0             # [-] fraction of icebank discharging capacity used (0 = no discharging, 1 = max. discharging)
        umin        = [0, 0, 0]     # [u_comp, u_charge, u_discharge]
        umax        = [1, 1, 1] 
        
        u_0 = min(umax[0],max(umin[0],u_comp))      # [-] compressorsturing
        u_1 = min(umax[1],max(umin[1],u_charge))    # [-] fraction of icebank charging capacity used (0 = no charging, 1 = max. charging)
        u_2 = min(umax[2],max(umin[2],u_discharge)) # [-] fraction of icebank discharging capacity used (0 = no discharging, 1 = max. discharging)

        u_  = np.array([u_0,u_1,u_2])

        return u_
    
def fcn_koudeopwekking(u_comp,Tbrine,Tamb,p):
    #   DESCRIPTION:
    #       submodel van koudeopwekking
    #   INPUTS:
    #         u_comp = used fraction of compressor capacity [-]
    #         Tamb = external temperature of ice bank [gC]
    #         ModelParameters = structure with variable no. of model parameter vectors
    #
    #   OUTPUTS:
    #         Pth_comp = thermal refrigeration capacity [kW]
    #         Pel_comp = electric power uptake by refrigeration plant (i.e. compressors) [kW el.]
    
    
    Pth_comp_max_nom = p['ko.pm(1)'] # [kW] max. nominaal thermisch koelvermogen van compressoren
    Tc_nom           = p['ko.pm(2)'] # [gC] condensatietemperatuur waarbij max. nominaal thermisch koelvermogen is gespecificeerd
    To_nom           = p['ko.pm(3)'] # [gC] verdampingstemperatuur waarbij max. nominaal thermisch koelvermogen is gespecificeerd
    eta_comp         = p['ko.pm(4)'] # [(kW th)/(kW el)] fixed compressor rendement
    
    if p['KlimaatCompuSettings.persdruk.RegelenOp'] == 'persdruk setpoint':
        Tc = p['KlimaatCompuSettings.persdruk.GewenstePersdruk'] # [gC] condensatietemperatuur
    else:
        Tc = max(p['KlimaatCompuSettings.persdruk.MinimaleAfschakeltemperatuur'],
               min(Tamb+p['KlimaatCompuSettings.persdruk.DifferentieOpBuitentemperatuur'], 
                   p['KlimaatCompuSettings.persdruk.MaximaleOpschakeltmeperatuur']))  # [gC] condensatietemperatuur

    if p['KlimaatCompuSettings.zuigdruk.regelfunctie'] == 'vast':
        To = Tbrine-3   # [gC] actual evaporation temperature
    else:
        To = 100        # [gC], in werkelijkheid: de koudste (desired temperature - suction temperature diff.) - VerliesLeidingwerk
    
    Pth_comp_max = max(0,Pth_comp_max_nom*(1+0.02*((Tc_nom-To_nom)-(Tc-To)))*
                       max(0,min(1,(Tc-Tamb)/p['KlimaatCompuSettings.persdruk.DifferentieOpBuitentemperatuur']))) # [kW]
    Pth_comp     = u_comp*Pth_comp_max              # [kW th]
    COP          = eta_comp*(To+273.15)/(Tc-To)     # [(kW th.)/(kW el.)]
    Pel_comp     = Pth_comp/COP                     # [kW el.]
    
    return [Pth_comp,Pel_comp]

def fcn_icebank(u,x,Tamb,p, Tbrine_in,Tc):
    #   DESCRIPTION:
    #       ice bank submodel
    
    #   INPUTS:
    #       u(1) = sturing laden voor ijsbank (1= max. laden, 0 = niet laden) [-]
    #       u(2) =  sturing ONTladen voor ijsbank (1= max. ONTladen, 0 = niet ONTladen) [-]
    #       x(1) = enthalpy content of ice buffer [kJ]
    #       Tamb = external temperature of ice bank [gC]
    #       ModelParameters =  structure with variable no. of model parameter vectors
    #       fc = structure met in het model gebruikte fysische constantes van water en ijs.
    #       Tbrine_in = glycoltemperatuur aan glycolingang van ijsbank
    #       Tc = condensation temperature 
    
    #   OUTPUTS:
    #       dHicebankdt = rate of change of enthalpy content of ice buffer [kJ/h]
    #       Ticebank = temperature of water/ice in ice bank [gC]
    #       Mice = mass of ice [kg]
    #       h_icebank = height of water leven in ice bank [kg]
    
    # MODEL PARAMETERS
    A_st                    = p['ib.pm(4)'] # [m2] external surface area of storage tank
    U_st                    = p['ib.pm(5)'] # [W/m2.gC] specific heat transmission coeff. of storage tank wall
    UA_ib_charge_melted     = p['ib.pm(6)'] # [W/gC] UA_ib_charge_melted, heat transfer coefficient from ice/water to glycol in case all content is water
    UA_ib_charge_iced       = p['ib.pm(7)'] # [W/gC] UA_ib_charge_iced, heat transfer coefficient from ice/water to glycol in case all content is ice
    UA_ib_discharge_melted  = p['ib.pm(8)'] # [W/gC] UA_ib_discharge_melted, heat transfer coefficient from ice/water to subcoolerglycol in case all content is water
    UA_ib_discharge_iced    = p['ib.pm(9)'] # [W/gC] UA_ib_discharge_iced, heat transfer coefficient from ice/water to sbucoolerglycol in case all content is ice
    Mib                     = p['ib.pm(1)'] # [kg] mass of water/ice in storage tank (9173 kg after initial filling by Jouke, Sept. '17)
    
    # ASSESS Ticebank AND m_ice
    Hicebank    = x # [kJ] enthalpy content of storage medium
    
    [Ticebank, Mice,h_icebank] = fcn_mTh_icebank(Hicebank,p)       # [gC] temperature of water and ice in icebank, [kg] mass of ice in icebank
    
    # CALCULATE ACTUAL HEAT FLOWS
    UA_st   = U_st*A_st                     # [W/gC] heat transmission coeff. of storage tank wall
    Q_amb   = UA_st *(Tamb-Ticebank)/1000   # [kW] heat flow from ambient into storage tank
    
    UA_ib_charge = (Mib-Mice)/Mib*UA_ib_charge_melted+Mice/Mib*UA_ib_charge_iced
    Q_ib_charge  = UA_ib_charge*u[0]*(Ticebank-Tbrine_in)/1000                   #[kW] heat flow from storage medium to brine
    
    UA_ib_discharge  = (Mib-Mice)/Mib*UA_ib_discharge_melted+Mice/Mib*UA_ib_discharge_iced
    Q_ib_discharge   = UA_ib_discharge*u[1]*(Tc-Ticebank)/1000                              #[kW] heat flow from hot gas refrigerant to storage medium
    dHicebankdt      = 3600*(Q_amb - Q_ib_charge+Q_ib_discharge)                            # [kJ/h] rate of change of enthalpy content of storage medium
    
    return [dHicebankdt, Ticebank, Mice, h_icebank, Q_ib_charge, Q_ib_discharge, UA_ib_discharge]




def fcn_mTh_icebank(Hicebank,p):
    #   DESCRIPTION:
    #       calculate ice mass and icebank temperature
    
    #   INPUTS:
    #         Hicebank = enthalpy content of ice buffer [kJ]
    #         ModelParameters = structure with variable no. of model parameter vectors
    #         fc = structure met in het model gebruikte fysische constantes van water en ijs.
    
    #   OUTPUTS:
    #         Mice = amount of ice in ice bank [kg]
    #         Ticebank = temperatur of water/ice in ice bank temperature [gC]
    #          h_icebank = height of water level in ice bank [m]
       
    
    # MODEL PARAMETERS
    m_ib = p['ib.pm(1)']                # [kg] mass of water/ice in storage tank
    r_ib = p['ib.pm(2)']                # [m] inner radius of storage tank
    
    if Hicebank<-m_ib*p['fc.h_melt']:   # storage medium is all ice
        
        m_ice    = m_ib                                                  # [kg] mass of ice in icebank
        Ticebank = (Hicebank+m_ice*p['fc.h_melt'])/(m_ice*p['fc.c_ice']) # [gC] temperature of water and ice in icebank
    
    elif Hicebank <= 0:                 # storage medium is mixture of ice and water
        
        m_ice = -Hicebank/p['fc.h_melt']            # [kg] mass of ice in icebank
        Ticebank = 0                                #[gC] temperature of water and ice in icebank
    else:                                           # storage medium is all water
        m_ice    = 0                                # [kg] mass of ice in icebank
        Ticebank = Hicebank/(m_ib*p['fc.c_water'])  # [gC] temperature of water and ice in icebank

    m_water    = m_ib-m_ice
    Vib        = m_water/p['fc.rho_water'] + m_ice/p['fc.rho_ice']  # [m3] ice bank volume
    h_icebank  = Vib/(np.pi*r_ib)                                   # [m] height of water level in ice bank
    
    return [Ticebank, m_ice,h_icebank]