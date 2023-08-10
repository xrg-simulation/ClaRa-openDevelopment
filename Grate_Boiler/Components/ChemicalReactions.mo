within Grate_Boiler.Components;
package ChemicalReactions
  extends ClaRa.Basics.Icons.PackageIcons.Components80;

  model Pyrolyis_Base "Base Model for Pyrolysis"
    extends ClaRa.Basics.Icons.ChemicalReactions;

  end Pyrolyis_Base;

  partial model PartialReactionZonePyro
    "Reaction Zone for solid fuels with volatiles to consider pyrolysis"
    extends Pyrolyis_Base;

  outer ClaRa.SimCenter simCenter;

  inner parameter TILMedia.GasTypes.BaseGas flueGas=simCenter.flueGasModel;

  parameter ClaRa.Basics.Media.FuelTypes.BaseFuel fuelModel=simCenter.fuelModel1;

  public
    ClaRa.Basics.Units.MassFlowRate m_flow_volatiles;
    ClaRa.Basics.Units.EnthalpyFlowRate H_pyrolysis_fuel;
    ClaRa.Basics.Units.EnthalpyFlowRate H_pyrolysis_fluegas;

  //________/ component properties \________
    ClaRa.Basics.Units.Mass m_fuel_vol;
    ClaRa.Basics.Units.Mass m_fuel_char;

  //__________/ masses of elements in fluegas created in pyrolysis \________
     ClaRa.Basics.Units.Mass m_vol_C;
     ClaRa.Basics.Units.Mass m_vol_H;
     ClaRa.Basics.Units.Mass m_vol_O;
  //
    ClaRa.Basics.Units.Mass m_fuel_C;
    ClaRa.Basics.Units.Mass m_fuel_H;
    ClaRa.Basics.Units.Mass m_fuel_O;
    ClaRa.Basics.Units.Mass m_fuel_N;

    ClaRa.Basics.Units.Mass m_vol_CO;
    ClaRa.Basics.Units.Mass m_vol_CO2;
    ClaRa.Basics.Units.Mass m_vol_H2O;
    ClaRa.Basics.Units.Mass m_vol_CH4;
    ClaRa.Basics.Units.Mass m_vol_H2;
    ClaRa.Basics.Units.Mass m_vol_N2;

  //
     ClaRa.Basics.Units.MassFlowRate m_flow_vol_C;
     ClaRa.Basics.Units.MassFlowRate m_flow_vol_H;
     ClaRa.Basics.Units.MassFlowRate m_flow_vol_O;
     ClaRa.Basics.Units.MassFlowRate m_flow_vol_N;

  //__________/ masses of fluegas components created in fuel \________
    ClaRa.Basics.Units.MassFraction prod_comp[flueGas.nc-1];

  //__________/ mass fractions in fluegas created in pyrolysis \________
  // protected
    ClaRa.Basics.Units.MassFraction xi_CO_fuel;
    ClaRa.Basics.Units.MassFraction xi_CH4_fuel;
    ClaRa.Basics.Units.MassFraction xi_H2_fuel;
    ClaRa.Basics.Units.MassFraction xi_H2O_fuel;
    ClaRa.Basics.Units.MassFraction xi_CO2_fuel;
    ClaRa.Basics.Units.MassFraction xi_N2_fuel;

    ClaRa.Basics.Units.EnthalpyMassSpecific h_fluegas "Specific enthalpy of fluegas after pyrolysis";
    ClaRa.Basics.Units.EnthalpyMassSpecific h_C "Specific enthalpy of the fuel carbon";
    ClaRa.Basics.Units.EnthalpyMassSpecific h_O "Specific enthalpy of the fuel Oxygen";
    ClaRa.Basics.Units.EnthalpyMassSpecific h_H  "Specific enthalpy of the fuel Hydrogen";

  //  ClaRa.Basics.Units.MassFraction xi_pyro "Composition of fluegas from pyrolysis ={N2,H2,CO,O2,H2O,CO2,CH4}";

  end PartialReactionZonePyro;

  model SolidFuelReactionZone "pyrolysis process - fixed ratios"

  extends Grate_Boiler.Components.ChemicalReactions.PartialReactionZonePyro;
    outer Bedsegment.Records.IComBedSegment iCom;

    //composition of pyrolysis gas will be constant

  // parameter Real xi_CO2_vol = 0.143 "Carbondioxide that is created from C and O of fuel" annotation (Dialog(group="Pyrolysis"));
  //
  // parameter Real xi_H2O_vol = 0.131  "Water that is created from H2 and O of fuel"   annotation (Dialog(group="Pyrolysis"));
  //
  // parameter Real xi_char_fuel= 0.45 "Fraction of carbon that is fixed as solid char";

  // additional pyrolysis products to be added (CH4, H2, CO)

  //parameter replacable machen und welche aus Rückert (Euler/Euler) zur Verfügung stellen
  parameter Real A_0(unit="1/s") = 650  "Coefficient of impact";
  parameter Real E(unit="J/mol") = 54044  "Reaction activating energy";
            Real k(unit="1/s") "Pyrolysis rate constant";
  parameter Real xi_CH4_conversion= 1/3 "Fraction of H2 that is reacting with C to Methane";
  parameter Real xi_CO_conversion= 1/3 "Fraction of H2 that is reacting with C to Methane";
  parameter Real xi_CO2_conversion= 1/3 "Fraction of H2 that is reacting with C to Methane";
  final parameter Real xi_H2O_conversion = 1-xi_CO2_conversion-xi_CO_conversion;

      /*
    In the literature source, there is also conversion to tar and the degradation
    of this as a additional step considered, which is left out for simplification
    */

  protected
    TILMedia.Gas ptr(gasType=flueGas, computeTransportProperties=true)
                                    annotation (Placement(transformation(extent={{30,10},{50,30}})));
  equation
    //total mass input
    m_fuel_C = iCom.m_fuel_waf*iCom.xi_vol[1];
    m_fuel_H = iCom.m_fuel_waf*iCom.xi_vol[2];
    m_fuel_O = iCom.m_fuel_waf*iCom.xi_vol[3];
    m_fuel_N = iCom.m_fuel_waf*iCom.xi_vol[4];
    //distribute educts to pyrolysis products
    m_vol_CO  = (xi_CO_conversion*m_fuel_O)*28/16;
    m_vol_CO2 = (xi_CO2_conversion*m_fuel_O)*44/32;
    m_vol_H2O = (xi_H2O_conversion*m_fuel_O)*18/16;

    m_vol_CH4 = (xi_CH4_conversion*m_fuel_H)*16/4;
    m_vol_H2  = (m_fuel_H - m_vol_CH4*4/16 - m_vol_H2O*2/18);
    m_vol_N2= m_fuel_N;

    //calculate total mass of volatiles
    m_fuel_vol = m_vol_CO + m_vol_CO2 + m_vol_H2O + m_vol_CH4 + m_vol_H2 + m_vol_N2;

    xi_CH4_fuel = m_vol_CH4/(max(m_fuel_vol,Modelica.Constants.eps));
    xi_CO_fuel  = m_vol_CO/(max(m_fuel_vol,Modelica.Constants.eps));
    xi_CO2_fuel = m_vol_CO2/(max(m_fuel_vol,Modelica.Constants.eps));
    xi_H2O_fuel = m_vol_H2O/(max(m_fuel_vol,Modelica.Constants.eps));
    xi_H2_fuel  = m_vol_H2/(max(m_fuel_vol,Modelica.Constants.eps));
    xi_N2_fuel  = m_vol_N2/(max(m_fuel_vol,Modelica.Constants.eps));

  //calculate pyrolysis rate
    k = A_0*exp(-E/(Modelica.Constants.R*iCom.T_fuel));

  //calculate masses
    m_fuel_char = m_fuel_C - m_vol_CO*12/28 - m_vol_CO2*12/44 - m_vol_CH4*12/16;

    m_vol_C = m_fuel_C - m_fuel_char;
    m_vol_H = m_fuel_H;
    m_vol_O = m_fuel_O;

  //calculate mass flows
    m_flow_volatiles = k * m_fuel_vol;
    m_flow_vol_C = m_vol_C * k;
    m_flow_vol_H = m_vol_H * k;
    m_flow_vol_O = m_vol_O * k;
    m_flow_vol_N = m_vol_N2 * k;

  //calculate gas composition of gas pyrolysis gas
     for i in 1:(flueGas.nc-1) loop
             if i==1 then prod_comp[1] = max(Modelica.Constants.eps,xi_CH4_fuel);
        else if i==2 then prod_comp[2] = max(Modelica.Constants.eps,xi_CO_fuel);
        else if i==3 then prod_comp[3] = max(Modelica.Constants.eps,xi_CO2_fuel);
        else if i==3 then prod_comp[5] = max(Modelica.Constants.eps,xi_N2_fuel);
        else if i==8 then prod_comp[8] = max(Modelica.Constants.eps,xi_H2O_fuel);
        else if i==9 then prod_comp[9] = max(Modelica.Constants.eps,xi_H2_fuel);
      else
          prod_comp[i] = Modelica.Constants.eps;
             end if; end if; end if; end if; end if; end if;
     end for;

  //enthalpy of the gaseous volatiles in the fuel
  //enthalpy of the gaseous volatiles in the fuel
    //h_fluegas = TILMedia.GasFunctions.specificEnthalpy_pTxi(flueGas, iCom.p_fluegas, iCom.T_fuel, prod_comp);
    h_fluegas = TILMedia.GasObjectFunctions.specificEnthalpy_pTxi(iCom.p_fluegas, iCom.T_fuel, prod_comp, ptr.gasPointer);
    h_C  = ClaRa.Basics.Media.FuelFunctions.enthalpy_pTxi(iCom.p_fluegas,iCom.T_fuel,{1,0,0,0,0,0},fuelModel);
    h_O  = ClaRa.Basics.Media.FuelFunctions.enthalpy_pTxi(iCom.p_fluegas,iCom.T_fuel,{0,0,1,0,0,0},fuelModel);
    h_H  = ClaRa.Basics.Media.FuelFunctions.enthalpy_pTxi(iCom.p_fluegas,iCom.T_fuel,{0,1,0,0,0,0},fuelModel);

  //_________/ Enthalpy of flue gas for global energy balances \______

    H_pyrolysis_fuel    = m_flow_volatiles*(h_C + h_O + h_H);//(h_C*m_vol_C/m_fuel_vol + h_O*m_vol_O/m_fuel_vol + h_H*m_vol_H/m_fuel_vol);//
    H_pyrolysis_fluegas = m_flow_volatiles*h_fluegas;

  end SolidFuelReactionZone;
end ChemicalReactions;
