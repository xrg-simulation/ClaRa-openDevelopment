within Grate_Boiler.Components;
package Gasphase "Combustion above fuel bed"
  extends ClaRa.Basics.Icons.PackageIcons.Components80;

  package Heat_Transfer "Models for radiation"
    extends ClaRa.Basics.Icons.PackageIcons.Basics50;

    model Radiation_gas2Wall_L2_Grate "All Geo || L2 || Radiation Between Gas and Wall"
      //___________________________________________________________________________//
      // Component of the ClaRa library, version: 1.4.1                            //
      //                                                                           //
      // Licensed by the DYNCAP/DYNSTART research team under Modelica License 2.   //
      // Copyright  2013-2019, DYNCAP/DYNSTART research team.                      //
      //___________________________________________________________________________//
      // DYNCAP and DYNSTART are research projects supported by the German Federal //
      // Ministry of Economic Affairs and Energy (FKZ 03ET2009/FKZ 03ET7060).      //
      // The research team consists of the following project partners:             //
      // Institute of Energy Systems (Hamburg University of Technology),           //
      // Institute of Thermo-Fluid Dynamics (Hamburg University of Technology),    //
      // TLK-Thermo GmbH (Braunschweig, Germany),                                  //
      // XRG Simulation GmbH (Hamburg, Germany).                                   //
      //___________________________________________________________________________//

      extends Grate_Boiler.Components.Gasphase.Heat_Transfer.Radiation_Base_L2_Grate;

      input Real CF_fouling=0.8 "Scaling factor accounting for the fouling of the wall"
                                                                                       annotation (Dialog(group="Heat Transfer"));
      parameter Real emissivity_wall=0.8 "Emissivity of the wall";
      parameter Real emissivity_flame=0.9 "Emissivity of the flame";
      parameter Real absorbance_flame=0.9 "Absorbance of the flame";
      parameter Integer heatSurfaceAlloc=1 "To be considered heat transfer area" annotation (dialog(enable=false, tab="Expert Setting"), choices(
          choice=1 "Lateral surface",
          choice=2 "Inner heat transfer surface",
          choice=3 "Selection to be extended"));

      outer ClaRa.Basics.ControlVolumes.Fundamentals.Geometry.HollowBlock geo;

      parameter String temperatureDifference="Outlet" "Temperature Difference" annotation (Dialog(group="Heat Transfer"), choices(
          choice="Arithmetic mean",
          choice="Inlet",
          choice="Outlet",
          choice = "Bulk"));

      ClaRa.Basics.Units.Temperature Delta_T_mean "Mean temperature";

    equation

      if temperatureDifference == "Arithmetic mean" then
        Delta_T_mean = (iCom.T_in + iCom.T_out)/2;
      elseif temperatureDifference == "Inlet" then
        Delta_T_mean = iCom.T_in;
      elseif temperatureDifference == "Outlet" then
        Delta_T_mean = iCom.T_out;
      elseif temperatureDifference == "Bulk" then
        Delta_T_mean = iCom.T_gas;
      else
        Delta_T_mean = -1;
        assert(true, "Unknown temperature difference option in HT model");
      end if;

      //According to VDI Waermeatlas for a wall surrounding a gas volume chapter Kc5.
      heat.Q_flow = geo.A_heat_CF[heatSurfaceAlloc]*CF_fouling*Modelica.Constants.sigma*emissivity_wall/(absorbance_flame + emissivity_wall - absorbance_flame*emissivity_wall)*(absorbance_flame*heat.T^4 - emissivity_flame*Delta_T_mean^4);

      annotation (Documentation(info="<html>
<p><b>Model description: </b>A simple correlation for radiant heat transfer between gas and wall inside furnaces</p>
<p><b>Contact:</b> Lasse Nielsen, TLK-Thermo GmbH</p>
<p><b>FEATURES</b> </p>
<p><ul>
<li>Emissivities of flue gas and wall are constant values</li>
<li>Heat transfer area is calculated from geometry</li>
</ul></p>
</html>"));
    end Radiation_gas2Wall_L2_Grate;

    model Radiation_gas2Gas_L2_Grate "All Geo || L2 || Radiation Between Gas Volumes"
      //___________________________________________________________________________//
      // Component of the ClaRa library, version: 1.4.1                            //
      //                                                                           //
      // Licensed by the DYNCAP/DYNSTART research team under Modelica License 2.   //
      // Copyright  2013-2019, DYNCAP/DYNSTART research team.                      //
      //___________________________________________________________________________//
      // DYNCAP and DYNSTART are research projects supported by the German Federal //
      // Ministry of Economic Affairs and Energy (FKZ 03ET2009/FKZ 03ET7060).      //
      // The research team consists of the following project partners:             //
      // Institute of Energy Systems (Hamburg University of Technology),           //
      // Institute of Thermo-Fluid Dynamics (Hamburg University of Technology),    //
      // TLK-Thermo GmbH (Braunschweig, Germany),                                  //
      // XRG Simulation GmbH (Hamburg, Germany).                                   //
      //___________________________________________________________________________//

      extends Grate_Boiler.Components.Gasphase.Heat_Transfer.Radiation_Base_L2_Grate;

      parameter Real emissivity_top=0.8 "Emissivity of the gas volume above";
      parameter Real emissivity_flame=0.9 "Emissivity of the flame";

      //ClaRa.Basics.Units.Area A_eff "Effective heat transfer area";

    equation
      //Nach VDI Waermeatlas fuer zwei sehr grosse parallele ebene Flaechen
      heat.Q_flow = iCom.A_front*Modelica.Constants.sigma/(1.0/emissivity_top + 1.0/emissivity_flame - 1.0)*(heat.T^4 - iCom.T_out^4);

      annotation (Documentation(info="<html>
<p><b>Model description: </b>A simple correlation for radiant heat transfer between gases inside furnaces</p>
<p><b>Contact:</b> Lasse Nielsen, TLK-Thermo GmbH</p>
<p><b>FEATURES</b> </p>
<p><ul>
<li>Emissivities of flue gases are constant values</li>
<li>Heat transfer area is calculated from geometry</li>
</ul></p>
</html>"));
    end Radiation_gas2Gas_L2_Grate;

    model Radiation_Base_L2_Grate "All Geo || L2 || Radiation Between Gas Volumes"
      //___________________________________________________________________________//
      // Component of the ClaRa library, version: 1.4.1                            //
      //                                                                           //
      // Licensed by the DYNCAP/DYNSTART research team under Modelica License 2.   //
      // Copyright  2013-2019, DYNCAP/DYNSTART research team.                      //
      //___________________________________________________________________________//
      // DYNCAP and DYNSTART are research projects supported by the German Federal //
      // Ministry of Economic Affairs and Energy (FKZ 03ET2009/FKZ 03ET7060).      //
      // The research team consists of the following project partners:             //
      // Institute of Energy Systems (Hamburg University of Technology),           //
      // Institute of Thermo-Fluid Dynamics (Hamburg University of Technology),    //
      // TLK-Thermo GmbH (Braunschweig, Germany),                                  //
      // XRG Simulation GmbH (Hamburg, Germany).                                   //
      //___________________________________________________________________________//

      extends ClaRa.Basics.ControlVolumes.Fundamentals.HeatTransport.Gas_HT.Radiation.HeatTransfer_L2;
      outer Grate_Boiler.Components.Bedsegment.Records.IComFlueGas iCom;
      extends ClaRa.Basics.Icons.Epsilon;

      //ClaRa.Basics.Units.Area A_eff "Effective heat transfer area";

    equation

      annotation (Documentation(info="<html>
<p><b>Model description: </b>A simple correlation for radiant heat transfer between gases inside furnaces</p>
<p><b>Contact:</b> Lasse Nielsen, TLK-Thermo GmbH</p>
<p><b>FEATURES</b> </p>
<p><ul>
<li>Emissivities of flue gases are constant values</li>
<li>Heat transfer area is calculated from geometry</li>
</ul></p>
</html>"));
    end Radiation_Base_L2_Grate;
  end Heat_Transfer;

  package Gascombustion "combustion of gaseous products"
    extends ClaRa.Basics.Icons.PackageIcons.Basics50;

    model GasCombustion_Base "base class"
      extends ClaRa.Basics.Icons.ChemicalReactions;
    equation

    end GasCombustion_Base;

    partial model PartialGasCombustion "model for combustion of gaseous components"
      extends GasCombustion_Base;

      outer ClaRa.SimCenter simCenter;
      inner parameter TILMedia.GasTypes.BaseGas flueGas=simCenter.flueGasModel;

      //Molar fractions
      Real n_CO;
      Real n_O;
      Real n_CH4;
      Real n_H;
      Real n_H2O;  // (unit="kmol/m3")

      //mass flows
      ClaRa.Basics.Units.MassFlowRate m_flow_CH4;
      ClaRa.Basics.Units.MassFlowRate m_flow_CO;
      ClaRa.Basics.Units.MassFlowRate m_flow_CO2;
      // ClaRa.Basics.Units.MassFlowRate m_flow_N2;
      ClaRa.Basics.Units.MassFlowRate m_flow_O2;
      // ClaRa.Basics.Units.MassFlowRate m_flow_NO;
      ClaRa.Basics.Units.MassFlowRate m_flow_H2O;
      ClaRa.Basics.Units.MassFlowRate m_flow_H2;
    //   ClaRa.Basics.Units.MassFlowRate m_flow_comb;

      //molar reaction enthalpies
      Real delta_h_CH4(unit="J/mol") = -74.8977e3;
      Real delta_h_CO( unit="J/mol") = -110.1177e3;
      Real delta_h_CO2(unit="J/mol") = -392.0652e3;
      Real delta_h_H2O(unit="J/mol") = -240.9363e3;
      //heat flows
      ClaRa.Basics.Units.EnthalpyFlowRate Q_reac;
      ClaRa.Basics.Units.EnthalpyFlowRate Q_reac_CO;
      ClaRa.Basics.Units.EnthalpyFlowRate Q_reac_CH4;
      ClaRa.Basics.Units.EnthalpyFlowRate Q_reac_H2;

    equation

    end PartialGasCombustion;

    model GasCombustion  "Model for combustion of volatiles"
      extends PartialGasCombustion;

      outer Grate_Boiler.Components.Bedsegment.Records.IComFlueGas iCom;
      parameter Real CF_T=1 "Correction Factor for Temperature Distribution";
      parameter Real CF_m=1 "Correction Factor for Mass Distribution";

      //__________/ Reaction rates \__________
      Real R_CO(  unit="mol/(m3.s)");
      Real R_CH4( unit="mol/(m3.s)");
      Real R_H2(  unit="mol/(m3.s)");

    equation
      //Mass fractions to concentrations
      n_CO = iCom.m_fluegas*iCom.xi_fluegas[2]/((ClaRa.Basics.Constants.M_C + ClaRa.Basics.Constants.M_O)*iCom.Volume);
      n_O  = iCom.m_fluegas*iCom.xi_fluegas[6]/((ClaRa.Basics.Constants.M_O)*iCom.Volume);
      n_CH4= iCom.m_fluegas*iCom.xi_fluegas[1]/((ClaRa.Basics.Constants.M_C + 4* ClaRa.Basics.Constants.M_H)*iCom.Volume);
      n_H  = iCom.m_fluegas*iCom.xi_fluegas[9]/((ClaRa.Basics.Constants.M_H)*iCom.Volume);
      n_H2O= iCom.m_fluegas*iCom.xi_fluegas[8]/(ClaRa.Basics.Constants.M_H2O*iCom.Volume);

    //reaction rates
      //(values taken from Karamarkovic "Modellierung einer Biomassefeuerung am Rost als Basis für die Entwicklung einer modellbasiertem, prädikativen Regelung"
      //reaction of CO and 1/2 O2 to CO2
      R_CO = 3.25e7 * exp(-15098/iCom.T_gas/CF_T)* max(n_CO,Modelica.Constants.eps) * (max(n_O,Modelica.Constants.eps))^0.5 * (max(n_H2O,Modelica.Constants.eps))^0.5 * CF_m;
      //reaction of CH4 and 3/2 O2 to CO and 2 H2O
      R_CH4 = 1.585e10 * exp(-24157/iCom.T_gas/CF_T)*(max(n_CH4,Modelica.Constants.eps))^0.7 * (max(n_O,Modelica.Constants.eps))^0.8 * CF_m;
      //reaction of 2 H2 and O2 to 2CO and 3 H2O
      R_H2 = 1.631e9 * (iCom.T_gas)^(-3/2) * exp(-3420/iCom.T_gas/CF_T)*(max(n_H,Modelica.Constants.eps))^1.5 * (max(n_O,Modelica.Constants.eps)) * CF_m;

    //massflows
      m_flow_CH4 =  -R_CH4*(ClaRa.Basics.Constants.M_C + 4*ClaRa.Basics.Constants.M_H)  * iCom.Volume;
      m_flow_CO  =  (R_CH4 - R_CO)*(ClaRa.Basics.Constants.M_C + ClaRa.Basics.Constants.M_O)   * iCom.Volume;
      m_flow_CO2 =  (R_CO)*(ClaRa.Basics.Constants.M_C + ClaRa.Basics.Constants.M_O2)   * iCom.Volume;
      m_flow_O2  =  (-3/2*R_CH4 - 1/2*R_H2 - 1/2*R_CO)*(ClaRa.Basics.Constants.M_O2)   * iCom.Volume;
      m_flow_H2O =  (2*R_CH4 + R_H2)*(ClaRa.Basics.Constants.M_H2O)  * iCom.Volume;
      m_flow_H2  =  (-R_H2)*(2* ClaRa.Basics.Constants.M_H)   * iCom.Volume;

    //   m_flow_comb= abs(m_flow_CH4) + abs(m_flow_CO) + abs(m_flow_CO2) + abs(m_flow_O2) + abs(m_flow_H2O) + abs(m_flow_H2);

    //Energy production
      Q_reac_CO = abs((delta_h_CO2 - delta_h_CO)*R_CO*iCom.Volume);
      Q_reac_CH4= abs((delta_h_CO + 2*delta_h_H2O - delta_h_CH4)*R_CH4*iCom.Volume);
      Q_reac_H2 = abs((delta_h_H2O) * R_H2*iCom.Volume);
      Q_reac = Q_reac_CO + Q_reac_CH4 + Q_reac_H2;

    end GasCombustion;
  end Gascombustion;

  model GasPhase_Tester "Tester for gasphase model"
    extends ClaRa.Basics.Icons.PackageIcons.ExecutableExample100;

    parameter Real xi1=0.1;
    parameter Real xi2=0.1;
    parameter Real xi3=0.1;
    parameter Real xi4=0.1;
    parameter Real xi5=0.1;
    parameter Real xi6=0.1;
    parameter Real xi7=0.1;
    parameter Real xi8=0.1;
    parameter Real xi9=0.1;

    final parameter Real xi_gas[9]={xi1,xi2,xi3,xi4,xi5,xi6,xi7,xi8,xi9};

    Gasphase_Radiation gasphase_Base(
      T_start=273.15 + 25,
      T_start_flueGas_out=273.15 + 25,
      xi_start_flueGas_out={0,0,0,0,0.78,0.21,0,0,0},
      length=1,
      height=1,
      width=1,
      redeclare model GasCombustion = Grate_Boiler.Components.Gasphase.Gascombustion.GasCombustion) annotation (Placement(transformation(extent={{-20,-10},{20,10}})));
    ClaRa.Components.BoundaryConditions.BoundaryGas_Txim_flow boundaryGas_Txim_flow(
      m_flow_const=0,
      variable_T=false,
      T_const=298.15)                                                                   annotation (Placement(transformation(extent={{-64,-10},{-44,10}})));
    ClaRa.Components.BoundaryConditions.BoundaryGas_Txim_flow CH4(m_flow_const=1, xi_const={1,0,0,0,0,0,0,0,0}) annotation (Placement(transformation(
          extent={{-10,-10},{10,10}},
          rotation=90,
          origin={-52,-46})));
    ClaRa.Components.BoundaryConditions.BoundaryGas_pTxi boundaryGas_pTxi annotation (Placement(transformation(
          extent={{-10,-10},{10,10}},
          rotation=-90,
          origin={0,58})));
    inner ClaRa.SimCenter simCenter(redeclare Media.FlueGasTILMedia_GrateBoiler flueGasModel, T_amb=298.15)                                                              annotation (Placement(transformation(extent={{60,-100},{100,-80}})));
    Modelica.Blocks.Sources.Ramp ramp(
      height=6,
      duration=0.8,
      offset=0,
      startTime=0.1)
                    annotation (Placement(transformation(extent={{-96,-90},{-76,-70}})));
    ClaRa.Components.BoundaryConditions.BoundaryGas_Txim_flow O2(
      variable_m_flow=true,
      m_flow_const=1,
      xi_const={0,0,0,0,0,1,0,0,0}) annotation (Placement(transformation(
          extent={{-10,-10},{10,10}},
          rotation=90,
          origin={-26,-48})));
    ClaRa.Components.BoundaryConditions.BoundaryGas_Txim_flow N2(
      variable_m_flow=true,
      m_flow_const=1,
      xi_const={0,0,0,0,1,0,0,0,0}) annotation (Placement(transformation(
          extent={{-10,-10},{10,10}},
          rotation=90,
          origin={0,-46})));
    ClaRa.Components.BoundaryConditions.BoundaryGas_Txim_flow H2(
      variable_T=false,
      m_flow_const=0.05,
      xi_const={0,0,0,0,0,0,0,0,1}) annotation (Placement(transformation(
          extent={{-10,-10},{10,10}},
          rotation=90,
          origin={26,-46})));
    Modelica.Blocks.Math.Gain gain(k=3.29) annotation (Placement(transformation(extent={{-26,-86},{-14,-74}})));
    ClaRa.Components.VolumesValvesFittings.Valves.GenericValveGas_L1 genericValveGas_L1_1(redeclare model PressureLoss = ClaRa.Components.VolumesValvesFittings.Valves.Fundamentals.LinearNominalPoint (Delta_p_nom=1000, m_flow_nom=10)) annotation (Placement(transformation(
          extent={{-10,-6},{10,6}},
          rotation=90,
          origin={0,28})));
    ClaRa.Components.VolumesValvesFittings.Fittings.JoinGas_L2_flex joinGas_L2_flex(N_ports_in=4) annotation (Placement(transformation(
          extent={{-7,-7},{7,7}},
          rotation=90,
          origin={0,-20})));
  equation

    connect(gasphase_Base.air_inlet, boundaryGas_Txim_flow.gas_a) annotation (Line(
        points={{-20,0},{-44,0}},
        color={118,106,98},
        thickness=0.5));
    connect(ramp.y, O2.m_flow) annotation (Line(points={{-75,-80},{-32,-80},{-32,-58}}, color={0,0,127}));
    connect(ramp.y, gain.u) annotation (Line(points={{-75,-80},{-27.2,-80}}, color={0,0,127}));
    connect(gain.y, N2.m_flow) annotation (Line(points={{-13.4,-80},{-6,-80},{-6,-56}}, color={0,0,127}));
    connect(gasphase_Base.fluegas_outlet, genericValveGas_L1_1.inlet) annotation (Line(
        points={{0,10},{0,18}},
        color={118,106,98},
        thickness=0.5));
    connect(genericValveGas_L1_1.outlet, boundaryGas_pTxi.gas_a) annotation (Line(
        points={{0,38},{0,48}},
        color={118,106,98},
        thickness=0.5));
    connect(gasphase_Base.fluegas_inlet, joinGas_L2_flex.outlet) annotation (Line(
        points={{0,-10},{0,-13},{1.33227e-15,-13}},
        color={118,106,98},
        thickness=0.5));
    connect(H2.gas_a, joinGas_L2_flex.inlet[4]) annotation (Line(
        points={{26,-36},{26,-32},{-0.2625,-32},{-0.2625,-27}},
        color={118,106,98},
        thickness=0.5));
    connect(N2.gas_a, joinGas_L2_flex.inlet[3]) annotation (Line(
        points={{4.44089e-16,-36},{0,-36},{0,-32},{-0.0875,-32},{-0.0875,-27}},
        color={118,106,98},
        thickness=0.5));
    connect(O2.gas_a, joinGas_L2_flex.inlet[2]) annotation (Line(
        points={{-26,-38},{-26,-32},{0.0875,-32},{0.0875,-27}},
        color={118,106,98},
        thickness=0.5));
    connect(CH4.gas_a, joinGas_L2_flex.inlet[1]) annotation (Line(
        points={{-52,-36},{-52,-32},{0.2625,-32},{0.2625,-27}},
        color={118,106,98},
        thickness=0.5));
    annotation (
      experiment(
        StopTime=1000,
        Tolerance=1e-05,
        __Dymola_Algorithm="Dassl"),
      __Dymola_experimentSetupOutput(equidistant=false),
      __Dymola_experimentFlags(
        Advanced(GenerateVariableDependencies=false, OutputModelicaCode=false),
        Evaluate=false,
        OutputCPUtime=true,
        OutputFlatModelica=false));
  end GasPhase_Tester;

  model Gasphase_Radiation "homogeneous phase, where combustion of gas components takes place"

     outer ClaRa.SimCenter simCenter;

    inner parameter Boolean useHomotopy=simCenter.useHomotopy "True, if homotopy method is used during initialisation";

    parameter ClaRa.Basics.Units.Length length = 1  "Length of gas volume"  annotation (Dialog(group="Geometry"));
    parameter ClaRa.Basics.Units.Length height = 1  "Hight of gas volume" annotation(Dialog(group="Geometry"));
    parameter ClaRa.Basics.Units.Length width = 1  "width of gas volume" annotation(Dialog(group="Geometry"));

    inner parameter TILMedia.GasTypes.BaseGas flueGas = simCenter.flueGasModel "Medium to be used in volume" annotation(choicesAllMatching, Dialog(group="Fundamental Definitions"));

    ClaRa.Basics.Units.Volume Volume;
    ClaRa.Basics.Units.Mass m_fluegas;
    ClaRa.Basics.Units.EnthalpyMassSpecific h_gas;
    ClaRa.Basics.Units.MassFraction xi_fluegas[flueGas.nc - 1] "Flue gas composition ";
    ClaRa.Basics.Units.MassFraction xi_fluegas_comp[flueGas.nc - 1] "Flue gas composition ";

    ClaRa.Basics.Units.Pressure p;
    ClaRa.Basics.Units.Pressure p_ps(start=p_start_flueGas_out);

  //    ClaRa.Basics.Units.HeatFlowRate Q_flow_bottom "Heat flow from bottom section";
  //    ClaRa.Basics.Units.HeatFlowRate Q_flow_top "Heat flow from top section";
  //    ClaRa.Basics.Units.HeatFlowRate Q_flow_wall "Heat flow from walls";

    Real drhodt  "Density derivative";
  //
    ClaRa.Basics.Units.MassFraction sum_xi_fluegas;
    //pressure loss
    parameter ClaRa.Basics.Units.MassFlowRate m_flow_nom=100  "Nominal mass flow rates at inlet for pressure loss";

   replaceable model PressureLoss = Bedsegment.Pressure_Loss.LinearPressureLoss_L2_Grate_GasPhase
          constrainedby ClaRa.Basics.ControlVolumes.Fundamentals.PressureLoss.PressureLossBaseGas_L2 annotation (Dialog(group="Combustion"), choicesAllMatching=true);

  //__________/ Radiation \___________
   ClaRa.Basics.Interfaces.HeatPort_a heatPort annotation (Placement(transformation(extent={{190,-10},{210,10}}), iconTransformation(extent={{189.5,-10.5},{210.5,10.5}})));
  //Gas2Gas
  //  parameter ClaRa.Basics.Units.Area A_front = length * width;
  //Gas2Wall
  //  parameter Integer N_heat=2 "No. of heat transfer areas";
  //  parameter  ClaRa.Basics.Units.Area A_heat[N_heat](each min=Modelica.Constants.eps) = ones(N_heat) "Heat transfer area: /1/ dedicated to lateral surface" annotation (Dialog(group="Essential Geometry Definition"));
  //  final parameter ClaRa.Basics.Units.Area A_heat_CF[N_heat](each min=Modelica.Constants.eps) = {A_heat[i]*CF_geo[i] for i in 1:N_heat} "Corrected heat transfer area: /1/ dedicated to lateral surface" annotation (Dialog(group="Essential Geometry Definition"));
  // parameter Real CF_geo[N_heat](each min=Modelica.Constants.eps) = ones(N_heat) "Correction factor for heat transfer area: /1/ dedicated to lateral surface" annotation(Dialog(group="Essential Geometry Definition"));

  //_____________/ Start values \_____________________________________________________________
    parameter ClaRa.Basics.Units.Temperature T_start = 298.15 "Start temperature for radiation" annotation (Dialog(tab="Initialisation"));
    parameter ClaRa.Basics.Units.Pressure p_start_flueGas_out=1.013e5 "Start pressure at outlet" annotation (Dialog(tab="Initialisation"));
    parameter ClaRa.Basics.Units.Temperature T_start_flueGas_out=298.15 "Start temperature at outlet" annotation (Dialog(tab="Initialisation"));
    parameter ClaRa.Basics.Units.MassFraction xi_start_flueGas_out[flueGas.nc - 1]={0.01,0,0.1,0,0.74,0.13,0,0.02,0} "Start composition of flue gas" annotation (Dialog(tab="Initialisation"));

    final parameter ClaRa.Basics.Units.EnthalpyMassSpecific h_start = TILMedia.GasFunctions.specificEnthalpy_pTxi(flueGas, p_start_flueGas_out, T_start_flueGas_out, xi_start_flueGas_out) "Start flue gas enthalpy" annotation (Dialog(tab="Initialisation"));

  //__________/ Combustion Model \__________
    replaceable model GasCombustion = Grate_Boiler.Components.Gasphase.Gascombustion.GasCombustion constrainedby Grate_Boiler.Components.Gasphase.Gascombustion.GasCombustion_Base
                                                                                         annotation (Dialog(group="Combustion"), choicesAllMatching=true);

  //__________/ Connectors \_____________
  //    ClaRa.Basics.Interfaces.HeatPort_a heatPort_bottom annotation (Placement(transformation(extent={{50,-110},{70,-90}})));
  //    ClaRa.Basics.Interfaces.HeatPort_a heatPort_top annotation (Placement(transformation(extent={{50,90},{70,110}})));
  //    ClaRa.Basics.Interfaces.HeatPort_a heatPort_wall annotation (Placement(transformation(extent={{150,-10},{170,10}})));
    ClaRa.Basics.Interfaces.GasPortOut fluegas_outlet(Medium=flueGas) annotation (Placement(transformation(extent={{-6,94},{6,106}}), iconTransformation(extent={{-10.75,89.25},{10.75,110.75}})));
    ClaRa.Basics.Interfaces.GasPortIn fluegas_inlet(Medium=flueGas) annotation (Placement(transformation(extent={{-6,-106},{6,-94}}), iconTransformation(extent={{-10.5,-110.5},{10.5,-89.5}})));
    ClaRa.Basics.Interfaces.GasPortIn air_inlet(Medium=flueGas) annotation (Placement(transformation(extent={{-206,-6},{-194,6}}), iconTransformation(extent={{-210.5,-10.5},{-189.5,10.5}})));

  //__________/ Heat transfer \____________
  //__________/ Media Objects \____________
     TILMedia.Gas_pT fluegas_in(
       p=fluegas_inlet.p,
       T=noEvent(actualStream(fluegas_inlet.T_outflow)),
       xi=noEvent(actualStream(fluegas_inlet.xi_outflow)),
       gasType=flueGas)                                      annotation (Placement(transformation(extent={{-22,-94},{-2,-74}})));//redeclare TILMedia.GasTypes.FlueGasTILMedia_GrateBoiler gasType

     TILMedia.Gas_pT fluegas_out(
       gasType=flueGas,
       p=fluegas_outlet.p,
       T=noEvent(actualStream(fluegas_outlet.T_outflow)),
       xi=noEvent(actualStream(fluegas_outlet.xi_outflow))) annotation (Placement(transformation(extent={{6,72},{26,92}})));

     TILMedia.Gas_ph gas(
       gasType=flueGas,
      p(displayUnit="Pa") = p,
       xi=xi_fluegas,
       h=h_gas) annotation (Placement(transformation(extent={{-10,-10},{10,10}})));

     TILMedia.Gas_pT air(
       gasType=flueGas,
       p=air_inlet.p,
       T=noEvent(actualStream(air_inlet.T_outflow)),
       xi=noEvent(actualStream(air_inlet.xi_outflow))) annotation (Placement(transformation(extent={{-180,-12},{-160,8}})));

  //__________/ iCom \__________
  public
   inner Bedsegment.Records.IComFlueGas iCom(
       m_fluegas=m_fluegas,
       Volume=Volume,
       T_gas=gas.T,
       xi_fluegas=xi_fluegas,
       m_flow_nom=m_flow_nom,
       m_flow_in=fluegas_inlet.m_flow + air_inlet.m_flow,
       T_in=fluegas_in.T,
       T_out=fluegas_out.T,
       p=p)
       annotation (Placement(transformation(extent={{178,-102},{198,-82}}))); //
       inner GasCombustion gasCombustion annotation(Placement(transformation(extent={{-60,-10},{-40,10}})));
    PressureLoss pressureLoss(Delta_p_nom=1000)
                                               annotation(Placement(transformation(extent={{30,-10},{50,10}})));

  initial equation
     h_gas = h_start;
     xi_fluegas_comp=xi_start_flueGas_out;
     p=p_ps;
  //   heatPort.T  = T_start;

  equation
  //______/ Geometry \____________
    Volume = length * height * width;  //calculates volume

     m_fluegas = gas.d  * Volume; // Calculates m_fluegas
     //steady state then 0 otherwise drhodt*Volume
     drhodt*Volume =  fluegas_inlet.m_flow + fluegas_outlet.m_flow + air_inlet.m_flow;  // Calculates fluegas_outlet.m_flow //
     drhodt = gas.drhodh_pxi*der(h_gas) + sum({gas.drhodxi_ph[i] * der(gas.xi[i]) for i in 1:flueGas.nc-1}) + gas.drhodp_hxi * der(p_ps);// Calculates drhodt  //suppose drhodp_hxi=0
      der(p) = 1000000*(p_ps-p);
  //Calculates xi_fluegas state
  //     der(xi_fluegas) = 1/max(m_fluegas,Modelica.Constants.eps)*(fluegas_inlet.m_flow*(fluegas_in.xi) + fluegas_outlet.m_flow*(fluegas_out.xi) + air_inlet.m_flow*(air.xi)
  //                     - xi_fluegas*drhodt*Volume + gasCombustion.m_flow_comb*gasCombustion.prod_comp);

    for i in 1:(flueGas.nc-1) loop
      if i==1 then der(xi_fluegas_comp[1]) = 1/max(m_fluegas,Modelica.Constants.eps)*(fluegas_inlet.m_flow*(fluegas_in.xi[1]) + fluegas_outlet.m_flow*(fluegas_out.xi[1]) + air_inlet.m_flow*(air.xi[1]) - xi_fluegas[1]*drhodt*Volume + gasCombustion.m_flow_CH4);//
      else if i==2 then der(xi_fluegas_comp[2]) = 1/max(m_fluegas,Modelica.Constants.eps)*(fluegas_inlet.m_flow  *(fluegas_in.xi[2])
                                                                                    + fluegas_outlet.m_flow*(fluegas_out.xi[2])
                                                                                    + air_inlet.m_flow*(air.xi[2])
                                                                                    - xi_fluegas[2]*drhodt*Volume + gasCombustion.m_flow_CO);//
      else if i==3 then der(xi_fluegas_comp[3]) = 1/max(m_fluegas,Modelica.Constants.eps)*(fluegas_inlet.m_flow*(fluegas_in.xi[3])
        + fluegas_outlet.m_flow*(fluegas_out.xi[3])
        + air_inlet.m_flow*(air.xi[3])
        - xi_fluegas[3]*drhodt*Volume
        + gasCombustion.m_flow_CO2);//
      else if i==6 then der(xi_fluegas_comp[6]) = 1/max(m_fluegas,Modelica.Constants.eps)*(fluegas_inlet.m_flow*(fluegas_in.xi[6])
        + fluegas_outlet.m_flow*(fluegas_out.xi[6])
        + air_inlet.m_flow*(air.xi[6])
        - xi_fluegas[6]*drhodt*Volume
        + gasCombustion.m_flow_O2);//
      else if i==8 then der(xi_fluegas_comp[8]) = 1/max(m_fluegas,Modelica.Constants.eps)*(fluegas_inlet.m_flow*(fluegas_in.xi[8])
        + fluegas_outlet.m_flow*(fluegas_out.xi[8])
        + air_inlet.m_flow*(air.xi[8])
        - xi_fluegas[8]*drhodt*Volume
        + gasCombustion.m_flow_H2O);//
      else if i==9 then der(xi_fluegas_comp[9]) = 1/max(m_fluegas,Modelica.Constants.eps)*(fluegas_inlet.m_flow*(fluegas_in.xi[9])
        + fluegas_outlet.m_flow*(fluegas_out.xi[9])
        + air_inlet.m_flow*(air.xi[9])
        - xi_fluegas[9]*drhodt*Volume + gasCombustion.m_flow_H2);//
      else
        der(xi_fluegas_comp[i]) = 1/max(m_fluegas,Modelica.Constants.eps)*(fluegas_inlet.m_flow*(fluegas_in.xi[i]) + fluegas_outlet.m_flow*(fluegas_out.xi[i]) + air_inlet.m_flow*(air.xi[i]) - xi_fluegas[i]*drhodt*Volume); //
      end if; end if; end if; end if; end if; end if;
    end for;
      for i in 1: (flueGas.nc-1) loop
        xi_fluegas[i]=max(Modelica.Constants.eps,xi_fluegas_comp[i]);
      end for;
      sum_xi_fluegas = sum(xi_fluegas);

  //energy balance
  //Calculates h_gas
      der(h_gas) = 1/max(m_fluegas,Modelica.Constants.eps)*(fluegas_inlet.m_flow*fluegas_in.h
                                                          + fluegas_outlet.m_flow*fluegas_out.h
                                                          + air_inlet.m_flow*air.h
                                                          + Volume*der(p)
                                                          - h_gas*drhodt*Volume
                                                          + gasCombustion.Q_reac
                                                          + heatPort.Q_flow);  //+ heatPort_bottom.Q_flow + heatPort_top.Q_flow + heatPort_wall.Q_flow

  /*
n_flow_C = (fluegas_inlet.m_flow*(fluegas_in.xi[1]) + air_inlet.m_flow*(air.xi[1]))/16 + ((fluegas_in.xi[2]) + air_inlet.m_flow*(air.xi[2]))/28 + ((fluegas_in.xi[3]) + air_inlet.m_flow*(air.xi[3]))/44;
n_flow_H = (fluegas_inlet.m_flow*(fluegas_in.xi[1]) + air_inlet.m_flow*(air.xi[1]))/16 + ((fluegas_in.xi[8]) + air_inlet.m_flow*(air.xi[8]))/18 + ((fluegas_in.xi[9]) + air_inlet.m_flow*(air.xi[9]))/1;
n_flow_O = (fluegas_inlet.m_flow*(fluegas_in.xi[2]) + air_inlet.m_flow*(air.xi[2]))/28 + ((fluegas_in.xi[3]) + air_inlet.m_flow*(air.xi[3]))/44 + ((fluegas_in.xi[6]) + air_inlet.m_flow*(air.xi[6]))/16 + ((fluegas_in.xi[8]) + air_inlet.m_flow*(air.xi[8]))/18;

m_flow_oxygen_req = (n_flow_C + n_flow_H/4 - n_flow_O/2)*ClaRa.Basics.Constants.M_O*2;
m_flow_air_req*max(Modelica.Constants.eps,air.xi[6]) = m_flow_oxygen_req;

lambdaComb = (fluegas_inlet.m_flow*instream(fluegas_inlet.xi_outflow[6]) + air_inlet.m_flow*instream(air_inlet.xi_outflow[6]))/(max(Modelica.Constants.eps, m_flow_oxygen_req));

lambdaCombActual = (fluegas_inlet.m_flow*instream(fluegas_inlet.xi_outflow[6]) + air_inlet.m_flow*instream(air_inlet.xi_outflow[6]))/(fluegas_inlet.m_flow*instream(fluegas_inlet.xi_outflow[6]) + air_inlet.m_flow*instream(air_inlet.xi_outflow[6]) + fluegas_outlet.m_flow*fluegas_out.xi[6]);
*/

  //______________/ pressure loss \___________
     fluegas_inlet.p = p + pressureLoss.Delta_p;  //Calculates p
     air_inlet.p  = p + pressureLoss.Delta_p;
     p = fluegas_outlet.p;  //Calculates fluegas_outlet.p

  // Calculates  fluegas_outlet.xi_outflow and fluegas_inlet.xi_outflow
      fluegas_outlet.xi_outflow = xi_fluegas;
      fluegas_inlet.xi_outflow = xi_fluegas;
      air_inlet.xi_outflow    = xi_fluegas;

  //Temperature dummy values
     fluegas_outlet.T_outflow  = gas.T;
     fluegas_inlet.T_outflow  = gas.T;
     air_inlet.T_outflow    = gas.T;
  //Heat ports temperatures
     heatPort.T  = gas.T;
  //      heatPort_bottom.T = fluegas_in.T;
  //      heatPort_wall.T  = gas.T;

  //____________/ Heat port temperatures and Q_flows \____________________________
  //     Q_flow_wall = heatPort_wall.Q_flow;
  //     Q_flow_top = heatPort_top.Q_flow;
  //     Q_flow_bottom = heatPort_bottom.Q_flow;

  //__________/ Connections \__________

                                                                annotation(Dialog(tab="Initialisation"),
      Diagram(coordinateSystem(extent={{-200,-100},{200,100}})),
      Icon(coordinateSystem(extent={{-200,-100},{200,100}}), graphics={Bitmap(
            extent={{-200,-100},{200,100}},
            imageSource="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"
                 + "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",
            fileName="modelica://ClaRa_Dev/Resources/Images/GasCombustion.png")}));
  end Gasphase_Radiation;
end Gasphase;
