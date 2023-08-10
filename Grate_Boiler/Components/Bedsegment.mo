within Grate_Boiler.Components;
package Bedsegment
extends ClaRa.Basics.Icons.PackageIcons.Components80;
  package Records
  extends ClaRa.Basics.Icons.PackageIcons.Basics50;
    record IComBedSegment
                   //used for exchanging data between model and replaceable models
      extends Grate_Boiler.Components.Bedsegment.Records.IComBedSegment_Base;

    //__________fuel_________________
    //declared in base class
    //__________volatiles____________
      ClaRa.Basics.Units.Mass m_fuel_waf;
      ClaRa.Basics.Units.MassFraction xi_vol[4];
    //__________char conversion______
      ClaRa.Basics.Units.MassFraction xi_fuel[6];
      ClaRa.Basics.Units.MassFraction xi_fluegas[9];
      ClaRa.Basics.Units.Pressure p_fluegas;
      ClaRa.Basics.Units.Mass m_fuel;
      ClaRa.Basics.Units.Mass m_fluegas;
      ClaRa.Basics.Units.Temperature T_fluegas;
    //__________pressure loss_____
      ClaRa.Basics.Units.MassFlowRate m_flow_in;
      ClaRa.Basics.Units.MassFlowRate m_flow_nom;

    end IComBedSegment;

    record IComBedSegment_Base "Record base"
                   //used for exchanging data between model and replaceable models
      extends ClaRa.Basics.Icons.IComIcon;

    //__________fuel_________________
      ClaRa.Basics.Units.Temperature T_fuel "Fuel Temperature";
      ClaRa.Basics.Units.Mass m_fuel_H2O "Water fraction in fuel";
      ClaRa.Basics.Units.Pressure p_fuel "Pressure in segment";

    end IComBedSegment_Base;

    record IComBedSegment_Pyro "Record for testing pyrolysis"
      extends IComBedSegment_Base;

    //__________volatiles____________
      ClaRa.Basics.Units.Mass m_fuel_waf;
      ClaRa.Basics.Units.MassFraction xi_vol[5];

    end IComBedSegment_Pyro;

    record IComFlueGas "Basic internal communication record for heat transfer"
        extends ClaRa.Basics.Icons.IComIcon;

    //____Variables for system description__________________________________________________
      ClaRa.Basics.Units.Mass m_fluegas;
      ClaRa.Basics.Units.MassFraction xi_fluegas[9];
      ClaRa.Basics.Units.Temperature T_in;
      ClaRa.Basics.Units.Temperature T_out;
      ClaRa.Basics.Units.Temperature T_gas;
      ClaRa.Basics.Units.Volume Volume;
      ClaRa.Basics.Units.Pressure p;
     // ClaRa.Basics.Units.Area A_heat_CF;
    //   ClaRa.Basics.Units.Area A_front;
    //__________pressure loss_____
      ClaRa.Basics.Units.MassFlowRate m_flow_in;
      ClaRa.Basics.Units.MassFlowRate m_flow_nom;

    end IComFlueGas;
  end Records;

  package Pressure_Loss "Pressure loss models from ClaRa"
  extends ClaRa.Basics.Icons.PackageIcons.Basics50;
    model LinearPressureLoss_L2_Grate "All geo || Linear pressure loss || Nominal pressure difference"
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
      outer Grate_Boiler.Components.Bedsegment.Records.IComBedSegment iCom;
      extends ClaRa.Basics.ControlVolumes.Fundamentals.PressureLoss.Generic_PL.PressureLoss_L2;
      extends ClaRa.Basics.ControlVolumes.Fundamentals.PressureLoss.TubeTypeVLE_L2;
      extends ClaRa.Basics.ControlVolumes.Fundamentals.PressureLoss.ShellTypeVLE_L2;

      parameter ClaRa.Basics.Units.Pressure Delta_p_nom=10 "Nominal pressure loss";

    equation
      Delta_p = Delta_p_nom/iCom.m_flow_nom*iCom.m_flow_in;

    end LinearPressureLoss_L2_Grate;

    model LinearPressureLoss_L2_Grate_GasPhase "All geo || Linear pressure loss || Nominal pressure difference"
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
      outer Grate_Boiler.Components.Bedsegment.Records.IComFlueGas iCom;
      extends ClaRa.Basics.ControlVolumes.Fundamentals.PressureLoss.Generic_PL.PressureLoss_L2;
      extends ClaRa.Basics.ControlVolumes.Fundamentals.PressureLoss.TubeTypeVLE_L2;
      extends ClaRa.Basics.ControlVolumes.Fundamentals.PressureLoss.ShellTypeVLE_L2;

      parameter ClaRa.Basics.Units.Pressure Delta_p_nom=10 "Nominal pressure loss";

    equation
      Delta_p = Delta_p_nom/iCom.m_flow_nom*iCom.m_flow_in;

    end LinearPressureLoss_L2_Grate_GasPhase;
  end Pressure_Loss;

  package Balance_Tests
    extends ClaRa.Basics.Icons.PackageIcons.Examples50;
    model Evaporation_Test "Test of function of replaceable models"
      extends ClaRa.Basics.Icons.PackageIcons.ExecutableExample100;
      ClaRa.Components.BoundaryConditions.BoundaryFuel_Txim_flow boundaryFuel_Txim_flow(                 variable_T=false, m_flow_const=0,
        xi_const={0.3,0.04,0.2,0.06,0.03,0.27})                                                                                annotation (Placement(transformation(extent={{-60,-10},{-40,10}})));
      ClaRa.Components.BoundaryConditions.BoundaryFuel_pTxi boundaryFuel_pTxi annotation (Placement(transformation(extent={{60,-10},{40,10}})));
      ClaRa.Components.BoundaryConditions.BoundaryGas_Txim_flow boundaryGas_Txim_flow(m_flow_const=1,
        xi_const={0,0,0,0,0.78,0.21,0,0,0},
        T_const=500,
        variable_T=true)                                                                                                    annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=90,
            origin={0,-50})));
      ClaRa.Components.BoundaryConditions.BoundaryGas_pTxi boundaryGas_pTxi annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=-90,
            origin={0,58})));
      inner ClaRa.SimCenter simCenter(redeclare Media.FlueGasTILMedia_GrateBoiler flueGasModel, redeclare Media.Waste_Warncke fuelModel1)
                                          annotation (Placement(transformation(extent={{-100,-100},{-60,-80}})));
      Modelica.Blocks.Sources.Ramp ramp(
        startTime=0,
        duration=0,
        offset=293,
        height=180)  annotation (Placement(transformation(extent={{-40,-80},{-20,-60}})));
      ClaRa.Components.VolumesValvesFittings.Valves.GenericValveGas_L1 adapter(redeclare model PressureLoss = ClaRa.Components.VolumesValvesFittings.Valves.Fundamentals.LinearNominalPoint (Delta_p_nom=10, m_flow_nom=10), checkValve=true) annotation (Placement(transformation(
            extent={{-7,-4},{7,4}},
            rotation=90,
            origin={-4.44089e-16,29})));
      Bedsegment_Base bedsegment_Base_ver_03_1(
        variable_v_grate=false,
        v_grate_const=0,
        m_segment_fuel_start=100,
        T_fuel_start=293.15,
        p_start_flueGas_out(displayUnit="bar") = 101315,
        T_start_flueGas_out=293.15,
        redeclare model Evaporation = Evaporation.Evaporation,
        redeclare model Pyrolysis = ChemicalReactions.SolidFuelReactionZone (A_0=0)) annotation (Placement(transformation(extent={{-20,-10},{20,10}})));
    equation
      connect(ramp.y, boundaryGas_Txim_flow.T) annotation (Line(points={{-19,-70},{0,-70},{0,-60}}, color={0,0,127}));
      connect(adapter.outlet, boundaryGas_pTxi.gas_a) annotation (Line(
          points={{0,36},{0,48}},
          color={118,106,98},
          thickness=0.5));
      connect(adapter.inlet, bedsegment_Base_ver_03_1.fluegas_outlet) annotation (Line(
          points={{0,22},{0,10}},
          color={118,106,98},
          thickness=0.5));
      connect(boundaryFuel_Txim_flow.fuel_a, bedsegment_Base_ver_03_1.fuel_inlet) annotation (Line(
          points={{-40,0},{-20,0}},
          color={27,36,42},
          pattern=LinePattern.Solid,
          thickness=0.5));
      connect(bedsegment_Base_ver_03_1.fuel_outlet, boundaryFuel_pTxi.fuel_a) annotation (Line(
          points={{20,0},{40,0}},
          color={27,36,42},
          pattern=LinePattern.Solid,
          thickness=0.5));
      connect(boundaryGas_Txim_flow.gas_a, bedsegment_Base_ver_03_1.fluegas_inlet) annotation (Line(
          points={{0,-40},{0,-10}},
          color={118,106,98},
          thickness=0.5));
      annotation (experiment(
          StopTime=2000,
          __Dymola_NumberOfIntervals=10000,
          Tolerance=1e-006,
          __Dymola_Algorithm="Dassl"),
        __Dymola_experimentSetupOutput(equidistant=false),
        __Dymola_experimentFlags(
          Advanced(GenerateVariableDependencies=false, OutputModelicaCode=false),
          Evaluate=false,
          OutputCPUtime=false,
          OutputFlatModelica=false));
    end Evaporation_Test;

    model CharCombustion_Test "Test of function of replaceable models"
      extends ClaRa.Basics.Icons.PackageIcons.ExecutableExample100;
      inner ClaRa.SimCenter simCenter(redeclare Media.FlueGasTILMedia_GrateBoiler flueGasModel, redeclare Media.Waste_Warncke fuelModel1)
                                          annotation (Placement(transformation(extent={{-100,-100},{-60,-80}})));
      Modelica.Blocks.Sources.Ramp ramp(
        startTime=0,
        duration=0,
        offset=0,
        height=773.15)
                     annotation (Placement(transformation(extent={{-40,-80},{-20,-60}})));
      ClaRa.Components.BoundaryConditions.BoundaryFuel_Txim_flow boundaryFuel_Txim_flow(
        variable_T=false,
        m_flow_const=0,
        xi_const={0.3,0.04,0.2,0.06,0.03,0.27})                                                                                annotation (Placement(transformation(extent={{-66,2},{-46,22}})));
      ClaRa.Components.BoundaryConditions.BoundaryFuel_pTxi boundaryFuel_pTxi annotation (Placement(transformation(extent={{54,2},{34,22}})));
      ClaRa.Components.BoundaryConditions.BoundaryGas_Txim_flow boundaryGas_Txim_flow(
        m_flow_const=1,
        xi_const={0,0,0,0,0.78,0.21,0,0,0},
        T_const=500,
        variable_T=true)                                                                                                    annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=90,
            origin={-6,-38})));
      ClaRa.Components.BoundaryConditions.BoundaryGas_pTxi boundaryGas_pTxi annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=-90,
            origin={-6,70})));
      ClaRa.Components.VolumesValvesFittings.Valves.GenericValveGas_L1 adapter(redeclare model PressureLoss = ClaRa.Components.VolumesValvesFittings.Valves.Fundamentals.LinearNominalPoint (Delta_p_nom=10, m_flow_nom=10), checkValve=true) annotation (Placement(transformation(
            extent={{-7,-4},{7,4}},
            rotation=90,
            origin={-6,41})));
      Bedsegment_Base bedsegment_Base_ver_03_1(
        variable_v_grate=false,
        v_grate_const=0,
        h_segment=1,
        m_segment_fuel_start=100,
        T_fuel_start=293.15,
        p_start_flueGas_out(displayUnit="bar") = 101315,
        T_start_flueGas_out=293.15,
        redeclare model Evaporation = Grate_Boiler.Components.Evaporation.Evaporation,
        redeclare model Pyrolysis = Grate_Boiler.Components.ChemicalReactions.SolidFuelReactionZone,
        redeclare model Charcombustion = Grate_Boiler.Components.CharCombustion.Charcombustion) annotation (Placement(transformation(extent={{-26,2},{14,22}})));
    equation
      connect(adapter.outlet,boundaryGas_pTxi. gas_a) annotation (Line(
          points={{-6,48},{-6,60}},
          color={118,106,98},
          thickness=0.5));
      connect(adapter.inlet,bedsegment_Base_ver_03_1. fluegas_outlet) annotation (Line(
          points={{-6,34},{-6,22}},
          color={118,106,98},
          thickness=0.5));
      connect(boundaryFuel_Txim_flow.fuel_a,bedsegment_Base_ver_03_1. fuel_inlet) annotation (Line(
          points={{-46,12},{-26,12}},
          color={27,36,42},
          pattern=LinePattern.Solid,
          thickness=0.5));
      connect(bedsegment_Base_ver_03_1.fuel_outlet,boundaryFuel_pTxi. fuel_a) annotation (Line(
          points={{14,12},{34,12}},
          color={27,36,42},
          pattern=LinePattern.Solid,
          thickness=0.5));
      connect(boundaryGas_Txim_flow.gas_a,bedsegment_Base_ver_03_1. fluegas_inlet) annotation (Line(
          points={{-6,-28},{-6,2}},
          color={118,106,98},
          thickness=0.5));
      connect(ramp.y, boundaryGas_Txim_flow.T) annotation (Line(points={{-19,-70},{-6,-70},{-6,-48}}, color={0,0,127}));
      annotation (experiment(
          StopTime=2000,
          __Dymola_NumberOfIntervals=1000,
          Tolerance=1e-006,
          __Dymola_Algorithm="Dassl"),
        __Dymola_experimentSetupOutput(equidistant=false),
        __Dymola_experimentFlags(
          Advanced(GenerateVariableDependencies=false, OutputModelicaCode=false),
          Evaluate=false,
          OutputCPUtime=false,
          OutputFlatModelica=false));
    end CharCombustion_Test;

    model Pyrolisys_Test "Test of function of replaceable models"
      extends ClaRa.Basics.Icons.PackageIcons.ExecutableExample100;
      ClaRa.Components.BoundaryConditions.BoundaryFuel_Txim_flow boundaryFuel_Txim_flow(                 variable_T=false, m_flow_const=0,
        xi_const={0.3,0.04,0.2,0.06,0.03,0.27})                                                                                annotation (Placement(transformation(extent={{-60,-10},{-40,10}})));
      ClaRa.Components.BoundaryConditions.BoundaryFuel_pTxi boundaryFuel_pTxi annotation (Placement(transformation(extent={{60,-10},{40,10}})));
      ClaRa.Components.BoundaryConditions.BoundaryGas_Txim_flow boundaryGas_Txim_flow(m_flow_const=1,
        xi_const={0,0,0,0,0.78,0.21,0,0,0},
        T_const=500,
        variable_T=true)                                                                                                    annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=90,
            origin={0,-50})));
      ClaRa.Components.BoundaryConditions.BoundaryGas_pTxi boundaryGas_pTxi annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=-90,
            origin={0,58})));
      inner ClaRa.SimCenter simCenter(redeclare Media.FlueGasTILMedia_GrateBoiler flueGasModel, redeclare Media.Waste_Warncke fuelModel1)
                                          annotation (Placement(transformation(extent={{-100,-100},{-60,-80}})));
      Modelica.Blocks.Sources.Ramp ramp(
        startTime=0,
        duration=0,
        offset=0,
        height=773.15)
                     annotation (Placement(transformation(extent={{-40,-80},{-20,-60}})));
      ClaRa.Components.VolumesValvesFittings.Valves.GenericValveGas_L1 adapter(redeclare model PressureLoss = ClaRa.Components.VolumesValvesFittings.Valves.Fundamentals.LinearNominalPoint (Delta_p_nom=10, m_flow_nom=10), checkValve=true) annotation (Placement(transformation(
            extent={{-7,-4},{7,4}},
            rotation=90,
            origin={-4.44089e-16,29})));
      Bedsegment_Base bedsegment_Base_ver_03_1(
        variable_v_grate=false,
        v_grate_const=0,
        h_segment=1,
        m_segment_fuel_start=100,
        T_fuel_start=293.15,
        p_start_flueGas_out(displayUnit="bar") = 101315,
        T_start_flueGas_out=293.15,
        redeclare model Evaporation = Grate_Boiler.Components.Evaporation.Evaporation,
        redeclare model Pyrolysis = ChemicalReactions.SolidFuelReactionZone,
        redeclare model Charcombustion = CharCombustion.Charcombustion (A_chem=0)) annotation (Placement(transformation(extent={{-20,-10},{20,10}})));
    equation
      connect(ramp.y, boundaryGas_Txim_flow.T) annotation (Line(points={{-19,-70},{0,-70},{0,-60}}, color={0,0,127}));
      connect(adapter.outlet, boundaryGas_pTxi.gas_a) annotation (Line(
          points={{0,36},{0,48}},
          color={118,106,98},
          thickness=0.5));
      connect(adapter.inlet, bedsegment_Base_ver_03_1.fluegas_outlet) annotation (Line(
          points={{0,22},{0,10}},
          color={118,106,98},
          thickness=0.5));
      connect(boundaryFuel_Txim_flow.fuel_a, bedsegment_Base_ver_03_1.fuel_inlet) annotation (Line(
          points={{-40,0},{-20,0}},
          color={27,36,42},
          pattern=LinePattern.Solid,
          thickness=0.5));
      connect(bedsegment_Base_ver_03_1.fuel_outlet, boundaryFuel_pTxi.fuel_a) annotation (Line(
          points={{20,0},{40,0}},
          color={27,36,42},
          pattern=LinePattern.Solid,
          thickness=0.5));
      connect(boundaryGas_Txim_flow.gas_a, bedsegment_Base_ver_03_1.fluegas_inlet) annotation (Line(
          points={{0,-40},{0,-10}},
          color={118,106,98},
          thickness=0.5));
      annotation (experiment(
          StopTime=2000,
          __Dymola_NumberOfIntervals=10000,
          Tolerance=1e-006,
          __Dymola_Algorithm="Dassl"),
        __Dymola_experimentSetupOutput(equidistant=false),
        __Dymola_experimentFlags(
          Advanced(GenerateVariableDependencies=false, OutputModelicaCode=false),
          Evaluate=false,
          OutputCPUtime=false,
          OutputFlatModelica=false));
    end Pyrolisys_Test;
  end Balance_Tests;

  model Bedsegment_Base "Segment of the discretized bed"
    import GrateBoiler_Code =
           Grate_Boiler;

    outer ClaRa.SimCenter simCenter;

     inner parameter Boolean useHomotopy=simCenter.useHomotopy "True, if homotopy method is used during initialisation";

    parameter ClaRa.Basics.Units.MassFraction epsilon=0.3 "Ratio of fluegas volume and volume of segment";  // minimum porosity is needed, so porosity of fuel is taken into account

    parameter Boolean variable_v_grate=true "True, if v_gate defined by variable input" annotation(Dialog(group="Geometry"));
    parameter ClaRa.Basics.Units.Velocity v_grate_const = 0.1  "Mean velocity of forward movement";
    Modelica.Blocks.Interfaces.RealInput v_grate_variable(value=v_grate) if (variable_v_grate) "Variable v_grate" annotation (Placement(transformation(extent={{-228,76},{-188,116}})));
    ClaRa.Basics.Units.Velocity v_grate "Velocity of the grate";  //Should be declared in global model of total bed - define as outer?
    //inner parameter ClaRa.Basics.Units.Velocity v_grate = 0.01 "Velocity of the grate";  //Should be declared in global model of total bed - define as outer?

    parameter ClaRa.Basics.Units.Length l_segment = 1  "Length of bed segment"  annotation (Dialog(group="Geometry"));
    parameter ClaRa.Basics.Units.Length h_segment = 0.1  "Hight of bed segment" annotation (Dialog(group="Geometry"));
    parameter ClaRa.Basics.Units.Length w_segment = 1  "width of bed segment" annotation (Dialog(group="Geometry"));
    ClaRa.Basics.Units.Length height_fuel;
    parameter Real A_specific( unit= "m2/kg") = 1.2  "Specific surface of fuel";  //value comes from Hamel and Krumm (Euler-Euler Source) maybe *V_fluegas/V_fuel

    inner parameter ClaRa.Basics.Media.FuelTypes.BaseFuel fuelModel=simCenter.fuelModel1 "Fuel elemental composition used for combustion" annotation (choices(choice=simCenter.fuelModel "Fuel model 1 as defined in simCenter"), Dialog(group="Media Definitions"));
    inner parameter TILMedia.GasTypes.BaseGas flueGas = simCenter.flueGasModel "Medium to be used as fluegas" annotation(choicesAllMatching, Dialog(group="Fundamental Definitions"));

  //________________/ Flue Gas Composition \_____________________
    ClaRa.Basics.Units.MassFraction xi_fluegas[flueGas.nc - 1](start=xi_start_flueGas_out) "Flue gas composition ";

  //________________/ Fuel Composition \_____________________
    ClaRa.Basics.Units.MassFraction xi_fuel[fuelModel.N_c-1];
  //________________/ Volatile Composition \__________
    ClaRa.Basics.Units.MassFraction xi_vol[4](start={1,1,1,1}/4) "internal vector for volatile calculation C,H,O";

  //Declaration
    ClaRa.Basics.Units.Volume V_segment(start=1);// = l_segment*w_segment*h_segment;
    ClaRa.Basics.Units.Volume V_fluegas;
    ClaRa.Basics.Units.Volume V_fuel;

    ClaRa.Basics.Units.HeatFlowRate Q_fuel_fluegas;
    parameter ClaRa.Basics.Units.CoefficientOfHeatTransfer alpha_fuel_fluegas=300;  //heat transfer coefficient between fluegas and fuel
    ClaRa.Basics.Units.HeatFlowRate Q_flow_top;
    ClaRa.Basics.Units.HeatFlowRate Q_flow_bottom;

    ClaRa.Basics.Units.Area A_fuel_fluegas;   //contact surface between fluegas and fuel - depends on fuel density

    ClaRa.Basics.Units.Mass m_segment_fuel;//(start=m_segment_fuel_start);
    parameter ClaRa.Basics.Units.Mass m_segment_fuel_start = 10  annotation (Dialog(tab="Initialisation"));//500*(1-epsilon)*V_segment;
    ClaRa.Basics.Units.Mass m_segment_fluegas;
    parameter ClaRa.Basics.Units.MassFraction xi_H2O_start=0.124 "Start composition of H2O" annotation (Dialog(tab="Initialisation"));
    ClaRa.Basics.Units.Mass m_fuel_H2O(start=1.24);

    ClaRa.Basics.Units.MassFlowRate m_flow_fuel_out(start=0.1);
    ClaRa.Basics.Units.Mass m_fuel_waf;    //Mass of volatiles is calculated from m_flow_fuel

    Real drhodt_fluegas(unit = "kg/(m3.s)");  //time derivative of density of fluegas
    Real drhodt_fuel(unit = "kg/(m3.s)");  //time derivative of density of fuel
    Real sum_xi_fg;
    Real sum_xi_fuel;

    ClaRa.Basics.Units.EnthalpyMassSpecific h_fuel(start = -5000);
    ClaRa.Basics.Units.EnthalpyMassSpecific h_gas(start = 2000);

    ClaRa.Basics.Units.Temperature T_ref = fuelModel.T_ref;
    parameter ClaRa.Basics.Units.Temperature T_fuel_start = 293.15 "Start temperature" annotation (Dialog(tab="Initialisation"));
    ClaRa.Basics.Units.Temperature T_fuel  "Temperature of fuel inside segment";
    ClaRa.Basics.Units.Pressure p(start=p_start_flueGas_out) "Pressure inside volume";

  parameter ClaRa.Basics.Units.MassFlowRate m_flow_nom=100;

    //_____________/ Start values \_____________________________________________________________
  public
    parameter ClaRa.Basics.Units.Pressure p_start_flueGas_out=1.013e5+200 "Start pressure at outlet" annotation (Dialog(tab="Initialisation"));
    parameter ClaRa.Basics.Units.Temperature T_start_flueGas_out=293.15 "Start temperature at outlet" annotation (Dialog(tab="Initialisation"));
    parameter ClaRa.Basics.Units.MassFraction xi_start_flueGas_out[flueGas.nc - 1]={0.001,0,0.1,0,0.74,0.13,0,0.02,0} "Start composition of flue gas" annotation (Dialog(tab="Initialisation"));

    final parameter ClaRa.Basics.Units.EnthalpyMassSpecific h_start = TILMedia.GasFunctions.specificEnthalpy_pTxi(flueGas, p_start_flueGas_out, T_start_flueGas_out, xi_start_flueGas_out) "Start flue gas enthalpy";

    parameter ClaRa.Basics.Units.MassFraction xi_fuel_start[fuelModel.N_c-1]={0.302,0.043,0.192,0.006,0.003,0.29}  "Start composition of fuel" annotation (Dialog(tab="Initialisation"));

  //_____________/ Media Poperties \_______________
  //   TILMedia.VLEFluid H2O_props_in(redeclare TILMedia.VLEFluidTypes.TILMedia_SplineWater vleFluidType);
  //   ClaRa.Basics.Units.DensityMassSpecific rho_bubble = TILMedia.VLEFluidObjectFunctions.bubbleDensity_pxi(p, {1}, H2O_props_in.vleFluidPointer); //bubble enthalpy of water
  //
  //_______________/ Pressure Loss \__________
    replaceable model PressureLoss = GrateBoiler_Code.Components.Bedsegment.Pressure_Loss.LinearPressureLoss_L2_Grate
        constrainedby ClaRa.Basics.ControlVolumes.Fundamentals.PressureLoss.PressureLossBaseGas_L2 annotation (Dialog(tab="Replaceable Models"), choicesAllMatching=true);

    //_____________/ Replaceable Models \___________
    replaceable model Evaporation = GrateBoiler_Code.Components.Evaporation.Evaporation
      constrainedby GrateBoiler_Code.Components.Evaporation.PartialEvaporationModel
      annotation (Dialog(tab="Replaceable Models"), choicesAllMatching=true);

    replaceable model Pyrolysis = GrateBoiler_Code.Components.ChemicalReactions.SolidFuelReactionZone
      constrainedby GrateBoiler_Code.Components.ChemicalReactions.PartialReactionZonePyro
      annotation (Dialog(tab="Replaceable Models"), choicesAllMatching=true);

    replaceable model Charcombustion = GrateBoiler_Code.Components.CharCombustion.Charcombustion
      constrainedby GrateBoiler_Code.Components.CharCombustion.Partial_Charcombustion
      annotation (Dialog(tab="Replaceable Models"), choicesAllMatching=true);

    //______________/ Connectors \____________
  ClaRa.Basics.Interfaces.Fuel_inlet fuel_inlet(fuelModel=fuelModel) annotation (Placement(transformation(extent={{-206,-6},{-194,6}}),
                          iconTransformation(extent={{-211,-11},{-189,11}})));
  ClaRa.Basics.Interfaces.Fuel_outlet fuel_outlet(fuelModel=fuelModel) annotation (Placement(transformation(extent={{194,-6},{206,6}}),
                         iconTransformation(extent={{190,-10},{210,10}})));

  ClaRa.Basics.Interfaces.GasPortOut fluegas_outlet(Medium=flueGas) annotation (Placement(transformation(extent={{-6,92},
              {6,104}}), iconTransformation(extent={{-11,89},{11,111}})));
  ClaRa.Basics.Interfaces.GasPortIn fluegas_inlet(Medium=flueGas) annotation (Placement(transformation(extent={{-6,-104},
              {6,-92}}), iconTransformation(extent={{-11,-111},{11,-89}})));

    ClaRa.Basics.Interfaces.HeatPort_b heatPort_bottom annotation (Placement(transformation(extent={{50,-110},{70,-90}})));
    ClaRa.Basics.Interfaces.HeatPort_a heatPort_top annotation (Placement(transformation(extent={{50,90},{70,110}})));

  //______________/ Media Objects \__________
  // protected
    TILMedia.Gas_ph gas(
      gasType=flueGas,
      p(displayUnit="Pa") = fluegas_inlet.p,
      xi=xi_fluegas,
      h=h_gas)         annotation (Placement(transformation(extent={{-10,-12},{10,
              8}})));

    TILMedia.Gas_pT gas_in(
      gasType=flueGas,
      p=fluegas_inlet.p,
      T=noEvent(actualStream(fluegas_inlet.T_outflow)),
      xi=noEvent(actualStream(fluegas_inlet.xi_outflow)))
      annotation (Placement(transformation(extent={{-10,-86},{10,-66}})));

    TILMedia.Gas_pT gas_out(
      gasType=flueGas,
      p=fluegas_outlet.p,
      T=noEvent(actualStream(fluegas_outlet.T_outflow)),
      xi=noEvent(actualStream(fluegas_outlet.xi_outflow)))
      annotation (Placement(transformation(extent={{-10,66},{10,86}})));

    ClaRa.Basics.Media.FuelObject       fuelObject_in(
      fuelModel=fuelModel,
      p=fuel_inlet.p,
      T=noEvent(actualStream(fuel_inlet.T_outflow)),
      xi_c=noEvent(actualStream(fuel_inlet.xi_outflow))) annotation (Placement(transformation(extent={{-190,-30},{-170,-10}})));

    ClaRa.Basics.Media.FuelObject       fuelObject_out(
      fuelModel=fuelModel,
      p=fuel_outlet.p,
      T=noEvent(actualStream(fuel_outlet.T_outflow)),
      xi_c=noEvent(actualStream(fuel_outlet.xi_outflow))) annotation (Placement(transformation(extent={{170,-30},{190,-10}})));
    ClaRa.Basics.Media.FuelObject       fuelObject_bulk(
      fuelModel=fuelModel,
      p=fuel_inlet.p,
      T=T_fuel,
      xi_c=xi_fuel) annotation (Placement(transformation(extent={{-10,-34},{10,-14}})));

  //________________/ iCom \____________
  public
    inner Records.IComBedSegment      iCom(
      T_fuel=T_fuel,
      T_fluegas=gas.T,
      m_fuel_H2O=m_fuel_H2O,
      m_fuel=m_segment_fuel,
      m_fluegas=m_segment_fluegas,
      p_fuel=fuel_inlet.p,
      p_fluegas=p,
      m_fuel_waf=m_fuel_waf,
      xi_vol=xi_vol,
      xi_fuel=xi_fuel,
      xi_fluegas=xi_fluegas,
      m_flow_nom=m_flow_nom,
      m_flow_in=fluegas_inlet.m_flow)
        annotation (Placement(transformation(extent={{176,-100},{196,-80}})));

  public
    inner Evaporation evaporation_model annotation(Placement(transformation(extent={{-60,40},{-40,60}})));
    inner Pyrolysis solidFuelReactionZone annotation (Placement(transformation(extent={{-60,-10},{-40,10}})));
    inner Charcombustion charcombustion(n=0.36)
                                        annotation (Placement(transformation(extent={{-60,-60},{-40,-40}})));
     PressureLoss pressureLoss(Delta_p_nom=10)
                               annotation (Placement(transformation(extent={{20,-10},{40,10}})));

  //::::::::::::::::::::::::/ Equations \::::::::::::::::::::::::
  initial equation
    h_gas = h_start;
    xi_fuel=xi_fuel_start;
    xi_fluegas=xi_start_flueGas_out;
   // m_fuel_waf = 1;
    T_fuel = T_fuel_start;
    m_segment_fuel=m_segment_fuel_start;

  equation
    if (not variable_v_grate) then
      v_grate=v_grate_const;
    end if;

    /*_______________/ Solid Phase \____________*/

  //Porosity of bed and mass
    V_fluegas = V_segment-V_fuel; //Calculates V_fluegas
    //V_fuel = (1-epsilon)*V_segment; //0.4 starting value

    height_fuel = V_fuel/(l_segment*w_segment);  //Calculates heigth_fuel
    epsilon = V_fluegas/V_segment;

    m_segment_fuel = fuelObject_bulk.rho*V_fuel; // Calculates V_fuel assuming known xi_fuel (below eq. for xi_fuel)
    m_fuel_H2O = m_segment_fuel*fuelObject_bulk.xi_h2o; // Calculates mass of water in fuel

  //information for pyrolysis
    //mass of volatiles, therefor only components C,H,O are needed
    m_fuel_waf = m_segment_fuel*(1  - xi_fuel[5] - xi_fuel[6] - (1-sum(xi_fuel)));   //substitute sulphur, ash and water from fuel, fuel for pyrolysis calculation - xi_fuel[4]
  //for volatile calculation
    for i in 1:4 loop
      xi_vol[i] = xi_fuel[i]*m_segment_fuel/m_fuel_waf;
    end for;

  //Calculation of surface area added here
    A_fuel_fluegas = m_segment_fuel*A_specific;  //adding porosity factors for each species //Adding V-ratio *V_fluegas/V_fuel
    m_flow_fuel_out = -v_grate*w_segment*height_fuel*fuelObject_out.rho;  //outgoing mass flow of fuel, due to the movement of the bed by the grate

    fuel_outlet.m_flow = m_flow_fuel_out;  //Give Connector the m_flow_fuel_out

  //mass balance
    der(m_segment_fuel) = fuel_inlet.m_flow + m_flow_fuel_out - evaporation_model.m_flow_fuel_evap - solidFuelReactionZone.m_flow_volatiles - charcombustion.m_flow_fuel; //Calculates m_segment_fuel
    drhodt_fuel = 1/V_fuel*(der(m_segment_fuel) - fuelObject_bulk.rho*der(V_fuel)); //Calculates drhodt_fuel (simplified from original: (V_fuel*(der(m_segment_fuel) - m_segment_fuel*der(V_fuel))/ V_fuel^2)

  //______________/ component balance \____________
  //Calculates xi_fuel state
  //   der(xi_fuel) = 1/max(m_segment_fuel,Modelica.Constants.eps)*(fuel_inlet.m_flow*fuelObject_in.xi_c + m_flow_fuel_out*fuelObject_out.xi_c - xi_fuel*der(m_segment_fuel));

    for i in 1:(fuelModel.N_c-1) loop
             if i==1 then der(xi_fuel[1]) = 1/max(m_segment_fuel,Modelica.Constants.eps)*(fuel_inlet.m_flow*fuelObject_in.xi_c[1] + fuel_outlet.m_flow*fuelObject_out.xi_c[1] - xi_fuel[1]*der(m_segment_fuel) - charcombustion.m_flow_comb_C - solidFuelReactionZone.m_flow_vol_C);  //*m_fuel_waf/max(m_segment_fuel,Modelica.Constants.eps)
        else if i==2 then der(xi_fuel[2]) = 1/max(m_segment_fuel,Modelica.Constants.eps)*(fuel_inlet.m_flow*fuelObject_in.xi_c[2] + fuel_outlet.m_flow*fuelObject_out.xi_c[2] - xi_fuel[2]*der(m_segment_fuel) - solidFuelReactionZone.m_flow_vol_H);  //*m_fuel_waf/max(m_segment_fuel,Modelica.Constants.eps)
        else if i==3 then der(xi_fuel[3]) = 1/max(m_segment_fuel,Modelica.Constants.eps)*(fuel_inlet.m_flow*fuelObject_in.xi_c[3] + fuel_outlet.m_flow*fuelObject_out.xi_c[3] - xi_fuel[3]*der(m_segment_fuel) - solidFuelReactionZone.m_flow_vol_O);  //*m_fuel_waf/max(m_segment_fuel,Modelica.Constants.eps)
        else if i==4 then der(xi_fuel[4]) = 1/max(m_segment_fuel,Modelica.Constants.eps)*(fuel_inlet.m_flow*fuelObject_in.xi_c[4] + fuel_outlet.m_flow*fuelObject_out.xi_c[4] - xi_fuel[4]*der(m_segment_fuel) - solidFuelReactionZone.m_flow_vol_N);   //- charcombustion.m_flow_comb_NO*14/30) //(1+ClaRa.Basics.Constants.M_O/ClaRa.Basics.Constants.M_N)
        else if i==5 then der(xi_fuel[5]) = 1/max(m_segment_fuel,Modelica.Constants.eps)*(fuel_inlet.m_flow*fuelObject_in.xi_c[5] + fuel_outlet.m_flow*fuelObject_out.xi_c[5] - xi_fuel[5]*der(m_segment_fuel) - charcombustion.m_flow_comb_SO2*32/64);  //(1+ClaRa.Basics.Constants.M_O2/ClaRa.Basics.Constants.M_S)
      else
        der(xi_fuel[i]) = 1/max(m_segment_fuel,Modelica.Constants.eps)*(fuel_inlet.m_flow*fuelObject_in.xi_c[i] + m_flow_fuel_out*fuelObject_out.xi_c[i] - xi_fuel[i]*der(m_segment_fuel));
        end if;
        end if;
       end if;
      end if;
     end if;
    end for;

  //energy balance
  //Calculates h_fuel
    der(h_fuel) = 1/max(m_segment_fuel,Modelica.Constants.eps)*(fuel_inlet.m_flow*fuelObject_in.h + m_flow_fuel_out*fuelObject_out.h
                 - evaporation_model.H_evap_fuel
                 - solidFuelReactionZone.H_pyrolysis_fuel
                 - charcombustion.Q_comb_fuel + Q_fuel_fluegas
                 - h_fuel*drhodt_fuel*V_fuel
                 - fuelObject_bulk.rho*h_fuel*der(V_fuel) + Q_flow_bottom + Q_flow_top);   // + p*der(V_fuel) + V_fuel*der(p)

   //Calculates T_fuel
    h_fuel = fuelObject_bulk.cp*(T_fuel - T_ref); //Suppose cp is independent of temperature
   //Calculates Q_fuel_fluegas
    Q_fuel_fluegas = alpha_fuel_fluegas*A_fuel_fluegas*(gas.T - T_fuel);

  //______________/ Fuel Objects \___________
     fuel_inlet.p = fuel_outlet.p;  //Calculates fuel_outlet.p //Assume no pressure loss since it is solid
     fuel_inlet.T_outflow = T_fuel; //Calculates fuel_inlet.T_outflow
     fuel_outlet.T_outflow = T_fuel; //Calculates fuel_outlet.T_outflow

  /*______________/ Gas Phase \______________*/
    m_segment_fluegas = gas.d * V_fluegas; // Calculates m_segment_fluegas

     drhodt_fluegas*V_fluegas = - gas.d*der(V_fluegas) + fluegas_inlet.m_flow + fluegas_outlet.m_flow + evaporation_model.m_flow_fuel_evap + solidFuelReactionZone.m_flow_volatiles + charcombustion.m_flow_fluegas; // Calculates fluegas_outlet.m_flow
     drhodt_fluegas = gas.drhodh_pxi*der(h_gas) + sum({gas.drhodxi_ph[i] * der(gas.xi[i]) for i in 1:flueGas.nc-1}) + gas.drhodp_hxi * der(p); // Calculates drhodt_fluegas + gas.drhodp_hxi * der(p)

     //Calculates xi_fluegas state
    der(xi_fluegas) = 1/max(m_segment_fluegas,Modelica.Constants.eps)*(fluegas_inlet.m_flow*(gas_in.xi)
                      + fluegas_outlet.m_flow*(gas_out.xi) + evaporation_model.m_flow_fuel_evap*({0,0,0,0,0,0,0,1,0})
                      - xi_fluegas*drhodt_fluegas*V_fluegas - xi_fluegas*gas.d*der(V_fluegas)
                      + solidFuelReactionZone.prod_comp*solidFuelReactionZone.m_flow_volatiles
                      + charcombustion.m_flow_fluegas*charcombustion.prod_comp);

  //energy balance
  //Calculates h_gas
    der(h_gas) = 1/max(m_segment_fluegas,Modelica.Constants.eps)*(fluegas_inlet.m_flow*gas_in.h + fluegas_outlet.m_flow*gas_out.h
                 + evaporation_model.H_evap
                 + solidFuelReactionZone.H_pyrolysis_fluegas
                 + charcombustion.Q_comb_fluegas - Q_fuel_fluegas
                 + p*der(V_fluegas) + V_fluegas*der(p) - h_gas*drhodt_fluegas*V_fluegas - gas.d*h_gas*der(V_fluegas));

  //definition of composition at connector outflow
  //Calculates xi according to size of xi vector in connector
    fuel_inlet.xi_outflow = xi_fuel;
    fuel_outlet.xi_outflow = xi_fuel;

  //______________/ pressure loss \___________
    fluegas_inlet.p = p + pressureLoss.Delta_p;  //Calculates p (p is meant to be pressure of fluegas) no pressure loss due to stream through porose material delta_p
    p = fluegas_outlet.p; //Calculates fluegas_outlet.p

  //______________/ component balance \_____________
  // Calculates  fluegas_outlet.xi_outflow and fluegas_inlet.xi_outflow
    fluegas_outlet.xi_outflow = xi_fluegas;
    fluegas_inlet.xi_outflow = xi_fluegas;

  // Calculates  fluegas_outlet.T_outflow and fluegas_inlet.T_outflow
    fluegas_inlet.T_outflow = gas.T;
    fluegas_outlet.T_outflow = gas.T;

  //heatPorts
    Q_flow_top = heatPort_top.Q_flow;
    Q_flow_bottom = heatPort_bottom.Q_flow;
    heatPort_bottom.T = T_fuel;
    heatPort_top.T = T_fuel;

    //check if sum of xi for fuel and fluegas is one
    sum_xi_fg   = sum(xi_fluegas);
    sum_xi_fuel = sum(xi_fuel);

    annotation (Diagram(coordinateSystem(extent={{-200,-100},{200,100}})), Icon(coordinateSystem(extent={{-200,-100},{200,100}}), graphics={Bitmap(
            extent={{-200,-100},{200,100}},
            imageSource="iVBORw0KGgoAAAANSUhEUgAABjkAAAMcCAYAAAACJMfeAAAABHNCSVQICAgIfAhkiAAAAAlwSFlzAAAkwQAAJMEB5obOCQAAABl0RVh0U29mdHdhcmUAd3d3Lmlua3NjYXBlLm9yZ5vuPBoAACAASURBVHic7N1nmFXluQbgd6giChYUKxYQa+wt0cSYGEti771EjdEkxt479l5QQVCqVKWIYBdBLEgVAZGqOCBNOgxM2+eHJ56TKFFgzV6zZu77H8PM876ZKwKzn/2tr2DLZrvkAgAAAAAAIGNqpL0AAAAAAADAmlByAAAAAAAAmaTkAAAAAAAAMknJAQAAAAAAZJKSAwAAAAAAyCQlBwAAAAAAkElKDgAAAAAAIJOUHAAAAAAAQCYpOQAAAAAAgExScgAAAAAAAJmk5AAAAAAAADJJyQEAAAAAAGSSkgMAAAAAAMgkJQcAAAAAAJBJSg4AAAAAACCTlBwAAAAAAEAmKTkAAAAAAIBMUnIAAAAAAACZpOQAAAAAAAAySckBAAAAAABkkpIDAAAAAADIJCUHAAAAAACQSUoOAAAAAAAgk5QcAAAAAABAJik5AAAAAACATKqVr0G77rxT7Ln7bvkaBwAAAAAApGDEqE9jwsRJeZmVt5Jjz913i0suvCBf4wAAAAAAgBQ8+exzeSs5PK4KAAAAAADIJCUHAAAAAACQSUoOAAAAAAAgk5QcAAAAAABAJik5AAAAAACATFJyAAAAAAAAmaTkAAAAAAAAMknJAQAAAAAAZJKSAwAAAAAAyCQlBwAAAAAAkElKDgAAAAAAIJOUHAAAAAAAQCYpOQAAAAAAgExScgAAAAAAAJmk5AAAAAAAADJJyQEAAAAAAGSSkgMAAAAAAMgkJQcAAAAAAJBJSg4AAAAAACCTlBwAAAAAAEAmKTkAAAAAAIBMUnIAAAAAAACZpOQAAAAAAAAySckBAAAAAABkkpIDAAAAAADIJCUHAAAAAACQSUoOAAAAAAAgk5QcAAAAAABAJik5AAAAAACATFJyAAAAAAAAmaTkAAAAAAAAMknJAQAAAAAAZJKSAwAAAAAAyCQlBwAAAAAAkElKDgAAAAAAIJOUHAAAAAAAQCYpOQAAAAAAgExScgAAAAAAAJmk5AAAAAAAADJJyQEAAAAAAGSSkgMAAAAAAMgkJQcAAAAAAJBJSg4AAAAAACCTlBwAAAAAAEAmKTkAAAAAAIBMUnIAAAAAAACZpOQAAAAAAAAySckBAAAAAABkkpIDAAAAAADIJCUHAAAAAACQSUoOAAAAAAAgk5QcAAAAAABAJik5AAAAAACATFJyAAAAAAAAmaTkAAAAAAAAMknJAQAAAAAAZJKSAwAAAAAAyCQlBwAAAAAAkElKDgAAAAAAIJOUHAAAAAAAQCYpOQAAAAAAgExScgAAAAAAAJmk5AAAAAAAADJJyQEAAAAAAGSSkgMAAAAAAMgkJQcAAAAAAJBJSg4AAAAAACCTlBwAAAAAAEAmKTkAAAAAAIBMUnIAAAAAAACZpOQAAAAAAAAySckBAAAAAABkkpIDAAAAAADIJCUHAAAAAACQSUoOAAAAAAAgk5QcAAAAAABAJik5AAAAAACATFJyAAAAAAAAmaTkAAAAAAAAMknJAQAAAAAAZJKSAwAAAAAAyCQlBwAAAAAAkElKDgAAAAAAIJOUHAAAAAAAQCYpOQAAAAAAgExScgAAAAAAAJmk5AAAAAAAADJJyQEAAAAAAGSSkgMAAAAAAMgkJQcAAAAAAJBJSg4AAAAAACCTlBwAAAAAAEAmKTkAAAAAAIBMUnIAAAAAAACZpOQAAAAAAAAySckBAAAAAABkkpIDAAAAAADIJCUHAAAAAACQSUoOAAAAAAAgk5QcAAAAAABAJik5AAAAAACATFJyAAAAAAAAmaTkAAAAAAAAMknJAQAAAAAAZJKSAwAAAAAAyCQlBwAAAAAAkElKDgAAAAAAIJOUHAAAAAAAQCYpOQAAAAAAgExScgAAAAAAAJmk5AAAAAAAADJJyQEAAAAAAGSSkgMAAAAAAMgkJQcAAAAAAJBJSg4AAAAAACCTlBwAAAAAAEAmKTkAAAAAAIBMUnIAAAAAAACZpOQAAAAAAAAySckBAAAAAABkkpIDAAAAAADIJCUHAAAAAACQSUoOAAAAAAAgk5QcAAAAAABAJik5AAAAAACATFJyAAAAAAAAmaTkAAAAAAAAMknJAQAAAAAAZJKSAwAAAAAAyCQlBwAAAAAAkElKDgAAAAAAIJOUHAAAAAAAQCYpOQAAAAAAgExScgAAAAAAAJmk5AAAAAAAADJJyQEAAAAAAGSSkgMAAAAAAMgkJQcAAAAAAJBJSg4AAAAAACCTlBwAAAAAAEAmKTkAAAAAAIBMUnIAAAAA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"
                 + "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"
                 + "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",
            fileName="modelica://ClaRa_Dev/Resources/Images/BedSegment_Icon.png")}),
      Documentation(info="<html>
<p>In this version of the bedsegment base, heatports and a pressure loss model are implementet</p>
</html>"));
  end Bedsegment_Base;
end Bedsegment;
