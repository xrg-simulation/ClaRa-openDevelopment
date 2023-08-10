within ClaRa_Solar.Examples.StaticCycles;
model Static_Cycle_SixBranchEvapThreeBranchSH_Turbine
  extends ClaRa.Basics.Icons.Init;
  import TILMedia.Internals.VLEFluidConfigurations.FullyMixtureCompatible.VLEFluidFunctions.*;
  import SI = ClaRa.Basics.Units;
  outer ClaRa.SimCenter simCenter;
  inner parameter Real P_target_= 1 "Value of load in p.u."    annotation(Dialog(group="Global parameter"));
  parameter TILMedia.VLEFluidTypes.BaseVLEFluid medium = simCenter.fluid1 "Medium in the component"
    annotation(choices(choice=simCenter.fluid1 "First fluid defined in global simCenter",
                       choice=simCenter.fluid2 "Second fluid defined in global simCenter",
                       choice=simCenter.fluid3 "Third fluid defined in global simCenter"),
                                                          Dialog(group="Fundamental Definitions"));

  parameter Modelica.Units.SI.SpecificEnthalpy h_evap=650e3  "Value of additional enthalpy in evaporator"    annotation(Dialog(group="Parrabolic Tube"));
  parameter Modelica.Units.SI.SpecificEnthalpy h_sh=400e3  "Value of additional enthalpy in superheater"    annotation(Dialog(group="Parrabolic Tube"));
  parameter Modelica.Units.SI.Length d_i=0.065;
  parameter Modelica.Units.SI.Length d_a=0.07;
  parameter Modelica.Units.SI.Length L=720;
  parameter Modelica.Units.SI.Length W=4.6;
  parameter Integer nNodes=36;
  parameter Integer nCells=3;
  parameter Integer nCellsHead=3;
  parameter Real zeta=8000;
  parameter Real deltaTime=lengthDistortion/velocityDistortion;
  parameter Real deltaIrr_rel=0.6;
  parameter Real deltaIrr=deltaIrr_rel*850;
  parameter Real Irr=850;
  parameter Modelica.Units.SI.SpecificEnthalpy h_in=600e3;
  parameter Modelica.Units.SI.Length d_iSH=0.065;
  parameter Modelica.Units.SI.Length d_aSH=0.07;
  parameter Modelica.Units.SI.Length LSH=240;
  parameter Integer nNodesSH=8;
  parameter Modelica.Units.SI.MassFlowRate m_flowCirc=3*3.2715;
  parameter Modelica.Units.SI.MassFlowRate deltam_flowCirc=3*3.2715;
  parameter Modelica.Units.SI.Pressure p_out=3000000;
  parameter Modelica.Units.SI.Pressure deltap_out=3000000;
  parameter Modelica.Units.SI.Time deltaTimeMass=500;
  parameter Modelica.Units.SI.Mass deltaMass=500;
  parameter Modelica.Units.SI.Time startingTimeDistortion=20000;
  parameter Modelica.Units.SI.Velocity velocityDistortion=1;
  parameter Modelica.Units.SI.Length lengthDistortion=320;
  parameter Modelica.Units.SI.Time deltaTimeDistortion=20;
  parameter Integer n_Qrad_Evap=36;
  parameter Integer n_Qrad_SH=8;
  parameter Integer nEvapLoop=6;
  parameter Integer nSHLoop=3;
  parameter Modelica.Units.SI.Pressure deltap_head=200000;

//Initialization
  parameter Modelica.Units.SI.Pressure p_startEvapIn=4450000;
  parameter Modelica.Units.SI.Pressure p_startDeltaValve=200000;
  parameter Modelica.Units.SI.Pressure p_startDeltaEvapCollInTotal=150000;
  parameter Modelica.Units.SI.Pressure p_startDeltaEvap=600000;

  ClaRa.StaticCycles.Storage.Separator separator annotation (Placement(transformation(extent={{-66,-74},{-46,-54}})));
  ClaRa.StaticCycles.ValvesConnects.Tube1 TubeIn(Delta_p_nom=0.99e5, frictionAtOutlet=true)
                                                 annotation (Placement(transformation(extent={{194,-68},{174,-60}})));
  ClaRa.StaticCycles.Fittings.Mixer1 join_Y annotation (Placement(transformation(extent={{167,-69},{157,-63}})));
  ClaRa.StaticCycles.ValvesConnects.Tube1 TubeRePumpOut(Delta_p_nom=1e4, frictionAtOutlet=true)
                                                        annotation (Placement(transformation(extent={{118,-100},{138,-92}})));
  ClaRa.StaticCycles.ValvesConnects.Tube1 TubeIn1(
    Delta_p_nom=0.99e5,
    frictionAtInlet=true,
    frictionAtOutlet=true)                        annotation (Placement(transformation(extent={{150,-68},{130,-60}})));
  ClaRa.StaticCycles.ValvesConnects.Tube1 TubeRePumpIn(Delta_p_nom=2e5, frictionAtInlet=true)
                                                       annotation (Placement(transformation(extent={{-50,-100},{-30,-92}})));
  ClaRa.StaticCycles.Fittings.Split3 split_Y1(splitRatio=1/6)
                                              annotation (Placement(transformation(
        extent={{-5,-3},{5,3}},
        rotation=180,
        origin={116,-62})));
  ClaRa.StaticCycles.ValvesConnects.Valve_cutPressure1 valveBranch annotation (Placement(transformation(extent={{95,-67},{85,-61}})));
  Components.Absorber.StaticCycle.Parapolic_Tube
                                          tube(Delta_p_nom=p_startDeltaEvap,
    Delta_h_nom=h_evap,                                                      frictionAtOutlet=true)
                                               annotation (Placement(transformation(extent={{66,-68},{46,-60}})));
  ClaRa.StaticCycles.Fittings.Mixer1 join_Y1 annotation (Placement(transformation(
        extent={{-5,-3},{5,3}},
        rotation=180,
        origin={20,-62})));
  ClaRa.StaticCycles.ValvesConnects.Tube1 TubeInSep(
    Delta_p_nom=1e4,
    frictionAtInlet=true,
    frictionAtOutlet=true)                          annotation (Placement(transformation(extent={{-10,-68},{-30,-60}})));
  ClaRa.StaticCycles.ValvesConnects.Tube1 TubeEvapCollIn1(
    Delta_p_nom=deltap_head/5,
    frictionAtInlet=true,
    frictionAtOutlet=true)                                annotation (Placement(transformation(
        extent={{10,-4},{-10,4}},
        rotation=270,
        origin={130,-34})));
  ClaRa.StaticCycles.Fittings.Split3 split_Y2(splitRatio=1/5)
                                              annotation (Placement(transformation(
        extent={{-5,-3},{5,3}},
        rotation=180,
        origin={114,-14})));
  ClaRa.StaticCycles.Fittings.Split3 split_Y3(splitRatio=1/4) annotation (Placement(transformation(
        extent={{-5,-3},{5,3}},
        rotation=180,
        origin={114,32})));
  ClaRa.StaticCycles.Fittings.Split3 split_Y6(splitRatio=1/2)
                                              annotation (Placement(transformation(
        extent={{-5,-3},{5,3}},
        rotation=180,
        origin={-132,-10})));
  ClaRa.StaticCycles.Fittings.Split3 split_Y45(splitRatio=1/3) annotation (Placement(transformation(
        extent={{-5,-3},{5,3}},
        rotation=180,
        origin={114,74})));
  ClaRa.StaticCycles.Fittings.Split3 split_Y5(splitRatio=1/2)
                                              annotation (Placement(transformation(
        extent={{-5,-3},{5,3}},
        rotation=180,
        origin={114,114})));
  ClaRa.StaticCycles.ValvesConnects.Tube1 TubeEvapCollIn2(
    Delta_p_nom=deltap_head/5,
    frictionAtInlet=true,
    frictionAtOutlet=true)                                annotation (Placement(transformation(
        extent={{10,-4},{-10,4}},
        rotation=270,
        origin={128,10})));
  ClaRa.StaticCycles.ValvesConnects.Tube1 TubeEvapCollIn3(
    Delta_p_nom=deltap_head/5,
    frictionAtInlet=true,
    frictionAtOutlet=true)                                annotation (Placement(transformation(
        extent={{10,-4},{-10,4}},
        rotation=270,
        origin={130,56})));
  ClaRa.StaticCycles.ValvesConnects.Tube1 TubeEvapCollIn4(
    Delta_p_nom=deltap_head/5,
    frictionAtInlet=true,
    frictionAtOutlet=true)                                annotation (Placement(transformation(
        extent={{10,-4},{-10,4}},
        rotation=270,
        origin={130,96})));
  ClaRa.StaticCycles.ValvesConnects.Tube1 TubeEvapCollIn5(
    Delta_p_nom=deltap_head/5,
    frictionAtInlet=true,
    frictionAtOutlet=true)                                annotation (Placement(transformation(
        extent={{10,-4},{-10,4}},
        rotation=270,
        origin={130,136})));
  ClaRa.StaticCycles.Fittings.Split3 split_Y(splitRatio=1/3)
                                             annotation (Placement(transformation(
        extent={{-5,-3},{5,3}},
        rotation=180,
        origin={-118,-54})));
  ClaRa.StaticCycles.ValvesConnects.Tube1 TubeInSH(
    Delta_p_nom=2e5,
    frictionAtInlet=true,
    frictionAtOutlet=true) annotation (Placement(transformation(extent={{-86,-60},{-106,-52}})));
  ClaRa.StaticCycles.ValvesConnects.Valve_cutPressure1 valveBranch1 annotation (Placement(transformation(extent={{89,-19},{79,-13}})));
  Components.Absorber.StaticCycle.Parapolic_Tube
                                          tube1(Delta_p_nom=p_startDeltaEvap,
    Delta_h_nom=h_evap,                                                       frictionAtOutlet=true)
                                                annotation (Placement(transformation(extent={{60,-20},{40,-12}})));
  ClaRa.StaticCycles.Fittings.Mixer1 join_Y2 annotation (Placement(transformation(
        extent={{-5,3},{5,-3}},
        rotation=270,
        origin={22,-16})));
  ClaRa.StaticCycles.ValvesConnects.Tube1 TubeEvapCollOut1(
    Delta_p_nom=deltap_head/5,
    frictionAtInlet=true,
    frictionAtOutlet=true)                                 annotation (Placement(transformation(
        extent={{10,-4},{-10,4}},
        rotation=90,
        origin={20,-40})));
  ClaRa.StaticCycles.Fittings.Mixer1 join_Y3 annotation (Placement(transformation(
        extent={{-5,3},{5,-3}},
        rotation=270,
        origin={22,30})));
  Components.Absorber.StaticCycle.Parapolic_Tube
                                          tube2(Delta_p_nom=p_startDeltaEvap,
    Delta_h_nom=h_evap,                                                       frictionAtOutlet=true)
                                                annotation (Placement(transformation(extent={{60,26},{40,34}})));
  ClaRa.StaticCycles.ValvesConnects.Valve_cutPressure1 valveBranch2 annotation (Placement(transformation(extent={{89,27},{79,33}})));
  ClaRa.StaticCycles.Fittings.Mixer1 join_Y4 annotation (Placement(transformation(
        extent={{-5,3},{5,-3}},
        rotation=270,
        origin={22,72})));
  Components.Absorber.StaticCycle.Parapolic_Tube
                                          tube3(Delta_p_nom=p_startDeltaEvap,
    Delta_h_nom=h_evap,                                                       frictionAtOutlet=true)
                                                annotation (Placement(transformation(extent={{60,68},{40,76}})));
  ClaRa.StaticCycles.ValvesConnects.Valve_cutPressure1 valveBranch3 annotation (Placement(transformation(extent={{89,69},{79,75}})));
  ClaRa.StaticCycles.Fittings.Mixer1 join_Y5 annotation (Placement(transformation(
        extent={{-5,3},{5,-3}},
        rotation=270,
        origin={22,112})));
  Components.Absorber.StaticCycle.Parapolic_Tube
                                          tube4(Delta_p_nom=p_startDeltaEvap,
    Delta_h_nom=h_evap,                                                       frictionAtOutlet=true)
                                                annotation (Placement(transformation(extent={{60,108},{40,116}})));
  ClaRa.StaticCycles.ValvesConnects.Valve_cutPressure1 valveBranch4 annotation (Placement(transformation(extent={{89,109},{79,115}})));
  Components.Absorber.StaticCycle.Parapolic_Tube
                                          tube5(Delta_p_nom=p_startDeltaEvap,
    Delta_h_nom=h_evap,                                                       frictionAtOutlet=true)
                                                annotation (Placement(transformation(extent={{58,148},{38,156}})));
  ClaRa.StaticCycles.ValvesConnects.Valve_dp_nom3 valveBranch5(Delta_p_nom=p_startDeltaValve)
                                                               annotation (Placement(transformation(extent={{89,149},{79,155}})));
  ClaRa.StaticCycles.ValvesConnects.Tube1 TubeEvapCollOut2(
    Delta_p_nom=deltap_head/5,
    frictionAtInlet=true,
    frictionAtOutlet=true)                                 annotation (Placement(transformation(
        extent={{10,-4},{-10,4}},
        rotation=90,
        origin={20,6})));
  ClaRa.StaticCycles.ValvesConnects.Tube1 TubeEvapCollOut3(
    Delta_p_nom=deltap_head/5,
    frictionAtInlet=true,
    frictionAtOutlet=true)                                 annotation (Placement(transformation(
        extent={{10,-4},{-10,4}},
        rotation=90,
        origin={20,52})));
  ClaRa.StaticCycles.ValvesConnects.Tube1 TubeEvapCollOut4(
    Delta_p_nom=deltap_head/5,
    frictionAtInlet=true,
    frictionAtOutlet=true)                                 annotation (Placement(transformation(
        extent={{10,-4},{-10,4}},
        rotation=90,
        origin={20,92})));
  ClaRa.StaticCycles.ValvesConnects.Tube1 TubeEvapCollOut5(Delta_p_nom=deltap_head/5, frictionAtOutlet=true)
                                                           annotation (Placement(transformation(
        extent={{10,-4},{-10,4}},
        rotation=90,
        origin={20,134})));
  ClaRa.StaticCycles.ValvesConnects.Tube1 TubeSHCollOut1(
    Delta_p_nom=0.1e5,
    frictionAtInlet=true,
    frictionAtOutlet=true)                               annotation (Placement(transformation(
        extent={{10,-4},{-10,4}},
        rotation=90,
        origin={-206,-36})));
  ClaRa.StaticCycles.Fittings.Mixer1 join_Y6 annotation (Placement(transformation(
        extent={{-5,-3},{5,3}},
        rotation=180,
        origin={-206,-54})));
  ClaRa.StaticCycles.Fittings.Mixer1 join_Y7 annotation (Placement(transformation(
        extent={{5,-3},{-5,3}},
        rotation=90,
        origin={-204,-12})));
  ClaRa.StaticCycles.ValvesConnects.Valve_cutPressure1 valveBranchSH annotation (Placement(transformation(extent={{-145,-59},{-155,-53}})));
  Components.Absorber.StaticCycle.Parapolic_Tube
                                          tubeSH(
    Delta_p_nom=4e5,
    Delta_h_nom=h_sh,
    frictionAtInlet=true,
    frictionAtOutlet=true)                       annotation (Placement(transformation(extent={{-166,-60},{-186,-52}})));
  ClaRa.StaticCycles.ValvesConnects.Valve_cutPressure1 valveBranchSH1 annotation (Placement(transformation(extent={{-145,-15},{-155,-9}})));
  Components.Absorber.StaticCycle.Parapolic_Tube
                                          tubeSH1(
    Delta_p_nom=4e5,
    Delta_h_nom=h_sh,
    frictionAtInlet=true,
    frictionAtOutlet=true)                        annotation (Placement(transformation(extent={{-166,-16},{-186,-8}})));
  ClaRa.StaticCycles.ValvesConnects.Tube1 TubeSHCollOut2(Delta_p_nom=0.1e5, frictionAtOutlet=true)
                                                         annotation (Placement(transformation(
        extent={{10,-4},{-10,4}},
        rotation=90,
        origin={-206,14})));
  ClaRa.StaticCycles.ValvesConnects.Tube1 TubeSHCollIn1(
    Delta_p_nom=0.1e5,
    frictionAtInlet=true,
    frictionAtOutlet=true)                              annotation (Placement(transformation(
        extent={{10,-4},{-10,4}},
        rotation=270,
        origin={-118,-36})));
  ClaRa.StaticCycles.ValvesConnects.Tube1 TubeSHCollIn2(Delta_p_nom=0.1e5, frictionAtInlet=true)
                                                        annotation (Placement(transformation(
        extent={{10,-4},{-10,4}},
        rotation=270,
        origin={-132,12})));
  Components.Absorber.StaticCycle.Parapolic_Tube
                                          tubeSH2(
    Delta_p_nom=4e5,
    Delta_h_nom=h_sh,
    frictionAtInlet=true,
    frictionAtOutlet=true)                        annotation (Placement(transformation(extent={{-170,24},{-190,32}})));
  ClaRa.StaticCycles.ValvesConnects.Valve_dp_nom3 valveBranch6(Delta_p_nom=1e3)
                                                               annotation (Placement(transformation(extent={{-145,25},{-155,31}})));
  ClaRa.StaticCycles.Quadruple quadruple annotation (Placement(transformation(extent={{-198,-82},{-178,-72}})));
  ClaRa.StaticCycles.Quadruple quadruple1 annotation (Placement(transformation(extent={{-188,40},{-168,50}})));
  ClaRa.StaticCycles.Quadruple quadruple2 annotation (Placement(transformation(extent={{-184,2},{-164,12}})));
  ClaRa.StaticCycles.Quadruple quadruple3 annotation (Placement(transformation(extent={{-188,-40},{-168,-30}})));
  ClaRa.StaticCycles.Quadruple quadruple4 annotation (Placement(transformation(extent={{-110,4},{-90,14}})));
  ClaRa.StaticCycles.Quadruple quadruple5 annotation (Placement(transformation(extent={{-132,32},{-112,42}})));
  ClaRa.StaticCycles.Quadruple quadruple6 annotation (Placement(transformation(extent={{-144,-36},{-124,-26}})));
  ClaRa.StaticCycles.Quadruple quadruple7 annotation (Placement(transformation(extent={{-60,-34},{-40,-24}})));
  ClaRa.StaticCycles.Quadruple quadruple8 annotation (Placement(transformation(extent={{-32,-56},{-12,-46}})));
  ClaRa.StaticCycles.Quadruple quadruple9 annotation (Placement(transformation(extent={{-40,-86},{-20,-76}})));
  ClaRa.StaticCycles.Quadruple quadruple10 annotation (Placement(transformation(extent={{-22,-122},{-2,-112}})));
  ClaRa.StaticCycles.Quadruple quadruple11 annotation (Placement(transformation(extent={{118,-122},{138,-112}})));
  ClaRa.StaticCycles.Quadruple quadruple13 annotation (Placement(transformation(extent={{158,-52},{178,-42}})));
  ClaRa.StaticCycles.Quadruple quadruple14 annotation (Placement(transformation(extent={{100,-48},{120,-38}})));
  ClaRa.StaticCycles.Quadruple quadruple15 annotation (Placement(transformation(extent={{72,-50},{92,-40}})));
  ClaRa.StaticCycles.Quadruple quadruple16 annotation (Placement(transformation(extent={{38,-48},{58,-38}})));
  ClaRa.StaticCycles.Quadruple quadruple17 annotation (Placement(transformation(extent={{-6,-56},{14,-46}})));
  ClaRa.StaticCycles.Quadruple quadruple18 annotation (Placement(transformation(extent={{98,-2},{118,8}})));
  ClaRa.StaticCycles.Quadruple quadruple19 annotation (Placement(transformation(extent={{66,-2},{86,8}})));
  ClaRa.StaticCycles.Quadruple quadruple20 annotation (Placement(transformation(extent={{32,-2},{52,8}})));
  ClaRa.StaticCycles.Quadruple quadruple21 annotation (Placement(transformation(extent={{34,40},{54,50}})));
  ClaRa.StaticCycles.Quadruple quadruple22 annotation (Placement(transformation(extent={{72,40},{92,50}})));
  ClaRa.StaticCycles.Quadruple quadruple23 annotation (Placement(transformation(extent={{102,40},{122,50}})));
  ClaRa.StaticCycles.Quadruple quadruple24 annotation (Placement(transformation(extent={{36,122},{56,132}})));
  ClaRa.StaticCycles.Quadruple quadruple25 annotation (Placement(transformation(extent={{70,120},{90,130}})));
  ClaRa.StaticCycles.Quadruple quadruple26 annotation (Placement(transformation(extent={{104,120},{124,130}})));
  ClaRa.StaticCycles.Quadruple quadruple27 annotation (Placement(transformation(extent={{34,162},{54,172}})));
  ClaRa.StaticCycles.Quadruple quadruple28 annotation (Placement(transformation(extent={{72,162},{92,172}})));
  ClaRa.StaticCycles.Quadruple quadruple29 annotation (Placement(transformation(extent={{118,162},{138,172}})));
  ClaRa.StaticCycles.Quadruple quadruple30 annotation (Placement(transformation(extent={{40,84},{60,94}})));
  ClaRa.StaticCycles.Quadruple quadruple31 annotation (Placement(transformation(extent={{68,84},{88,94}})));
  ClaRa.StaticCycles.Quadruple quadruple32 annotation (Placement(transformation(extent={{98,82},{118,92}})));
  ClaRa.StaticCycles.Machines.Pump1 pumpRe annotation (Placement(transformation(extent={{36,-106},{56,-86}})));
  ClaRa.StaticCycles.ValvesConnects.PressureAnchor_constFlow1 pressureAnchor_constFlow1_1(p_nom=34.41e5) annotation (Placement(transformation(extent={{-5,-99},{5,-93}})));
  ClaRa.StaticCycles.Quadruple quadruple33 annotation (Placement(transformation(extent={{18,-122},{38,-112}})));
  ClaRa.StaticCycles.ValvesConnects.Tube1 TubeMainValveIn(Delta_p_nom=1e4, frictionAtOutlet=false) annotation (Placement(transformation(extent={{254,-68},{234,-60}})));
  ClaRa.StaticCycles.ValvesConnects.Tube2 TubeMainPumpIn(Delta_p_nom=1e4,
    frictionAtInlet=true,                                                 frictionAtOutlet=false) annotation (Placement(transformation(
        extent={{10,-4},{-10,4}},
        rotation=180,
        origin={194,-172})));
  ClaRa.StaticCycles.ValvesConnects.Valve_dp_nom3 valveMain(Delta_p_nom=20000) annotation (Placement(transformation(extent={{217,-67},{207,-61}})));
  ClaRa.StaticCycles.Machines.Pump1 pumpMain annotation (Placement(transformation(
        extent={{10,-10},{-10,10}},
        rotation=270,
        origin={264,-144})));
  ClaRa.StaticCycles.Quadruple quadruple12 annotation (Placement(transformation(extent={{200,-50},{220,-40}})));
  ClaRa.StaticCycles.Quadruple quadruple34 annotation (Placement(transformation(extent={{228,-50},{248,-40}})));
  ClaRa.StaticCycles.Quadruple quadruple35 annotation (Placement(transformation(extent={{272,-48},{292,-38}})));
  ClaRa.StaticCycles.Quadruple quadruple36 annotation (Placement(transformation(extent={{250,-190},{270,-180}})));
  ClaRa.StaticCycles.Quadruple quadruple37 annotation (Placement(transformation(extent={{180,-196},{200,-186}})));
  ClaRa.StaticCycles.HeatExchanger.Condenser condenser annotation (Placement(transformation(extent={{-180,-160},{-160,-140}})));
  ClaRa.StaticCycles.Machines.Turbine turbine annotation (Placement(transformation(extent={{-216,-128},{-204,-108}})));
  ClaRa.StaticCycles.ValvesConnects.PressureAnchor_constFlow1 pressureAnchor_constFlow1_2(p_nom=p_out) annotation (Placement(transformation(
        extent={{-5,-2},{5,2}},
        rotation=270,
        origin={-256,-68})));
  ClaRa.StaticCycles.Storage.Feedwatertank3 feedwatertank(p_FWT_nom=3.7e5, m_flow_nom=6) annotation (Placement(transformation(extent={{-26,-176},{-6,-164}})));
  ClaRa.StaticCycles.Machines.Pump1 Pump_cond annotation (Placement(transformation(extent={{-120,-182},{-100,-162}})));
  ClaRa.StaticCycles.Quadruple quadruple38 annotation (Placement(transformation(extent={{-70,-192},{-50,-182}})));
  ClaRa.StaticCycles.Quadruple quadruple39 annotation (Placement(transformation(extent={{-152,-192},{-132,-182}})));
  ClaRa.StaticCycles.Quadruple quadruple40 annotation (Placement(transformation(extent={{-214,-150},{-194,-140}})));
  ClaRa.StaticCycles.Quadruple quadruple41 annotation (Placement(transformation(extent={{-250,-138},{-230,-128}})));
  ClaRa.StaticCycles.Fittings.Split5 split_turbine annotation (Placement(transformation(extent={{-249,-103},{-239,-97}})));
  ClaRa.StaticCycles.ValvesConnects.Valve_cutPressure2 valve_SteamFWT annotation (Placement(transformation(extent={{-155,-100},{-145,-96}})));
  ClaRa.StaticCycles.Quadruple quadruple42 annotation (Placement(transformation(extent={{-178,-114},{-158,-104}})));
  ClaRa.StaticCycles.Quadruple quadruple43 annotation (Placement(transformation(extent={{-64,-132},{-44,-122}})));
  ClaRa.StaticCycles.Quadruple quadruple44 annotation (Placement(transformation(extent={{-252,-88},{-232,-78}})));
equation
  connect(TubeRePumpOut.outlet, join_Y.inlet_2) annotation (Line(points={{138.5,-96},{162,-96},{162,-69.5}}, color={0,131,169}));
  connect(join_Y.inlet_1, TubeIn.outlet) annotation (Line(points={{167.5,-64},{173.5,-64}}, color={0,131,169}));
  connect(TubeIn1.inlet, join_Y.outlet) annotation (Line(points={{150.5,-64},{156.5,-64}}, color={0,131,169}));
  connect(separator.outlet1, TubeRePumpIn.inlet) annotation (Line(points={{-56,-74},{-56,-96},{-50.5,-96}}, color={0,131,169}));
  connect(split_Y1.inlet, TubeIn1.outlet) annotation (Line(points={{121.5,-64},{129.5,-64}}, color={0,131,169}));
  connect(valveBranch.inlet, split_Y1.outlet_1) annotation (Line(points={{95.5,-64},{110.5,-64}}, color={0,131,169}));
  connect(tube.inlet, valveBranch.outlet) annotation (Line(points={{66.5,-64},{84.5,-64}}, color={0,131,169}));
  connect(join_Y1.inlet_1, tube.outlet) annotation (Line(points={{25.5,-64},{45.5,-64}}, color={0,131,169}));
  connect(separator.inlet, TubeInSep.outlet) annotation (Line(points={{-46,-64},{-30.5,-64}}, color={0,131,169}));
  connect(TubeInSep.inlet, join_Y1.outlet) annotation (Line(points={{-9.5,-64},{14.5,-64}}, color={0,131,169}));
  connect(split_Y1.outlet_2, TubeEvapCollIn1.inlet) annotation (Line(points={{116,-58.5},{116,-52},{130,-52},{130,-44.5}}, color={0,131,169}));
  connect(TubeEvapCollIn1.outlet,split_Y2. inlet) annotation (Line(points={{130,-23.5},{130,-16},{119.5,-16}}, color={0,131,169}));
  connect(split_Y2.outlet_2, TubeEvapCollIn2.inlet) annotation (Line(points={{114,-10.5},{114,-8},{128,-8},{128,-0.5}}, color={0,131,169}));
  connect(TubeEvapCollIn2.outlet, split_Y3.inlet) annotation (Line(points={{128,20.5},{128,30},{119.5,30}}, color={0,131,169}));
  connect(TubeEvapCollIn3.inlet, split_Y3.outlet_2) annotation (Line(points={{130,45.5},{130,40},{114,40},{114,35.5}}, color={0,131,169}));
  connect(TubeEvapCollIn3.outlet, split_Y45.inlet) annotation (Line(points={{130,66.5},{130,72},{119.5,72}}, color={0,131,169}));
  connect(split_Y45.outlet_2, TubeEvapCollIn4.inlet) annotation (Line(points={{114,77.5},{114,80},{130,80},{130,85.5}}, color={0,131,169}));
  connect(TubeEvapCollIn4.outlet,split_Y5. inlet) annotation (Line(points={{130,106.5},{130,112},{119.5,112}}, color={0,131,169}));
  connect(split_Y5.outlet_2, TubeEvapCollIn5.inlet) annotation (Line(points={{114,117.5},{114,120},{130,120},{130,125.5}}, color={0,131,169}));
  connect(TubeInSH.inlet, separator.outlet2) annotation (Line(points={{-85.5,-56},{-72,-56},{-72,-48},{-56,-48},{-56,-54}}, color={0,131,169}));
  connect(split_Y.inlet, TubeInSH.outlet) annotation (Line(points={{-112.5,-56},{-106.5,-56}}, color={0,131,169}));
  connect(tube1.inlet, valveBranch1.outlet) annotation (Line(points={{60.5,-16},{78.5,-16}}, color={0,131,169}));
  connect(valveBranch1.inlet,split_Y2. outlet_1) annotation (Line(points={{89.5,-16},{108.5,-16}}, color={0,131,169}));
  connect(join_Y1.inlet_2, TubeEvapCollOut1.outlet) annotation (Line(points={{20,-58.5},{20,-50.5}}, color={0,131,169}));
  connect(TubeEvapCollOut1.inlet, join_Y2.outlet) annotation (Line(points={{20,-29.5},{20,-21.5}}, color={0,131,169}));
  connect(join_Y2.inlet_2, tube1.outlet) annotation (Line(points={{25.5,-16},{39.5,-16}}, color={0,131,169}));
  connect(join_Y3.inlet_2, tube2.outlet) annotation (Line(points={{25.5,30},{39.5,30}}, color={0,131,169}));
  connect(tube2.inlet, valveBranch2.outlet) annotation (Line(points={{60.5,30},{78.5,30}}, color={0,131,169}));
  connect(join_Y4.inlet_2, tube3.outlet) annotation (Line(points={{25.5,72},{39.5,72}}, color={0,131,169}));
  connect(tube3.inlet, valveBranch3.outlet) annotation (Line(points={{60.5,72},{78.5,72}}, color={0,131,169}));
  connect(join_Y5.inlet_2, tube4.outlet) annotation (Line(points={{25.5,112},{39.5,112}}, color={0,131,169}));
  connect(tube4.inlet, valveBranch4.outlet) annotation (Line(points={{60.5,112},{78.5,112}}, color={0,131,169}));
  connect(tube5.inlet, valveBranch5.outlet) annotation (Line(points={{58.5,152},{78.5,152}}, color={0,131,169}));
  connect(valveBranch4.inlet,split_Y5. outlet_1) annotation (Line(points={{89.5,112},{108.5,112}}, color={0,131,169}));
  connect(valveBranch3.inlet, split_Y45.outlet_1) annotation (Line(points={{89.5,72},{108.5,72}}, color={0,131,169}));
  connect(valveBranch2.inlet, split_Y3.outlet_1) annotation (Line(points={{89.5,30},{108.5,30}}, color={0,131,169}));
  connect(valveBranch5.inlet, TubeEvapCollIn5.outlet) annotation (Line(points={{89.5,152},{130,152},{130,146.5}}, color={0,131,169}));
  connect(tube5.outlet, TubeEvapCollOut5.inlet) annotation (Line(points={{37.5,152},{20,152},{20,144.5}}, color={0,131,169}));
  connect(TubeEvapCollOut5.outlet, join_Y5.inlet_1) annotation (Line(points={{20,123.5},{20,117.5}}, color={0,131,169}));
  connect(join_Y5.outlet, TubeEvapCollOut4.inlet) annotation (Line(points={{20,106.5},{20,102.5}}, color={0,131,169}));
  connect(TubeEvapCollOut4.outlet, join_Y4.inlet_1) annotation (Line(points={{20,81.5},{20,77.5}}, color={0,131,169}));
  connect(join_Y4.outlet, TubeEvapCollOut3.inlet) annotation (Line(points={{20,66.5},{20,62.5}}, color={0,131,169}));
  connect(TubeEvapCollOut3.outlet,join_Y3. inlet_1) annotation (Line(points={{20,41.5},{20,35.5}}, color={0,131,169}));
  connect(join_Y3.outlet, TubeEvapCollOut2.inlet) annotation (Line(points={{20,24.5},{20,16.5}}, color={0,131,169}));
  connect(TubeEvapCollOut2.outlet, join_Y2.inlet_1) annotation (Line(points={{20,-4.5},{20,-10.5}}, color={0,131,169}));
  connect(tubeSH.inlet, valveBranchSH.outlet) annotation (Line(points={{-165.5,-56},{-155.5,-56}}, color={0,131,169}));
  connect(join_Y6.inlet_1, tubeSH.outlet) annotation (Line(points={{-200.5,-56},{-186.5,-56}}, color={0,131,169}));
  connect(valveBranchSH.inlet, split_Y.outlet_1) annotation (Line(points={{-144.5,-56},{-123.5,-56}}, color={0,131,169}));
  connect(tubeSH1.inlet, valveBranchSH1.outlet) annotation (Line(points={{-165.5,-12},{-155.5,-12}}, color={0,131,169}));
  connect(valveBranchSH1.inlet,split_Y6. outlet_1) annotation (Line(points={{-144.5,-12},{-137.5,-12}}, color={0,131,169}));
  connect(TubeSHCollOut1.inlet,join_Y7. outlet) annotation (Line(points={{-206,-25.5},{-206,-17.5}}, color={0,131,169}));
  connect(join_Y7.inlet_2, tubeSH1.outlet) annotation (Line(points={{-200.5,-12},{-186.5,-12}}, color={0,131,169}));
  connect(TubeSHCollOut1.outlet,join_Y6. inlet_2) annotation (Line(points={{-206,-46.5},{-206,-50.5}}, color={0,131,169}));
  connect(split_Y.outlet_2, TubeSHCollIn1.inlet) annotation (Line(points={{-118,-50.5},{-118,-46.5}}, color={0,131,169}));
  connect(split_Y6.inlet, TubeSHCollIn1.outlet) annotation (Line(points={{-126.5,-12},{-118,-12},{-118,-25.5}}, color={0,131,169}));
  connect(split_Y6.outlet_2, TubeSHCollIn2.inlet) annotation (Line(points={{-132,-6.5},{-132,1.5}}, color={0,131,169}));
  connect(join_Y7.inlet_1, TubeSHCollOut2.outlet) annotation (Line(points={{-206,-6.5},{-206,3.5}}, color={0,131,169}));
  connect(TubeSHCollOut2.inlet, tubeSH2.outlet) annotation (Line(points={{-206,24.5},{-206,28},{-190.5,28}}, color={0,131,169}));
  connect(tubeSH2.inlet, valveBranch6.outlet) annotation (Line(points={{-169.5,28},{-155.5,28}}, color={0,131,169}));
  connect(valveBranch6.inlet, TubeSHCollIn2.outlet) annotation (Line(points={{-144.5,28},{-132,28},{-132,22.5}}, color={0,131,169}));
  connect(quadruple.steamSignal, join_Y6.outlet) annotation (Line(points={{-198,-76.8},{-214,-76.8},{-214,-56},{-211.5,-56}}, color={0,131,169}));
  connect(quadruple3.steamSignal, tubeSH.outlet) annotation (Line(points={{-188,-34.8},{-192,-34.8},{-192,-56},{-186.5,-56}}, color={0,131,169}));
  connect(quadruple2.steamSignal, tubeSH1.outlet) annotation (Line(points={{-184,7.2},{-190,7.2},{-190,-12},{-186.5,-12}}, color={0,131,169}));
  connect(quadruple1.steamSignal, tubeSH2.outlet) annotation (Line(points={{-188,45.2},{-188,28},{-190.5,28}}, color={0,131,169}));
  connect(quadruple6.steamSignal, split_Y.outlet_1) annotation (Line(points={{-144,-30.8},{-140,-30.8},{-140,-56},{-123.5,-56}}, color={0,131,169}));
  connect(quadruple4.steamSignal, split_Y6.outlet_1) annotation (Line(points={{-110,9.2},{-126,9.2},{-126,-12},{-137.5,-12}}, color={0,131,169}));
  connect(quadruple5.steamSignal, TubeSHCollIn2.outlet) annotation (Line(points={{-132,37.2},{-136,37.2},{-136,28},{-132,28},{-132,22.5}}, color={0,131,169}));
  connect(quadruple7.steamSignal, separator.outlet2) annotation (Line(points={{-60,-28.8},{-64,-28.8},{-64,-48},{-56,-48},{-56,-54}}, color={0,131,169}));
  connect(quadruple9.steamSignal, TubeRePumpIn.inlet) annotation (Line(points={{-40,-80.8},{-56,-80.8},{-56,-96},{-50.5,-96}}, color={0,131,169}));
  connect(quadruple8.steamSignal, TubeInSep.outlet) annotation (Line(points={{-32,-50.8},{-36,-50.8},{-36,-64},{-30.5,-64}}, color={0,131,169}));
  connect(quadruple11.steamSignal, TubeRePumpOut.inlet) annotation (Line(points={{118,-116.8},{114,-116.8},{114,-96},{117.5,-96}}, color={0,131,169}));
  connect(quadruple13.steamSignal, join_Y.outlet) annotation (Line(points={{158,-46.8},{156,-46.8},{156,-64},{156.5,-64}}, color={0,131,169}));
  connect(quadruple17.steamSignal, join_Y1.outlet) annotation (Line(points={{-6,-50.8},{-6,-64},{14.5,-64}}, color={0,131,169}));
  connect(quadruple16.steamSignal, tube.outlet) annotation (Line(points={{38,-42.8},{38,-64},{45.5,-64}}, color={0,131,169}));
  connect(quadruple15.steamSignal, valveBranch.outlet) annotation (Line(points={{72,-44.8},{70,-44.8},{70,-64},{84.5,-64}}, color={0,131,169}));
  connect(quadruple14.steamSignal, split_Y1.outlet_1) annotation (Line(points={{100,-42.8},{100,-64},{110.5,-64}}, color={0,131,169}));
  connect(quadruple29.steamSignal, TubeEvapCollIn5.outlet) annotation (Line(points={{118,167.2},{114,167.2},{114,152},{130,152},{130,146.5}}, color={0,131,169}));
  connect(quadruple28.steamSignal, valveBranch5.outlet) annotation (Line(points={{72,167.2},{70,167.2},{70,152},{78.5,152}}, color={0,131,169}));
  connect(quadruple27.steamSignal, TubeEvapCollOut5.inlet) annotation (Line(points={{34,167.2},{30,167.2},{30,152},{20,152},{20,144.5}}, color={0,131,169}));
  connect(quadruple24.steamSignal, tube4.outlet) annotation (Line(points={{36,127.2},{36,112},{39.5,112}}, color={0,131,169}));
  connect(quadruple25.steamSignal, valveBranch4.outlet) annotation (Line(points={{70,125.2},{70,112},{78.5,112}}, color={0,131,169}));
  connect(quadruple26.steamSignal, split_Y5.outlet_1) annotation (Line(points={{104,125.2},{102,125.2},{102,112},{108.5,112}}, color={0,131,169}));
  connect(quadruple30.steamSignal, tube3.outlet) annotation (Line(points={{40,89.2},{36,89.2},{36,72},{39.5,72}}, color={0,131,169}));
  connect(quadruple31.steamSignal, valveBranch3.outlet) annotation (Line(points={{68,89.2},{68,72},{78.5,72}}, color={0,131,169}));
  connect(quadruple32.steamSignal, split_Y45.outlet_1) annotation (Line(points={{98,87.2},{98,72},{108.5,72}}, color={0,131,169}));
  connect(quadruple21.steamSignal, tube2.outlet) annotation (Line(points={{34,45.2},{34,30},{39.5,30}}, color={0,131,169}));
  connect(quadruple22.steamSignal, valveBranch2.outlet) annotation (Line(points={{72,45.2},{72,30},{78.5,30}}, color={0,131,169}));
  connect(quadruple23.steamSignal, split_Y3.outlet_1) annotation (Line(points={{102,45.2},{102,30},{108.5,30}}, color={0,131,169}));
  connect(quadruple20.steamSignal, tube1.outlet) annotation (Line(points={{32,3.2},{34,3.2},{34,-16},{39.5,-16}}, color={0,131,169}));
  connect(quadruple19.steamSignal, valveBranch1.outlet) annotation (Line(points={{66,3.2},{66,-16},{78.5,-16}}, color={0,131,169}));
  connect(quadruple18.steamSignal, split_Y2.outlet_1) annotation (Line(points={{98,3.2},{98,-16},{108.5,-16}}, color={0,131,169}));
  connect(pumpRe.outlet, TubeRePumpOut.inlet) annotation (Line(points={{56.5,-96},{117.5,-96}}, color={0,131,169}));
  connect(TubeRePumpIn.outlet, pressureAnchor_constFlow1_1.inlet) annotation (Line(points={{-29.5,-96},{-5.5,-96}}, color={0,131,169}));
  connect(pressureAnchor_constFlow1_1.outlet, pumpRe.inlet) annotation (Line(points={{5.5,-96},{35.5,-96}}, color={0,131,169}));
  connect(quadruple10.steamSignal, pressureAnchor_constFlow1_1.inlet) annotation (Line(points={{-22,-116.8},{-26,-116.8},{-26,-96},{-5.5,-96}}, color={0,131,169}));
  connect(quadruple33.steamSignal, pumpRe.inlet) annotation (Line(points={{18,-116.8},{16,-116.8},{16,-96},{35.5,-96}}, color={0,131,169}));
  connect(valveMain.outlet, TubeIn.inlet) annotation (Line(points={{206.5,-64},{194.5,-64}}, color={0,131,169}));
  connect(valveMain.inlet, TubeMainValveIn.outlet) annotation (Line(points={{217.5,-64},{233.5,-64}}, color={0,131,169}));
  connect(TubeMainValveIn.inlet, pumpMain.outlet) annotation (Line(points={{254.5,-64},{264,-64},{264,-133.5}},color={0,131,169}));
  connect(TubeMainPumpIn.outlet, pumpMain.inlet) annotation (Line(points={{204.5,-172},{264,-172},{264,-154.5}},
                                                                                                          color={0,131,169}));
  connect(quadruple12.steamSignal, TubeIn.inlet) annotation (Line(points={{200,-44.8},{200,-64},{194.5,-64}}, color={0,131,169}));
  connect(quadruple34.steamSignal, TubeMainValveIn.outlet) annotation (Line(points={{228,-44.8},{228,-64},{233.5,-64}}, color={0,131,169}));
  connect(quadruple35.steamSignal, pumpMain.outlet) annotation (Line(points={{272,-42.8},{264,-42.8},{264,-133.5}},        color={0,131,169}));
  connect(quadruple36.steamSignal, pumpMain.inlet) annotation (Line(points={{250,-184.8},{244,-184.8},{244,-172},{264,-172},{264,-154.5}},
                                                                                                                                color={0,131,169}));
  connect(quadruple37.steamSignal, TubeMainPumpIn.inlet) annotation (Line(points={{180,-190.8},{174,-190.8},{174,-172},{183.5,-172}},
                                                                                                                   color={0,131,169}));
  connect(turbine.outlet, condenser.inlet) annotation (Line(points={{-203.5,-126},{-203.5,-134},{-170,-134},{-170,-139.5}}, color={0,131,169}));
  connect(pressureAnchor_constFlow1_2.inlet, join_Y6.outlet) annotation (Line(points={{-256,-62.5},{-256,-56},{-211.5,-56}}, color={0,131,169}));
  connect(feedwatertank.cond_out, TubeMainPumpIn.inlet) annotation (Line(points={{-5.5,-172},{183.5,-172}}, color={0,131,169}));
  connect(condenser.outlet, Pump_cond.inlet) annotation (Line(points={{-170,-160.5},{-170,-172},{-120.5,-172}}, color={0,131,169}));
  connect(Pump_cond.outlet, feedwatertank.cond_in) annotation (Line(points={{-99.5,-172},{-26.5,-172}}, color={0,131,169}));
  connect(quadruple41.steamSignal, turbine.inlet) annotation (Line(points={{-250,-132.8},{-220,-132.8},{-220,-114},{-216.5,-114}}, color={0,131,169}));
  connect(quadruple40.steamSignal, condenser.inlet) annotation (Line(points={{-214,-144.8},{-190,-144.8},{-190,-134},{-170,-134},{-170,-139.5}}, color={0,131,169}));
  connect(quadruple39.steamSignal, Pump_cond.inlet) annotation (Line(points={{-152,-186.8},{-158,-186.8},{-158,-172},{-120.5,-172}}, color={0,131,169}));
  connect(quadruple38.steamSignal, feedwatertank.cond_in) annotation (Line(points={{-70,-186.8},{-76,-186.8},{-76,-172},{-26.5,-172}}, color={0,131,169}));
  connect(pressureAnchor_constFlow1_2.outlet, split_turbine.inlet) annotation (Line(points={{-256,-73.5},{-256,-98},{-249.5,-98}}, color={0,131,169}));
  connect(split_turbine.outlet_2, turbine.inlet) annotation (Line(points={{-244,-103.5},{-243.875,-103.5},{-243.875,-114},{-216.5,-114}}, color={0,131,169}));
  connect(split_turbine.outlet_1, valve_SteamFWT.inlet) annotation (Line(points={{-238.5,-98},{-155.5,-98}}, color={0,131,169}));
  connect(valve_SteamFWT.outlet, feedwatertank.tap_in) annotation (Line(points={{-144.5,-98},{-112,-98},{-112,-136},{-16,-136},{-16,-163.5}}, color={0,131,169}));
  connect(quadruple42.steamSignal, valve_SteamFWT.inlet) annotation (Line(points={{-178,-108.8},{-188,-108.8},{-188,-98},{-155.5,-98}}, color={0,131,169}));
  connect(quadruple43.steamSignal, feedwatertank.tap_in) annotation (Line(points={{-64,-126.8},{-68,-126.8},{-68,-136},{-16,-136},{-16,-163.5}}, color={0,131,169}));
  connect(quadruple44.steamSignal, split_turbine.inlet) annotation (Line(points={{-252,-82.8},{-252,-84},{-256,-84},{-256,-98},{-249.5,-98}}, color={0,131,169}));
  annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-300,-200},{300,200}})),
              Diagram(coordinateSystem(extent={{-300,-200},{300,200}})));
end Static_Cycle_SixBranchEvapThreeBranchSH_Turbine;
