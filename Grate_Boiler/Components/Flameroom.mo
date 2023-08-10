within Grate_Boiler.Components;
package Flameroom "Model of the complete Flameroom"
 extends ClaRa.Basics.Icons.PackageIcons.Subsystems80;
  model Grate_Boiler_Example "Comparison between flameroom model and grate boiler of Nielsen et al."
    extends ClaRa.Basics.Icons.PackageIcons.ExecutableExample100;
    inner ClaRa.SimCenter simCenter(
      redeclare Media.FlueGasTILMedia_GrateBoiler flueGasModel,
      useHomotopy=true,
      redeclare Media.Waste_NielsenEtAl fuelModel1)                                                                                                                                   annotation (Placement(transformation(extent={{100,-140},{140,-120}})));
    Grate.Grate grate(
      T_start=(273.15 + 25)*ones(grate.N_cv_pipe),
      Grate_L=bed_Model.length,
      Grate_W=bed_Model.width,
      Grate_H=0.1,
      elements=bed_Model.n_segments) annotation (Placement(transformation(extent={{-30,-92},{30,-82}})));
    Radiation.RadiationModel_Nielsen_pseudoStates radiation(
      height=16,
      m=2,
      A_down_scaling=bed_Model.length/radiation.length,
      k=1,
      length=6,
      width=4,
      n=bed_Model.n_segments,
      timeConst_ps=10)        annotation (Placement(transformation(extent={{49,50},{71,70}})));
    Gasphase.Gasphase_Radiation gasphase2(
      m_flow_nom=10,
      redeclare model PressureLoss = Bedsegment.Pressure_Loss.LinearPressureLoss_L2_Grate_GasPhase (Delta_p_nom=1000),
      T_start=273.15 + 260,
      T_start_flueGas_out=273.15 + 260,
      xi_start_flueGas_out={0.0,0,0.0,0,0.78,0.21,0,0.0,0},
      length=radiation.length,
      width=radiation.width,
      height=radiation.height/2,
      redeclare model GasCombustion = Grate_Boiler.Components.Gasphase.Gascombustion.GasCombustion (CF_T=6)) annotation (Placement(transformation(extent={{-36,76},{4,96}})));
    FlameRoomWalls.FlameroomWall_3Layers flameroomWall_3L_front_back(
      elements=1,
      d_ceramics=0.03,
      d_concrete=0.03,
      d_steel=0.005,
      T_start=(533.15)*ones(flameroomWall_3L_front_back.elements),
      length=radiation.width + radiation.length,
      width=radiation.height) annotation (Placement(transformation(extent={{86,24},{92,56}})));
    FlameRoomWalls.FlameroomWall_2Layers flameroomWall_2L(
      elements=radiation.k,
      d_ceramics=0.002,
      d_steel=0.005,
      T_start=(533.15)*ones(flameroomWall_2L.elements),
      length=radiation.length,
      width=radiation.width,
      redeclare model ceramic = ClaRa.Basics.ControlVolumes.SolidVolumes.ThinPlateWall_L4 (
          redeclare model Material = Grate_Boiler.Components.FlameRoomWalls.Media_Inconel,
          T_start=(273.15 + 250)*ones(flameroomWall_2L.elements),
          stateLocation=2),
      redeclare model steel = ClaRa.Basics.ControlVolumes.SolidVolumes.ThinPlateWall_L4 (redeclare model Material = Grate_Boiler.Components.FlameRoomWalls.Media_Steel_lambda25, stateLocation=2)) annotation (Placement(transformation(
          extent={{-2.75,-19.5},{2.75,19.5}},
          rotation=90,
          origin={60.5,100})));
    ClaRa.Components.BoundaryConditions.BoundaryFuel_Txim_flow boundaryFuel_Txim_flow(
      variable_m_flow=true,
      variable_xi=false,
      m_flow_const=4.2,
      T_const(displayUnit="degC") = 298.15,
      xi_const={0.3,0.04,0.2,0.06,0.03,0.27})                                             annotation (Placement(transformation(extent={{-80,-66},{-60,-46}})));
    ClaRa.Components.BoundaryConditions.BoundaryFuel_pTxi boundaryFuel_pTxi annotation (Placement(transformation(extent={{112,-66},{92,-46}})));
    ClaRa.Components.BoundaryConditions.BoundaryGas_Txim_flow boundaryGas_Txim_flow(
      variable_m_flow=true,
      T_const=273.15 + 150,
      m_flow_const=5,
      xi_const={0.0,0,0.0,0,0.78,0.21,0,0.0,0})                                         annotation (Placement(transformation(
          extent={{-10,-10},{10,10}},
          rotation=0,
          origin={-70,-120})));
    ClaRa.Components.BoundaryConditions.BoundaryGas_pTxi boundaryGas_pTxi annotation (Placement(transformation(extent={{10,-10},{-10,10}},
          rotation=90,
          origin={-16,130})));
    ClaRa.Components.BoundaryConditions.BoundaryGas_Txim_flow boundaryGas_Txim_flow2(
      variable_m_flow=true,
      m_flow_const=5,
      T_const=273.15 + 150,
      xi_const={0.0,0,0.0,0,0.78,0.21,0,0.0,0})                                          annotation (Placement(transformation(extent={{-82,76},{-62,96}})));
    Gasphase.Gasphase_Radiation gasphase1(
      m_flow_nom=10,
      redeclare model PressureLoss = Bedsegment.Pressure_Loss.LinearPressureLoss_L2_Grate_GasPhase (Delta_p_nom=1000),
      T_start=260 + 273.15,
      T_start_flueGas_out=260 + 273.15,
      length=radiation.length,
      width=radiation.width,
      height=radiation.height/2,
      xi_start_flueGas_out={0.0,0,0.0,0,0.78,0.21,0,0.0,0},
      redeclare model GasCombustion = Gasphase.Gascombustion.GasCombustion) annotation (Placement(transformation(extent={{-36,32},{4,52}})));
    ClaRa.Components.BoundaryConditions.BoundaryGas_Txim_flow boundaryGas_Txim_flow1(
      T_const=273 + 600,
      variable_m_flow=false,
      xi_const={0.0,0,0.0,0,0.78,0.21,0,0.0,0},
      m_flow_const=0*10)                                                                 annotation (Placement(transformation(extent={{-80,32},{-60,52}})));
    Bed.Bed_Model bed_Model(
      length=8,
      width=4,
      m_segment_fuel_start={1100,990,880,770,660,550,440,330,220,110},
      T_fuel_start=273.15 + 25,
      T_start_flueGas_out=273.15 + 25,
      n_segments=10,
      xi_fuel_start={0.357,0.0429,0.307,0.0,0.0,0.00714},
      height=0.1) annotation (Placement(transformation(extent={{-30,-66},{30,-46}})));
    ClaRa.Components.VolumesValvesFittings.Fittings.JoinGas_L2_flex joinGas_L2_flex(N_ports_in=bed_Model.n_segments, T_start=1000)
                                                                                                                  annotation (Placement(transformation(
          extent={{6,-6},{-6,6}},
          rotation=-90,
          origin={-16,20})));
    Modelica.Thermal.HeatTransfer.Sources.FixedTemperature fixedTemperature1(T(displayUnit="K") = 533.15)
                                                                                       annotation (Placement(transformation(
          extent={{-10,-10},{10,10}},
          rotation=180,
          origin={130,130})));
    ClaRa.Components.VolumesValvesFittings.Valves.GenericValveGas_L1 adapter[bed_Model.n_segments](redeclare model PressureLoss = ClaRa.Components.VolumesValvesFittings.Valves.Fundamentals.LinearNominalPoint (Delta_p_nom=1000, m_flow_nom=10))
                                                                                                                                                                                                    annotation (Placement(transformation(
          extent={{-8,-4},{8,4}},
          rotation=90,
          origin={-16,-4})));
    ClaRa.Components.VolumesValvesFittings.Valves.GenericValveGas_L1 genericValveGas_L1_2(redeclare model PressureLoss = ClaRa.Components.VolumesValvesFittings.Valves.Fundamentals.LinearNominalPoint (Delta_p_nom=1000, m_flow_nom=10))
                                                                                                                                                                                                    annotation (Placement(transformation(
          extent={{-8,-4},{8,4}},
          rotation=90,
          origin={-16,64})));
    ClaRa.Components.VolumesValvesFittings.Valves.GenericValveGas_L1 genericValveGas_L1_3(redeclare model PressureLoss = ClaRa.Components.VolumesValvesFittings.Valves.Fundamentals.LinearNominalPoint (Delta_p_nom=1000, m_flow_nom=10))
                                                                                                                                                                                                    annotation (Placement(transformation(
          extent={{-8,-4},{8,4}},
          rotation=90,
          origin={-16,108})));
    Modelica.Blocks.Sources.Ramp Fuel_ramp(height=6.5, duration=3600) annotation (Placement(transformation(extent={{-140,-60},{-120,-40}})));
    FlameRoomWalls.FlameroomWall_3Layers flameroomWall_3L_side(
      elements=1,
      d_ceramics=0.03,
      d_concrete=0.03,
      d_steel=0.005,
      T_start=(533.15)*ones(flameroomWall_3L_side.elements),
      length=radiation.length + radiation.width,
      width=radiation.height) annotation (Placement(transformation(extent={{86,65},{92,95}})));
    ClaRa.Components.VolumesValvesFittings.Valves.GenericValveGas_L1 adapter1   [bed_Model.n_segments](redeclare model PressureLoss = ClaRa.Components.VolumesValvesFittings.Valves.Fundamentals.LinearNominalPoint (Delta_p_nom=1000, m_flow_nom=10))
                                                                                                                                                                                                    annotation (Placement(transformation(
          extent={{-8,-4},{8,4}},
          rotation=90,
          origin={-16,-106})));
    VolumesValvesFittings.Fittings.SplitGas_L1_flex split(N_ports_out=bed_Model.n_segments, K_split={0.03,0.03,0.06,0.06,0.09,0.13,0.17,0.21,0.17})   annotation (Placement(transformation(extent={{-50,-130},{-30,-110}})));
    ClaRa.Components.BoundaryConditions.PrescribedHeatFlowScalar Aux_burner_1 annotation (Placement(transformation(extent={{64,-88},{44,-68}})));
    ClaRa.Components.BoundaryConditions.PrescribedHeatFlowScalar Aux_burner_2 annotation (Placement(transformation(extent={{64,-102},{44,-82}})));
    ClaRa.Components.BoundaryConditions.PrescribedHeatFlowScalar Aux_burner_3 annotation (Placement(transformation(extent={{64,-116},{44,-96}})));
    Modelica.Blocks.Math.Gain gain(k=0.5) annotation (Placement(transformation(extent={{85,-83},{75,-73}})));
    Modelica.Blocks.Math.Gain gain1(k=0.35) annotation (Placement(transformation(extent={{85,-97},{75,-87}})));
    Modelica.Blocks.Math.Gain gain2(k=0.15) annotation (Placement(transformation(extent={{85,-111},{75,-101}})));
    Modelica.Blocks.Sources.Ramp Aux_burner_ramp(
      height=-3.7e6,
      duration=600,
      offset=3.7e6,
      startTime=3000) annotation (Placement(transformation(extent={{140,-102},{120,-82}})));
    Modelica.Blocks.Sources.Ramp SA_ramp(
      height=-5,
      duration=0,
      offset=6.5*0.6*4.58*0.6,
      startTime=10800) annotation (Placement(transformation(extent={{-140,82},{-120,102}})));
    Modelica.Blocks.Sources.Ramp PA_ramp(
      height=-10,
      duration=0,
      offset=6.5*0.6*4.58*1,
      startTime=9000)  annotation (Placement(transformation(extent={{-140,-124},{-120,-104}})));
    Modelica.Blocks.Sources.Ramp Fuel_ramp1(height=0.005, duration=3600) annotation (Placement(transformation(extent={{-140,-24},{-120,-4}})));
    Modelica.Thermal.HeatTransfer.Components.ThermalConductor thermalConductor(G=2.8e4) annotation (Placement(transformation(extent={{115,125},{105,135}})));
  equation

    connect(gasphase2.air_inlet, boundaryGas_Txim_flow2.gas_a) annotation (Line(
        points={{-36,86},{-62,86}},
        color={118,106,98},
        thickness=0.5));
    connect(boundaryGas_Txim_flow1.gas_a,gasphase1. air_inlet) annotation (Line(
        points={{-60,42},{-36,42}},
        color={118,106,98},
        thickness=0.5));
    connect(boundaryFuel_Txim_flow.fuel_a,bed_Model. fuel_inlet) annotation (Line(
        points={{-60,-56},{-30,-56}},
        color={27,36,42},
        pattern=LinePattern.Solid,
        thickness=0.5));
    connect(joinGas_L2_flex.outlet,gasphase1. fluegas_inlet) annotation (Line(
        points={{-16,26},{-16,32}},
        color={118,106,98},
        thickness=0.5));
    connect(gasphase1.fluegas_outlet,genericValveGas_L1_2. inlet) annotation (Line(
        points={{-16,52},{-16,56}},
        color={118,106,98},
        thickness=0.5));
    connect(genericValveGas_L1_2.outlet, gasphase2.fluegas_inlet) annotation (Line(
        points={{-16,72},{-16,76}},
        color={118,106,98},
        thickness=0.5));
    connect(gasphase2.fluegas_outlet, genericValveGas_L1_3.inlet) annotation (Line(
        points={{-16,96},{-16,100}},
        color={118,106,98},
        thickness=0.5));
    connect(genericValveGas_L1_3.outlet,boundaryGas_pTxi. gas_a) annotation (Line(
        points={{-16,116},{-16,120}},
        color={118,106,98},
        thickness=0.5));
    connect(Fuel_ramp.y, boundaryFuel_Txim_flow.m_flow) annotation (Line(points={{-119,-50},{-80,-50}}, color={0,0,127}));
    connect(bed_Model.fluegas_outlet,adapter. inlet) annotation (Line(
        points={{-16,-46},{-16,-12}},
        color={118,106,98},
        thickness=0.5));
    connect(adapter.outlet,joinGas_L2_flex. inlet) annotation (Line(
        points={{-16,4},{-16,14}},
        color={118,106,98},
        thickness=0.5));
    connect(bed_Model.fuel_outlet,boundaryFuel_pTxi. fuel_a) annotation (Line(
        points={{30,-56},{92,-56}},
        color={27,36,42},
        pattern=LinePattern.Solid,
        thickness=0.5));
    connect(boundaryGas_Txim_flow.gas_a,split. inlet) annotation (Line(
        points={{-60,-120},{-50,-120}},
        color={118,106,98},
        thickness=0.5));
    connect(split.outlet,adapter1. inlet) annotation (Line(
        points={{-30,-120},{-16,-120},{-16,-114}},
        color={118,106,98},
        thickness=0.5));
    connect(bed_Model.heatPort_top,radiation. heatPorts_bottom) annotation (Line(
        points={{2,-46},{2,20},{60,20},{60,50}},
        color={167,25,48},
        thickness=0.5));
    connect(adapter1.outlet, grate.gasPortIn) annotation (Line(
        points={{-16,-98},{-16,-92}},
        color={118,106,98},
        thickness=0.5));
    connect(grate.gasPortOut, bed_Model.fluegas_inlet) annotation (Line(
        points={{-16,-82},{-16,-66}},
        color={118,106,98},
        thickness=0.5));
    connect(grate.heatPort_fuel, bed_Model.heatPort_bottom) annotation (Line(
        points={{2,-82},{2,-66}},
        color={167,25,48},
        thickness=0.5));
    connect(Aux_burner_1.Q_flow, gain.y) annotation (Line(points={{64,-78},{74.5,-78}},             color={0,0,127}));
    connect(Aux_burner_2.Q_flow, gain1.y) annotation (Line(points={{64,-92},{74.5,-92}},                   color={0,0,127}));
    connect(Aux_burner_3.Q_flow, gain2.y) annotation (Line(points={{64,-106},{74.5,-106}},                 color={0,0,127}));
    connect(Aux_burner_ramp.y, gain.u) annotation (Line(points={{119,-92},{104,-92},{104,-78},{86,-78}},
                                                                                                   color={0,0,127}));
    connect(Aux_burner_ramp.y, gain1.u) annotation (Line(points={{119,-92},{86,-92}},                   color={0,0,127}));
    connect(Aux_burner_ramp.y, gain2.u) annotation (Line(points={{119,-92},{104,-92},{104,-106},{86,-106}},
                                                                                                        color={0,0,127}));
    connect(SA_ramp.y, boundaryGas_Txim_flow2.m_flow) annotation (Line(points={{-119,92},{-82,92}}, color={0,0,127}));
    connect(PA_ramp.y, boundaryGas_Txim_flow.m_flow) annotation (Line(points={{-119,-114},{-80,-114}},  color={0,0,127}));
    connect(Fuel_ramp1.y, bed_Model.v_grate_input) annotation (Line(points={{-119,-14},{-46,-14},{-46,-46.2},{-31,-46.2}}, color={0,0,127}));
    connect(flameroomWall_3L_front_back.heatPort_Inner[1], radiation.heatPorts_walls[2]) annotation (Line(
        points={{86,40},{78,40},{78,60.25},{71,60.25}},
        color={167,25,48},
        thickness=0.5));
    connect(flameroomWall_3L_side.heatPort_Inner[1], radiation.heatPorts_walls[1]) annotation (Line(
        points={{86,80},{78,80},{78,59.75},{71,59.75}},
        color={167,25,48},
        thickness=0.5));
    connect(gasphase2.heatPort, radiation.heatPorts_gasPhase[1]) annotation (Line(
        points={{4,86},{26,86},{26,59.75},{49,59.75}},
        color={167,25,48},
        thickness=0.5));
    connect(gasphase1.heatPort, radiation.heatPorts_gasPhase[2]) annotation (Line(
        points={{4,42},{26,42},{26,60.25},{49,60.25}},
        color={167,25,48},
        thickness=0.5));
    connect(radiation.heatPorts_top[1], flameroomWall_2L.heatPort_Inner[1]) annotation (Line(
        points={{60,70},{60.5,70},{60.5,97.25}},
        color={167,25,48},
        thickness=0.5));
    connect(Aux_burner_1.port, bed_Model.heatPort_bottom[8]) annotation (Line(
        points={{44,-78},{2,-78},{2,-65.75}},
        color={167,25,48},
        thickness=0.5));
    connect(Aux_burner_2.port, bed_Model.heatPort_bottom[9]) annotation (Line(
        points={{44,-92},{38,-92},{38,-78},{2,-78},{2,-65.65}},
        color={167,25,48},
        thickness=0.5));
    connect(Aux_burner_3.port, bed_Model.heatPort_bottom[10]) annotation (Line(
        points={{44,-106},{38,-106},{38,-78},{2,-78},{2,-65.55}},
        color={167,25,48},
        thickness=0.5));
    connect(flameroomWall_2L.heatPort_Outer[1], thermalConductor.port_b) annotation (Line(
        points={{60.5,102.75},{60.5,130},{105,130}},
        color={167,25,48},
        thickness=0.5));
    connect(flameroomWall_3L_side.heatPort_Outer[1], thermalConductor.port_b) annotation (Line(
        points={{92,80},{100,80},{100,130},{105,130}},
        color={167,25,48},
        thickness=0.5));
    connect(flameroomWall_3L_front_back.heatPort_Outer[1], thermalConductor.port_b) annotation (Line(
        points={{92,40},{100,40},{100,130},{105,130}},
        color={167,25,48},
        thickness=0.5));
    connect(thermalConductor.port_a, fixedTemperature1.port) annotation (Line(points={{115,130},{120,130}},                     color={191,0,0}));
    annotation (Diagram(coordinateSystem(extent={{-140,-140},{140,140}})),
      experiment(
        StopTime=15000,
        Tolerance=1e-05,
        __Dymola_Algorithm="Sdirk34hw"),
      __Dymola_experimentSetupOutput(equidistant=false),
      __Dymola_experimentFlags(
        Advanced(GenerateVariableDependencies=false, OutputModelicaCode=false),
        Evaluate=false,
        OutputCPUtime=true,
        OutputFlatModelica=false));
  end Grate_Boiler_Example;
end Flameroom;
