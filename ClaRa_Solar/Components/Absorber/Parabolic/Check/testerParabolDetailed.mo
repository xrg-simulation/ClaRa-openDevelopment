within ClaRa_Solar.Components.Absorber.Parabolic.Check;
model testerParabolDetailed
  extends ClaRa.Basics.Icons.PackageIcons.ExecutableExample100;
  ParabolSolarliteDetailled parabolSolarliteDetailled(
    Zs=0,
    n_ax=pipe.N_cv,
    L=pipe.length,
    dx=pipe.Delta_x) annotation (Placement(transformation(extent={{-20,12},{0,32}})));
  SunModels.Sun sun(
    dayOfYear=1,
    startTimeSimulation=8,
    SL=-105,
    L=15,
    LL=-99)                                               annotation (Placement(transformation(extent={{-80,52},{-60,72}})));
  Modelica.Blocks.Sources.RealExpression realExpression(y=850)
    annotation (Placement(transformation(extent={{-80,12},{-60,32}})));
  Modelica.Blocks.Sources.RealExpression realExpression1(y=0)
    annotation (Placement(transformation(extent={{-80,30},{-60,50}})));
  ClaRa.Components.VolumesValvesFittings.Pipes.PipeFlowVLE_L4_Simple
                                                                  pipe(
    redeclare model HeatTransfer = ClaRa.Basics.ControlVolumes.Fundamentals.HeatTransport.VLE_HT.NusseltPipe_L4,
    h_start=ones(3)*100e3,                                             length=100) annotation (Placement(transformation(extent={{-4,-17},{24,-7}})));
  ClaRa.Components.BoundaryConditions.BoundaryVLE_phxi massFlowSink(
    variable_p=true,
    h_const=100e3,
    m_flow_nom=100,
    p_const=1000000,
    Delta_p=100000) annotation (Placement(transformation(
        extent={{10,-10},{-10,10}},
        rotation=180,
        origin={-38,-12})));
  Modelica.Blocks.Sources.Step inlet_pressure(
    offset=1e5,
    startTime=100,
    height=1e4)
    annotation (Placement(transformation(extent={{-80,-28},{-60,-8}})));
  ClaRa.Components.BoundaryConditions.BoundaryVLE_hxim_flow massFlowSource(
    m_flow_const=0.1,
    variable_m_flow=true,
    h_const=200e3,
    m_flow_nom=0,
    p_nom=1000,
    variable_h=false)
                annotation (Placement(transformation(extent={{60,-22},{40,-2}})));
  Modelica.Blocks.Sources.Step outlet_massFlow(
    height=0,
    startTime=0,
    offset=-1)
    annotation (Placement(transformation(extent={{96,-16},{76,4}})));
  inner ClaRa.SimCenter simCenter annotation (Placement(transformation(extent={{-80,-80},{-40,-60}})));
  ClaRa.Basics.ControlVolumes.SolidVolumes.CylindricalThinWall_L4 cylindricalThinWall(
    redeclare model Material = TILMedia.SolidTypes.TILMedia_Steel,
    N_ax=pipe.N_cv,
    diameter_o=pipe.diameter_i + 0.1,
    diameter_i=pipe.diameter_i,
    length=pipe.length) annotation (Placement(transformation(extent={{-4,-2},{24,8}})));
equation
  connect(realExpression.y, parabolSolarliteDetailled.Q_radiation[1]) annotation (Line(
      points={{-59,22},{-20,22}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(realExpression1.y, parabolSolarliteDetailled.Alpha) annotation (Line(
      points={{-59,40},{-4,40},{-4,32}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(sun.y, parabolSolarliteDetailled.SunIn) annotation (Line(
      points={{-59,56},{-15,56},{-15,32}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(inlet_pressure.y, massFlowSink.p) annotation (Line(points={{-59,-18},{-48,-18}},           color={0,0,127}));
  connect(massFlowSink.steam_a, pipe.inlet) annotation (Line(
      points={{-28,-12},{-4,-12}},
      color={0,131,169},
      thickness=0.5));
  connect(massFlowSource.steam_a, pipe.outlet) annotation (Line(
      points={{40,-12},{24,-12}},
      color={0,131,169},
      thickness=0.5));
  connect(massFlowSource.m_flow, outlet_massFlow.y) annotation (Line(points={{62,-6},{75,-6}},                   color={0,0,127}));
  connect(pipe.heat, cylindricalThinWall.innerPhase) annotation (Line(
      points={{10,-8},{10,-2}},
      color={167,25,48},
      thickness=0.5));
  connect(cylindricalThinWall.outerPhase, parabolSolarliteDetailled.port) annotation (Line(
      points={{10,8},{10,22},{0,22}},
      color={167,25,48},
      thickness=0.5));
  annotation (Diagram(coordinateSystem(preserveAspectRatio=true,  extent={{-100,
            -100},{100,100}})),
    experiment(StopTime=28800, __Dymola_Algorithm="Dassl"),
    __Dymola_experimentSetupOutput);
end testerParabolDetailed;
