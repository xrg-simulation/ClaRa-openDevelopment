within ClaRa_Solar.Examples;
model SixBranchesEvapThreeBranchsSH_Turbine_Detailed

  extends ClaRa.Basics.Icons.PackageIcons.ExecutableExample100;

  import Modelica.Constants.pi;

  parameter Modelica.Units.SI.Length d_i=0.065 annotation(Dialog(tab="Parameter"));
  parameter Modelica.Units.SI.Length d_a=0.07 annotation(Dialog(tab="Parameter"));
  parameter Modelica.Units.SI.Length L=720 annotation(Dialog(tab="Parameter"));
  parameter Modelica.Units.SI.Length W=4.6 annotation(Dialog(tab="Parameter"));
  parameter Modelica.Units.SI.Area A_cross=d_i^2*pi annotation(Dialog(tab="Parameter"));
  parameter Integer nNodes=36 annotation(Dialog(tab="Parameter"));
parameter Integer nCells=3 annotation(Dialog(tab="Parameter"));
parameter Integer nCellsHead=3 annotation(Dialog(tab="Parameter"));
parameter Real zeta=8000 annotation(Dialog(tab="Parameter"));
parameter Real deltaTime=lengthDistortion/velocityDistortion annotation(Dialog(tab="Parameter"));
parameter Real deltaIrr_rel=0.6 annotation(Dialog(tab="Parameter"));
parameter Real deltaIrr=deltaIrr_rel*850 annotation(Dialog(tab="Parameter"));
parameter Real Irr=850 annotation(Dialog(tab="Parameter"));
  parameter Modelica.Units.SI.SpecificEnthalpy h_in=600e3 annotation(Dialog(tab="Parameter"));
  parameter Modelica.Units.SI.Length d_iSH=0.065 annotation(Dialog(tab="Parameter"));
  parameter Modelica.Units.SI.Length d_aSH=0.07 annotation(Dialog(tab="Parameter"));
  parameter Modelica.Units.SI.Length LSH=240 annotation(Dialog(tab="Parameter"));
  parameter Integer nNodesSH=8 annotation(Dialog(tab="Parameter"));
  parameter Modelica.Units.SI.MassFlowRate m_flowCirc=3*3.2715 annotation(Dialog(tab="Parameter"));
  parameter Modelica.Units.SI.MassFlowRate deltam_flowCirc=3*3.2715 annotation(Dialog(tab="Parameter"));
  parameter Modelica.Units.SI.Pressure p_out=3000000 annotation(Dialog(tab="Parameter"));
  parameter Modelica.Units.SI.Pressure deltap_out=3000000 annotation(Dialog(tab="Parameter"));
  parameter Modelica.Units.SI.Time deltaTimeMass=500 annotation(Dialog(tab="Parameter"));
  parameter Modelica.Units.SI.Mass deltaMass=500 annotation(Dialog(tab="Parameter"));
  parameter Modelica.Units.SI.Time startingTimeDistortion=20000 annotation(Dialog(tab="Parameter"));
  parameter Modelica.Units.SI.Velocity velocityDistortion=1 annotation(Dialog(tab="Parameter"));
  parameter Modelica.Units.SI.Length lengthDistortion=320 annotation(Dialog(tab="Parameter"));
  parameter Modelica.Units.SI.Time deltaTimeDistortion=20 annotation(Dialog(tab="Parameter"));
  parameter Integer n_Qrad_Evap=36 annotation(Dialog(tab="Parameter"));
parameter Integer n_Qrad_SH=8 annotation(Dialog(tab="Parameter"));
parameter Integer nEvapLoop=6 annotation(Dialog(tab="Parameter"));
parameter Integer nSHLoop=3 annotation(Dialog(tab="Parameter"));
  parameter Modelica.Units.SI.Pressure deltap_head=200000 annotation(Dialog(tab="Parameter"));

//Initialization
  parameter Modelica.Units.SI.Temperature t_startWallEvap=473.15 annotation(Dialog(tab="Initialisation"));
  parameter Modelica.Units.SI.Temperature t_startWallSH=473.15 + 80
                                                                  annotation(Dialog(tab="Initialisation"));
  parameter Modelica.Units.SI.Pressure p_startEvapIn=4450000 annotation(Dialog(tab="Initialisation"));
  parameter Modelica.Units.SI.Pressure p_startDeltaValve=200000 annotation(Dialog(tab="Initialisation"));
  parameter Modelica.Units.SI.Pressure p_startDeltaEvapCollInTotal=150000 annotation(Dialog(tab="Initialisation"));
  parameter Modelica.Units.SI.Pressure p_startDeltaEvap=600000 annotation(Dialog(tab="Initialisation"));

//Sun
  parameter Real dayOfYear=1  "Day of Year" annotation(Dialog(group="Time and Date"));
  parameter Real startTimeSimulation=7  "Local Time at Simulation Start in Hours" annotation(Dialog(group="Time and Date"));
  parameter Real DS=0 "Daylight Saving (if daylight saving (summertime) 60 else 0)"
                                                                                   annotation(Dialog(group="Time and Date"));
  parameter Real SL=-105 "Standard Longitude in Degree (location east of Greenwich is negative)" annotation(Dialog(group="Location"));
  parameter Real LLat=15   "Local Latitude in Degree" annotation(Dialog(group="Location"));
  parameter Real LLong=-99.0 "Local Longitude in Degree (location east of Greenwich is negative)"
                                                                                                 annotation(Dialog(group="Location"));
  parameter Boolean Green=true "True if Location East of Greenwich" annotation(Dialog(group="Location"));

  inner ClaRa.SimCenter          modelProperties(   useHomotopy=false, redeclare replaceable TILMedia.VLEFluidTypes.TILMedia_SplineWater fluid1)
    annotation (Placement(transformation(extent={{380,360},{420,380}})));

  ClaRa.Basics.ControlVolumes.SolidVolumes.CylindricalThinWall_L4
                                 wall(
    redeclare replaceable model Material = TILMedia.SolidTypes.TILMedia_Steel,
    N_ax=tube.N_cv,
    Delta_x=tube.Delta_x,
    diameter_o=d_a,
    diameter_i=d_i,
    length=L,
    T_start=ones(nNodes)*(t_startWallEvap))
         annotation (Placement(transformation(extent={{54,-22},{74,-14}})));
  Components.Absorber.Parabolic.ParabolSolarliteDetailled parabolSolarliteDetailled5(
    L=L,
    n_ax=nNodes,
    W=W,
    dx=tube1.Delta_x,
    n_Qrad_ax=1,
    Zs=0) annotation (Placement(transformation(extent={{14,40},{34,60}})));
  Modelica.Blocks.Sources.Constant const(k=0)
    annotation (Placement(transformation(extent={{-6,60},{14,80}})));
  ClaRa.Components.MechanicalSeparation.SteamSeparatorVLE_L3
                                                           seperator(
    yps_start=0,
    initOption=0,
    m_flow_nom=NOM.separator.m_flow_1,
    redeclare model PressureLoss = ClaRa.Basics.ControlVolumes.Fundamentals.PressureLoss.Generic_PL.LinearParallelZones_L3,
    p_start=INIT.separator.p,
    p_nom=NOM.separator.p,
    levelOutput=true)
                     annotation (Placement(transformation(extent={{-90,-40},{-110,-20}})));

  ClaRa.Basics.ControlVolumes.SolidVolumes.CylindricalThinWall_L4
                                                  wallSH(
    redeclare replaceable model Material = TILMedia.SolidTypes.TILMedia_Steel,
    N_ax=nNodesSH,
    Delta_x=tubeSH.Delta_x,
    diameter_o=d_aSH,
    diameter_i=d_iSH,
    length=LSH,
    T_start=ones(nNodesSH)*(t_startWallSH))
         annotation (Placement(transformation(extent={{-240,-22},{-220,-14}})));
  Components.Absorber.Parabolic.ParabolSolarliteDetailled parabolSolarliteDetailled1(
    L=L,
    n_ax=nNodes,
    W=W,
    dx=tube.Delta_x,
    n_Qrad_ax=1,
    Zs=0) annotation (Placement(transformation(extent={{14,-20},{34,0}})));
  Modelica.Blocks.Sources.Constant const1(
                                         k=0)
    annotation (Placement(transformation(extent={{-6,0},{14,20}})));
  Components.Absorber.Parabolic.ParabolSolarliteDetailled parabolSolarliteDetailled7(
    W=W,
    L=LSH,
    n_ax=nNodesSH,
    dx=tubeSH.Delta_x,
    n_Qrad_ax=1,
    Zs=0) annotation (Placement(transformation(extent={{-268,-20},{-248,0}})));
  Modelica.Blocks.Sources.Constant const2(
                                         k=0)
    annotation (Placement(transformation(extent={{-292,10},{-272,30}})));

  Modelica.Blocks.Sources.RealExpression realExpression20(y=850)
    annotation (Placement(transformation(extent={{-300,-20},{-280,0}})));
  Modelica.Blocks.Sources.RealExpression realExpression22(y=850)
    annotation (Placement(transformation(extent={{-26,40},{-6,60}})));
  Modelica.Blocks.Sources.RealExpression realExpression21(y=850)
    annotation (Placement(transformation(extent={{-26,-20},{-6,0}})));

  Components.Absorber.Parabolic.ParabolSolarliteDetailled parabolSolarliteDetailled4(
    L=L,
    n_ax=nNodes,
    W=W,
    dx=tube1.Delta_x,
    n_Qrad_ax=1,
    Zs=0) annotation (Placement(transformation(extent={{14,100},{34,120}})));
  Modelica.Blocks.Sources.Constant const3(
                                         k=0)
    annotation (Placement(transformation(extent={{-6,120},{14,140}})));
  Modelica.Blocks.Sources.RealExpression realExpression23(y=850)
    annotation (Placement(transformation(extent={{-26,100},{-6,120}})));

  Components.Absorber.Parabolic.ParabolSolarliteDetailled parabolSolarliteDetailled3(
    L=L,
    n_ax=nNodes,
    W=W,
    dx=tube1.Delta_x,
    n_Qrad_ax=1,
    Zs=0) annotation (Placement(transformation(extent={{14,160},{34,180}})));
  Modelica.Blocks.Sources.Constant const4(
                                         k=0)
    annotation (Placement(transformation(extent={{-6,180},{14,200}})));
  Modelica.Blocks.Sources.RealExpression realExpression24(y=850)
    annotation (Placement(transformation(extent={{-26,160},{-6,180}})));

  Modelica.Blocks.Sources.Constant const5(
                                         k=0)
    annotation (Placement(transformation(extent={{-290,80},{-270,100}})));
  Modelica.Blocks.Sources.RealExpression realExpression19(y=850)
    annotation (Placement(transformation(extent={{-300,50},{-280,70}})));
  Modelica.Blocks.Sources.RealExpression realExpression2(y=-seperator.outlet2.m_flow
        /0.4/nEvapLoop)
    annotation (Placement(transformation(extent={{254,0},{234,20}})));
  Modelica.Blocks.Sources.RealExpression realExpression6(y=tube.inlet.m_flow)
    annotation (Placement(transformation(extent={{254,-22},{234,-2}})));
  ClaRa.Components.Utilities.Blocks.LimPID         PIDValveBranch1(
    y_max=1,
    y_min=0.02,
    Tau_i=10,
    Tau_d=500,
    initOption=503,
    limitsAtInit=true,
    controllerType=Modelica.Blocks.Types.SimpleController.PI,
    k=0.01,
    y_start=1,
    sign=1)    annotation (Placement(transformation(extent={{214,0},{194,20}})));
  ClaRa.Components.VolumesValvesFittings.Valves.GenericValveVLE_L1   valveBranch(
    redeclare model PressureLoss = ClaRa.Components.VolumesValvesFittings.Valves.Fundamentals.LinearNominalPoint (Delta_p_nom=2e4, m_flow_nom=2.7),
    openingInputIsActive=true,
    checkValve=false)
    annotation (Placement(transformation(extent={{134,-36},{114,-24}})));
  Modelica.Blocks.Sources.RealExpression realExpression3(y=-seperator.outlet2.m_flow
        /0.4/nEvapLoop)
    annotation (Placement(transformation(extent={{254,60},{234,80}})));
  Modelica.Blocks.Sources.RealExpression realExpression5(y=tube1.inlet.m_flow)
    annotation (Placement(transformation(extent={{254,38},{234,58}})));
  Modelica.Blocks.Sources.RealExpression realExpression7(y=-seperator.outlet2.m_flow
        /0.4/nEvapLoop)
    annotation (Placement(transformation(extent={{254,120},{234,140}})));
  Modelica.Blocks.Sources.RealExpression realExpression8(y=tube2.inlet.m_flow)
    annotation (Placement(transformation(extent={{254,98},{234,118}})));
  Modelica.Blocks.Sources.RealExpression realExpression9(y=-seperator.outlet2.m_flow
        /0.4/nEvapLoop)
    annotation (Placement(transformation(extent={{254,180},{234,200}})));
  Modelica.Blocks.Sources.RealExpression realExpression10(y=tube3.inlet.m_flow)
    annotation (Placement(transformation(extent={{254,158},{234,178}})));

  ClaRa.Components.VolumesValvesFittings.Fittings.JoinVLE_L2_Y
                                                            join_Y(
    volume=0.001,
    p_nom=NOM.join_Y.p,
    h_nom=NOM.join_Y.h1,
    h_start=INIT.join_Y.h1,
    p_start=INIT.join_Y.p)
    annotation (Placement(transformation(extent={{250,-20},{230,-40}})));
  ClaRa.Components.VolumesValvesFittings.Fittings.SplitVLE_L2_Y
                                          split_Y(
    redeclare model PressureLossIn = ClaRa.Components.VolumesValvesFittings.Fittings.Fundamentals.Linear (m_flow_nom=2.2*nSHLoop, dp_nom=1e4),
    volume=0.001,
    p_nom=NOM.split_Y.p,
    h_nom=NOM.split_Y.h1,
    h_start=INIT.split_Y.h1,
    p_start=INIT.split_Y.p) annotation (Placement(transformation(extent={{-160,-20},{-180,-40}})));
  ClaRa.Components.VolumesValvesFittings.Fittings.SplitVLE_L2_Y
                                                             split_Y1(
    volume=0.001,
    m_flow_out_nom={1,5},
    p_nom=NOM.split_Y1.p,
    h_nom=NOM.split_Y1.h1,
    p_start=INIT.split_Y1.p,
    h_start=INIT.split_Y1.h1)
    annotation (Placement(transformation(extent={{176,-20},{156,-40}})));
  ClaRa.Components.VolumesValvesFittings.Fittings.SplitVLE_L2_Y split_Y2(
    volume=0.001,
    m_flow_out_nom={1,4},
    p_nom=NOM.split_Y2.p,
    h_nom=NOM.split_Y2.h1,
    h_start=INIT.split_Y2.h1,
    p_start=INIT.split_Y2.p)                               annotation (Placement(transformation(extent={{176,40},{156,20}})));
  ClaRa.Components.VolumesValvesFittings.Fittings.SplitVLE_L2_Y
                                                             split_Y3(
    volume=0.001,
    m_flow_out_nom={1,3},
    p_nom=NOM.split_Y3.p,
    h_nom=NOM.split_Y3.h1,
    h_start=INIT.split_Y3.h1,
    p_start=INIT.split_Y3.p)
    annotation (Placement(transformation(extent={{176,100},{156,80}})));
  ClaRa.Components.VolumesValvesFittings.Fittings.JoinVLE_L2_Y
                                                            join_Y3(
    volume=0.001,
    p_nom=NOM.join_Y3.p,
    h_nom=NOM.join_Y3.h1,
    h_start=INIT.join_Y3.h1,
    p_start=INIT.join_Y3.p)                                       annotation (
      Placement(transformation(
        extent={{10,10},{-10,-10}},
        rotation=90,
        origin={-36,90})));
  ClaRa.Components.VolumesValvesFittings.Fittings.JoinVLE_L2_Y
                                                            join_Y2(
    volume=0.001,
    p_nom=NOM.join_Y2.p,
    h_nom=NOM.join_Y2.h1,
    h_start=INIT.join_Y2.h1,
    p_start=INIT.join_Y2.p)                                       annotation (
      Placement(transformation(
        extent={{10,10},{-10,-10}},
        rotation=90,
        origin={-36,30})));
  ClaRa.Components.VolumesValvesFittings.Fittings.JoinVLE_L2_Y join_Y1(
    volume=0.001,
    p_nom=NOM.join_Y1.p,
    h_nom=NOM.join_Y1.h1,
    h_start=NOM.join_Y1.h1,
    p_start=NOM.join_Y1.p)                                        annotation (Placement(transformation(
        extent={{-10,10},{10,-10}},
        rotation=180,
        origin={-36,-30})));
  ClaRa.Components.VolumesValvesFittings.Fittings.JoinVLE_L2_Y join_Y6(
    volume=0.001,
    p_nom=NOM.join_Y6.p,
    h_nom=NOM.join_Y6.h1,
    h_start=INIT.join_Y6.h1,
    p_start=INIT.join_Y6.p)
                          annotation (Placement(transformation(
        extent={{-10,10},{10,-10}},
        rotation=180,
        origin={-320,-30})));
  ClaRa.Components.VolumesValvesFittings.Fittings.JoinVLE_L2_Y
                                                            join_Y7(
    volume=0.001,
    p_nom=NOM.join_Y7.p,
    h_nom=NOM.join_Y7.h1,
    h_start=INIT.join_Y7.h1,
    p_start=INIT.join_Y7.p)                                       annotation (
      Placement(transformation(
        extent={{10,10},{-10,-10}},
        rotation=90,
        origin={-320,40})));
  ClaRa.Components.VolumesValvesFittings.Fittings.SplitVLE_L2_Y
                                          split_Y6(
    redeclare model PressureLossIn = ClaRa.Components.VolumesValvesFittings.Fittings.Fundamentals.Linear (m_flow_nom=2.2*(nSHLoop - 1), dp_nom=1e4),
    volume=0.001,
    p_nom=NOM.split_Y6.p,
    h_nom=NOM.split_Y6.h1,
    h_start=INIT.split_Y6.h1,
    p_start=INIT.split_Y6.p)      annotation (Placement(transformation(extent={{-180,50},{-200,30}})));

  Components.Absorber.Parabolic.ParabolSolarliteDetailled parabolSolarliteDetailled6(
    W=W,
    L=LSH,
    n_ax=nNodesSH,
    dx=tubeSH.Delta_x,
    n_Qrad_ax=1,
    Zs=0) annotation (Placement(transformation(extent={{-268,120},{-248,140}})));
  Modelica.Blocks.Sources.Constant const6(
                                         k=0)
    annotation (Placement(transformation(extent={{-290,150},{-270,170}})));
  Modelica.Blocks.Sources.RealExpression realExpression(y=850)
    annotation (Placement(transformation(extent={{-300,120},{-280,140}})));
  Components.Absorber.Parabolic.ParabolSolarliteDetailled parabolSolarliteDetailled2(
    L=L,
    n_ax=nNodes,
    W=W,
    dx=tube1.Delta_x,
    n_Qrad_ax=1,
    Zs=0) annotation (Placement(transformation(extent={{14,220},{34,240}})));
  Modelica.Blocks.Sources.Constant const7(
                                         k=0)
    annotation (Placement(transformation(extent={{-6,240},{14,260}})));
  Modelica.Blocks.Sources.RealExpression realExpression28(y=850)
    annotation (Placement(transformation(extent={{-26,220},{-6,240}})));
  Components.Absorber.Parabolic.ParabolSolarliteDetailled parabolSolarliteDetailled(
    L=L,
    n_ax=nNodes,
    W=W,
    dx=tube1.Delta_x,
    Zs=0) annotation (Placement(transformation(extent={{14,280},{34,300}})));
  Modelica.Blocks.Sources.Constant const8(
                                         k=0)
    annotation (Placement(transformation(extent={{-6,300},{14,320}})));
  Modelica.Blocks.Sources.RealExpression realExpression18(y=850)
    annotation (Placement(transformation(extent={{-26,280},{-6,300}})));
  Modelica.Blocks.Sources.RealExpression realExpression11(
                                                         y=-seperator.outlet2.m_flow
        /0.4/nEvapLoop)
    annotation (Placement(transformation(extent={{254,240},{234,260}})));
  Modelica.Blocks.Sources.RealExpression realExpression12(y=tube4.inlet.m_flow)
    annotation (Placement(transformation(extent={{254,218},{234,238}})));
  Modelica.Blocks.Sources.RealExpression realExpression13(
                                                         y=-seperator.outlet2.m_flow
        /0.4/nEvapLoop)
    annotation (Placement(transformation(extent={{254,300},{234,320}})));
  Modelica.Blocks.Sources.RealExpression realExpression14(y=tube5.inlet.m_flow)
    annotation (Placement(transformation(extent={{254,278},{234,298}})));
  ClaRa.Components.VolumesValvesFittings.Fittings.JoinVLE_L2_Y
                                                            join_Y4(
    volume=0.001,
    p_nom=NOM.join_Y4.p,
    h_nom=NOM.join_Y4.h1,
    h_start=INIT.join_Y4.h1,
    p_start=INIT.join_Y4.p)                                       annotation (
      Placement(transformation(
        extent={{10,10},{-10,-10}},
        rotation=90,
        origin={-36,150})));
  ClaRa.Components.VolumesValvesFittings.Fittings.JoinVLE_L2_Y join_Y5(
    volume=0.001,
    p_nom=NOM.join_Y5.p,
    h_nom=NOM.join_Y5.h1,
    h_start=INIT.join_Y5.h1,
    p_start=INIT.join_Y5.p)                                       annotation (Placement(transformation(
        extent={{10,10},{-10,-10}},
        rotation=90,
        origin={-36,210})));
  ClaRa.Components.VolumesValvesFittings.Fittings.SplitVLE_L2_Y split_Y4(
    volume=0.001,
    m_flow_out_nom={1,2},
    p_nom=NOM.split_Y45.p,
    h_nom=NOM.split_Y45.h1,
    h_start=INIT.split_Y45.h1,
    p_start=INIT.split_Y45.p)                                annotation (Placement(transformation(extent={{174,160},{154,140}})));
  ClaRa.Components.VolumesValvesFittings.Fittings.SplitVLE_L2_Y
                                                             split_Y5(
    volume=0.001,
    m_flow_out_nom={1,1},
    p_nom=NOM.split_Y5.p,
    h_nom=NOM.split_Y5.h1,
    h_start=INIT.split_Y5.h1,
    p_start=INIT.split_Y5.p)
    annotation (Placement(transformation(extent={{174,220},{154,200}})));
  ClaRa.Components.VolumesValvesFittings.Pipes.PipeFlowVLE_L4_Simple TubeIn(
    frictionAtOutlet=true,
    redeclare model PressureLoss = ClaRa.Basics.ControlVolumes.Fundamentals.PressureLoss.Generic_PL.LinearPressureLoss_L4,
    m_flow_nom=NOM.TubeIn.m_flow,
    Delta_p_nom=NOM.TubeIn.Delta_p_nom,
    h_start=ones(nCells)*INIT.TubeIn.h_in,
    p_start=linspace(
        INIT.TubeIn.p_in,
        INIT.TubeIn.p_out,
        nCells),
    showData=true,
    N_cv=nCells) annotation (Placement(transformation(extent={{294,-35},{266,-25}})));
  ClaRa.Components.VolumesValvesFittings.Pipes.PipeFlowVLE_L4_Simple TubeRePumpOut(
    frictionAtOutlet=true,
    redeclare model PressureLoss = ClaRa.Basics.ControlVolumes.Fundamentals.PressureLoss.Generic_PL.LinearPressureLoss_L4,
    m_flow_nom=NOM.TubeRePumpOut.m_flow,
    Delta_p_nom=NOM.TubeRePumpOut.Delta_p_nom,
    h_start=ones(nCells)*INIT.TubeRePumpOut.h_in,
    p_start=linspace(
        INIT.TubeRePumpOut.p_in,
        INIT.TubeRePumpOut.p_out,
        nCells),
    showData=true,
    N_cv=nCells) annotation (Placement(transformation(extent={{186,-73},{214,-63}})));
  ClaRa.Components.VolumesValvesFittings.Pipes.PipeFlowVLE_L4_Simple TubeIn1(
    frictionAtInlet=true,
    frictionAtOutlet=true,
    redeclare model PressureLoss = ClaRa.Basics.ControlVolumes.Fundamentals.PressureLoss.Generic_PL.LinearPressureLoss_L4,
    m_flow_nom=NOM.TubeIn1.m_flow,
    Delta_p_nom=NOM.TubeIn1.Delta_p_nom,
    h_start=ones(nCells)*INIT.TubeIn1.h_in,
    p_start=linspace(
        INIT.TubeIn1.p_in,
        INIT.TubeIn1.p_out,
        nCells),
    N_cv=nCells) annotation (Placement(transformation(extent={{214,-35},{186,-25}})));
  ClaRa.Components.Utilities.Blocks.LimPID         PIDValveBranch2(
    y_max=1,
    y_min=0.02,
    Tau_i=10,
    Tau_d=500,
    initOption=503,
    limitsAtInit=true,
    controllerType=Modelica.Blocks.Types.SimpleController.PI,
    k=0.01,
    y_start=1,
    sign=1)    annotation (Placement(transformation(extent={{214,60},{194,80}})));
  ClaRa.Components.Utilities.Blocks.LimPID         PIDValveBranch3(
    y_max=1,
    y_min=0.02,
    Tau_i=10,
    Tau_d=500,
    initOption=503,
    limitsAtInit=true,
    controllerType=Modelica.Blocks.Types.SimpleController.PI,
    k=0.01,
    y_start=1,
    sign=1)    annotation (Placement(transformation(extent={{214,120},{194,140}})));
  ClaRa.Components.Utilities.Blocks.LimPID         PIDValveBranch4(
    y_max=1,
    y_min=0.02,
    Tau_i=10,
    Tau_d=500,
    initOption=503,
    limitsAtInit=true,
    controllerType=Modelica.Blocks.Types.SimpleController.PI,
    k=0.01,
    y_start=1,
    sign=1)    annotation (Placement(transformation(extent={{214,180},{194,200}})));
  ClaRa.Components.Utilities.Blocks.LimPID PIDValveBranch5(
    y_max=1,
    y_min=0.02,
    Tau_i=10,
    Tau_d=500,
    initOption=503,
    limitsAtInit=true,
    controllerType=Modelica.Blocks.Types.SimpleController.PI,
    k=0.01,
    y_start=1,
    sign=1) annotation (Placement(transformation(extent={{214,240},{194,260}})));
  ClaRa.Components.Utilities.Blocks.LimPID PIDValveBranch6(
    y_max=1,
    y_min=0.02,
    Tau_i=10,
    Tau_d=500,
    initOption=503,
    limitsAtInit=true,
    controllerType=Modelica.Blocks.Types.SimpleController.PI,
    k=0.01,
    y_start=1,
    sign=1) annotation (Placement(transformation(extent={{214,300},{194,320}})));
  ClaRa.Components.VolumesValvesFittings.Valves.GenericValveVLE_L1   valveBranch1(redeclare model PressureLoss = ClaRa.Components.VolumesValvesFittings.Valves.Fundamentals.LinearNominalPoint (Delta_p_nom=2e4, m_flow_nom=2.7), openingInputIsActive=true,
    checkValve=false)
    annotation (Placement(transformation(extent={{134,24},{114,36}})));
  ClaRa.Components.VolumesValvesFittings.Valves.GenericValveVLE_L1   valveBranch2(redeclare model PressureLoss = ClaRa.Components.VolumesValvesFittings.Valves.Fundamentals.LinearNominalPoint (Delta_p_nom=2e4, m_flow_nom=2.7), openingInputIsActive=true,
    checkValve=false)
    annotation (Placement(transformation(extent={{134,84},{114,96}})));
  ClaRa.Components.VolumesValvesFittings.Valves.GenericValveVLE_L1   valveBranch3(redeclare model PressureLoss = ClaRa.Components.VolumesValvesFittings.Valves.Fundamentals.LinearNominalPoint (Delta_p_nom=2e4, m_flow_nom=2.7), openingInputIsActive=true,
    checkValve=false)
    annotation (Placement(transformation(extent={{134,144},{114,156}})));
  ClaRa.Components.VolumesValvesFittings.Valves.GenericValveVLE_L1   valveBranch4(redeclare model PressureLoss = ClaRa.Components.VolumesValvesFittings.Valves.Fundamentals.LinearNominalPoint (Delta_p_nom=2e4, m_flow_nom=2.7), openingInputIsActive=true,
    checkValve=false)
    annotation (Placement(transformation(extent={{134,204},{114,216}})));
  ClaRa.Components.VolumesValvesFittings.Valves.GenericValveVLE_L1   valveBranch5(redeclare model PressureLoss = ClaRa.Components.VolumesValvesFittings.Valves.Fundamentals.LinearNominalPoint (Delta_p_nom=2e4, m_flow_nom=2.7), openingInputIsActive=true,
    checkValve=false)
    annotation (Placement(transformation(extent={{134,264},{114,276}})));
  ClaRa.Components.VolumesValvesFittings.Pipes.PipeFlowVLE_L4_Simple TubeEvapCollIn1(
    frictionAtInlet=true,
    frictionAtOutlet=true,
    Delta_p_nom=NOM.TubeEvapCollIn1.Delta_p_nom,
    length=52.8,
    diameter_i=d_i,
    N_cv=nCellsHead,
    p_nom=ones(nCellsHead)*1e5,
    h_nom=ones(nCellsHead)*1e5,
    h_start=ones(nCellsHead)*INIT.TubeEvapCollIn1.h_in,
    p_start=linspace(
        INIT.TubeEvapCollIn1.p_in,
        INIT.TubeEvapCollIn1.p_out,
        nCellsHead),
    redeclare model PressureLoss = ClaRa.Basics.ControlVolumes.Fundamentals.PressureLoss.Generic_PL.LinearPressureLoss_L4,
    m_flow_nom=NOM.TubeEvapCollIn1.m_flow)
                       annotation (Placement(transformation(
        extent={{14,-5},{-14,5}},
        rotation=270,
        origin={166,2})));
  ClaRa.Components.VolumesValvesFittings.Pipes.PipeFlowVLE_L4_Simple TubeEvapCollIn2(
    frictionAtInlet=true,
    frictionAtOutlet=true,
    Delta_p_nom=NOM.TubeEvapCollIn2.Delta_p_nom,
    length=52.8,
    diameter_i=d_i,
    N_cv=nCellsHead,
    p_nom=ones(nCellsHead)*1e5,
    h_nom=ones(nCellsHead)*1e5,
    h_start=ones(nCellsHead)*INIT.TubeEvapCollIn2.h_in,
    p_start=linspace(
        INIT.TubeEvapCollIn2.p_in,
        INIT.TubeEvapCollIn2.p_out,
        nCellsHead),
    redeclare model PressureLoss = ClaRa.Basics.ControlVolumes.Fundamentals.PressureLoss.Generic_PL.LinearPressureLoss_L4,
    m_flow_nom=NOM.TubeEvapCollIn2.m_flow)
                       annotation (Placement(transformation(
        extent={{14,-5},{-14,5}},
        rotation=270,
        origin={166,62})));
  ClaRa.Components.VolumesValvesFittings.Pipes.PipeFlowVLE_L4_Simple TubeEvapCollIn3(
    frictionAtInlet=true,
    frictionAtOutlet=true,
    Delta_p_nom=NOM.TubeEvapCollIn3.Delta_p_nom,
    length=52.8,
    diameter_i=d_i,
    N_cv=nCellsHead,
    p_nom=ones(nCellsHead)*1e5,
    h_nom=ones(nCellsHead)*1e5,
    h_start=ones(nCellsHead)*INIT.TubeEvapCollIn3.h_in,
    p_start=linspace(
        INIT.TubeEvapCollIn3.p_in,
        INIT.TubeEvapCollIn3.p_out,
        nCellsHead),
    redeclare model PressureLoss = ClaRa.Basics.ControlVolumes.Fundamentals.PressureLoss.Generic_PL.LinearPressureLoss_L4,
    m_flow_nom=NOM.TubeEvapCollIn3.m_flow)
                       annotation (Placement(transformation(
        extent={{14,-5},{-14,5}},
        rotation=270,
        origin={166,122})));
  ClaRa.Components.VolumesValvesFittings.Pipes.PipeFlowVLE_L4_Simple TubeEvapCollIn4(
    frictionAtInlet=true,
    frictionAtOutlet=true,
    Delta_p_nom=NOM.TubeEvapCollIn4.Delta_p_nom,
    length=52.8,
    diameter_i=d_i,
    N_cv=nCellsHead,
    p_nom=ones(nCellsHead)*1e5,
    h_nom=ones(nCellsHead)*1e5,
    h_start=ones(nCellsHead)*INIT.TubeEvapCollIn4.h_in,
    p_start=linspace(
        INIT.TubeEvapCollIn4.p_in,
        INIT.TubeEvapCollIn4.p_out,
        nCellsHead),
    redeclare model PressureLoss = ClaRa.Basics.ControlVolumes.Fundamentals.PressureLoss.Generic_PL.LinearPressureLoss_L4,
    m_flow_nom=NOM.TubeEvapCollIn4.m_flow)
                       annotation (Placement(transformation(
        extent={{14,-5},{-14,5}},
        rotation=270,
        origin={164,182})));
  ClaRa.Components.VolumesValvesFittings.Pipes.PipeFlowVLE_L4_Simple TubeEvapCollIn5(
    frictionAtInlet=true,
    frictionAtOutlet=false,
    Delta_p_nom=NOM.TubeEvapCollIn5.Delta_p_nom,
    showData=true,
    length=52.8,
    diameter_i=d_i,
    N_cv=nCellsHead,
    p_nom=ones(nCellsHead)*1e5,
    h_nom=ones(nCellsHead)*1e5,
    h_start=ones(nCellsHead)*INIT.TubeEvapCollIn5.h_in,
    p_start=linspace(
        INIT.TubeEvapCollIn5.p_in,
        INIT.TubeEvapCollIn5.p_out,
        nCellsHead),
    redeclare model PressureLoss = ClaRa.Basics.ControlVolumes.Fundamentals.PressureLoss.Generic_PL.LinearPressureLoss_L4,
    m_flow_nom=NOM.TubeEvapCollIn5.m_flow)
                       annotation (Placement(transformation(
        extent={{14,-5},{-14,5}},
        rotation=270,
        origin={164,240})));
  ClaRa.Components.VolumesValvesFittings.Pipes.PipeFlowVLE_L4_Simple tube(
    frictionAtOutlet=true,
    Delta_p_nom=NOM.tube.Delta_p_nom,
    showData=true,
    redeclare model HeatTransfer = ClaRa.Basics.ControlVolumes.Fundamentals.HeatTransport.VLE_HT.NusseltPipe_L4,
    m_flow_nom=NOM.tube.m_flow,
    h_start=linspace(
        INIT.tube.h_in,
        INIT.tube.h_out,
        nNodes),
    p_start=linspace(
        INIT.tube.p_in,
        INIT.tube.p_out,
        nNodes),
    redeclare model PressureLoss = ClaRa.Basics.ControlVolumes.Fundamentals.PressureLoss.Generic_PL.LinearPressureLoss_L4,
    length=L,
    diameter_i=d_i,
    N_cv=nNodes) annotation (Placement(transformation(
        extent={{14,-5},{-14,5}},
        rotation=0,
        origin={64,-30})));
  ClaRa.Basics.ControlVolumes.SolidVolumes.CylindricalThinWall_L4
                                 wall1(redeclare replaceable model Material = TILMedia.SolidTypes.TILMedia_Steel,
    N_ax=tube1.N_cv,
    Delta_x=tube1.Delta_x,
    diameter_o=d_a,
    diameter_i=d_i,
    length=L,
    T_start=ones(nNodes)*(t_startWallEvap))
         annotation (Placement(transformation(extent={{54,38},{74,46}})));
  ClaRa.Components.VolumesValvesFittings.Pipes.PipeFlowVLE_L4_Simple tube1(
    frictionAtOutlet=true,
    Delta_p_nom=NOM.tube1.Delta_p_nom,
    showData=true,
    redeclare model HeatTransfer = ClaRa.Basics.ControlVolumes.Fundamentals.HeatTransport.VLE_HT.NusseltPipe_L4,
    m_flow_nom=NOM.tube1.m_flow,
    h_start=linspace(
        INIT.tube1.h_in,
        INIT.tube1.h_out,
        nNodes),
    p_start=linspace(
        INIT.tube1.p_in,
        INIT.tube1.p_out,
        nNodes),
    redeclare model PressureLoss = ClaRa.Basics.ControlVolumes.Fundamentals.PressureLoss.Generic_PL.LinearPressureLoss_L4,
    length=L,
    diameter_i=d_i,
    N_cv=nNodes) annotation (Placement(transformation(
        extent={{14,-5},{-14,5}},
        rotation=0,
        origin={64,30})));
  ClaRa.Basics.ControlVolumes.SolidVolumes.CylindricalThinWall_L4
                                 wall2(redeclare replaceable model Material = TILMedia.SolidTypes.TILMedia_Steel,
    N_ax=tube2.N_cv,
    Delta_x=tube2.Delta_x,
    diameter_o=d_a,
    diameter_i=d_i,
    length=L,
    T_start=ones(nNodes)*(t_startWallEvap))
         annotation (Placement(transformation(extent={{54,98},{74,106}})));
  ClaRa.Components.VolumesValvesFittings.Pipes.PipeFlowVLE_L4_Simple tube2(
    frictionAtOutlet=true,
    Delta_p_nom=NOM.tube2.Delta_p_nom,
    showData=true,
    redeclare model HeatTransfer = ClaRa.Basics.ControlVolumes.Fundamentals.HeatTransport.VLE_HT.NusseltPipe_L4,
    m_flow_nom=NOM.tube2.m_flow,
    h_start=linspace(
        INIT.tube2.h_in,
        INIT.tube2.h_out,
        nNodes),
    p_start=linspace(
        INIT.tube2.p_in,
        INIT.tube2.p_out,
        nNodes),
    redeclare model PressureLoss = ClaRa.Basics.ControlVolumes.Fundamentals.PressureLoss.Generic_PL.LinearPressureLoss_L4,
    length=L,
    diameter_i=d_i,
    N_cv=nNodes) annotation (Placement(transformation(
        extent={{14,-5},{-14,5}},
        rotation=0,
        origin={64,90})));
  ClaRa.Basics.ControlVolumes.SolidVolumes.CylindricalThinWall_L4
                                 wall3(redeclare replaceable model Material = TILMedia.SolidTypes.TILMedia_Steel,
    N_ax=tube3.N_cv,
    Delta_x=tube3.Delta_x,
    diameter_o=d_a,
    diameter_i=d_i,
    length=L,
    T_start=ones(nNodes)*(t_startWallEvap))
         annotation (Placement(transformation(extent={{54,158},{74,166}})));
  ClaRa.Components.VolumesValvesFittings.Pipes.PipeFlowVLE_L4_Simple tube3(
    frictionAtOutlet=true,
    Delta_p_nom=NOM.tube3.Delta_p_nom,
    showData=true,
    redeclare model HeatTransfer = ClaRa.Basics.ControlVolumes.Fundamentals.HeatTransport.VLE_HT.NusseltPipe_L4,
    m_flow_nom=NOM.tube3.m_flow,
    h_start=linspace(
        INIT.tube3.h_in,
        INIT.tube3.h_out,
        nNodes),
    p_start=linspace(
        INIT.tube3.p_in,
        INIT.tube3.p_out,
        nNodes),
    redeclare model PressureLoss = ClaRa.Basics.ControlVolumes.Fundamentals.PressureLoss.Generic_PL.LinearPressureLoss_L4,
    length=L,
    diameter_i=d_i,
    N_cv=nNodes) annotation (Placement(transformation(
        extent={{14,-5},{-14,5}},
        rotation=0,
        origin={64,150})));
  ClaRa.Basics.ControlVolumes.SolidVolumes.CylindricalThinWall_L4
                                 wall4(redeclare replaceable model Material = TILMedia.SolidTypes.TILMedia_Steel,
    N_ax=tube4.N_cv,
    Delta_x=tube4.Delta_x,
    diameter_o=d_a,
    diameter_i=d_i,
    length=L,
    T_start=ones(nNodes)*(t_startWallEvap))
         annotation (Placement(transformation(extent={{54,218},{74,226}})));
  ClaRa.Components.VolumesValvesFittings.Pipes.PipeFlowVLE_L4_Simple tube4(
    frictionAtOutlet=true,
    Delta_p_nom=NOM.tube4.Delta_p_nom,
    showData=true,
    redeclare model HeatTransfer = ClaRa.Basics.ControlVolumes.Fundamentals.HeatTransport.VLE_HT.NusseltPipe_L4,
    m_flow_nom=NOM.tube4.m_flow,
    h_start=linspace(
        INIT.tube4.h_in,
        INIT.tube4.h_out,
        nNodes),
    p_start=linspace(
        INIT.tube4.p_in,
        INIT.tube4.p_out,
        nNodes),
    redeclare model PressureLoss = ClaRa.Basics.ControlVolumes.Fundamentals.PressureLoss.Generic_PL.LinearPressureLoss_L4,
    length=L,
    diameter_i=d_i,
    N_cv=nNodes) annotation (Placement(transformation(
        extent={{14,-5},{-14,5}},
        rotation=0,
        origin={64,210})));
  ClaRa.Basics.ControlVolumes.SolidVolumes.CylindricalThinWall_L4
                                 wall5(redeclare replaceable model Material = TILMedia.SolidTypes.TILMedia_Steel,
    N_ax=tube5.N_cv,
    Delta_x=tube5.Delta_x,
    diameter_o=d_a,
    diameter_i=d_i,
    length=L,
    T_start=ones(nNodes)*(t_startWallEvap))
         annotation (Placement(transformation(extent={{54,278},{74,286}})));
  ClaRa.Components.VolumesValvesFittings.Pipes.PipeFlowVLE_L4_Simple tube5(
    frictionAtOutlet=true,
    Delta_p_nom=NOM.tube5.Delta_p_nom,
    showData=true,
    redeclare model HeatTransfer = ClaRa.Basics.ControlVolumes.Fundamentals.HeatTransport.VLE_HT.NusseltPipe_L4,
    m_flow_nom=NOM.tube5.m_flow,
    h_start=linspace(
        INIT.tube5.h_in,
        INIT.tube5.h_out,
        nNodes),
    p_start=linspace(
        INIT.tube5.p_in,
        INIT.tube5.p_out,
        nNodes),
    redeclare model PressureLoss = ClaRa.Basics.ControlVolumes.Fundamentals.PressureLoss.Generic_PL.LinearPressureLoss_L4,
    length=L,
    diameter_i=d_i,
    N_cv=nNodes) annotation (Placement(transformation(
        extent={{14,-5},{-14,5}},
        rotation=0,
        origin={64,270})));
  ClaRa.Components.VolumesValvesFittings.Pipes.PipeFlowVLE_L4_Simple TubeEvapCollOut1(
    frictionAtInlet=true,
    frictionAtOutlet=true,
    Delta_p_nom=NOM.TubeEvapCollOut1.Delta_p_nom,
    length=52.8,
    diameter_i=d_i,
    N_cv=nCellsHead,
    p_nom=ones(nCellsHead)*1e5,
    h_nom=ones(nCellsHead)*1e5,
    h_start=ones(nCellsHead)*INIT.TubeEvapCollOut1.h_in,
    p_start=linspace(
        INIT.TubeEvapCollOut1.p_in,
        INIT.TubeEvapCollOut1.p_out,
        nCellsHead),
    redeclare model PressureLoss = ClaRa.Basics.ControlVolumes.Fundamentals.PressureLoss.Generic_PL.LinearPressureLoss_L4,
    m_flow_nom=NOM.TubeEvapCollOut1.m_flow)
                       annotation (Placement(transformation(
        extent={{14,-5},{-14,5}},
        rotation=90,
        origin={-36,0})));
  ClaRa.Components.VolumesValvesFittings.Pipes.PipeFlowVLE_L4_Simple TubeEvapCollOut2(
    frictionAtInlet=true,
    frictionAtOutlet=true,
    Delta_p_nom=NOM.TubeEvapCollOut2.Delta_p_nom,
    length=52.8,
    diameter_i=d_i,
    N_cv=nCellsHead,
    p_nom=ones(nCellsHead)*1e5,
    h_nom=ones(nCellsHead)*1e5,
    h_start=ones(nCellsHead)*INIT.TubeEvapCollOut2.h_in,
    p_start=linspace(
        INIT.TubeEvapCollOut2.p_in,
        INIT.TubeEvapCollOut2.p_out,
        nCellsHead),
    redeclare model PressureLoss = ClaRa.Basics.ControlVolumes.Fundamentals.PressureLoss.Generic_PL.LinearPressureLoss_L4,
    m_flow_nom=NOM.TubeEvapCollOut2.m_flow)
                       annotation (Placement(transformation(
        extent={{14,-5},{-14,5}},
        rotation=90,
        origin={-36,60})));
  ClaRa.Components.VolumesValvesFittings.Pipes.PipeFlowVLE_L4_Simple TubeEvapCollOut3(
    frictionAtInlet=true,
    frictionAtOutlet=true,
    Delta_p_nom=NOM.TubeEvapCollOut3.Delta_p_nom,
    length=52.8,
    diameter_i=d_i,
    N_cv=nCellsHead,
    p_nom=ones(nCellsHead)*1e5,
    h_nom=ones(nCellsHead)*1e5,
    h_start=ones(nCellsHead)*INIT.TubeEvapCollOut3.h_in,
    p_start=linspace(
        INIT.TubeEvapCollOut3.p_in,
        INIT.TubeEvapCollOut3.p_out,
        nCellsHead),
    redeclare model PressureLoss = ClaRa.Basics.ControlVolumes.Fundamentals.PressureLoss.Generic_PL.LinearPressureLoss_L4,
    m_flow_nom=NOM.TubeEvapCollOut3.m_flow)
                       annotation (Placement(transformation(
        extent={{14,-5},{-14,5}},
        rotation=90,
        origin={-36,120})));
  ClaRa.Components.VolumesValvesFittings.Pipes.PipeFlowVLE_L4_Simple TubeEvapCollOut4(
    frictionAtInlet=true,
    frictionAtOutlet=true,
    Delta_p_nom=NOM.TubeEvapCollOut4.Delta_p_nom,
    length=52.8,
    diameter_i=d_i,
    N_cv=nCellsHead,
    p_nom=ones(nCellsHead)*1e5,
    h_nom=ones(nCellsHead)*1e5,
    h_start=ones(nCellsHead)*INIT.TubeEvapCollOut4.h_in,
    p_start=ones(nCellsHead)*(p_startEvapIn - p_startDeltaEvapCollInTotal/5*4),
    redeclare model PressureLoss = ClaRa.Basics.ControlVolumes.Fundamentals.PressureLoss.Generic_PL.LinearPressureLoss_L4,
    m_flow_nom=NOM.TubeEvapCollOut4.m_flow)
                       annotation (Placement(transformation(
        extent={{14,-5},{-14,5}},
        rotation=90,
        origin={-36,180})));
  ClaRa.Components.VolumesValvesFittings.Pipes.PipeFlowVLE_L4_Simple TubeEvapCollOut5(
    frictionAtOutlet=true,
    Delta_p_nom=NOM.TubeEvapCollOut5.Delta_p_nom,
    length=52.8,
    diameter_i=d_i,
    N_cv=nCellsHead,
    p_nom=ones(nCellsHead)*1e5,
    h_nom=ones(nCellsHead)*1e5,
    h_start=ones(nCellsHead)*INIT.TubeEvapCollOut5.h_in,
    p_start=linspace(
        INIT.TubeEvapCollOut5.p_in,
        INIT.TubeEvapCollOut5.p_out,
        nCellsHead),
    redeclare model PressureLoss = ClaRa.Basics.ControlVolumes.Fundamentals.PressureLoss.Generic_PL.LinearPressureLoss_L4,
    m_flow_nom=NOM.TubeEvapCollOut5.m_flow)
                       annotation (Placement(transformation(
        extent={{14,-5},{-14,5}},
        rotation=90,
        origin={-36,240})));
  ClaRa.Components.VolumesValvesFittings.Pipes.PipeFlowVLE_L4_Simple TubeInSep(
    frictionAtInlet=true,
    frictionAtOutlet=false,
    redeclare model PressureLoss = ClaRa.Basics.ControlVolumes.Fundamentals.PressureLoss.Generic_PL.LinearPressureLoss_L4,
    m_flow_nom=NOM.TubeInSep.m_flow,
    Delta_p_nom=NOM.TubeInSep.Delta_p_nom,
    h_start=ones(nCells)*INIT.TubeInSep.h_in,
    p_start=linspace(
        INIT.TubeInSep.p_in,
        INIT.TubeInSep.p_out,
        nCells),
    showData=true,
    diameter_i=d_i,
    N_cv=nCells) annotation (Placement(transformation(extent={{-50,-35},{-78,-25}})));
  ClaRa.Components.VolumesValvesFittings.Pipes.PipeFlowVLE_L4_Simple TubeSHIN(
    frictionAtInlet=true,
    frictionAtOutlet=false,
    redeclare model PressureLoss = ClaRa.Basics.ControlVolumes.Fundamentals.PressureLoss.Generic_PL.LinearPressureLoss_L4,
    m_flow_nom=NOM.TubeInSH.m_flow,
    Delta_p_nom=NOM.TubeInSH.Delta_p_nom,
    h_start=ones(6)*INIT.TubeInSH.h_in,
    p_start=linspace(
        INIT.TubeInSH.p_in,
        INIT.TubeInSH.p_out,
        6),
    length=240,
    diameter_i=d_i,
    N_cv=6) annotation (Placement(transformation(extent={{-120,-35},{-148,-25}})));
  ClaRa.Visualisation.DynamicBar dynamicBar(provideInputConnectors=true) annotation (Placement(transformation(extent={{-92,-38},{-82,-18}})));
  ClaRa.Components.VolumesValvesFittings.Pipes.PipeFlowVLE_L4_Simple TubeRePumpIn(
    frictionAtInlet=true,
    frictionAtOutlet=false,
    redeclare model PressureLoss = ClaRa.Basics.ControlVolumes.Fundamentals.PressureLoss.Generic_PL.LinearPressureLoss_L4,
    m_flow_nom=NOM.TubeRePumpIn.m_flow,
    Delta_p_nom=NOM.TubeRePumpIn.Delta_p_nom,
    h_start=linspace(
        INIT.TubeRePumpIn.h_in,
        INIT.TubeRePumpIn.h_in,
        nCells),
    p_start=linspace(
        INIT.TubeRePumpIn.p_in,
        INIT.TubeRePumpIn.p_out,
        nCells),
    showData=true,
    N_cv=nCells) annotation (Placement(transformation(extent={{14,5},{-14,-5}},
        rotation=180,
        origin={-20,-68})));
  ClaRa.Components.VolumesValvesFittings.Pipes.PipeFlowVLE_L4_Simple TubeSHCollIn1(
    frictionAtInlet=true,
    Delta_p_nom=NOM.TubeSHCollIn1.Delta_p_nom,
    length=77/2,
    diameter_i=d_i,
    N_cv=nCellsHead,
    p_nom=ones(nCellsHead)*1e5,
    h_nom=ones(nCellsHead)*1e5,
    h_start=ones(nCellsHead)*INIT.TubeSHCollIn1.h_in,
    p_start=linspace(
        INIT.TubeSHCollIn1.p_in,
        INIT.TubeSHCollIn1.p_out,
        nCellsHead),
    redeclare model PressureLoss = ClaRa.Basics.ControlVolumes.Fundamentals.PressureLoss.Generic_PL.LinearPressureLoss_L4,
    m_flow_nom=NOM.TubeSHCollIn1.m_flow)
                      annotation (Placement(transformation(
        extent={{14,-5},{-14,5}},
        rotation=270,
        origin={-170,6})));
  ClaRa.Components.VolumesValvesFittings.Pipes.PipeFlowVLE_L4_Simple TubeSHCollIn2(
    frictionAtInlet=true,
    Delta_p_nom=NOM.TubeSHCollIn2.Delta_p_nom,
    showData=true,
    length=77/2,
    diameter_i=d_i,
    N_cv=nCellsHead,
    p_nom=ones(nCellsHead)*1e5,
    h_nom=ones(nCellsHead)*1e5,
    h_start=ones(nCellsHead)*INIT.TubeSHCollIn2.h_in,
    p_start=linspace(
        INIT.TubeSHCollIn2.p_in,
        INIT.TubeSHCollIn2.p_out,
        nCellsHead),
    redeclare model PressureLoss = ClaRa.Basics.ControlVolumes.Fundamentals.PressureLoss.Generic_PL.LinearPressureLoss_L4,
    m_flow_nom=NOM.TubeSHCollIn2.m_flow)
                    annotation (Placement(transformation(
        extent={{14,-5},{-14,5}},
        rotation=270,
        origin={-170,86})));
  ClaRa.Components.VolumesValvesFittings.Pipes.PipeFlowVLE_L4_Simple TubeSHCollOut1(
    frictionAtInlet=true,
    frictionAtOutlet=true,
    Delta_p_nom=NOM.TubeSHCollOut1.Delta_p_nom,
    length=77/2,
    diameter_i=d_i,
    N_cv=nCellsHead,
    p_nom=ones(nCellsHead)*1e5,
    h_nom=ones(nCellsHead)*1e5,
    h_start=ones(nCellsHead)*INIT.TubeSHCollOut1.h_in,
    p_start=linspace(
        INIT.TubeSHCollOut1.p_in,
        INIT.TubeSHCollOut1.p_out,
        nCellsHead),
    redeclare model PressureLoss = ClaRa.Basics.ControlVolumes.Fundamentals.PressureLoss.Generic_PL.LinearPressureLoss_L4,
    m_flow_nom=NOM.TubeSHCollOut1.m_flow)
                      annotation (Placement(transformation(
        extent={{14,-5},{-14,5}},
        rotation=90,
        origin={-320,6})));
  ClaRa.Components.VolumesValvesFittings.Pipes.PipeFlowVLE_L4_Simple TubeSHCollOut2(
    frictionAtOutlet=true,
    Delta_p_nom=NOM.TubeSHCollOut2.Delta_p_nom,
    length=77/2,
    diameter_i=d_i,
    N_cv=nCellsHead,
    p_nom=ones(nCellsHead)*1e5,
    h_nom=ones(nCellsHead)*1e5,
    h_start=ones(nCellsHead)*INIT.TubeSHCollOut2.h_in,
    p_start=linspace(
        INIT.TubeSHCollOut2.p_in,
        INIT.TubeSHCollOut2.p_out,
        nCellsHead),
    redeclare model PressureLoss = ClaRa.Basics.ControlVolumes.Fundamentals.PressureLoss.Generic_PL.LinearPressureLoss_L4,
    m_flow_nom=NOM.TubeSHCollOut2.m_flow)
                    annotation (Placement(transformation(
        extent={{14,-5},{-14,5}},
        rotation=90,
        origin={-320,86})));
  ClaRa.Components.VolumesValvesFittings.Pipes.PipeFlowVLE_L4_Simple tubeSH(
    frictionAtInlet=true,
    frictionAtOutlet=true,
    Delta_p_nom=NOM.tubeSH.Delta_p_nom,
    showData=true,
    redeclare model HeatTransfer = ClaRa.Basics.ControlVolumes.Fundamentals.HeatTransport.VLE_HT.NusseltPipe_L4,
    m_flow_nom=NOM.tubeSH.m_flow,
    h_start=linspace(
        INIT.tubeSH.h_in,
        INIT.tubeSH.h_out,
        nNodesSH),
    p_start=linspace(
        INIT.tubeSH.p_in,
        INIT.tubeSH.p_out,
        nNodesSH),
    redeclare model PressureLoss = ClaRa.Basics.ControlVolumes.Fundamentals.PressureLoss.Generic_PL.LinearPressureLoss_L4,
    length=LSH,
    diameter_i=d_iSH,
    N_cv=nNodesSH) annotation (Placement(transformation(
        extent={{14,-5},{-14,5}},
        rotation=0,
        origin={-230,-30})));
  ClaRa.Basics.ControlVolumes.SolidVolumes.CylindricalThinWall_L4
                                                  wallSH1(
    redeclare replaceable model Material = TILMedia.SolidTypes.TILMedia_Steel,
    N_ax=nNodesSH,
    Delta_x=tubeSH.Delta_x,
    diameter_o=d_aSH,
    diameter_i=d_iSH,
    length=LSH,
    T_start=ones(nNodesSH)*(t_startWallSH))
         annotation (Placement(transformation(extent={{-240,48},{-220,56}})));
  ClaRa.Components.VolumesValvesFittings.Pipes.PipeFlowVLE_L4_Simple tubeSH1(
    frictionAtInlet=true,
    frictionAtOutlet=true,
    Delta_p_nom=NOM.tubeSH1.Delta_p_nom,
    showData=true,
    redeclare model HeatTransfer = ClaRa.Basics.ControlVolumes.Fundamentals.HeatTransport.VLE_HT.NusseltPipe_L4,
    m_flow_nom=NOM.tubeSH1.m_flow,
    h_start=linspace(
        INIT.tubeSH1.h_in,
        INIT.tubeSH1.h_out,
        nNodesSH),
    p_start=linspace(
        INIT.tubeSH1.p_in,
        INIT.tubeSH1.p_out,
        nNodesSH),
    redeclare model PressureLoss = ClaRa.Basics.ControlVolumes.Fundamentals.PressureLoss.Generic_PL.LinearPressureLoss_L4,
    length=LSH,
    diameter_i=d_iSH,
    N_cv=nNodesSH) annotation (Placement(transformation(
        extent={{14,-5},{-14,5}},
        rotation=0,
        origin={-230,40})));
  Components.Absorber.Parabolic.ParabolSolarliteDetailled parabolSolarliteDetailled8(
    W=W,
    L=LSH,
    n_ax=nNodesSH,
    dx=tubeSH.Delta_x,
    n_Qrad_ax=1,
    Zs=0) annotation (Placement(transformation(extent={{-268,50},{-248,70}})));
  ClaRa.Basics.ControlVolumes.SolidVolumes.CylindricalThinWall_L4
                                                  wallSH2(
    redeclare replaceable model Material = TILMedia.SolidTypes.TILMedia_Steel,
    N_ax=nNodesSH,
    Delta_x=tubeSH.Delta_x,
    diameter_o=d_aSH,
    diameter_i=d_iSH,
    length=LSH,
    T_start=ones(nNodesSH)*(t_startWallSH))
         annotation (Placement(transformation(extent={{-240,118},{-220,126}})));
  ClaRa.Components.VolumesValvesFittings.Pipes.PipeFlowVLE_L4_Simple tubeSH2(
    frictionAtInlet=true,
    frictionAtOutlet=true,
    Delta_p_nom=NOM.tubeSH2.Delta_p_nom,
    showData=true,
    redeclare model HeatTransfer = ClaRa.Basics.ControlVolumes.Fundamentals.HeatTransport.VLE_HT.NusseltPipe_L4,
    m_flow_nom=NOM.tubeSH2.m_flow,
    h_start=linspace(
        INIT.tubeSH2.h_in,
        INIT.tubeSH2.h_out,
        nNodesSH),
    p_start=linspace(
        INIT.tubeSH2.p_in,
        INIT.tubeSH2.p_out,
        nNodesSH),
    redeclare model PressureLoss = ClaRa.Basics.ControlVolumes.Fundamentals.PressureLoss.Generic_PL.LinearPressureLoss_L4,
    length=LSH,
    diameter_i=d_iSH,
    N_cv=nNodesSH) annotation (Placement(transformation(
        extent={{14,-5},{-14,5}},
        rotation=0,
        origin={-230,110})));
  StaticCycles.Static_Cycle_SixBranchEvapThreeBranchSH_Turbine
                                                       INIT annotation (Placement(transformation(extent={{280,360},{300,380}})));
  StaticCycles.Static_Cycle_SixBranchEvapThreeBranchSH_Turbine
                                                       NOM annotation (Placement(transformation(extent={{320,360},{340,380}})));
  ClaRa.Visualisation.Quadruple quadruple1 annotation (Placement(transformation(extent={{-242,92},{-218,104}})));
  ClaRa.Visualisation.Quadruple quadruple2 annotation (Placement(transformation(extent={{-242,22},{-218,34}})));
  ClaRa.Visualisation.Quadruple quadruple3 annotation (Placement(transformation(extent={{-244,-50},{-220,-38}})));
  ClaRa.Visualisation.Quadruple quadruple4 annotation (Placement(transformation(extent={{-216,-50},{-192,-38}})));
  ClaRa.Visualisation.Quadruple quadruple5 annotation (Placement(transformation(extent={{-202,20},{-178,32}})));
  ClaRa.Visualisation.Quadruple quadruple6 annotation (Placement(transformation(extent={{-204,92},{-180,104}})));
  ClaRa.Visualisation.Quadruple quadruple7 annotation (Placement(transformation(extent={{-108,-14},{-84,-2}})));
  ClaRa.Visualisation.Quadruple quadruple8 annotation (Placement(transformation(extent={{-90,-58},{-66,-46}})));
  ClaRa.Visualisation.Quadruple quadruple9 annotation (Placement(transformation(extent={{-374,36},{-350,48}})));
  ClaRa.Visualisation.Quadruple quadruple10 annotation (Placement(transformation(extent={{-74,-14},{-50,-2}})));
  ClaRa.Visualisation.Quadruple quadruple12 annotation (Placement(transformation(extent={{122,-86},{146,-74}})));
  ClaRa.Visualisation.Quadruple quadruple14 annotation (Placement(transformation(extent={{198,-50},{222,-38}})));
  ClaRa.Visualisation.Quadruple quadruple15 annotation (Placement(transformation(extent={{102,-58},{126,-46}})));
  ClaRa.Visualisation.Quadruple quadruple16 annotation (Placement(transformation(extent={{94,2},{118,14}})));
  ClaRa.Visualisation.Quadruple quadruple17 annotation (Placement(transformation(extent={{94,62},{118,74}})));
  ClaRa.Visualisation.Quadruple quadruple18 annotation (Placement(transformation(extent={{96,120},{120,132}})));
  ClaRa.Visualisation.Quadruple quadruple19 annotation (Placement(transformation(extent={{96,182},{120,194}})));
  ClaRa.Visualisation.Quadruple quadruple20 annotation (Placement(transformation(extent={{96,242},{120,254}})));
  ClaRa.Visualisation.Quadruple quadruple21 annotation (Placement(transformation(extent={{140,274},{164,286}})));
  ClaRa.Visualisation.Quadruple quadruple22 annotation (Placement(transformation(extent={{130,162},{154,174}})));
  ClaRa.Visualisation.Quadruple quadruple23 annotation (Placement(transformation(extent={{132,222},{156,234}})));
  ClaRa.Visualisation.Quadruple quadruple24 annotation (Placement(transformation(extent={{132,104},{156,116}})));
  ClaRa.Visualisation.Quadruple quadruple25 annotation (Placement(transformation(extent={{132,42},{156,54}})));
  ClaRa.Visualisation.Quadruple quadruple26 annotation (Placement(transformation(extent={{134,-12},{158,0}})));
  ClaRa.Visualisation.Quadruple quadruple27 annotation (Placement(transformation(extent={{50,-50},{74,-38}})));
  ClaRa.Visualisation.Quadruple quadruple28 annotation (Placement(transformation(extent={{52,10},{76,22}})));
  ClaRa.Visualisation.Quadruple quadruple29 annotation (Placement(transformation(extent={{52,68},{76,80}})));
  ClaRa.Visualisation.Quadruple quadruple30 annotation (Placement(transformation(extent={{52,130},{76,142}})));
  ClaRa.Visualisation.Quadruple quadruple31 annotation (Placement(transformation(extent={{50,190},{74,202}})));
  ClaRa.Visualisation.Quadruple quadruple32 annotation (Placement(transformation(extent={{52,252},{76,264}})));
  ClaRa.Visualisation.Quadruple quadruple33 annotation (Placement(transformation(extent={{-76,206},{-52,218}})));
  ClaRa.Visualisation.Quadruple quadruple34 annotation (Placement(transformation(extent={{-78,148},{-54,160}})));
  ClaRa.Visualisation.Quadruple quadruple35 annotation (Placement(transformation(extent={{-76,90},{-52,102}})));
  ClaRa.Visualisation.Quadruple quadruple36 annotation (Placement(transformation(extent={{-78,26},{-54,38}})));
  ClaRa.Visualisation.Quadruple quadruple37 annotation (Placement(transformation(extent={{-42,-54},{-18,-42}})));
  ClaRa.Components.TurboMachines.Pumps.PumpVLE_L1_simple pumpRe(eta_mech=NOM.pumpRe.efficiency, m_flow_start=INIT.pumpRe.m_flow)
                    annotation (Placement(transformation(extent={{80,-78},{100,-58}})));
  Modelica.Blocks.Sources.RealExpression realExpression15(y=min(valveBranch.pressureLoss.Delta_p, min(valveBranch1.pressureLoss.Delta_p, min(valveBranch2.pressureLoss.Delta_p, min(valveBranch3.pressureLoss.Delta_p, min(valveBranch5.pressureLoss.Delta_p, valveBranch4.pressureLoss.Delta_p))))))
    annotation (Placement(transformation(extent={{-182,-86},{-162,-66}})));
  Modelica.Blocks.Sources.RealExpression realExpression16(y=0.5e5)
    annotation (Placement(transformation(extent={{-182,-66},{-162,-46}})));
  ClaRa.Components.Utilities.Blocks.LimPID PIDRecircPump(
    sign=1,
    y_max=50e3,
    y_min=500,
    k=0.005,
    Tau_i=10,
    Tau_d=1e-3,
    y_start=5000,
    controllerType=Modelica.Blocks.Types.SimpleController.PI) annotation (Placement(transformation(extent={{-136,-66},{-116,-46}})));
  ClaRa.Visualisation.Quadruple quadruple11 annotation (Placement(transformation(extent={{4,-86},{28,-74}})));
  ClaRa.Visualisation.Quadruple quadruple38 annotation (Placement(transformation(extent={{224,-86},{248,-74}})));
  ClaRa.Components.TurboMachines.Pumps.PumpVLE_L1_simple pumpMain(eta_mech=NOM.pumpMain.efficiency, m_flow_start=INIT.pumpMain.m_flow)
                    annotation (Placement(transformation(extent={{334,-154},{354,-134}})));
  Modelica.Blocks.Sources.RealExpression realExpression1(y=0.2e5)
    annotation (Placement(transformation(extent={{268,-110},{288,-90}})));
  Modelica.Blocks.Sources.RealExpression realExpression27(y=valveMain.pressureLoss.Delta_p)
    annotation (Placement(transformation(extent={{268,-130},{288,-110}})));
  ClaRa.Components.Utilities.Blocks.LimPID PIDMainPump(
    controllerType=Modelica.Blocks.Types.SimpleController.PI,
    y_max=1000e3,
    y_min=100,
    Tau_i=15,
    Tau_d=1e-3,
    y_start=50e3,
    sign=1)  annotation (Placement(transformation(extent={{308,-110},{328,-90}})));
  Modelica.Blocks.Sources.RealExpression realExpression26(y=seperator.level - (TubeIn.outlet.m_flow - join_Y6.outlet.m_flow))
    annotation (Placement(transformation(extent={{280,18},{300,38}})));
  Modelica.Blocks.Sources.RealExpression realExpression25(y=0.5)
    annotation (Placement(transformation(extent={{280,40},{300,60}})));
  ClaRa.Components.Utilities.Blocks.LimPID PIDMainValve(
    sign=1,
    y_max=1,
    y_min=0.01,
    Tau_i=2000,
    Tau_d=500,
    limitsAtInit=true,
    controllerType=Modelica.Blocks.Types.SimpleController.PI,
    y_start=0.8,
    k=2) annotation (Placement(transformation(extent={{320,40},{340,60}})));
  ClaRa.Components.VolumesValvesFittings.Valves.GenericValveVLE_L1   valveMain(redeclare model PressureLoss = ClaRa.Components.VolumesValvesFittings.Valves.Fundamentals.LinearNominalPoint (Delta_p_nom=NOM.valveMain.Delta_p_nom, m_flow_nom=NOM.valveMain.m_flow), openingInputIsActive=true)
    annotation (Placement(transformation(extent={{356,-36},{336,-24}})));
  ClaRa.Components.VolumesValvesFittings.Pipes.PipeFlowVLE_L4_Simple TubeMainValveIn(
    frictionAtOutlet=false,
    redeclare model PressureLoss = ClaRa.Basics.ControlVolumes.Fundamentals.PressureLoss.Generic_PL.LinearPressureLoss_L4,
    m_flow_nom=NOM.TubeMainValveIn.m_flow,
    Delta_p_nom=NOM.TubeMainValveIn.Delta_p_nom,
    h_start=ones(nCells)*INIT.TubeMainValveIn.h_in,
    p_start=linspace(
        INIT.TubeMainValveIn.p_in,
        INIT.TubeMainValveIn.p_out,
        nCells),
    showData=true,
    length=10,
    N_cv=nCells) annotation (Placement(transformation(extent={{414,-35},{386,-25}})));
  ClaRa.Components.VolumesValvesFittings.Pipes.PipeFlowVLE_L4_Simple TubeMainPumpIn(
    frictionAtInlet=true,
    frictionAtOutlet=false,
    redeclare model PressureLoss = ClaRa.Basics.ControlVolumes.Fundamentals.PressureLoss.Generic_PL.LinearPressureLoss_L4,
    m_flow_nom=NOM.TubeMainPumpIn.m_flow,
    Delta_p_nom=NOM.TubeMainPumpIn.Delta_p_nom,
    h_start=ones(nCells)*INIT.TubeMainPumpIn.h_in,
    p_start=linspace(
        INIT.TubeMainPumpIn.p_in,
        INIT.TubeMainPumpIn.p_out,
        nCells),
    showData=true,
    length=10,
    N_cv=nCells) annotation (Placement(transformation(extent={{166,-149},{194,-139}})));
  ClaRa.Visualisation.Quadruple quadruple13 annotation (Placement(transformation(extent={{268,-60},{292,-48}})));
  ClaRa.Visualisation.Quadruple quadruple39 annotation (Placement(transformation(extent={{338,-60},{362,-48}})));
  ClaRa.Visualisation.Quadruple quadruple40 annotation (Placement(transformation(extent={{388,-60},{412,-48}})));
  ClaRa.Visualisation.Quadruple quadruple42 annotation (Placement(transformation(extent={{316,-172},{340,-160}})));
  ClaRa.Visualisation.Quadruple quadruple43 annotation (Placement(transformation(extent={{366,-172},{390,-160}})));
  ClaRa.Components.VolumesValvesFittings.Valves.GenericValveVLE_L1 valve_turbine(redeclare model PressureLoss = ClaRa.Components.VolumesValvesFittings.Valves.Fundamentals.LinearNominalPoint (Delta_p_nom=1e3), openingInputIsActive=true) annotation (Placement(transformation(extent={{-344,-36},{-364,-24}})));
  ClaRa.Components.TurboMachines.Turbines.SteamTurbineVLE_L1 turbine(
    p_nom=NOM.turbine.p_in,
    m_flow_nom=NOM.turbine.m_flow,
    Pi=NOM.turbine.p_out/NOM.turbine.p_in,
    rho_nom=NOM.turbine.rho_in) annotation (Placement(transformation(extent={{-400,-100},{-390,-80}})));
  ClaRa.Components.TurboMachines.Pumps.PumpVLE_L1_simple Pump_cond(
    showExpertSummary=true,
    eta_mech=NOM.Pump_cond.efficiency,
    m_flow_start=INIT.Pump_cond.summary.inlet.m_flow)                                                  annotation (Placement(transformation(extent={{-134,-112},{-114,-132}})));
  ClaRa.Components.Utilities.Blocks.LimPID PI_Pump_cond(
    sign=-1,
    y_ref=1e6,
    Tau_d=30,
    controllerType=Modelica.Blocks.Types.SimpleController.PI,
    y_max=NOM.Pump_cond.P_pump*10,
    y_min=NOM.Pump_cond.P_pump/200,
    y_start=INIT.Pump_cond.P_pump,
    k=10,
    Tau_i=100,
    initOption=796) annotation (Placement(transformation(extent={{-180,-160},{-160,-140}})));
  Modelica.Blocks.Sources.RealExpression setPoint_condenser(y=0.5/6) annotation (Placement(transformation(extent={{-216,-156},{-202,-144}})));
  Modelica.Blocks.Continuous.FirstOrder measurement(
    initType=Modelica.Blocks.Types.Init.InitialOutput,
    y_start=0.1,
    T=10)
    annotation (Placement(transformation(extent={{-214,-174},{-206,-166}})));
  ClaRa.Components.BoundaryConditions.BoundaryVLE_Txim_flow
                                                      boundaryVLE_Txim_flow(T_const=273.15 + 15, m_flow_const=25000) annotation (Placement(transformation(extent={{-360,-160},{-340,-140}})));
  ClaRa.Visualisation.DynamicBar
                           fillingLevel_condenser(
    u_set=0.5/6,
    u_high=0.5/3,
    u_low=0.5/12,
    provideInputConnectors=true)
                  annotation (Placement(transformation(extent={{-278,-144},{-288,-124}})));
  ClaRa.Components.HeatExchangers.HEXvle2vle_L3_2ph_BU_simple
                                             condenser(
    height=5,
    width=5,
    redeclare model PressureLossShell = ClaRa.Basics.ControlVolumes.Fundamentals.PressureLoss.Generic_PL.LinearParallelZones_L3 (Delta_p_nom={100,100,100}),
    z_in_shell=4.9,
    z_out_shell=0.1,
    level_rel_start=0.5/6,
    m_flow_nom_shell=NOM.condenser.m_flow_in,
    p_nom_shell=NOM.condenser.p_condenser,
    p_start_shell=INIT.condenser.p_condenser,
    initOptionShell=204,
    levelOutput=true,
    z_in_aux1=4.9,
    z_in_aux2=4.9,
    height_hotwell=2,
    width_hotwell=1,
    length_hotwell=10,
    diameter_i=0.008,
    diameter_o=0.01,
    m_flow_nom_tubes=10000,
    p_nom_tubes=2e5,
    h_start_tubes=85e3,
    p_start_tubes=2e5,
    redeclare model HeatTransferTubes = ClaRa.Basics.ControlVolumes.Fundamentals.HeatTransport.VLE_HT.NusseltPipe1ph_L2,
    length=12,
    redeclare model HeatTransfer_Shell = ClaRa.Basics.ControlVolumes.Fundamentals.HeatTransport.VLE_HT.Constant_L3_ypsDependent (alpha_nom={3000,12000}),
    N_tubes=15000)       annotation (Placement(transformation(extent={{-286,-146},{-306,-126}})));
  ClaRa.Components.VolumesValvesFittings.Valves.GenericValveVLE_L1
                                                             valveControl_preheater_LP1(redeclare model PressureLoss = ClaRa.Components.VolumesValvesFittings.Valves.Fundamentals.LinearNominalPoint (
        CL_valve=[0,0; 1,1],
        Delta_p_nom=1000,
        m_flow_nom=25000)) annotation (Placement(transformation(
        extent={{10,-6},{-10,6}},
        rotation=0,
        origin={-322,-130})));
  ClaRa.Components.BoundaryConditions.BoundaryVLE_phxi
                                                 boundaryVLE_phxi(p_const=2e5) annotation (Placement(transformation(extent={{-360,-140},{-340,-120}})));
  ClaRa.Visualisation.Quadruple quadruple41
                                           annotation (Placement(transformation(extent={{-328,-52},{-304,-40}})));
  ClaRa.Visualisation.Quadruple quadruple44
                                           annotation (Placement(transformation(extent={{-392,-52},{-368,-40}})));
  ClaRa.Visualisation.Quadruple quadruple45
                                           annotation (Placement(transformation(extent={{-376,-96},{-352,-84}})));
  ClaRa.Visualisation.Quadruple quadruple46
                                           annotation (Placement(transformation(extent={{-288,-112},{-264,-100}})));
  ClaRa.Visualisation.Quadruple quadruple47
                                           annotation (Placement(transformation(extent={{-108,-118},{-84,-106}})));
  ClaRa.Visualisation.Quadruple quadruple48
                                           annotation (Placement(transformation(extent={{-20,-112},{4,-100}})));
  ClaRa.Components.Utilities.Blocks.LimPID PID_TurbineValve(sign=1, y_min=0.9) annotation (Placement(transformation(extent={{-380,0},{-360,20}})));
  Modelica.Blocks.Sources.RealExpression realExpression4(y=30e5)
    annotation (Placement(transformation(extent={{-420,0},{-400,20}})));
  Modelica.Blocks.Sources.RealExpression realExpression17(y=turbine.p_in)
    annotation (Placement(transformation(extent={{-420,-22},{-400,-2}})));
  ClaRa.Components.VolumesValvesFittings.Fittings.SplitVLE_L2_Y splitTurbine(
    m_flow_out_nom={NOM.turbine.m_flow,NOM.valve_SteamFWT.m_flow},
    p_nom=NOM.split_turbine.p,
    h_nom=NOM.split_turbine.h1,
    h_start=INIT.split_turbine.h1,
    p_start=INIT.split_turbine.p) annotation (Placement(transformation(
        extent={{-10,10},{10,-10}},
        rotation=270,
        origin={-400,-60})));
  ClaRa.Components.VolumesValvesFittings.Valves.GenericValveVLE_L1 ValveSteamFWT(redeclare model PressureLoss = ClaRa.Components.VolumesValvesFittings.Valves.Fundamentals.LinearNominalPoint (Delta_p_nom=NOM.valve_SteamFWT.Delta_p, m_flow_nom=NOM.valve_SteamFWT.m_flow)) annotation (Placement(transformation(extent={{-230,-102},{-210,-90}})));
  ClaRa.Visualisation.Quadruple quadruple49
                                           annotation (Placement(transformation(extent={{-200,-112},{-176,-100}})));
  ClaRa.Visualisation.Quadruple quadruple50
                                           annotation (Placement(transformation(extent={{-352,-82},{-328,-70}})));
  ClaRa.Visualisation.Quadruple quadruple51
                                           annotation (Placement(transformation(extent={{-384,-82},{-360,-70}})));
  ClaRa.Components.MechanicalSeparation.FeedWaterTank_L3 feedWaterTank(
    level_rel_start=0.5,
    diameter=5,
    orientation=ClaRa.Basics.Choices.GeometryOrientation.horizontal,
    p_start(displayUnit="bar") = INIT.feedwatertank.p_FWT,
    z_tapping=4.5,
    z_vent=4.5,
    z_condensate=4.5,
    redeclare model PressureLoss = ClaRa.Basics.ControlVolumes.Fundamentals.PressureLoss.Generic_PL.LinearParallelZones_L3 (Delta_p_nom={1000,1000,1000}),
    m_flow_cond_nom=NOM.feedwatertank.m_flow_cond,
    p_nom=NOM.feedwatertank.p_FWT,
    h_nom=NOM.feedwatertank.h_cond_in,
    m_flow_heat_nom=NOM.valve_SteamFWT.m_flow,
    initOption=204,
    T_wall_start=ones(feedWaterTank.wall.N_rad)*(120 + 273.15),
    showLevel=true,
    length=12,
    z_aux=1,
    equalPressures=false,
    absorbInflow=0.6) annotation (Placement(transformation(extent={{-20,-138},{-80,-118}})));
  Components.SunModels.Sun sun(
    dayOfYear=dayOfYear,
    startTimeSimulation=startTimeSimulation,
    DS=DS,
    SL=SL,
    L=LLat,
    LL=LLong,
    Green=Green)
                annotation (Placement(transformation(extent={{44,304},{24,324}})));
  Components.SunModels.Sun sun1(
    dayOfYear=dayOfYear,
    startTimeSimulation=startTimeSimulation,
    DS=DS,
    SL=SL,
    L=LLat,
    LL=LLong,
    Green=Green)
                annotation (Placement(transformation(extent={{44,244},{24,264}})));
  Components.SunModels.Sun sun2(
    dayOfYear=dayOfYear,
    startTimeSimulation=startTimeSimulation,
    DS=DS,
    SL=SL,
    L=LLat,
    LL=LLong,
    Green=Green)
                annotation (Placement(transformation(extent={{44,184},{24,204}})));
  Components.SunModels.Sun sun3(
    dayOfYear=dayOfYear,
    startTimeSimulation=startTimeSimulation,
    DS=DS,
    SL=SL,
    L=LLat,
    LL=LLong,
    Green=Green)
                annotation (Placement(transformation(extent={{44,124},{24,144}})));
  Components.SunModels.Sun sun4(
    dayOfYear=dayOfYear,
    startTimeSimulation=startTimeSimulation,
    DS=DS,
    SL=SL,
    L=LLat,
    LL=LLong,
    Green=Green)
                annotation (Placement(transformation(extent={{44,64},{24,84}})));
  Components.SunModels.Sun sun5(
    dayOfYear=dayOfYear,
    startTimeSimulation=startTimeSimulation,
    DS=DS,
    SL=SL,
    L=LLat,
    LL=LLong,
    Green=Green)
                annotation (Placement(transformation(extent={{44,4},{24,24}})));
  Components.SunModels.Sun sun6(
    dayOfYear=dayOfYear,
    startTimeSimulation=startTimeSimulation,
    DS=DS,
    SL=SL,
    L=LLat,
    LL=LLong,
    Green=Green)
                annotation (Placement(transformation(extent={{-226,142},{-246,162}})));
  Components.SunModels.Sun sun7(
    dayOfYear=dayOfYear,
    startTimeSimulation=startTimeSimulation,
    DS=DS,
    SL=SL,
    L=LLat,
    LL=LLong,
    Green=Green)
                annotation (Placement(transformation(extent={{-228,72},{-248,92}})));
  Components.SunModels.Sun sun8(
    dayOfYear=dayOfYear,
    startTimeSimulation=startTimeSimulation,
    DS=DS,
    SL=SL,
    L=LLat,
    LL=LLong,
    Green=Green)
                annotation (Placement(transformation(extent={{-228,2},{-248,22}})));
equation
  connect(const.y, parabolSolarliteDetailled5.Alpha) annotation (Line(
      points={{15,70},{30,70},{30,60}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(const1.y, parabolSolarliteDetailled1.Alpha) annotation (Line(
      points={{15,10},{30,10},{30,0}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(const2.y, parabolSolarliteDetailled7.Alpha) annotation (Line(
      points={{-271,20},{-252,20},{-252,0}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(const3.y, parabolSolarliteDetailled4.Alpha) annotation (Line(
      points={{15,130},{30,130},{30,120}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(const4.y, parabolSolarliteDetailled3.Alpha) annotation (Line(
      points={{15,190},{30,190},{30,180}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(realExpression6.y, PIDValveBranch1.u_m) annotation (Line(
      points={{233,-12},{204,-12},{204,-2},{203.9,-2}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(realExpression2.y, PIDValveBranch1.u_s) annotation (Line(
      points={{233,10},{216,10}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(PIDValveBranch1.y, valveBranch.opening_in) annotation (Line(
      points={{193,10},{124,10},{124,-21}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(const6.y, parabolSolarliteDetailled6.Alpha) annotation (Line(
      points={{-269,160},{-252,160},{-252,140}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(const7.y, parabolSolarliteDetailled2.Alpha) annotation (Line(
      points={{15,250},{30,250},{30,240}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(const8.y,parabolSolarliteDetailled. Alpha) annotation (Line(
      points={{15,310},{30,310},{30,300}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(TubeRePumpOut.outlet, join_Y.inlet2) annotation (Line(
      points={{214,-68},{240,-68},{240,-40}},
      color={0,131,169},
      pattern=LinePattern.Solid,
      thickness=0.5));
  connect(join_Y.inlet1, TubeIn.outlet) annotation (Line(
      points={{250,-30},{266,-30}},
      color={0,131,169},
      pattern=LinePattern.Solid,
      thickness=0.5));
  connect(TubeIn1.inlet, join_Y.outlet) annotation (Line(
      points={{214,-30},{230,-30}},
      color={0,131,169},
      pattern=LinePattern.Solid,
      thickness=0.5));
  connect(split_Y1.inlet, TubeIn1.outlet) annotation (Line(
      points={{176,-30},{186,-30}},
      color={0,131,169},
      pattern=LinePattern.Solid,
      thickness=0.5));
  connect(valveBranch.inlet, split_Y1.outlet1) annotation (Line(
      points={{134,-30},{156,-30}},
      color={0,131,169},
      pattern=LinePattern.Solid,
      thickness=0.5));
  connect(PIDValveBranch2.u_s, realExpression3.y) annotation (Line(points={{216,70},{233,70}}, color={0,0,127}));
  connect(PIDValveBranch2.u_m, realExpression5.y) annotation (Line(points={{203.9,58},{204,58},{204,48},{233,48}}, color={0,0,127}));
  connect(PIDValveBranch3.u_s, realExpression7.y) annotation (Line(points={{216,130},{233,130}}, color={0,0,127}));
  connect(PIDValveBranch3.u_m, realExpression8.y) annotation (Line(points={{203.9,118},{203.9,108},{233,108}}, color={0,0,127}));
  connect(PIDValveBranch4.u_s, realExpression9.y) annotation (Line(points={{216,190},{233,190}}, color={0,0,127}));
  connect(PIDValveBranch4.u_m, realExpression10.y) annotation (Line(points={{203.9,178},{203.9,168},{233,168}}, color={0,0,127}));
  connect(PIDValveBranch5.u_s, realExpression11.y) annotation (Line(points={{216,250},{233,250}}, color={0,0,127}));
  connect(PIDValveBranch5.u_m, realExpression12.y) annotation (Line(points={{203.9,238},{206,238},{206,228},{233,228}}, color={0,0,127}));
  connect(PIDValveBranch6.u_s, realExpression13.y) annotation (Line(points={{216,310},{233,310}}, color={0,0,127}));
  connect(PIDValveBranch6.u_m, realExpression14.y) annotation (Line(points={{203.9,298},{203.9,288},{233,288}}, color={0,0,127}));
  connect(TubeEvapCollIn1.inlet, split_Y1.outlet2) annotation (Line(
      points={{166,-12},{166,-20}},
      color={0,131,169},
      pattern=LinePattern.Solid,
      thickness=0.5));
  connect(valveBranch5.inlet, TubeEvapCollIn5.outlet) annotation (Line(
      points={{134,270},{164,270},{164,254}},
      color={0,131,169},
      pattern=LinePattern.Solid,
      thickness=0.5));
  connect(TubeEvapCollIn1.outlet, split_Y2.inlet) annotation (Line(
      points={{166,16},{180,16},{180,30},{176,30}},
      color={0,131,169},
      pattern=LinePattern.Solid,
      thickness=0.5));
  connect(split_Y2.outlet2, TubeEvapCollIn2.inlet) annotation (Line(
      points={{166,40},{166,48}},
      color={0,131,169},
      pattern=LinePattern.Solid,
      thickness=0.5));
  connect(split_Y2.outlet1, valveBranch1.inlet) annotation (Line(
      points={{156,30},{134,30}},
      color={0,131,169},
      pattern=LinePattern.Solid,
      thickness=0.5));
  connect(PIDValveBranch2.y, valveBranch1.opening_in) annotation (Line(points={{193,70},{124,70},{124,39}}, color={0,0,127}));
  connect(valveBranch2.inlet,split_Y3. outlet1) annotation (Line(
      points={{134,90},{156,90}},
      color={0,131,169},
      pattern=LinePattern.Solid,
      thickness=0.5));
  connect(TubeEvapCollIn2.outlet,split_Y3. inlet) annotation (Line(
      points={{166,76},{180,76},{180,90},{176,90}},
      color={0,131,169},
      pattern=LinePattern.Solid,
      thickness=0.5));
  connect(split_Y3.outlet2, TubeEvapCollIn3.inlet) annotation (Line(
      points={{166,100},{166,108}},
      color={0,131,169},
      pattern=LinePattern.Solid,
      thickness=0.5));
  connect(TubeEvapCollIn3.outlet, split_Y4.inlet) annotation (Line(
      points={{166,136},{178,136},{178,150},{174,150}},
      color={0,131,169},
      pattern=LinePattern.Solid,
      thickness=0.5));
  connect(valveBranch3.inlet, split_Y4.outlet1) annotation (Line(
      points={{134,150},{154,150}},
      color={0,131,169},
      pattern=LinePattern.Solid,
      thickness=0.5));
  connect(TubeEvapCollIn4.inlet, split_Y4.outlet2) annotation (Line(
      points={{164,168},{164,160}},
      color={0,131,169},
      pattern=LinePattern.Solid,
      thickness=0.5));
  connect(TubeEvapCollIn4.outlet,split_Y5. inlet) annotation (Line(
      points={{164,196},{178,196},{178,210},{174,210}},
      color={0,131,169},
      pattern=LinePattern.Solid,
      thickness=0.5));
  connect(split_Y5.outlet2, TubeEvapCollIn5.inlet) annotation (Line(
      points={{164,220},{164,226}},
      color={0,131,169},
      pattern=LinePattern.Solid,
      thickness=0.5));
  connect(valveBranch4.inlet,split_Y5. outlet1) annotation (Line(
      points={{134,210},{154,210}},
      color={0,131,169},
      pattern=LinePattern.Solid,
      thickness=0.5));
  connect(PIDValveBranch5.y, valveBranch4.opening_in) annotation (Line(points={{193,250},{124,250},{124,219}}, color={0,0,127}));
  connect(PIDValveBranch6.y, valveBranch5.opening_in) annotation (Line(points={{193,310},{124,310},{124,279}}, color={0,0,127}));
  connect(PIDValveBranch4.y, valveBranch3.opening_in) annotation (Line(points={{193,190},{124,190},{124,159}}, color={0,0,127}));
  connect(PIDValveBranch3.y, valveBranch2.opening_in) annotation (Line(points={{193,130},{124,130},{124,99}}, color={0,0,127}));
  connect(tube.inlet, valveBranch.outlet) annotation (Line(
      points={{78,-30},{114,-30}},
      color={0,131,169},
      pattern=LinePattern.Solid,
      thickness=0.5));
  connect(join_Y1.inlet1, tube.outlet) annotation (Line(
      points={{-26,-30},{50,-30}},
      color={0,131,169},
      pattern=LinePattern.Solid,
      thickness=0.5));
  connect(wall.innerPhase, tube.heat) annotation (Line(
      points={{64,-22},{64,-26}},
      color={167,25,48},
      thickness=0.5));
  connect(wall.outerPhase, parabolSolarliteDetailled1.port) annotation (Line(
      points={{64,-14},{64,-10},{34,-10}},
      color={167,25,48},
      thickness=0.5));
  connect(wall1.innerPhase, tube1.heat) annotation (Line(
      points={{64,38},{64,34}},
      color={167,25,48},
      thickness=0.5));
  connect(tube1.inlet, valveBranch1.outlet) annotation (Line(
      points={{78,30},{114,30}},
      color={0,131,169},
      pattern=LinePattern.Solid,
      thickness=0.5));
  connect(wall1.outerPhase, parabolSolarliteDetailled5.port) annotation (Line(
      points={{64,46},{64,50},{34,50}},
      color={167,25,48},
      thickness=0.5));
  connect(wall2.innerPhase, tube2.heat) annotation (Line(
      points={{64,98},{64,94}},
      color={167,25,48},
      thickness=0.5));
  connect(valveBranch2.outlet, tube2.inlet) annotation (Line(
      points={{114,90},{78,90}},
      color={0,131,169},
      pattern=LinePattern.Solid,
      thickness=0.5));
  connect(wall2.outerPhase, parabolSolarliteDetailled4.port) annotation (Line(
      points={{64,106},{64,110},{34,110}},
      color={167,25,48},
      thickness=0.5));
  connect(wall3.innerPhase, tube3.heat) annotation (Line(
      points={{64,158},{64,154}},
      color={167,25,48},
      thickness=0.5));
  connect(valveBranch3.outlet, tube3.inlet) annotation (Line(
      points={{114,150},{78,150}},
      color={0,131,169},
      pattern=LinePattern.Solid,
      thickness=0.5));
  connect(wall3.outerPhase, parabolSolarliteDetailled3.port) annotation (Line(
      points={{64,166},{64,170},{34,170}},
      color={167,25,48},
      thickness=0.5));
  connect(wall4.innerPhase, tube4.heat) annotation (Line(
      points={{64,218},{64,214}},
      color={167,25,48},
      thickness=0.5));
  connect(tube4.inlet, valveBranch4.outlet) annotation (Line(
      points={{78,210},{114,210}},
      color={0,131,169},
      pattern=LinePattern.Solid,
      thickness=0.5));
  connect(wall4.outerPhase, parabolSolarliteDetailled2.port) annotation (Line(
      points={{64,226},{64,230},{34,230}},
      color={167,25,48},
      thickness=0.5));
  connect(wall5.innerPhase, tube5.heat) annotation (Line(
      points={{64,278},{64,274}},
      color={167,25,48},
      thickness=0.5));
  connect(tube5.inlet, valveBranch5.outlet) annotation (Line(
      points={{78,270},{114,270}},
      color={0,131,169},
      pattern=LinePattern.Solid,
      thickness=0.5));
  connect(wall5.outerPhase,parabolSolarliteDetailled. port) annotation (Line(
      points={{64,286},{64,290},{34,290}},
      color={167,25,48},
      thickness=0.5));
  connect(join_Y2.outlet, TubeEvapCollOut1.inlet) annotation (Line(
      points={{-36,20},{-36,14}},
      color={0,131,169},
      pattern=LinePattern.Solid,
      thickness=0.5));
  connect(TubeEvapCollOut1.outlet, join_Y1.inlet2) annotation (Line(
      points={{-36,-14},{-36,-20}},
      color={0,131,169},
      pattern=LinePattern.Solid,
      thickness=0.5));
  connect(tube5.outlet, TubeEvapCollOut5.inlet) annotation (Line(
      points={{50,270},{-36,270},{-36,254}},
      color={0,131,169},
      pattern=LinePattern.Solid,
      thickness=0.5));
  connect(TubeEvapCollOut5.outlet, join_Y5.inlet1) annotation (Line(
      points={{-36,226},{-36,220}},
      color={0,131,169},
      pattern=LinePattern.Solid,
      thickness=0.5));
  connect(join_Y5.inlet2, tube4.outlet) annotation (Line(
      points={{-26,210},{50,210}},
      color={0,131,169},
      pattern=LinePattern.Solid,
      thickness=0.5));
  connect(join_Y5.outlet, TubeEvapCollOut4.inlet) annotation (Line(
      points={{-36,200},{-36,194}},
      color={0,131,169},
      pattern=LinePattern.Solid,
      thickness=0.5));
  connect(TubeEvapCollOut3.outlet,join_Y3. inlet1) annotation (Line(
      points={{-36,106},{-36,100}},
      color={0,131,169},
      pattern=LinePattern.Solid,
      thickness=0.5));
  connect(join_Y4.outlet, TubeEvapCollOut3.inlet) annotation (Line(
      points={{-36,140},{-36,134}},
      color={0,131,169},
      pattern=LinePattern.Solid,
      thickness=0.5));
  connect(join_Y4.inlet2, tube3.outlet) annotation (Line(
      points={{-26,150},{50,150}},
      color={0,131,169},
      pattern=LinePattern.Solid,
      thickness=0.5));
  connect(TubeEvapCollOut4.outlet,join_Y4. inlet1) annotation (Line(
      points={{-36,166},{-36,160}},
      color={0,131,169},
      pattern=LinePattern.Solid,
      thickness=0.5));
  connect(join_Y3.outlet, TubeEvapCollOut2.inlet) annotation (Line(
      points={{-36,80},{-36,74}},
      color={0,131,169},
      pattern=LinePattern.Solid,
      thickness=0.5));
  connect(join_Y3.inlet2, tube2.outlet) annotation (Line(
      points={{-26,90},{50,90}},
      color={0,131,169},
      pattern=LinePattern.Solid,
      thickness=0.5));
  connect(join_Y2.inlet2, tube1.outlet) annotation (Line(
      points={{-26,30},{50,30}},
      color={0,131,169},
      pattern=LinePattern.Solid,
      thickness=0.5));
  connect(TubeEvapCollOut2.outlet, join_Y2.inlet1) annotation (Line(
      points={{-36,46},{-36,40}},
      color={0,131,169},
      pattern=LinePattern.Solid,
      thickness=0.5));
  connect(seperator.inlet, TubeInSep.outlet) annotation (Line(
      points={{-90,-30},{-78,-30}},
      color={0,131,169},
      pattern=LinePattern.Solid,
      thickness=0.5));
  connect(TubeInSep.inlet, join_Y1.outlet) annotation (Line(
      points={{-50,-30},{-46,-30}},
      color={0,131,169},
      pattern=LinePattern.Solid,
      thickness=0.5));
  connect(seperator.outlet2, TubeSHIN.inlet) annotation (Line(
      points={{-100,-20},{-112,-20},{-112,-30},{-120,-30}},
      color={0,131,169},
      pattern=LinePattern.Solid,
      thickness=0.5));
  connect(seperator.level, dynamicBar.u_in) annotation (Line(points={{-111,-38},{-93,-38}},color={0,0,127}));
  connect(seperator.outlet1, TubeRePumpIn.inlet) annotation (Line(
      points={{-100,-40},{-100,-68},{-34,-68}},
      color={0,131,169},
      pattern=LinePattern.Solid,
      thickness=0.5));
  connect(split_Y.inlet, TubeSHIN.outlet) annotation (Line(
      points={{-160,-30},{-148,-30}},
      color={0,131,169},
      pattern=LinePattern.Solid,
      thickness=0.5));
  connect(split_Y.outlet2, TubeSHCollIn1.inlet) annotation (Line(
      points={{-170,-20},{-170,-8}},
      color={0,131,169},
      pattern=LinePattern.Solid,
      thickness=0.5));
  connect(TubeSHCollIn1.outlet,split_Y6. inlet) annotation (Line(
      points={{-170,20},{-170,40},{-180,40}},
      color={0,131,169},
      pattern=LinePattern.Solid,
      thickness=0.5));
  connect(split_Y6.outlet2, TubeSHCollIn2.inlet) annotation (Line(
      points={{-190,50},{-190,60},{-170,60},{-170,72}},
      color={0,131,169},
      pattern=LinePattern.Solid,
      thickness=0.5));
  connect(join_Y7.outlet, TubeSHCollOut1.inlet) annotation (Line(
      points={{-320,30},{-320,20}},
      color={0,131,169},
      pattern=LinePattern.Solid,
      thickness=0.5));
  connect(join_Y6.inlet2, TubeSHCollOut1.outlet) annotation (Line(
      points={{-320,-20},{-320,-8}},
      color={0,131,169},
      pattern=LinePattern.Solid,
      thickness=0.5));
  connect(join_Y7.inlet1, TubeSHCollOut2.outlet) annotation (Line(
      points={{-320,50},{-320,72}},
      color={0,131,169},
      pattern=LinePattern.Solid,
      thickness=0.5));
  connect(tubeSH.inlet, split_Y.outlet1) annotation (Line(
      points={{-216,-30},{-180,-30}},
      color={0,131,169},
      pattern=LinePattern.Solid,
      thickness=0.5));
  connect(tubeSH.outlet, join_Y6.inlet1) annotation (Line(
      points={{-244,-30},{-310,-30}},
      color={0,131,169},
      pattern=LinePattern.Solid,
      thickness=0.5));
  connect(tubeSH.heat, wallSH.innerPhase) annotation (Line(
      points={{-230,-26},{-230,-22}},
      color={167,25,48},
      thickness=0.5));
  connect(wallSH.outerPhase, parabolSolarliteDetailled7.port) annotation (Line(
      points={{-230,-14},{-230,-10},{-248,-10}},
      color={167,25,48},
      thickness=0.5));
  connect(tubeSH1.heat, wallSH1.innerPhase) annotation (Line(
      points={{-230,44},{-230,48}},
      color={167,25,48},
      thickness=0.5));
  connect(tubeSH1.inlet,split_Y6. outlet1) annotation (Line(
      points={{-216,40},{-200,40}},
      color={0,131,169},
      pattern=LinePattern.Solid,
      thickness=0.5));
  connect(join_Y7.inlet2, tubeSH1.outlet) annotation (Line(
      points={{-310,40},{-244,40}},
      color={0,131,169},
      pattern=LinePattern.Solid,
      thickness=0.5));
  connect(const5.y, parabolSolarliteDetailled8.Alpha) annotation (Line(points={{-269,90},{-252,90},{-252,70}}, color={0,0,127}));
  connect(wallSH1.outerPhase, parabolSolarliteDetailled8.port) annotation (Line(
      points={{-230,56},{-230,60},{-248,60}},
      color={167,25,48},
      thickness=0.5));
  connect(tubeSH2.heat, wallSH2.innerPhase) annotation (Line(
      points={{-230,114},{-230,118}},
      color={167,25,48},
      thickness=0.5));
  connect(wallSH2.outerPhase, parabolSolarliteDetailled6.port) annotation (Line(
      points={{-230,126},{-230,130},{-248,130}},
      color={167,25,48},
      thickness=0.5));
  connect(tubeSH2.inlet, TubeSHCollIn2.outlet) annotation (Line(
      points={{-216,110},{-170,110},{-170,100}},
      color={0,131,169},
      pattern=LinePattern.Solid,
      thickness=0.5));
  connect(TubeSHCollOut2.inlet, tubeSH2.outlet) annotation (Line(
      points={{-320,100},{-320,110},{-244,110}},
      color={0,131,169},
      pattern=LinePattern.Solid,
      thickness=0.5));
  connect(tubeSH2.eye, quadruple1.eye) annotation (Line(points={{-244.6,106.6},{-244.6,98.35},{-242,98.35},{-242,98}}, color={190,190,190}));
  connect(quadruple2.eye, tubeSH1.eye) annotation (Line(points={{-242,28},{-244.6,28},{-244.6,36.6}}, color={190,190,190}));
  connect(quadruple3.eye, tubeSH.eye) annotation (Line(points={{-244,-44},{-244,-40},{-244.6,-40},{-244.6,-33.4}}, color={190,190,190}));
  connect(TubeSHCollIn2.eye, quadruple6.eye) annotation (Line(points={{-173.4,100.6},{-178.7,100.6},{-178.7,98},{-204,98}}, color={190,190,190}));
  connect(quadruple5.eye, split_Y6.eye[1]) annotation (Line(points={{-202,26},{-202,46},{-200,46}}, color={190,190,190}));
  connect(quadruple4.eye, split_Y.eye[1]) annotation (Line(points={{-216,-44},{-216,-56},{-188,-56},{-188,-24},{-180,-24}}, color={190,190,190}));
  connect(quadruple8.eye, seperator.eye_out1) annotation (Line(points={{-90,-52},{-96,-52},{-96,-41},{-104,-41}}, color={190,190,190}));
  connect(quadruple7.eye, seperator.eye_out2) annotation (Line(points={{-108,-8},{-108,-19},{-104,-19}}, color={190,190,190}));
  connect(quadruple9.eye, join_Y7.eye) annotation (Line(points={{-374,42},{-384,42},{-384,30},{-328,30}}, color={190,190,190}));
  connect(quadruple10.eye, TubeInSep.eye) annotation (Line(points={{-74,-8},{-78.6,-8},{-78.6,-33.4}}, color={190,190,190}));
  connect(quadruple14.eye, join_Y.eye) annotation (Line(points={{198,-44},{226,-44},{226,-22},{230,-22}}, color={190,190,190}));
  connect(quadruple15.eye, valveBranch.eye) annotation (Line(points={{102,-52},{100,-52},{100,-34},{114,-34}}, color={190,190,190}));
  connect(quadruple16.eye, valveBranch1.eye) annotation (Line(points={{94,8},{92,8},{92,26},{114,26}}, color={190,190,190}));
  connect(quadruple17.eye, valveBranch2.eye) annotation (Line(points={{94,68},{92,68},{92,86},{114,86}}, color={190,190,190}));
  connect(quadruple18.eye, valveBranch3.eye) annotation (Line(points={{96,126},{94,126},{94,146},{114,146}}, color={190,190,190}));
  connect(quadruple19.eye, valveBranch4.eye) annotation (Line(points={{96,188},{92,188},{92,206},{114,206}}, color={190,190,190}));
  connect(quadruple20.eye, valveBranch5.eye) annotation (Line(points={{96,248},{92,248},{92,266},{114,266}}, color={190,190,190}));
  connect(quadruple21.eye, TubeEvapCollIn5.eye) annotation (Line(points={{140,280},{140,296},{172,296},{172,268},{160.6,268},{160.6,254.6}}, color={190,190,190}));
  connect(quadruple22.eye, split_Y4.eye[1]) annotation (Line(points={{130,168},{130,156},{154,156}}, color={190,190,190}));
  connect(quadruple23.eye, split_Y5.eye[1]) annotation (Line(points={{132,228},{132,216},{154,216}}, color={190,190,190}));
  connect(quadruple24.eye, split_Y3.eye[1]) annotation (Line(points={{132,110},{132,96},{156,96}}, color={190,190,190}));
  connect(quadruple25.eye, split_Y2.eye[1]) annotation (Line(points={{132,48},{132,36},{156,36}}, color={190,190,190}));
  connect(quadruple26.eye, split_Y1.eye[1]) annotation (Line(points={{134,-6},{132,-6},{132,-24},{156,-24}}, color={190,190,190}));
  connect(quadruple27.eye, tube.eye) annotation (Line(points={{50,-44},{49.4,-44},{49.4,-33.4}}, color={190,190,190}));
  connect(quadruple28.eye, tube1.eye) annotation (Line(points={{52,16},{49.4,16},{49.4,26.6}}, color={190,190,190}));
  connect(quadruple29.eye, tube2.eye) annotation (Line(points={{52,74},{49.4,74},{49.4,86.6}}, color={190,190,190}));
  connect(quadruple30.eye, tube3.eye) annotation (Line(points={{52,136},{49.4,136},{49.4,146.6}}, color={190,190,190}));
  connect(quadruple31.eye, tube4.eye) annotation (Line(points={{50,196},{49.4,196},{49.4,206.6}}, color={190,190,190}));
  connect(quadruple32.eye, tube5.eye) annotation (Line(points={{52,258},{49.4,258},{49.4,266.6}}, color={190,190,190}));
  connect(quadruple33.eye, join_Y5.eye) annotation (Line(points={{-76,212},{-84,212},{-84,200},{-44,200}}, color={190,190,190}));
  connect(quadruple34.eye, join_Y4.eye) annotation (Line(points={{-78,154},{-84,154},{-84,140},{-44,140}}, color={190,190,190}));
  connect(quadruple35.eye, join_Y3.eye) annotation (Line(points={{-76,96},{-80,96},{-80,80},{-44,80}}, color={190,190,190}));
  connect(quadruple36.eye, join_Y2.eye) annotation (Line(points={{-78,32},{-82,32},{-82,20},{-44,20}}, color={190,190,190}));
  connect(quadruple37.eye, join_Y1.eye) annotation (Line(points={{-42,-48},{-42,-47},{-46,-47},{-46,-38}}, color={190,190,190}));
  connect(PIDRecircPump.y,pumpRe. P_drive) annotation (Line(points={{-115,-56},{90,-56}},  color={0,0,127}));
  connect(realExpression15.y,PIDRecircPump. u_m) annotation (Line(points={{-161,-76},{-126,-76},{-126,-68},{-125.9,-68}}, color={0,0,127}));
  connect(realExpression16.y,PIDRecircPump. u_s) annotation (Line(points={{-161,-56},{-138,-56}}, color={0,0,127}));
  connect(TubeRePumpIn.outlet, pumpRe.inlet) annotation (Line(
      points={{-6,-68},{80,-68}},
      color={0,131,169},
      pattern=LinePattern.Solid,
      thickness=0.5));
  connect(pumpRe.outlet, TubeRePumpOut.inlet) annotation (Line(
      points={{100,-68},{186,-68}},
      color={0,131,169},
      pattern=LinePattern.Solid,
      thickness=0.5));
  connect(quadruple12.eye, pumpRe.eye) annotation (Line(points={{122,-80},{101,-80},{101,-74}}, color={190,190,190}));
  connect(TubeRePumpIn.eye, quadruple11.eye) annotation (Line(points={{-5.4,-71.4},{-5.4,-80},{4,-80}}, color={190,190,190}));
  connect(quadruple38.eye, TubeRePumpOut.eye) annotation (Line(points={{224,-80},{214.6,-80},{214.6,-71.4}}, color={190,190,190}));
  connect(realExpression1.y,PIDMainPump. u_s) annotation (Line(
      points={{289,-100},{306,-100}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(realExpression27.y,PIDMainPump. u_m) annotation (Line(
      points={{289,-120},{318,-120},{318,-112},{318.1,-112}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(PIDMainPump.y,pumpMain. P_drive) annotation (Line(
      points={{329,-100},{344,-100},{344,-132}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(realExpression26.y,PIDMainValve. u_m)
                                       annotation (Line(
      points={{301,28},{330,28},{330,38},{330.1,38}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(PIDMainValve.y,valveMain. opening_in)
                                             annotation (Line(
      points={{341,50},{346,50},{346,-21}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(valveMain.outlet, TubeIn.inlet) annotation (Line(
      points={{336,-30},{294,-30}},
      color={0,131,169},
      pattern=LinePattern.Solid,
      thickness=0.5));
  connect(pumpMain.outlet, TubeMainValveIn.inlet) annotation (Line(
      points={{354,-144},{424,-144},{424,-30},{414,-30}},
      color={0,131,169},
      pattern=LinePattern.Solid,
      thickness=0.5));
  connect(valveMain.inlet, TubeMainValveIn.outlet) annotation (Line(
      points={{356,-30},{386,-30}},
      color={0,131,169},
      pattern=LinePattern.Solid,
      thickness=0.5));
  connect(realExpression25.y, PIDMainValve.u_s) annotation (Line(points={{301,50},{318,50}}, color={0,0,127}));
  connect(quadruple13.eye, TubeIn.eye) annotation (Line(points={{268,-54},{265.4,-54},{265.4,-33.4}}, color={190,190,190}));
  connect(quadruple39.eye, valveMain.eye) annotation (Line(points={{338,-54},{336,-54},{336,-34}}, color={190,190,190}));
  connect(quadruple40.eye, TubeMainValveIn.eye) annotation (Line(points={{388,-54},{385.4,-54},{385.4,-33.4}}, color={190,190,190}));
  connect(TubeMainPumpIn.outlet, pumpMain.inlet) annotation (Line(
      points={{194,-144},{334,-144}},
      color={0,131,169},
      pattern=LinePattern.Solid,
      thickness=0.5));
  connect(TubeMainPumpIn.eye, quadruple42.eye) annotation (Line(points={{194.6,-147.4},{194.6,-156},{304,-156},{304,-166},{316,-166}},
                                                                                                                    color={190,190,190}));
  connect(quadruple43.eye, pumpMain.eye) annotation (Line(points={{366,-166},{360,-166},{360,-150},{355,-150}},
                                                                                                        color={190,190,190}));
  connect(join_Y6.outlet, valve_turbine.inlet) annotation (Line(
      points={{-330,-30},{-344,-30}},
      color={0,131,169},
      pattern=LinePattern.Solid,
      thickness=0.5));
  connect(PI_Pump_cond.y,Pump_cond. P_drive) annotation (Line(
      points={{-159,-150},{-159,-152},{-124,-152},{-124,-134}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(measurement.y,PI_Pump_cond. u_m) annotation (Line(
      points={{-205.6,-170},{-170,-170},{-170,-162},{-169.9,-162}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(setPoint_condenser.y,PI_Pump_cond. u_s) annotation (Line(points={{-201.3,-150},{-182,-150}},
                                                                                                     color={0,0,127}));
  connect(condenser.level,measurement. u) annotation (Line(points={{-304,-147},{-304,-170},{-214.8,-170}},
                                                                                                     color={0,0,127}));
  connect(boundaryVLE_Txim_flow.steam_a,condenser. In2) annotation (Line(
      points={{-340,-150},{-319,-150},{-319,-140},{-306,-140}},
      color={0,131,169},
      pattern=LinePattern.Solid,
      thickness=0.5));
  connect(condenser.Out2,valveControl_preheater_LP1. inlet) annotation (Line(
      points={{-306,-130},{-312,-130}},
      color={0,131,169},
      pattern=LinePattern.Solid,
      thickness=0.5));
  connect(boundaryVLE_phxi.steam_a, valveControl_preheater_LP1.outlet) annotation (Line(
      points={{-340,-130},{-332,-130}},
      color={0,131,169},
      pattern=LinePattern.Solid,
      thickness=0.5));
  connect(condenser.Out1, Pump_cond.inlet) annotation (Line(
      points={{-296,-146},{-296,-162},{-252,-162},{-252,-122},{-134,-122}},
      color={0,131,169},
      pattern=LinePattern.Solid,
      thickness=0.5));
  connect(turbine.outlet, condenser.In1) annotation (Line(
      points={{-390,-100},{-296,-100},{-296,-126.2}},
      color={0,131,169},
      pattern=LinePattern.Solid,
      thickness=0.5));
  connect(join_Y6.eye, quadruple41.eye) annotation (Line(points={{-330,-38},{-330,-42},{-328,-42},{-328,-46}}, color={190,190,190}));
  connect(valve_turbine.eye, quadruple44.eye) annotation (Line(points={{-364,-34},{-392,-34},{-392,-46}}, color={190,190,190}));
  connect(turbine.eye, quadruple45.eye) annotation (Line(points={{-389,-96},{-384,-96},{-384,-90},{-376,-90}}, color={190,190,190}));
  connect(quadruple46.eye, condenser.eye1) annotation (Line(points={{-288,-106},{-290,-106},{-290,-120},{-266,-120},{-266,-154},{-300,-154},{-300,-147}}, color={190,190,190}));
  connect(Pump_cond.eye, quadruple47.eye) annotation (Line(points={{-113,-116},{-113,-112},{-108,-112}}, color={190,190,190}));
  connect(PID_TurbineValve.y, valve_turbine.opening_in) annotation (Line(points={{-359,10},{-354,10},{-354,-21}}, color={0,0,127}));
  connect(realExpression4.y, PID_TurbineValve.u_s) annotation (Line(points={{-399,10},{-382,10}}, color={0,0,127}));
  connect(realExpression17.y, PID_TurbineValve.u_m) annotation (Line(points={{-399,-12},{-369.9,-12},{-369.9,-2}}, color={0,0,127}));
  connect(turbine.inlet, splitTurbine.outlet1) annotation (Line(
      points={{-400,-84},{-400,-70}},
      color={0,131,169},
      pattern=LinePattern.Solid,
      thickness=0.5));
  connect(splitTurbine.inlet, valve_turbine.outlet) annotation (Line(
      points={{-400,-50},{-400,-30},{-364,-30}},
      color={0,131,169},
      pattern=LinePattern.Solid,
      thickness=0.5));
  connect(splitTurbine.outlet2, ValveSteamFWT.inlet) annotation (Line(
      points={{-390,-60},{-244,-60},{-244,-96},{-230,-96}},
      color={0,131,169},
      pattern=LinePattern.Solid,
      thickness=0.5));
  connect(splitTurbine.eye[2], quadruple50.eye) annotation (Line(points={{-394,-70},{-396,-70},{-396,-76},{-388,-76},{-388,-64},{-352,-64},{-352,-76}}, color={190,190,190}));
  connect(splitTurbine.eye[1], quadruple51.eye) annotation (Line(points={{-394,-70},{-394,-76},{-384,-76}}, color={190,190,190}));
  connect(ValveSteamFWT.eye, quadruple49.eye) annotation (Line(points={{-210,-100},{-204,-100},{-204,-106},{-200,-106}}, color={190,190,190}));
  connect(Pump_cond.outlet, feedWaterTank.condensate) annotation (Line(
      points={{-114,-122},{-70,-122}},
      color={0,131,169},
      pattern=LinePattern.Solid,
      thickness=0.5));
  connect(ValveSteamFWT.outlet, feedWaterTank.heatingSteam) annotation (Line(
      points={{-210,-96},{-30,-96},{-30,-120}},
      color={0,131,169},
      pattern=LinePattern.Solid,
      thickness=0.5));
  connect(TubeMainPumpIn.inlet, feedWaterTank.feedwater) annotation (Line(
      points={{166,-144},{-24,-144},{-24,-138}},
      color={0,131,169},
      pattern=LinePattern.Solid,
      thickness=0.5));
  connect(quadruple48.eye, feedWaterTank.eye) annotation (Line(points={{-20,-106},{-20,-148},{-28,-148},{-28,-139}}, color={190,190,190}));
  connect(fillingLevel_condenser.u_in, condenser.level) annotation (Line(points={{-277,-144},{-277,-147},{-304,-147}}, color={0,0,127}));
  connect(parabolSolarliteDetailled.SunIn, sun.y) annotation (Line(points={{19,300},{19,308},{23,308}}, color={0,0,127}));
  connect(parabolSolarliteDetailled2.SunIn, sun1.y) annotation (Line(points={{19,240},{19,248},{23,248}}, color={0,0,127}));
  connect(sun2.y, parabolSolarliteDetailled3.SunIn) annotation (Line(points={{23,188},{19,188},{19,180}}, color={0,0,127}));
  connect(sun3.y, parabolSolarliteDetailled4.SunIn) annotation (Line(points={{23,128},{19,128},{19,120}}, color={0,0,127}));
  connect(sun4.y, parabolSolarliteDetailled5.SunIn) annotation (Line(points={{23,68},{19,68},{19,60}}, color={0,0,127}));
  connect(sun5.y, parabolSolarliteDetailled1.SunIn) annotation (Line(points={{23,8},{19,8},{19,0}}, color={0,0,127}));
  connect(sun6.y, parabolSolarliteDetailled6.SunIn) annotation (Line(points={{-247,146},{-264,146},{-264,140},{-263,140}}, color={0,0,127}));
  connect(parabolSolarliteDetailled8.SunIn, sun7.y) annotation (Line(points={{-263,70},{-264,70},{-264,76},{-249,76}}, color={0,0,127}));
  connect(parabolSolarliteDetailled7.SunIn, sun8.y) annotation (Line(points={{-263,0},{-263,6},{-249,6}}, color={0,0,127}));
  connect(realExpression18.y, parabolSolarliteDetailled.Q_radiation[1]) annotation (Line(points={{-5,290},{14,290}}, color={0,0,127}));
  connect(realExpression28.y, parabolSolarliteDetailled2.Q_radiation[1]) annotation (Line(points={{-5,230},{14,230}}, color={0,0,127}));
  connect(realExpression24.y, parabolSolarliteDetailled3.Q_radiation[1]) annotation (Line(points={{-5,170},{14,170}}, color={0,0,127}));
  connect(realExpression23.y, parabolSolarliteDetailled4.Q_radiation[1]) annotation (Line(points={{-5,110},{14,110}}, color={0,0,127}));
  connect(realExpression22.y, parabolSolarliteDetailled5.Q_radiation[1]) annotation (Line(points={{-5,50},{14,50}}, color={0,0,127}));
  connect(realExpression21.y, parabolSolarliteDetailled1.Q_radiation[1]) annotation (Line(points={{-5,-10},{14,-10}}, color={0,0,127}));
  connect(realExpression20.y, parabolSolarliteDetailled7.Q_radiation[1]) annotation (Line(points={{-279,-10},{-268,-10}}, color={0,0,127}));
  connect(realExpression19.y, parabolSolarliteDetailled8.Q_radiation[1]) annotation (Line(points={{-279,60},{-268,60}}, color={0,0,127}));
  connect(realExpression.y, parabolSolarliteDetailled6.Q_radiation[1]) annotation (Line(points={{-279,130},{-268,130}}, color={0,0,127}));
  annotation (Diagram(coordinateSystem(preserveAspectRatio=true, extent={{-440,-180},{440,400}}), graphics={
                                                                         Text(
          extent={{-420,412},{-222,326}},
          textColor={28,108,200},
          textStyle={TextStyle.Bold,TextStyle.UnderLine},
          textString="Thermal-Solar-Power-Plant"), Rectangle(
          extent={{270,390},{430,350}},
          lineColor={0,0,0},
          fillColor={215,215,215},
          fillPattern=FillPattern.Solid)}),
    experiment(
      StopTime=32400,
      Tolerance=1e-05,
      __Dymola_Algorithm="Dassl"),
    __Dymola_experimentSetupOutput(equdistant=false),
    Icon(coordinateSystem(preserveAspectRatio=true)));
end SixBranchesEvapThreeBranchsSH_Turbine_Detailed;
