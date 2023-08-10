within ClaRa_Solar.Components.Absorber.Parabolic.BaseClasses;
model Parabol_L1 "Calculating Heat Flow of Absorber Tube from Radiation with Simple Absorption Model, difference angle alpha for control of heat input"
  parameter Modelica.Units.SI.Length L "Length of cylinder" annotation (Dialog(group="Geometry"));
  parameter Integer n_ax = 3 "Number of axial elements" annotation(Dialog(group="Discretisation"));
  parameter Modelica.Units.SI.Length dx[n_ax]=ClaRa.Basics.Functions.GenerateGrid(
      {1,-1},
      L,
      n_ax) "Discretisation scheme" annotation (Dialog(group="Discretisation"));
  parameter Modelica.Units.SI.Length W "Width of parabol";
  parameter Integer n_Qrad_ax = 1 "Number of Zones with different Radiation"
                                                                           annotation(Dialog(group="Radiation"));
  parameter Real facOpt "Optic Losses of Parabol" annotation(Dialog(group="Radiation"));

  Modelica.Blocks.Interfaces.RealInput Q_radiation[n_Qrad_ax]
        annotation (Placement(transformation(
        origin={-100,0},
        extent={{20,-20},{-20,20}},
        rotation=180)));
  Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_b port[n_ax]
                             annotation (Placement(transformation(extent={{90,
            -10},{110,10}}, rotation=0)));
protected
parameter Integer n_Qrad_disc=integer(floor(n_ax/n_Qrad_ax));

    Real facAlpha; // modelling the change of heat flow due to angular movement of the parabol
    Real HeaFlowAbsorber[n_Qrad_ax];

public
  Modelica.Blocks.Interfaces.RealInput Alpha(start= 0)
        annotation (Placement(transformation(
        origin={0,100},
        extent={{20,-20},{-20,20}},
        rotation=90), iconTransformation(
        extent={{20,-20},{-20,20}},
        rotation=90,
        origin={-54,100})));

equation
  facAlpha=0.9*exp(-4*(Alpha)^2)+0.1; //simple function of angular movement of parabol, representing only an estimation of the dynamic behaviour
  for j in 1:n_Qrad_ax loop
      HeaFlowAbsorber[j]=facOpt*facAlpha*Q_radiation[j];
for i in 1:n_Qrad_disc loop
  port[i+(j-1)*n_Qrad_disc].Q_flow = -HeaFlowAbsorber[j]*dx[i+(j-1)*n_Qrad_disc]*W;

  end for;
  end for;
  annotation (
    Icon(coordinateSystem(preserveAspectRatio=true, extent={{-100,-100},{100,
            100}}),     graphics={
        Rectangle(
          extent={{-100,100},{100,-100}},
          lineColor={255,255,255},
          fillColor={255,255,255},
          fillPattern=FillPattern.Solid),
        Text(
          extent={{-146,-88},{154,-128}},
          textString="%name",
          lineColor={0,0,255}),
        Ellipse(
          extent={{-92,94},{92,-92}},
          lineColor={0,0,0},
          lineThickness=1),
        Rectangle(
          extent={{-100,100},{40,-48}},
          lineColor={255,255,255},
          fillColor={255,255,255},
          fillPattern=FillPattern.Solid),
        Ellipse(
          extent={{-16,16},{14,-14}},
          lineColor={255,255,255},
          fillColor={0,0,0},
          fillPattern=FillPattern.Solid),
        Line(
          points={{-74,92},{-20,20}},
          color={255,128,0},
          smooth=Smooth.None,
          arrow={Arrow.None,Arrow.Filled},
          thickness=0.5),
        Line(
          points={{-118,66},{-60,-8}},
          color={255,128,0},
          smooth=Smooth.None,
          arrow={Arrow.None,Arrow.Filled},
          thickness=0.5),
        Line(
          points={{-30,118},{28,44}},
          color={255,128,0},
          smooth=Smooth.None,
          arrow={Arrow.None,Arrow.Filled},
          thickness=0.5)}),
    Documentation(info="<HTML>
<p>
This model allows a specified amount of heat flow rate to be \"injected\"
into a thermal system at a given port.  The amount of heat
is given by the input signal Q_flow into the model. The heat flows into the
component to which the component PrescribedHeatFlow is connected,
if the input signal is positive.
</p>
<p>
If parameter alpha is > 0, the heat flow is mulitplied by (1 + alpha*(port.T - T_ref))
in order to simulate temperature dependent losses (which are given an reference temperature T_ref).
</p>
</HTML>
"), Diagram(coordinateSystem(preserveAspectRatio=true, extent={{-100,-100},{100,
            100}}),      graphics={
        Line(
          points={{-60,-20},{68,-20}},
          color={191,0,0},
          thickness=0.5),
        Line(
          points={{-60,20},{68,20}},
          color={191,0,0},
          thickness=0.5),
        Line(
          points={{-80,0},{-60,-20}},
          color={191,0,0},
          thickness=0.5),
        Line(
          points={{-80,0},{-60,20}},
          color={191,0,0},
          thickness=0.5),
        Polygon(
          points={{60,0},{60,40},{90,20},{60,0}},
          lineColor={191,0,0},
          fillColor={191,0,0},
          fillPattern=FillPattern.Solid),
        Polygon(
          points={{60,-40},{60,0},{90,-20},{60,-40}},
          lineColor={191,0,0},
          fillColor={191,0,0},
          fillPattern=FillPattern.Solid)}));
end Parabol_L1;
