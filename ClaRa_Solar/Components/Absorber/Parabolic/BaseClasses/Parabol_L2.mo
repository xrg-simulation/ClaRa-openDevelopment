within ClaRa_Solar.Components.Absorber.Parabolic.BaseClasses;
model Parabol_L2 "Calculating Heat Flow with respect to incidence angle"
  import Modelica.Constants.pi;
  import Modelica.Constants.eps;
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

parameter Real Zs "parabol azimuth angle, angle between normal to parabol from true south, westward is positive"
                                                                                                    annotation(Dialog(group="Geometry"));
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

    Real facAlpha "modelling the change of heat flow due to angular movement of the parabol";
    Real HeaFlowAbsorber[n_Qrad_ax];
    Real theta "incidence angle";
    Real beta_dummy "parabol tilt angle from horizontal without control angle";
    Real beta "parabol tilt angle from horizontal";
    Real hourAngle=SunIn[1] "hour angle";
    Real decl=SunIn[2] "declination";
    Real Lat=SunIn[3] "latitude";
    Real alphaSun=max(SunIn[4],eps) "altitude angle sun";
    Real z=SunIn[5] "azimuth angle sun";
    Real facCosinus "cosinus losses of parabol due to horizontal part of irradiance vector";
    Real IAM "incidence angle modifier, losses due to bad reflection of beams with incidence angle > 0";
public
  Modelica.Blocks.Interfaces.RealInput SunIn[5]
        annotation (Placement(transformation(
        origin={-50,100},
        extent={{20,-20},{-20,20}},
        rotation=90)));

public
  Modelica.Blocks.Interfaces.RealInput Alpha(start=0)
        annotation (Placement(transformation(
        origin={60,100},
        extent={{20,-20},{-20,20}},
        rotation=90)));
equation
  //calculation of tracking angle, source Solar Energy Engineering, Kalogirou 2009
tan(beta_dummy*2*pi/360)=sin((z-(Zs+90))*2*pi/360)/tan(alphaSun*2*pi/360);
  //calculation of real parabol tilt angle
  beta=-beta_dummy+Alpha;
  //calculation of incidence angle, source Solar Energy Engineering, Kalogirou 2009
  theta*2*pi/360=acos(cos(beta*2*pi/360)*cos((90-alphaSun)*2*pi/360)+sin(beta*2*pi/360)*sin((90-alphaSun)*2*pi/360)*cos((z-Zs)*2*pi/360));
 //cos(theta*2*pi/360)=sin(Lat*2*pi/360)*sin(decl*2*pi/360)*cos((90-beta)*2*pi/360)-cos(Lat*2*pi/360)*sin(decl*2*pi/360)*sin((90-beta)*2*pi/360)*cos(Zs*2*pi/360)+cos(Lat*2*pi/360)*cos(decl*2*pi/360)*cos((90-beta)*2*pi/360)*cos(hourAngle*2*pi/360)+sin(Lat*2*pi/360)*cos(decl*2*pi/360)*sin((90-beta)*2*pi/360)*cos(hourAngle*2*pi/360)*cos(Zs*2*pi/360)+sin(Zs*2*pi/360)*cos(decl*2*pi/360)*sin((90-beta)*2*pi/360)*sin(hourAngle*2*pi/360);
  //cosinus losses
  facCosinus=cos(theta*2*pi/360);
  // incidence angle modifier
  IAM= cos(theta*2*pi/360)*(1+(sin(theta*2*pi/360))^3);
  facAlpha=0.9*exp(-4*(Alpha)^2)+0.1; //simple function of angular movement of parabol, representing only an estimation of the dynamic behaviour
  for j in 1:n_Qrad_ax loop
      HeaFlowAbsorber[j]=facOpt*facCosinus*facAlpha*IAM*Q_radiation[j];
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
        Rectangle(
          extent={{-88,70},{72,-36}},
          lineColor={255,255,255},
          fillColor={255,255,255},
          fillPattern=FillPattern.Solid),
        Text(
          extent={{-150,-88},{150,-128}},
          textString="%name",
          lineColor={0,0,255}),
        Ellipse(
          extent={{-92,94},{92,-92}},
          lineColor={0,0,0},
          lineThickness=1),
        Rectangle(
          extent={{-76,100},{64,-48}},
          lineColor={255,255,255},
          fillColor={255,255,255},
          fillPattern=FillPattern.Solid),
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
          points={{-118,66},{-60,-8}},
          color={255,128,0},
          thickness=0.5,
          smooth=Smooth.None,
          arrow={Arrow.None,Arrow.Filled}),
        Line(
          points={{-74,92},{-20,20}},
          color={255,128,0},
          thickness=0.5,
          smooth=Smooth.None,
          arrow={Arrow.None,Arrow.Filled}),
        Line(
          points={{-30,118},{28,44}},
          color={255,128,0},
          thickness=0.5,
          smooth=Smooth.None,
          arrow={Arrow.None,Arrow.Filled})}),
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
end Parabol_L2;
