within ClaRa_Solar.Components.SunModels;
model Sun
          //Quelle Solar Energy Engineering, Kalogirou(2009)
  import Modelica.Constants.pi;
parameter Real dayOfYear=69 "Day of Year" annotation(Dialog(group="Time and Date"));
parameter Real startTimeSimulation=12 "Local Time at Simulation Start in Hours" annotation(Dialog(group="Time and Date"));
parameter Real DS=0 "Daylight Saving (if daylight saving (summertime) 60 else 0)"
                                                                                 annotation(Dialog(group="Time and Date"));
parameter Real SL=-105 "Standard Longitude in Degree (location east of Greenwich is negative)" annotation(Dialog(group="Location"));
parameter Real L=14.7 "Local Latitude in Degree" annotation(Dialog(group="Location"));
parameter Real LL=-99.0 "Local Longitude in Degree (location east of Greenwich is negative)"
                                                                                            annotation(Dialog(group="Location"));
parameter Boolean Green=true "True if Location East of Greenwich" annotation(Dialog(group="Location"));

Real decl "angle of declination";
Real alpha "solar altitude angle";
Real z "solar azimuth angel";
Real hourAngle "Hour Angle";
Real AST "apparent solar time";
Real LST "local standard time";
Real B;
Real ET "equation of time";
Real z_dummy;
Real hourAngle_dummy;
Real c1;
Real c2;
Real c3;

  Modelica.Blocks.Interfaces.RealOutput y[5]
    annotation (Placement(transformation(extent={{100,-70},{120,-50}})));
equation
decl=23.45*sin(360/365*(284+dayOfYear)/360*2*pi);
B=(dayOfYear-81)*360/364;
ET=9.87*sin(2*B/360*2*pi)-7.53*cos(B/360*2*pi)-1.5*sin(B/360*2*pi);
LST=startTimeSimulation*60+time/60;
if Green==true then
 AST=LST+ET-4*(SL-LL)-DS;
else
  AST=LST+ET+4*(SL-LL)-DS;
end if
  annotation (experiment(StopTime=7200), __Dymola_experimentSetupOutput);
hourAngle=(AST/60-12)*15;
sin(alpha/360*2*pi)=sin(decl/360*2*pi)*sin(L/360*2*pi)+cos(decl/360*2*pi)*cos(L/360*2*pi)*cos(hourAngle/360*2*pi);
//sin(z/360*2*pi)=cos(decl/360*2*pi)*sin(hourAngle/360*2*pi)/cos(alpha/360*2*pi); //only valid if the sun is not behind east_west line

//Calculation of Solar Azimuth Angle, Source Solar Engineering of Thermal Processes, Duffie, Beckmann 1991
z=c1*c2*z_dummy+c3*((1-c1*c2)/2)*180;

sin(z_dummy/360*2*pi)=cos(decl/360*2*pi)*sin(hourAngle/360*2*pi)/sin((90-alpha)/360*2*pi);
if tan(decl/360*2*pi)/tan(L/360*2*pi)>1 then
  cos(hourAngle_dummy/360*2*pi)=1;
elseif tan(decl/360*2*pi)/tan(L/360*2*pi)<(-1) then
  cos(hourAngle_dummy/360*2*pi)=-1;
  else
cos(hourAngle_dummy/360*2*pi)=tan(decl/360*2*pi)/tan(L/360*2*pi);
end if;

if abs(hourAngle)>hourAngle_dummy then
  c1=-1;
else
  c1=1;
end if;

if L-decl>=0 then
  c2=1;
else
  c2=-1;
end if;

if hourAngle>=0 then
  c3=1;
else
  c3=-1;
end if;

//  if cos(hourAngle/360*2*pi)>tan(decl/360*2*pi)/tan(L/360*2*pi) then
//  sin(z/360*2*pi)=cos(decl/360*2*pi)*sin(hourAngle/360*2*pi)/cos(alpha/360*2*pi); //only valid if the sun is not behind east_west line
//  else
//   pi-(z/360*2*pi)=asin(cos(decl/360*2*pi)*sin(hourAngle/360*2*pi)/cos(alpha/360*2*pi));
//  end if;
// dummy=cos(hourAngle/360*2*pi);
// dummy2=tan(decl/360*2*pi)/tan(L/360*2*pi);
y[1]=hourAngle;
y[2]=decl;
y[3]=L;
y[4]=alpha;
y[5]=z;
  annotation (experiment(StopTime=3000),
    Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{100,100}}),
        graphics={
        Rectangle(
          extent={{-100,100},{100,-100}},
          lineColor={0,0,0},
          fillPattern=FillPattern.Solid,
          lineThickness=1,
          fillColor={215,215,215}),
        Ellipse(
          extent={{-60,62},{60,-54}},
          lineColor={0,0,0},
          fillColor={255,255,0},
          fillPattern=FillPattern.Solid),
        Polygon(
          points={{-12,-58},{12,-58},{0,-94},{-12,-58}},
          lineColor={0,0,0},
          smooth=Smooth.None,
          fillColor={255,255,0},
          fillPattern=FillPattern.Solid),
        Polygon(
          points={{-12,15},{12,15},{-1,-20},{-12,15}},
          lineColor={0,0,0},
          smooth=Smooth.None,
          fillColor={255,255,0},
          fillPattern=FillPattern.Solid,
          origin={80,1},
          rotation=90),
        Polygon(
          points={{12,-15},{-12,-15},{-1,20},{12,-15}},
          lineColor={0,0,0},
          smooth=Smooth.None,
          fillColor={255,255,0},
          fillPattern=FillPattern.Solid,
          origin={-80,1},
          rotation=90),
        Polygon(
          points={{-12,66},{12,66},{0,100},{-12,66}},
          lineColor={0,0,0},
          smooth=Smooth.None,
          fillColor={255,255,0},
          fillPattern=FillPattern.Solid),
        Polygon(
          points={{-56,38},{-38,56},{-74,72},{-56,38}},
          lineColor={0,0,0},
          smooth=Smooth.None,
          fillColor={255,255,0},
          fillPattern=FillPattern.Solid),
        Polygon(
          points={{-3,-15},{15,3},{-19,19},{-3,-15}},
          lineColor={0,0,0},
          smooth=Smooth.None,
          fillColor={255,255,0},
          fillPattern=FillPattern.Solid,
          origin={-53,-47},
          rotation=90),
        Polygon(
          points={{-3,-15},{15,3},{-19,19},{-3,-15}},
          lineColor={0,0,0},
          smooth=Smooth.None,
          fillColor={255,255,0},
          fillPattern=FillPattern.Solid,
          origin={55,-45},
          rotation=180),
        Polygon(
          points={{-3,-15},{15,3},{-17,21},{-3,-15}},
          lineColor={0,0,0},
          smooth=Smooth.None,
          fillColor={255,255,0},
          fillPattern=FillPattern.Solid,
          origin={57,53},
          rotation=270)}),
    Diagram(coordinateSystem(preserveAspectRatio=true,  extent={{-100,-100},{
            100,100}}),
                    graphics));
end Sun;
