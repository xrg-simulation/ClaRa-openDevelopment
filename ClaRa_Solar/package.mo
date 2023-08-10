within ;
package ClaRa_Solar "ClaRa Solar Library"

annotation (                       Icon(graphics={
        Rectangle(
          extent={{-100,100},{100,-100}},
          lineColor={0,0,0},
          fillPattern=FillPattern.Solid,
          lineThickness=1,
          fillColor={255,255,255}),
        Ellipse(
          extent={{-30,30},{30,-30}},
          lineColor={0,0,0},
          fillColor={175,175,175},
          fillPattern=FillPattern.Solid),
        Polygon(
          points={{-12,-42},{12,-42},{0,-78},{-12,-42}},
          lineColor={0,0,0},
          smooth=Smooth.None,
          fillColor={175,175,175},
          fillPattern=FillPattern.Solid),
        Polygon(
          points={{-12,15},{12,15},{-1,-20},{-12,15}},
          lineColor={0,0,0},
          smooth=Smooth.None,
          fillColor={175,175,175},
          fillPattern=FillPattern.Solid,
          origin={58,1},
          rotation=90),
        Polygon(
          points={{12,-15},{-12,-15},{-1,20},{12,-15}},
          lineColor={0,0,0},
          smooth=Smooth.None,
          fillColor={175,175,175},
          fillPattern=FillPattern.Solid,
          origin={-58,1},
          rotation=90),
        Polygon(
          points={{-12,42},{12,42},{0,76},{-12,42}},
          lineColor={0,0,0},
          smooth=Smooth.None,
          fillColor={175,175,175},
          fillPattern=FillPattern.Solid),
        Polygon(
          points={{-40,22},{-22,40},{-58,56},{-40,22}},
          lineColor={0,0,0},
          smooth=Smooth.None,
          fillColor={175,175,175},
          fillPattern=FillPattern.Solid),
        Polygon(
          points={{-3,-15},{15,3},{-19,19},{-3,-15}},
          lineColor={0,0,0},
          smooth=Smooth.None,
          fillColor={175,175,175},
          fillPattern=FillPattern.Solid,
          origin={-37,-37},
          rotation=90),
        Polygon(
          points={{-3,-15},{15,3},{-19,19},{-3,-15}},
          lineColor={0,0,0},
          smooth=Smooth.None,
          fillColor={175,175,175},
          fillPattern=FillPattern.Solid,
          origin={37,-37},
          rotation=180),
        Polygon(
          points={{-3,-15},{15,3},{-17,21},{-3,-15}},
          lineColor={0,0,0},
          smooth=Smooth.None,
          fillColor={175,175,175},
          fillPattern=FillPattern.Solid,
          origin={37,37},
          rotation=270)}), uses(
      ClaRa(version="1.8.1"),
      Modelica(version="4.0.0"),
      TILMedia(version="1.8.1 ClaRa")));
end ClaRa_Solar;
