within ClaRa_Solar.Examples.StaticCycles;
model Static_Cycle_Execute
  extends ClaRa.Basics.Icons.PackageIcons.ExecutableExample100;

  inner ClaRa.SimCenter simCenter(redeclare TILMedia.VLEFluidTypes.TILMedia_SplineWater fluid1)
                                  annotation (Placement(transformation(extent={{-20,60},{20,80}})));
  Static_Cycle_SixBranchEvapThreeBranchSH_Turbine static_Cycle_SixBranchEvapThreeBranchSH_Turbine annotation (Placement(transformation(extent={{-30,-30},{30,30}})));
  annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{100,100}})));
end Static_Cycle_Execute;
