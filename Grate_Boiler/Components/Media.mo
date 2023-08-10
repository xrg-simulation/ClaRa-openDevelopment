within Grate_Boiler.Components;
package Media "Media for Fluegas and Fuel"
extends ClaRa.Basics.Icons.PackageIcons.Basics80;
  record Waste_NielsenEtAl "Waste composition  as described in GrateBoilerBeta - Beschreibung des Transports in der Brenstoffaufgabe"
    extends ClaRa.Basics.Media.FuelTypes.BaseFuel(
      N_c=7,
      N_e=7,
      C_LHV={33907e3,142324e3 - 2512e3*9,-142324e3/8,0,10465e3,0,-2512e3},
      C_rho={600,450,450,450,450,721,1000},
      C_cp={710,1300,1300,1300,1300,1300,4190},
      waterIndex=7,
      ashIndex=6,
      defaultComposition={0.357,0.0429,0.307,0.0,0.0,0.00714},
      xi_e_waf=zeros(0, 0),
      T_ref=273.15);

  end Waste_NielsenEtAl;

  record Waste_Warncke "Waste composition  as described in Warnecke - Beschreibung des Transports in der Brenstoffaufgabe"
    extends ClaRa.Basics.Media.FuelTypes.BaseFuel(
      N_c=7,
      N_e=7,
      C_LHV={33907e3,142324e3 - 2512e3*9,-142324e3/8,0,10465e3,0,-2512e3},
      C_rho={1400,1400,1400,1400,1400,1400,1000},
      C_cp={1266.67,1266.67,1266.67,1266.67,1266.67,1000,4190},
      waterIndex=7,
      ashIndex=6,
      defaultComposition={0.302,0.043,0.192,0.006,0.003,0.29},
      xi_e_waf=zeros(0, 0),
      T_ref=273.15);
      //composition from waste of GKS Schweinfurt
  end Waste_Warncke;

  record FlueGasTILMedia_GrateBoiler "Flue gas TILMedia (CH4,CO,CO2,SO2,N2,O2,NO,H2O,H2,Ar)"
    extends TILMedia.GasTypes.BaseGas(
      final fixedMixingRatio=false,
      final nc_propertyCalculation=10,
      final gasNames={"TILMediaXTR.methane","TILMediaXTR.Carbon_Monoxide",
          "TILMediaXTR.Carbon_Dioxide","TILMediaXTR.Sulfur_Dioxide",
          "TILMediaXTR.Nitrogen","TILMediaXTR.Oxygen",
          "TILMediaXTR.Nitrous_Oxide","TILMediaXTR.Water",
          "TILMediaXTR.Hydrogen","TILMediaXTR.Argon"},
      final condensingIndex=0,
      final mixingRatio_propertyCalculation={0.00,0.00,0.001,0.00,0.78,0.21,
          0.0,0.01,0.0,0.009});
  end FlueGasTILMedia_GrateBoiler;
end Media;
