within Grate_Boiler.Components;
package VolumesValvesFittings "All kinds of tube-shaped components"
//___________________________________________________________________________//
// Package of the ClaRaPlus library, version: 1.2.1                        //
//                                                                           //
// Copyright © 2018-2019, XRG Simulation GmbH and TLK Thermo GmbH.           //
//___________________________________________________________________________//
  extends ClaRa.Basics.Icons.PackageIcons.Components80;

  package Fittings "Bends, T-connectors,..."
  //___________________________________________________________________________//
  // Package of the ClaRaPlus library, version: 1.2.1                        //
  //                                                                           //
  // Copyright © 2018-2019, XRG Simulation GmbH and TLK Thermo GmbH.           //
  //___________________________________________________________________________//
    extends ClaRa.Basics.Icons.PackageIcons.Components60;

    model SplitGas_L1_flex "A voluminous split for an arbitrary number of inputs NOT CAPABLE FOR PHASE-CHANGE"
    //___________________________________________________________________________//
    // Component of the ClaRa library, version: 1.4.1                            //
    //                                                                           //
    // Licensed by the DYNCAP/DYNSTART research team under Modelica License 2.   //
    // Copyright  2013-2019, DYNCAP/DYNSTART research team.                      //
    //___________________________________________________________________________//
    // DYNCAP and DYNSTART are research projects supported by the German Federal //
    // Ministry of Economic Affairs and Energy (FKZ 03ET2009/FKZ 03ET7060).      //
    // The research team consists of the following project partners:             //
    // Institute of Energy Systems (Hamburg University of Technology),           //
    // Institute of Thermo-Fluid Dynamics (Hamburg University of Technology),    //
    // TLK-Thermo GmbH (Braunschweig, Germany),                                  //
    // XRG Simulation GmbH (Hamburg, Germany).                                   //
    //___________________________________________________________________________//

    //  extends ClaRa.Basics.Interfaces.DataInterfaceVector(N_sets=N_ports_out);
      extends ClaRa.Basics.Icons.Adapter5_fw;

      extends ClaRa.Basics.Icons.ComplexityLevel(complexity="L1");

      outer ClaRa.SimCenter simCenter;

      parameter TILMedia.GasTypes.BaseGas medium = simCenter.flueGasModel "Medium to be used in tubes" annotation (choicesAllMatching, Dialog(group="Fundamental Definitions"));

      parameter Integer N_ports_out(min=1)=1 "Number of outlet  ports" annotation(Evaluate=true, Dialog(tab="General",group="Fundamental Definitions"));//connectorSizing=true,
      parameter Real K_split[N_ports_out-1] = fill(0, N_ports_out-1) "fixed split ratio" annotation(Dialog(tab="General",group="Fundamental Definitions"));

      ClaRa.Basics.Interfaces.GasPortIn inlet(Medium=medium) "Inlet port" annotation (Placement(transformation(extent={{-110,-10},{-90,10}})));
      ClaRa.Basics.Interfaces.GasPortOut outlet[N_ports_out](each Medium=medium) "Outlet port" annotation (Placement(transformation(extent={{90,-10},{110,10}})));
      TILMedia.Gas_pT gas(gasType = medium,p=inlet.p, T=noEvent(actualStream(inlet.T_outflow)),xi=noEvent(actualStream(inlet.xi_outflow))) annotation (Placement(transformation(extent={{-80,-10},{-60,10}})));

    equation
    //~~~~~~~~~~~~~~~~~~~~~~~~~
    // Boundary conditions ~~~~
      for i in 1:N_ports_out-1 loop
        outlet[i].T_outflow = inStream(inlet.T_outflow);
        outlet[i].m_flow = -K_split[i]*inlet.m_flow;
        outlet[i].xi_outflow = inStream(inlet.xi_outflow);

      end for;
        outlet[N_ports_out].T_outflow = inStream(inlet.T_outflow);
        outlet[N_ports_out].m_flow = -(1-sum( K_split))*inlet.m_flow;
         outlet[N_ports_out].xi_outflow = inStream(inlet.xi_outflow);

         inlet.p = outlet[1].p;
        inlet.T_outflow=1000; // dummy, backflow is not supported
         inlet.xi_outflow = medium.xi_default; // dummy, backflow is not supported;

      annotation (Icon(graphics), Diagram(graphics));
    end SplitGas_L1_flex;
  end Fittings;
end VolumesValvesFittings;
