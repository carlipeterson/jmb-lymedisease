DRIVERInfection.m - run this file to find fixed points and 2-cycle points of the tick-host infection system and determine local stability. This file generates Figure 3 our paper, which depicts nymphal infection prevalence and the density of infected nymphs with respect to the questing parameter c.

R0.m - run this file to calculate the basic reproductive number R_0 on the fixed point and 2-cycle. This file generates Figures 1 and 2 in our paper, which depict R_0 with respect to (1) the demographic reproductive number and (2) the questing parameter c.

makeJacobian.m - this is called by the "R0.m" file to assess local stability; it numerically finds the Jacobian of the tick-host system without infection.

makeJacobianInfection.m - this is called by the "DRIVERInfection.m" file to assess local stability; it numerically finds the Jacobian of the tick-host system with infection.

mySystem.m - this is called by the "R0.m" files during the system solving process; it iterates the tickMap function once to find fixed points.

mySystem2.m - this is called by the "R0.m" file during the system solving process; it iterates the tickMap function twice to find 2-cycle points.

mySystemInfection.m - this is called by the "DRVERInfection.m" file during the system solving process; it iterates the tickMapInfection function once to find fixed points.

mySystem2Infection.m - this is called by the "DRIVERInfection.m" file during the system solving process; it iterates the tickMapInfection function twice to find 2-cycle points.

tickMap.m - this is called by the "mySystem.m" and "mySystem2.m" files during the system solving process; describes the tick-host model without infection.

tickMapInfection.m - this is called by the "mySystemInfection.m" and "mySystem2Infection.m" files during the system solving process; describes the tick-host model with infection.

hostsimulations.m - this script produces Figure 4 from our manuscript, which depicts the infected hosts populations as the infection scaling parameters beta and alpha vary.
