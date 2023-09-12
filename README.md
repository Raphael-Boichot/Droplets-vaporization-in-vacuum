## Numerical modeling of the droplet vaporization for design and operation of liquid pulsed CVD

**Supporting data for paper https://doi.org/10.1002/cvde.201507191**

The repo contains the Matlab codes and the theory that explains what happens when you send a droplet of a liquid in vaccum, in particular a droplet of binary liquid. This leading question, surprisingly, had no answer until I decided to tackle the set of coupled non linear differential equations equations solving the problem. The output of the codes is a vaporization time that helps to size the reactors used in Pressure Pulsed MOCVD process. The vaporization is made in two steps: a flash vaporization where all the sensible heat of the liquid is drained to sustain boiling. The mass transport in this regime is limited by Hertz-Knudsen flux (other said, the maximum escape velocity of the atoms from a finite surface at a given temperature). After sufficient cooling down, the droplet approches a radiative thermal equilibrium with its environment and the vaporization becomes limited by the radiative thermal flux that can reach it (we are in vacuum so no convective flux). In this state, the droplet behaves more or less as a comet in space and can survive for quite a long time, more or less frozen (not billions of years but seconds, which is yet very detrimental for the process and the reaction rate uniformity). The fact that the liquid can be a binary just complicates a bit the phenomena depending on the vapor pressure of the two components.

The code is of course lagrangian and particles are considered as isolated, even if the total pressure in reactor is also calculated based on the total quantity of droplets introduced and their shrinking radius.

## Some scheme to explain what happens in brief

![Flash_Vaporization_Steps](Documentation/Flash_Vaporization.png)

The set of non linear differential equations can be solved with your prefered numerical method but Runge-Kutta derived ones do the job with just some minor instabilities at the transition between the two regimes that do not change the overall results. The main parameter to explain the vaporization time is (of course) the ambient radiating temperature, knowing that the limited step is the pseudo-equilibrium vaporization. This was more or less empirically known by people using Direct Liquid Injection systems in MOCVD but to my knowledge they never document any calculus before this study.

## (Not) funfact

The initial paper was published with some unit errors (kCal and kJ, the classical one) that change a bit the paper conclusions so the paper comes with an Erratum explaining this. The quite specialized Journal of CVD in which the paper was published was discontinued after this issue and became the more broad and coveted Advanced Materials.


