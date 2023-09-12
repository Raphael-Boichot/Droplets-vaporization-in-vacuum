% Simulation of droplet evaporation
function y=droplet(t,x);

global pumping_rate
global initial_mass
global reactor_volume
global reactor_temperature
global alpha_evap
global rho_sol
global rho_prec
global M_sol
global M_prec
global Cp_sol
global Cp_prec
global delta_H_sol
global delta_H_prec
global R
global A_sol
global A_prec
global B_sol
global B_prec
global C_sol
global C_prec
global N_droplets
global eps_prec
global eps_sol
global sigma_prec
global sigma_sol
global reactor_diameter
global NA
global droplet_initial_mass

seuil_coupure_masse=1e-6;
droplet_prec_mass=x(1,1);
droplet_sol_mass=x(2,1);
prec_pressure=x(3,1);
sol_pressure=x(4,1);
T_droplet=x(5,1);
z_droplet=x(6,1);
dz_droplet=x(7,1);

%calculation of the molar fraction of solvent
moles_sol=droplet_sol_mass/M_sol;
moles_prec=droplet_prec_mass/M_prec;
x_mole_sol=moles_sol/(moles_sol+moles_prec);
x_mole_prec=1-x_mole_sol;

if (droplet_prec_mass+droplet_sol_mass)>(droplet_initial_mass*seuil_coupure_masse);
sol_volume=droplet_sol_mass/rho_sol;
prec_volume=droplet_prec_mass/rho_prec;    
V_droplet=sol_volume+prec_volume;
droplet_surface=4*pi*(3*V_droplet/(4*pi))^(2/3);

%calculation of mass of solvent evaporated per droplet;
Psat_sol_mmHg=10.^(A_sol-B_sol./(C_sol+(T_droplet-273.15)));
vapor_pressure_sol=133.322*Psat_sol_mmHg;
phi_sol_vaporated=x_mole_sol*(droplet_surface*M_sol*(2/(2-alpha_evap)*(1/(2*pi*M_sol*R)^0.5)*(alpha_evap*vapor_pressure_sol/T_droplet^0.5-alpha_evap*sol_pressure/reactor_temperature^0.5)));

%calculation of mass of solvent evaporated per droplet;
Psat_prec_mmHg=10.^(A_prec-B_prec./(C_prec+(T_droplet-273.15)));
vapor_pressure_prec=133.322*Psat_prec_mmHg;
phi_prec_vaporated=x_mole_prec*(droplet_surface*M_prec*(2/(2-alpha_evap)*(1/(2*pi*M_prec*R)^0.5)*(alpha_evap*vapor_pressure_prec/T_droplet^0.5-alpha_evap*prec_pressure/reactor_temperature^0.5)));

%mass balance per droplet
d_droplet_sol_mass=-phi_sol_vaporated;
d_droplet_prec_mass=-phi_prec_vaporated;

%Computation of the pressure change
d_sol_pressure=N_droplets*phi_sol_vaporated*R*reactor_temperature/(M_sol*reactor_volume)-pumping_rate*sol_pressure/reactor_volume;
d_prec_pressure=N_droplets*phi_prec_vaporated*R*reactor_temperature/(M_prec*reactor_volume)-pumping_rate*prec_pressure/reactor_volume;

%Computation of the droplet temperature change
Cp_moyen=(droplet_prec_mass*Cp_prec+droplet_sol_mass*Cp_sol)/(droplet_prec_mass+droplet_sol_mass);
rho_droplet=(droplet_prec_mass+droplet_sol_mass)/V_droplet;
thermal_prec_flux=phi_prec_vaporated*delta_H_prec;
thermal_sol_flux=phi_sol_vaporated*delta_H_sol;
thermal_flux=thermal_prec_flux+thermal_sol_flux;
dT=(-thermal_flux-droplet_surface*1*4*5.67e-8*(T_droplet^4-reactor_temperature^4))/(V_droplet*rho_droplet*Cp_moyen);

%computation of the droplet fall with time
droplet_mass=droplet_prec_mass+droplet_sol_mass;
droplet_diameter=2.*(3.*V_droplet./(4.*pi)).^(1/3);
t_star=reactor_temperature/eps_sol;
omega_viscosity=1.16145/(t_star^0.14874)+0.52487/(exp(0.7732*t_star))+2.16178/(exp(2.43787*t_star));
gas_viscosity=2.6693e-6*(1000*M_sol*reactor_temperature)^0.5/(sigma_sol^2*omega_viscosity);
reactor_gas_velocity=-pumping_rate/(pi*reactor_diameter^2/4);
reactor_pressure=prec_pressure+sol_pressure;
c_kn=NA.*reactor_pressure./(R.*reactor_temperature);
c=reactor_pressure./(R.*reactor_temperature);
rho_gas=c*M_sol;
lpm=1./(2^(0.5).*pi.*(sigma_sol*1e-10).^2.*c_kn);
Cc=1+(2*lpm/droplet_diameter)*(1.257+0.4*exp(-0.55*lpm/droplet_diameter));
reynolds_droplet=droplet_diameter*rho_gas*abs(reactor_gas_velocity-dz_droplet)/gas_viscosity;
d2z_droplet=-9.81+((3*pi*gas_viscosity*droplet_diameter/Cc)*(reactor_gas_velocity-dz_droplet)*(1+0.15*reynolds_droplet^0.687))/droplet_mass;
% relaxation_time=Cc*rho_sol*droplet_diameter^2/(18*gas_viscosity);
%diffusion_coefficient=1.38064852e-23*reactor_temperature*Cc/(3*pi*gas_viscosity*droplet_diameter)
y(1,1)=d_droplet_prec_mass;
y(2,1)=d_droplet_sol_mass;
y(3,1)=d_prec_pressure;
y(4,1)=d_sol_pressure;
y(5,1)=dT;
y(6,1)=dz_droplet;
y(7,1)=d2z_droplet;

end;

if (droplet_prec_mass+droplet_sol_mass)<=(droplet_initial_mass*seuil_coupure_masse);
%Computation of the pressure change
d_sol_pressure=-pumping_rate*sol_pressure/reactor_volume;
d_prec_pressure=-pumping_rate*prec_pressure/reactor_volume;

y(1,1)=0;
y(2,1)=0;
y(3,1)=d_prec_pressure;
y(4,1)=d_sol_pressure;
y(5,1)=0;
y(6,1)=0;
y(7,1)=0;
end;


