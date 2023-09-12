%Initialisation du problème
clear;
clc;
format long
%data=load('data.txt');
tic
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
global droplet_initial_mass
global eps_prec
global eps_sol
global sigma_prec
global sigma_sol
global reactor_diameter
global NA

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%PROCESS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
initial_x_solvent=0.9;                  %mass ratio solvent/(solvent+precursor)
initial_droplet_diameter=18e-6;         %initial radius of droplets m
pas_temps_flash=1e-7;                   %integration step time for the flash vaporisation in s
pas_temps_pseudo_equi=1e-4;             %integration step time for the pseudo-equilibrium in s

%process parameter
alpha_evap=0.1;                         %evaporation coefficient SA for solvent and precursor
initial_prec_pressure=0;                %initial pressure of precursor in Pa
initial_sol_pressure=0;                 %initial pressure of solvent in Pa
initial_droplet_temperature=300;        %initial liquid temperature in K
pumping_rate=15/3600;                   %volume puming rate in m3/s
initial_mass=0.0001;                    %mixture mass introduced in kg
reactor_volume=0.01;                    %reactor volume in m3
reactor_temperature=300;                %reactor wall temperature in K
pulse_duration=10;                      %time between pulses    
reactor_diameter=0.1;                   %reactor diameter in m
susceptor_distance=0.25;                %distance nozzle-susceptor in m

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%PRECURSORS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%TTIP properties in USI
A_prec=8.2157;                          %Antoine equation PsatmmHg=10^(A-B/(C+T°C))
B_prec=2736.8;                          %Antoine equation PsatmmHg=10^(A-B/(C+T°C))
C_prec=273.15;                          %Antoine equation PsatmmHg=10^(A-B/(C+T°C))
rho_prec=960;                           %density of solvent kg.m-3
M_prec=0.284215;                        %molar mass of solvent kg.mol-1
Cp_prec = 2207;                         %took the one of dodecane J.kg-1.K-1
delta_H_prec=219.2e3;                    %heat for phase change boiling j.kg-1
melting_point_prec=17;                  %melting temperature in °C
sigma_prec=8.145;                       %Lennard Jones parameter in A  
eps_prec=556.8;                         %Lennard Jones parameter in K

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%PRECURSORS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Aluminum sec butoxide properties in USI
% A_prec=24.86;                          %Antoine equation PsatmmHg=10^(A-B/(C+T°C))
% B_prec=11083;                          %Antoine equation PsatmmHg=10^(A-B/(C+T°C))
% C_prec=273.15;                          %Antoine equation PsatmmHg=10^(A-B/(C+T°C))
% rho_prec=960;                           %density of solvent kg.m-3
% M_prec=0.24632;                        %molar mass of solvent kg.mol-1
% Cp_prec = 2411;                         %took 3* the one of 2 butanol J.kg-1.K-1
% delta_H_prec=330.87e3;                    %heat for phase change boiling j.kg-1
% melting_point_prec=-60;                  %melting temperature in °C
% sigma_prec=7.76;                       %Lennard Jones parameter in A  
% eps_prec=409;                         %Lennard Jones parameter in K

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%SOLVENTS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Hexane properties in USI
% A_sol=7.01051;                        %Antoine equation PsatmmHg=10^(A-B/(C+T°C))
% B_sol=1246.33;                        %Antoine equation PsatmmHg=10^(A-B/(C+T°C))
% C_sol=232.988;                        %Antoine equation PsatmmHg=10^(A-B/(C+T°C))
% rho_sol=659.4;                        %density of solvent kg.m-3
% M_sol=0.08617;                        %molar mass of solvent kg.mol-1
% Cp_sol = 2293;                        %heat capacity for solvent J.kg-1.K-1
% delta_H_sol=334e3;                    %heat for phase change boiling j.kg-1
% melting_point_sol=-93.5;              %melting temperature in °C
% sigma_sol=5.949;                      %Lennard Jones parameter in A  
% eps_sol=399.3;                        %Lennard Jones parameter in K

%Toluene properties in USI
A_sol=6.95464;                        %Antoine equation PsatmmHg=10^(A-B/(C+T°C))
B_sol=1344.8;                         %Antoine equation PsatmmHg=10^(A-B/(C+T°C))
C_sol=219.482;                        %Antoine equation PsatmmHg=10^(A-B/(C+T°C))
rho_sol=870;                          %density of solvent kg.m-3
M_sol=0.09214;                        %molar mass of solvent kg.mol-1
Cp_sol=1703;                         %heat capacity for solvent J.kg-1.K-1
delta_H_sol=401.6e3;                  %heat for phase change boiling j.kg-1
melting_point_sol=-95.0;              %melting temperature in °C
sigma_sol=5.92;                       %Lennard Jones parameter in A  
eps_sol=410;                          %Lennard Jones parameter in K

%%%%%%%%%%%%%%%%%%%%%DO NOT MODIFY ANYTHING AFTER THIS LINE%%%%%%%%%%%%%%%%

R=8.314;                                %gas constant in J.mol-1.K-1
NA=6.02214129e23;                       %Avogadro number mol-1.
initial_mass_prec=(1-initial_x_solvent)*initial_mass;
initial_mass_sol=initial_x_solvent*initial_mass;
rho_droplet=(1-initial_x_solvent)*rho_prec+initial_x_solvent*rho_sol;
droplet_initial_volume=(4/3)*pi*(initial_droplet_diameter/2)^3;
droplet_initial_mass=droplet_initial_volume*rho_droplet;
N_droplets=initial_mass/droplet_initial_mass;
maximal_pressure=(initial_mass_prec/M_prec+initial_mass_sol/M_sol)*R*reactor_temperature/reactor_volume;
maximal_prec_pressure=(initial_mass_prec/M_prec)*R*reactor_temperature/reactor_volume;
maximal_sol_pressure=(initial_mass_sol/M_sol)*R*reactor_temperature/reactor_volume;
gas_residence_time=reactor_volume/pumping_rate;
droplet_ini_prec_mass=initial_mass_prec/N_droplets;
droplet_ini_sol_mass=initial_mass_sol/N_droplets;

disp(['Maximal total pressure if 100% flash vaporization (Pa)= ',num2str(maximal_pressure)]);
disp(['Maximal precursor pressure if 100% flash vaporization (Pa)= ',num2str(maximal_prec_pressure)]);

%calculation of the semi-equilibrium flash evaporation state (temperature,
%pressure and amount evaprated X 

weigth_fraction_flash_vaporized=1;
for k=1:1:100;
pseudo_equilibrium_pressure=(weigth_fraction_flash_vaporized*initial_mass/M_sol)*R*reactor_temperature/reactor_volume;
log10pmmHg=log10(pseudo_equilibrium_pressure/133.322);
pseudo_equilibrium_temp=-B_sol/(log10pmmHg-A_sol)-C_sol;
pseudo_equilibrium_temp_K=pseudo_equilibrium_temp+273.15;
weigth_fraction_flash_vaporized=Cp_sol*(initial_droplet_temperature-pseudo_equilibrium_temp_K)/delta_H_sol;
end;

disp('---------------------------------------------------------------------');
disp('State after the flash vaporization (relevant only if X<50%) :');
disp(['Pseudo-equilibrium mass fraction vaporized % = ',num2str(weigth_fraction_flash_vaporized*100)]);
disp(['Pseudo-equilibrium pressure (Pa) = ',num2str(pseudo_equilibrium_pressure)]);
disp(['Pseudo-equilibrium temperature (°C) = ',num2str(pseudo_equilibrium_temp)]);

t_ini=0;
%we compute what appends at a droplet scale (mass) and at the reactor scale
%(pressures)
y0=[droplet_ini_prec_mass,droplet_ini_sol_mass,initial_prec_pressure,initial_sol_pressure,initial_droplet_temperature,0,0];
time=[0 pulse_duration];

[t,y]=ode45('droplet',time,y0);

temps=t;

droplet_prec_mass=y(:,1);
droplet_sol_mass=y(:,2);
droplet_mass=droplet_prec_mass+droplet_sol_mass;
total_pressure=y(:,3)+y(:,4);
prec_pressure=y(:,3);
sol_pressure=y(:,4);
temperature=y(:,5);
droplet_position=y(:,6);
droplet_velocity=y(:,7);

sol_mass_fraction=droplet_sol_mass./(droplet_prec_mass+droplet_sol_mass);

    time_evap=0;
    flag=0;
    indice_evap=0;
    for i=2:1:size(droplet_mass);
        if droplet_mass(i)/droplet_mass(1)<=1e-6 && flag==0;
            time_evap=temps(i);
            flag=1;
            indice_evap=i;
        end;
    end;
    
    time_reach=0;
    flag=0;
    indice_reach=0;
    for i=2:1:size(droplet_position);
        if abs(droplet_position(i))-susceptor_distance>0 && flag==0;
            time_reach=temps(i);
            flag=1;
            indice_reach=i;
        end;
    end;
    
    if indice_evap==0;
        indice_evap=length(droplet_mass);
        time_evap=pulse_duration;
        warning('WARNING: Liquid phase still present at the end of pulse !!!')
    end;
        
    %calculation of the droplet diameter and surface
    sol_volume=droplet_sol_mass./rho_sol;
    prec_volume=droplet_prec_mass./rho_prec;    
    V_droplet=sol_volume+prec_volume;
    droplet_diameter=2.*(3.*V_droplet./(4.*pi)).^(1/3);
    droplet_surface=4.*pi.*(3.*V_droplet./(4.*pi)).^(2/3);
    
    %calculation of the molar fraction of solvent
    moles_sol=droplet_sol_mass./M_sol;
    moles_prec=droplet_prec_mass./M_prec;
    x_mole_sol=moles_sol./(moles_sol+moles_prec);
    x_mole_prec=1-x_mole_sol;
    
    %calculation of the flux of liquid solvant to the gas phase
    Psat_sol_mmHg=10.^(A_sol-B_sol./(C_sol+(temperature-273.15)));
    vapor_pressure_sol=133.322*Psat_sol_mmHg;
    phi_sol_vaporated=x_mole_sol.*abs(droplet_surface.*M_sol.*(2./(2-alpha_evap).*(1./(2.*pi.*M_sol.*R).^0.5).*(alpha_evap.*sol_pressure./reactor_temperature.^0.5-alpha_evap.*vapor_pressure_sol./temperature.^0.5)));

    %calculation of the flux of liquid precursor to the gas phase
    Psat_prec_mmHg=10.^(A_prec-B_prec./(C_prec+(temperature-273.15)));
    vapor_pressure_prec=133.322*Psat_prec_mmHg;
    phi_prec_vaporated=x_mole_prec.*abs(droplet_surface.*M_prec.*(2./(2-alpha_evap).*(1./(2.*pi.*M_prec.*R).^0.5).*(alpha_evap.*prec_pressure./reactor_temperature.^0.5-alpha_evap.*vapor_pressure_prec./temperature.^0.5)));
    
    figure(1)
    plot(temps(1:indice_evap)/time_evap,droplet_diameter(1:indice_evap)/initial_droplet_diameter,'-k.');
    %title('Droplet diameter','FontSize',16);
    xlabel('Time/vaporization time','FontSize',16);
    ylabel('Droplet diameter/Initial droplet diameter','FontSize',16);
    axis([0,1,0,1]);
    set(gca,'FontSize',16)
%    text(time_evap,0,['   100% vaporization time=',num2str(time_evap,6),' s'],'VerticalAlignment','top','Rotation',90);
    saveas(gcf,'Figure_1.png');
    
    figure(2)
    plot(temps(1:indice_evap)/time_evap,(temperature(1:indice_evap)-273.15),'-k.');
%     hold on
%     plot([0 time_evap],[melting_point_prec melting_point_prec],'k--*');
%     plot([0 time_evap],[melting_point_sol melting_point_sol],'k--°');
%     plot([0 time_evap],[reactor_temperature-273 reactor_temperature-273],'k--.');
%     hold off
    %title('Droplet temperature','FontSize',16);
    xlabel('Time/vaporization time','FontSize',16);
    ylabel('Droplet temperature in °C','FontSize',16);
    set(gca,'FontSize',16)
%    legend('Droplet','Precursor freezing','Solvent freezing','Reactor wall');
    saveas(gcf,'Figure_2.png');
    
    figure(3)
    semilogy(temps/pulse_duration,total_pressure/maximal_pressure,'-k.',temps/pulse_duration,prec_pressure/maximal_pressure,'-r.');
%     hold on
%     plot([time_evap time_evap],[0 max(total_pressure)],'k--');
%     hold off
    %title('Reactor pressure','FontSize',16);
    xlabel('Time/pulse duration','FontSize',16);
    ylabel('Reactor pressure/maximal pressure','FontSize',16);
    set(gca,'FontSize',16)
%    legend('Total','Droplet vaporization','Precursor');
    saveas(gcf,'Figure_3.png');
    
    figure(4)
    semilogy(temps(1:indice_evap)/time_evap,phi_sol_vaporated(1:indice_evap),'k.');
    hold on
    semilogy(temps(1:indice_evap)/time_evap,phi_prec_vaporated(1:indice_evap),'r.');
    hold off
    %title('Droplet evaporation flux');
    xlabel('Time/vaporization time');
    ylabel('Flux in kg.s');
    legend('Solvent','Precursor');
    set(gca,'FontSize',16)
    axis([0 inf 1e-16 inf]);
    saveas(gcf,'Figure_4.png');
        
    figure(5)
    semilogy(temperature-273.15,vapor_pressure_sol,'k.');
    hold on
    semilogy(temperature-273.15,vapor_pressure_prec,'r.');
    hold off
    %title('Vapor pressures');
    xlabel('Temperature in °C');
    ylabel('Pressure in Pa');
    legend('Solvent','Precursor');
    set(gca,'FontSize',16)
    saveas(gcf,'Figure_5.png');
    
    figure(6)
    plot(temps(1:indice_evap)/time_evap,N_droplets.*droplet_mass(1:indice_evap),'k.');
    hold on
    plot(temps(1:indice_evap)/time_evap,N_droplets.*droplet_prec_mass(1:indice_evap),'r.');
    hold off
    %title('Mass vs time');
    xlabel('Time/vaporization time');
    ylabel('Mass in kg');
    legend('Solvent + precursor','Precursor');
    set(gca,'FontSize',16)
    saveas(gcf,'Figure_6.png');
    
    figure(7)
    plot(temps(1:indice_evap)/time_evap,droplet_position(1:indice_evap),'k.'); 
    hold on
    plot([0 time_evap],[-susceptor_distance -susceptor_distance],'k--');
    hold off
    %title('Droplet position vs time');
    xlabel('Time/vaporization time');
    ylabel('Droplet vertical position in m');
    legend('Particle position','Susceptor position');
    set(gca,'FontSize',16)
    saveas(gcf,'Figure_7.png');
    
    reactor_gas_velocity=-pumping_rate/(pi*reactor_diameter^2/4);
    
    figure(8)
    plot(temps(1:indice_evap)/time_evap,-droplet_velocity(1:indice_evap)/reactor_gas_velocity,'k.'); 
    hold on
    plot([0 time_evap/time_evap],[-reactor_gas_velocity/reactor_gas_velocity  -reactor_gas_velocity/reactor_gas_velocity],'k--');
    hold off
    %title('Droplet velocity vs time');
    xlabel('Time/vaporization time');
    ylabel('Droplet velocity/reactor gas velocity');
    set(gca,'FontSize',16)
    legend('Particle velocity','Gas velocity');
    saveas(gcf,'Figure_8.png');

disp('--------------------------------------------------------------------- ');
disp('Case of a patatoidal reactor :');
disp(['Gas residence time (s) = ',num2str(gas_residence_time)]);
disp(['Droplet 100 % vaporization time (s) = ',num2str(time_evap)]);
disp(['Droplet 100 % vaporization time / Gas residence time = ',num2str(time_evap/gas_residence_time)]);
disp('--------------------------------------------------------------------- ');
disp('Case of a plug flow reactor :');
disp(['Distance to travel to susceptor (m) = ',num2str(susceptor_distance)]);
if not(indice_reach==0);
disp(['Time to reach susceptor (s) = ',num2str(temps(indice_reach))]);
end;
disp(['Falling distance before 100% vaporization (m) = ',num2str(abs(droplet_position(indice_evap)))]);
disp(['Falling distance / Distance to susceptor = ',num2str(abs(droplet_position(indice_evap))/susceptor_distance)]);
disp('--------------------------------------------------------------------- ');
disp(['DeltaH/Cp solvent : ',num2str(delta_H_sol/Cp_sol)]);
disp(['DeltaH/Cp precursor : ',num2str(delta_H_prec/Cp_prec)]);
disp('--------------------------------------------------------------------- ');
disp(['Matlab calculation time : ',num2str(toc)]);

    Psat_sol_mmHg=10.^(A_sol-B_sol./(C_sol+(300-273.15)));
    sol_vap_300=133.322*Psat_sol_mmHg;
disp(['Solvent vapor pressure at 300K : ', num2str(sol_vap_300)]);

    Psat_prec_mmHg=10.^(A_prec-B_prec./(C_prec+(300-273.15)));
    prec_vap_300=133.322*Psat_prec_mmHg;
disp(['Precursor vapor pressure at 300K : ', num2str(prec_vap_300)]);
disp('--------------------------------------------------------------------- ');
