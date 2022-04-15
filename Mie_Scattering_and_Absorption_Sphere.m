% Baldi Lab 15/04/2022
% This script reads the relative permittivity data of a material and uses
% Mie theory to compute the scattering, absorption, and extinction
% cross-sections for a spherical particle of that material embedded in a
% lossless medium. The permittivity file has to be a tab-delimited text
% file with three columns: energy (in eV), epsilon1, epsilon2
clear all
close all

% MANUAL INPUT
prompt = {'Sphere radius (nm)','Refractive index of the surrounding medium (air = 1, water = 1.333)'};
dlg_title = 'Input parameters';
num_lines = 1;
def = {'',''};
answer = inputdlg(prompt,dlg_title,num_lines,def);
r=str2double(answer{1});
index=str2double(answer{2});

% DEFAULT SETTINGS
Esteps = 1000; % Number of energy steps for the data interpolation
nmax = 10;     % Maximum order of the Bessel function 

% CALCULATIONS
read=dlmread(uigetfile('*','Select a dielectric function file')); % Load the permittivity file (energy in eV, eps1, eps2)
e = 1.60217646e-19; % Elementary charge in SI units
h = 6.626068e-34; % Planck's constant in SI units
hbar = 1.05457148e-34; % hbar in SI units
me = 9.10938215e-31; % Electron rest mass in SI units
c = 2.99792458e8; % Light speed in SI units

%Load Experimental data
eV_read=read(:,1);
e1_read=read(:,2);
e2_read=read(:,3);
Emin=eV_read(1);
Emax=eV_read(size(eV_read,1));

%Interpolates the tabulated experimental data to make the data set smooth
energy = (Emin:(Emax-Emin)/Esteps:Emax)';
e1_read = interp1(eV_read, e1_read, energy, 'spline');
e2_read = interp1(eV_read, e2_read, energy, 'spline');

%Converts the frequency data to wavelength in meters
lambda = h*c./(e*energy);

%Creates the wavenumber k
k = 2*pi*index./lambda;

%Creates the total permittivity values of the particle and the medium
n_read = (((e1_read.^2 + e2_read.^2).^(1/2) + e1_read)./2).^(1/2);
k_read = (((e1_read.^2 + e2_read.^2).^(1/2) - e1_read)./2).^(1/2);
etot_read = n_read + 1i*k_read;
emed= (index^2)*ones(size(etot_read));
m=etot_read./index;

radius = r*1e-9; % Convertion to meters
x = k.*radius; % Size parameter
mx = m.*x;
scaele = 0; % initialise scattering matrix element
extele = 0; % initialise extinction matrix element
for n=1:nmax
    jnx = sqrt(pi./(2.*x)).*besselj(n+0.5,x);
    jnminx = sqrt(pi./(2.*x)).*besselj(n-0.5,x);
    jnmx = sqrt(pi./(2.*mx)).*besselj(n+0.5,mx);
    jnminmx = sqrt(pi./(2.*mx)).*besselj(n-0.5,mx);
    hnx = sqrt(pi./(2.*x)).*besselh(n+0.5,x);
    hnminx = sqrt(pi./(2.*x)).*besselh(n-0.5,x);
    xjnxdiff = x.*jnminx - n.*jnx;
    mxjnmxdiff = mx.*jnminmx - n.*jnmx;
    xhnxdiff = x.*hnminx - n.*hnx;
    an = (m.^2.*jnmx.*xjnxdiff-jnx.*mxjnmxdiff)./(m.^2.*jnmx.*xhnxdiff-hnx.*mxjnmxdiff);
    bn = (jnmx.*xjnxdiff-jnx.*mxjnmxdiff)./(jnmx.*xhnxdiff-hnx.*mxjnmxdiff);
    scaele = scaele + (2*n+1).*(an.*conj(an)+bn.*conj(bn));
    extele = extele + (2*n+1).*real(an+bn);
end

% Calculates scattering, extinction, and absorption cross sections
Csca = 2*pi./(k.^2).*scaele;
Cext = 2*pi./(k.^2).*extele;
Cabs = Cext-Csca;

% Calculates scattering, extinction, and absorption efficiencies
Qsca = Csca./(pi*radius^2);
Qext = Cext./(pi*radius^2);
Qabs = Cabs./(pi*radius^2);

% Plots
% Relative permittivity
figure(1)
hold on
plot(read(:,1),read(:,2),'ok')
plot(read(:,1),read(:,3),'or')
plot(energy, e1_read,'k','linewidth',2);
plot(energy, e2_read,'r','linewidth',2);
xlabel ('energy (eV)', 'FontSize',10);
title ('relative permittivity, \epsilon', 'FontSize',10);
leg1=legend('epsilon1 - data','epsilon2 - data','epsilon1 - interpolated','epsilon2 - interpolated','Location','SouthEast');

% Scattering
figure (2)
hold on
plot(energy,Csca,'k','linewidth',2);
yline(pi*radius^2,'--');
xlabel('energy (eV)', 'FontSize',10);
title(['scattering cross-section of a ', num2str(radius*1E9), ' nm radius sphere'], 'FontSize',10);
leg1=legend('scattering cross-section','geometric cross-section');

% Absorption
figure (3)
hold on
plot(energy,Cabs,'k','linewidth',2);
yline(pi*radius^2,'--');
xlabel('energy (eV)', 'FontSize',10);
title(['absorption cross-section of a ', num2str(radius*1E9), ' nm radius sphere'], 'FontSize',10);
leg1=legend('absorption cross-section','geometric cross-section');

% Extinction
figure (4)
hold on
plot(energy,Cext,'k','linewidth',2);
yline(pi*radius^2,'--');
xlabel('energy (eV)', 'FontSize',10);
title(['extinction cross-section of a ', num2str(radius*1E9), ' nm radius sphere'], 'FontSize',10);
leg1=legend('extinction cross-section','geometric cross-section');
