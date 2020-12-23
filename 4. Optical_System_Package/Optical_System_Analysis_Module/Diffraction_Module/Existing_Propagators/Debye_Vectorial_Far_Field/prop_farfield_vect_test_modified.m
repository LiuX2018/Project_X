%
%  Test Mansuripur's implmentation
%
%   27.07.08  H.Gross
%   15.07.13  G.Fente
%  
% modified by Xin Liu
% Email: liuxin2018@zju.edu.cn
% Dec.18, 2020
clear; clc;
% close all;


pupilRes = 128;  % resolution of pupil
wavelength = 1;  % mm
NA = 0.95;  % Numerical Aperture
n = 1;

scalx = 20;
scaly = scalx ;
R0 = 1.0;  % radius of pupil
f0 = R0/(NA/n);  % focal length
RE = wavelength /NA/NA;
dz = .0*RE; % Defocus
ipol = 5;  % 1=lin x, 2=lin y, 3=rad, 4=azi, 5=circ

% pupil function
E0 = ones(pupilRes,pupilRes).*exp(1i*(vortexPhasePlate(pupilRes,1)-pi/2));
dx = 2*R0/(pupilRes-2);  % pixel size of pupil
xp = zeros(pupilRes,1);
for j=1:pupilRes
    xp(j) = -R0-dx+dx*(j-1);
end

dy = 2*R0/(pupilRes-2);
yp = zeros(pupilRes,1);
for j=1:pupilRes
    yp(j) = -R0-dy+dy*(j-1);
end

[xpm,ypm] = meshgrid(xp,yp);
rp = sqrt(xpm.^2+ypm.^2);
% if rcurv > 0 ; 
% efd = efd .* exp( 2*1i*pi/wl*(rcurv-sqrt(rcurv^2-rpm.^2)) ); end 
% //The if statement isn't nececessary for an aplanatic system.
% 
ind = find( rp > R0) ; E0(ind) = 0 ;
phi = atan2(ypm,xpm);
cophi = cos(phi);   siphi = sin(phi);

%% Polarization definition 
% Different polarization can be realized by varying the polarization number.
if ipol ==1  %x-polarization
   px = E0; py = zeros(pupilRes,pupilRes,1);
elseif ipol == 2  %y-polarization
   px = zeros(pupilRes,pupilRes,1); py = E0;
elseif ipol == 3  % radial polarization
   px = E0.*cophi; py = E0.*siphi;
elseif ipol == 4  % azithmutal polarization
   px = -E0.*siphi; py = E0.*cophi;
elseif ipol == 5  % left circular polarization
   px = E0; py = 1i*E0;
end

%% field calculation
xs = linspace(-5,5,200);
ys = xs;

tic;
[Ex,Ey,Ez,xs,ys] = prop_farfield_vect_modified(px,py,wavelength,xp,yp,...
    f0,dz,scalx,scaly,xs,ys);
toc
phsEx = mod(angle(Ex),2*pi);
phsEy = mod(angle(Ey),2*pi);
phsEz = mod(angle(Ez),2*pi);

Ix=abs(Ex).^2; Iy=abs(Ey).^2; Iz=abs(Ez).^2;
I = Ix + Iy +Iz;
Imax = max(I(:));

% Normalization by the peak intensity
Ix=Ix./Imax; Iy=Iy./Imax; Iz=Iz./Imax; I = I./Imax;
Ixmax=max(max(Ix)); Iymax=max(max(Iy)); Izmax=max(max(Iz));

%% plot
figure;
subplot(2,2,1);
imagesc(xs/wavelength,ys/wavelength,Ix);
axis image xy;
title(['Ix ',num2str(Ixmax)])
set(gcf,'Color',[1,1,1])
set( gca, 'FontSize' , 12, 'fontweight','bold' );

subplot(2,2,2);
imagesc(xs/wavelength,ys/wavelength,Iy);
axis image xy;
title(['Iy  ',num2str(Iymax)])
set(gcf,'Color',[1,1,1])
set( gca, 'FontSize' , 12, 'fontweight','bold' );

subplot(2,2,3);
imagesc(xs/wavelength,ys/wavelength,Iz);
axis image xy;
title(['Iz  ',num2str(Izmax)])
set(gcf,'Color',[1,1,1])
set( gca, 'FontSize' , 12, 'fontweight','bold' );

subplot(2,2,4);
imagesc(xs/wavelength,ys/wavelength,I);
axis image xy;
title(['I  ',num2str(max(max(I)))])
set(gcf,'Color',[1,1,1])
set( gca, 'FontSize' , 12, 'fontweight','bold' );

figure;
subplot(1,3,1);
imagesc(xs/wavelength,ys/wavelength,phsEx);
axis image xy;
title('phase Ex')
set(gcf,'Color',[1,1,1])
set( gca, 'FontSize' , 12, 'fontweight','bold' );

subplot(1,3,2);
imagesc(xs/wavelength,ys/wavelength,phsEy);
axis image xy;
title('phase Ey')
set(gcf,'Color',[1,1,1])
set( gca, 'FontSize' , 12, 'fontweight','bold' );

subplot(1,3,3);
imagesc(xs/wavelength,ys/wavelength,phsEz);
axis image xy;
title('phase Ez')
set(gcf,'Color',[1,1,1])
set( gca, 'FontSize' , 12, 'fontweight','bold' );
