%% monte carlo simulation of light propagation in tissue from EO-Flex
close all;
NA=0.28;
r0=1.9;

% Pinput=0.208; %mW
% Pinput=0.189; %mW
% Pinput=0.152; %mW
% Pinput=0.061; %mW
% Pinput=0.029; %mW
% Pinput=0.0057; %mW

% Pinput=0.0570; %mW
% [B, Irr0, B_irr ]=Intensity_profile(NA, r0,Pinput);

Pinput=[0.208; 0.189; 0.152; 0.061; 0.029; 0.0057;];
% Pinput=[0.208];
rmax =2500; % 2*r0; % Transverse plot dimensions(micrometers)
zmax =rmax*2; % 6*r0; % Longitudinal plot dimensions(micrometers)
dr =1; 
dz =1;

Irradiances=cell(1,2);
for i=1:length(Pinput)
   
    [B, R, A Irr0, B_irr, Bnorm]=Intensity_profile(NA, r0,Pinput(i), rmax, zmax, dr, dz);
    PlotIntensity(B, Irr0, B_irr, Bnorm,Pinput(i), rmax, zmax,dr,dz)
    Irradiances{1,1}=[Irradiances{1,1}; Pinput(i)];
    Irradiances{1,2}=[Irradiances{1,2}; Irr0];
    
    % simulate and plot a zoomed in region of the fiber
    rmax =50; % 2*r0; % Transverse plot dimensions(micrometers)
    zmax =100; % 6*r0; % Longitudinal plot dimensions(micrometers)
    dz = 0.05; % z increment (micrometers)
    dr = 0.05; % r increment (micrometers)
    [Bzoom, Rzoom, Azoom, Irr0, B_irr, Bnorm]=Intensity_profile(NA, r0,Pinput(i), rmax, zmax, dr, dz);
    PlotIntensity(Bzoom, Irr0, B_irr, Bnorm,Pinput(i), rmax, zmax, dr, dz)
end
save(sprintf('Irradiances_r0_%f.mat',r0), 'Irradiances');

%%
function [B, R, A, Irr0, B_irr, Bnorm]= Intensity_profile(NA, r0,Pinput, rmax, zmax, dr, dz)
n0 = 1.36; % Tissue refractive index
% dz increment (micrometers)
% dr increment (micrometers)

% from Foutz et al
K = 0.000125; % absorption coefficient (um^-1)
S = 0.00737; % Scattering coefficient (um^-1)
i=0;
j=0; %indices
theta = asin(NA/n0); % Divergence angle
for z = 0:dz:zmax
    z0(i+1)=i*dz;
    i=i+1;
    for r = -rmax:dr:rmax
        j=j+1;
        geometric = r0/(r0 + z*tan(theta)); %Geometric 
        gaussian = (1/(2*pi))*exp(-2*(r/(r0 + z*tan(theta)))^2); %Gaussian 
        a = 1+(K/S);
        b = ((a^2)-1)^0.5;
        transmission = b/((a*sinh(b*S*((r^2+z^2)^0.5)))+(b*cosh(b*S*((r^2+z^2)^0.5)))); %Transmission after Absorption and scattering
        reflection = sinh(b*S*((r^2+z^2)^0.5))*((a*sinh(b*S*((r^2+z^2)^0.5)))+(b*cosh(b*S*((r^2+z^2)^0.5)))); 
        B(j,i)=geometric*gaussian*transmission;
        R(j,i)=geometric*gaussian*reflection;
        A(j,i)=geometric*gaussian*(1-(reflection+transmission));
    end
    j=0;
end
%display image with probe pointing down
B=B';
R=R';
A=A';

Irr0=Pinput/(pi*(r0*10^-3)^2); %mW/mm^2
ChR2_thr=1; % mW/mm^2
Bnorm=(B-min(min(B)))/(max(max(B))-min(min(B)));
% irradiance for a given Pinput
B_irr=Bnorm*Irr0;
% map=map(80:(end), :);
% input_intensity=
end

function PlotIntensity(B, Irr0, B_irr, Bnorm, Pinput, rmax, zmax,dr ,dz);
Fontsize=44;
FontsizeAxis=30;
ChR2_thr=1; % mW/mm^2

convert_factor=1000; % scale microns to mm

map=parula;

% figure;
% imagesc((-rmax:rmax)/convert_factor,(0:zmax)/convert_factor,B);
% %colormap(hot);
% colormap(map);
% ax=ancestor(gca, 'axes');
% xrule=ax.XAxis;
% yrule=ax.YAxis;
% 
% xL=xlabel('Axial distance (mm)');
% % xticks([PCA_time(1) round(PCA_time(end))]);
% yL=ylabel('Transverse distance (mm)');
% 
% xrule.FontSize =FontsizeAxis;
% yrule.FontSize =FontsizeAxis;
% xL.FontSize=Fontsize;
% yL.FontSize=Fontsize;
% 
% title(sprintf('Irradiance profile %i mW', Pinput));
% colorbar;
% 
%Normalized irradiance
figure;

imagesc((-rmax:rmax)/convert_factor,(0:zmax)/convert_factor,Bnorm);
hold on;

[~, cthr]=contour((-rmax:dr:rmax)/convert_factor,(0:dz:zmax)/convert_factor,Bnorm, [(1/Irr0) (1/Irr0)], 'k--');
cthr.LineWidth=2;
hold off;

colormap(map);
ax=ancestor(gca, 'axes');
xrule=ax.XAxis;
yrule=ax.YAxis;

xL=xlabel(['Axial distance (mm)'], 'Interpreter', 'tex');

% xticks([PCA_time(1) round(PCA_time(end))]);
yL=ylabel(['Tissue Depth (mm)'], 'Interpreter', 'tex');


xrule.FontSize =FontsizeAxis;
yrule.FontSize =FontsizeAxis;
xL.FontSize=Fontsize;
yL.FontSize=Fontsize;

title('Irradiance profile');
c=colorbar;
c.FontSize=26;
c.Label.String='Normalized Irradiance';
% set(gca, 'FontSize', 26);
height=7;
width=9;

set(gcf, 'Position', [100 100 width*100, height*100]); %<- Set size
set(gcf,'color','w');
print('-painters','-bestfit',sprintf('%f_mW_normalized_%i_HalfWindow.pdf',Pinput,rmax),'-dpdf');
print('-painters',sprintf('%f_mW_normalized_%i_HalfWindow.svg',Pinput,rmax),'-dsvg');


% figure;
% 
% imagesc((-rmax:rmax)/convert_factor,(0:zmax)/convert_factor,B_irr);
% hold on;
% [~, cthr]=contour((-rmax:dr:rmax)/convert_factor,(0:dz:zmax)/convert_factor,B_irr, [ChR2_thr ChR2_thr], 'k-');
% % [~, cthr_max]=contour((-rmax:dr:rmax),(0:dz:zmax),B_irr, [300 300], 'k');
% 
% cthr.LineWidth=2;
% cthr_max.LineWidth=2;
% % caxis([0 19345]);
% hold off;
% colormap(map);
% 
% ax=ancestor(gca, 'axes');
% xrule=ax.XAxis;
% yrule=ax.YAxis;
% 
% % xL=xlabel(['Axial distance (',char(181),'m)'], 'Interpreter', 'tex');
% % 
% % % xticks([PCA_time(1) round(PCA_time(end))]);
% % yL=ylabel(['Distance (',char(181),'m)'], 'Interpreter', 'tex');
% xL=xlabel(['Axial distance (mm)'], 'Interpreter', 'tex');
% yL=ylabel(['Distance (mm)'], 'Interpreter', 'tex');
% xrule.FontSize =FontsizeAxis;
% yrule.FontSize =FontsizeAxis;
% xL.FontSize=Fontsize;
% yL.FontSize=Fontsize;
% 
% c=colorbar;
% c.FontSize=26;
% c.Label.String='Irradiance (mW mm^{-2})';
% % ylabel(c,'Irradiance (mW/mm^{2})');
% set(gca, 'FontSize', 26);
% height=7;
% width=9;
% set(gcf, 'Position', [100 100 width*100, height*100]); %<- Set size
% set(gcf,'color','w');
% % xlim([-rmax rmax]);
% % xticks([-rmax -rmax/2 0 rmax/2 rmax])
% % xticklabels({-rmax,-rmax/2,0,rmax/2,rmax})
% print('-painters','-bestfit',sprintf('%f_mW_profile_%i_cutoff.pdf',Pinput,ChR2_thr),'-dpdf');
% print('-painters',sprintf('%f_mW_profile_%i_cutoff.svg',Pinput,ChR2_thr),'-dsvg');

ChR2_thr1=1; % mW/mm^2
% irradiance for a given Pinput
% figure;
% B_irr=Bnorm*Irr0;
% imagesc((-rmax:rmax)/convert_factor,(0:zmax)/convert_factor,B_irr);
% hold on;
% [~, cthr]=contour((-rmax:dr:rmax)/convert_factor,(0:dz:zmax)/convert_factor,B_irr, [ChR2_thr ChR2_thr], 'k-');
% % [~, cthr_max]=contour((-rmax:dr:rmax),(0:dz:zmax),B_irr, [300 300], 'k');
% 
% cthr.LineWidth=2;
% cthr_max.LineWidth=2;
% % caxis([0 19345]);
% hold off;
% colormap(map);
% 
% % c=colorbar;
% % c.FontSize=26;
% % c.Label.String='Irradiance (mW mm^{-2})';
% % ylabel(c,'Irradiance (mW/mm^{2})');
% % set(gca, 'FontSize', 26);
% ax=ancestor(gca, 'axes');
% xrule=ax.XAxis;
% yrule=ax.YAxis;
% 
% xL=xlabel(['Axial distance (mm)'], 'Interpreter', 'tex');
% 
% % xticks([PCA_time(1) round(PCA_time(end))]);
% yL=ylabel(['Tissue Depth (mm)'], 'Interpreter', 'tex');
% 
% 
% xrule.FontSize =FontsizeAxis;
% yrule.FontSize =FontsizeAxis;
% xL.FontSize=Fontsize;
% yL.FontSize=Fontsize;
% height=7;
% width=7;
% set(gcf, 'Position', [100 100 width*100, height*100]); %<- Set size
% set(gcf,'color','w');
% % xlim([-rmax rmax]);
% % xticks([-rmax -rmax/2 0 rmax/2 rmax])
% % xticklabels({-rmax,-rmax/2,0,rmax/2,rmax})
% print('-painters','-bestfit',sprintf('%f_mW_profile_%i_cutoff_no_legend.pdf',Pinput,ChR2_thr),'-dpdf');
% print('-painters',sprintf('%f_mW_profile_%i_cutoff_no_legend.svg',Pinput,ChR2_thr),'-dsvg');

% ChR2_thr=10; % mW/mm^2
% % irradiance for a given Pinput
% figure;
% B_irr=Bnorm*Irr0;
% imagesc((-rmax:rmax)/convert_factor,(0:zmax)/convert_factor,B_irr);
% hold on;
% [~, cthr]=contour((-rmax:dr:rmax)/convert_factor,(0:dz:zmax)/convert_factor,B_irr, [ChR2_thr ChR2_thr], 'k-');
% % [~, cthr_max]=contour((-rmax:dr:rmax),(0:dz:zmax),B_irr, [300 300], 'k');
% 
% cthr.LineWidth=2;
% cthr_max.LineWidth=2;
% % caxis([0 19345]);
% hold off;
% colormap(map);
% ax=ancestor(gca, 'axes');
% xrule=ax.XAxis;
% yrule=ax.YAxis;
% 
% xL=xlabel(['Axial distance (mm)'], 'Interpreter', 'tex');
% yL=ylabel(['Tissue Depth (mm)'], 'Interpreter', 'tex');
% xrule.FontSize =FontsizeAxis;
% yrule.FontSize =FontsizeAxis;
% xL.FontSize=Fontsize;
% yL.FontSize=Fontsize;
% 
% c=colorbar;
% c.FontSize=26;
% c.Label.String='Irradiance (mW mm^{-2})';
% % ylabel(c,'Irradiance (mW/mm^{2})');
% % set(gca, 'FontSize', 26);
% height=7;
% width=9;
% set(gcf, 'Position', [100 100 width*100, height*100]); %<- Set size
% set(gcf,'color','w');
% % xlim([-rmax rmax]);
% % xticks([-rmax -rmax/2 0 rmax/2 rmax])
% % xticklabels({-rmax,-rmax/2,0,rmax/2,rmax})
% print('-painters','-bestfit',sprintf('%f_mW_profile_%i_cutoff.pdf',Pinput,ChR2_thr),'-dpdf');
% print('-painters',sprintf('%f_mW_profile_%i_cutoff.svg',Pinput,ChR2_thr),'-dsvg');

ChR2_thr2=10; % mW/mm^2
% irradiance for a given Pinput with both dashed lines for thr
figure;
B_irr=Bnorm*Irr0;
imagesc((-rmax:rmax)/convert_factor,(0:zmax)/convert_factor,B_irr);
hold on;
if min(B_irr) < ChR2_thr2
    [~, cthr]=contour((-rmax:dr:rmax)/convert_factor,(0:dz:zmax)/convert_factor,B_irr, [ChR2_thr1 ChR2_thr2], 'k-');
    [~, cthr_max]=contour((-rmax:dr:rmax)/convert_factor,(0:dz:zmax)/convert_factor,B_irr, [1 1], 'k--');
    cthr.LineWidth=2;
    cthr_max.LineWidth=2;

else
    [~, cthr]=contour((-rmax:dr:rmax)/convert_factor,(0:dz:zmax)/convert_factor,B_irr, [1 1], 'k--');
    cthr.LineWidth=2;
end

% caxis([0 19345]);
hold off;
colormap(map);
ax=ancestor(gca, 'axes');
xrule=ax.XAxis;
yrule=ax.YAxis;

xL=xlabel(['Axial distance (mm)'], 'Interpreter', 'tex');
yL=ylabel(['Tissue Depth (mm)'], 'Interpreter', 'tex');
xrule.FontSize =FontsizeAxis;
yrule.FontSize =FontsizeAxis;
xL.FontSize=Fontsize;
yL.FontSize=Fontsize;

c=colorbar;
c.FontSize=26;
c.Label.String='Irradiance (mW mm^{-2})';
% ylabel(c,'Irradiance (mW/mm^{2})');
% set(gca, 'FontSize', 26);
height=7;
width=9;
set(gcf, 'Position', [100 100 width*100, height*100]); %<- Set size
set(gcf,'color','w');
% xlim([-rmax rmax]);
% xticks([-rmax -rmax/2 0 rmax/2 rmax])
% xticklabels({-rmax,-rmax/2,0,rmax/2,rmax})
print('-painters','-bestfit',sprintf('%f_mW_profile_both_%i_cutoff_%i_dr_%i_HalfWindow.pdf',Pinput,ChR2_thr,dr,rmax),'-dpdf');
print('-painters',sprintf('%f_mW_profile_both_%i_cutoff_%i_dr_%i_HalfWindow.svg',Pinput,ChR2_thr,dr,rmax),'-dsvg');


% plot tilted profile
ChR2_thr=1; % mW/mm^2
% irradiance for a given Pinput
figure;
B_irr=Bnorm*Irr0;

B_irr_rot=imrotate(B_irr,45);

[val,pos]=max(max(B_irr_rot,[],1))
xmax=pos;
[val,pos]=max(max(B_irr_rot,[],2))
ymax=pos;
% B_irr_rot=B_irr_rot(xmax:(xmax+length(B_irr)-1),1:(1+length(B_irr)-1));
B_irr_rot=B_irr_rot(1:(1+length(B_irr)-1),1:(1+length(B_irr)-1));


imagesc((-rmax:rmax)/convert_factor,(0:zmax)/convert_factor,B_irr_rot);
hold on;
%if min(B_irr) < ChR2_thr
    [~, cthr]=contour((-rmax:dr:rmax)/convert_factor,(0:dz:zmax)/convert_factor,B_irr_rot, [ChR2_thr ChR2_thr], 'k-');

    cthr.LineWidth=2;
    cthr_max.LineWidth=2;
%end
% caxis([0 19345]);
hold off;
colormap(map);

ax=ancestor(gca, 'axes');
xrule=ax.XAxis;
yrule=ax.YAxis;

xL=xlabel(['Axial distance (mm)'], 'Interpreter', 'tex');
yL=ylabel(['Tissue Depth (mm)'], 'Interpreter', 'tex');
xrule.FontSize =FontsizeAxis;
yrule.FontSize =FontsizeAxis;
xL.FontSize=Fontsize;
yL.FontSize=Fontsize;
c=colorbar;
c.FontSize=26;
c.Label.String='Irradiance (mW mm^{-2})';
% c=colorbar;
% c.FontSize=26;
% c.Label.String='Irradiance (mW mm^{-2})';
% ylabel(c,'Irradiance (mW/mm^{2})');
% set(gca, 'FontSize', 26);
height=7;
width=7;
set(gcf, 'Position', [100 100 width*100, height*100]); %<- Set size
set(gcf,'color','w');
xlim([-rmax/convert_factor rmax/convert_factor]);
% xticks([-rmax/convert_factor -rmax/2/convert_factor 0 rmax/2/convert_factor rmax/convert_factor])
% xticklabels(round([0,2*rmax/4/convert_factor,4*rmax/4/convert_factor,6*rmax/4/convert_factor, 2*rmax/convert_factor],1))
if dr >= 1
    xticks([-rmax/convert_factor -0.75 -0.25 0.25 0.75 rmax/convert_factor])
    xticklabels(round([0,0.5,1,1.5, 2, 2.5],1))
end
print('-painters','-bestfit',sprintf('%f_mW_profile_%i_tilted_%i_dr_%i_HalfWindow.pdf',Pinput,ChR2_thr,dr,rmax),'-dpdf');
print('-painters',sprintf('%f_mW_profile_%i_tilted_%i_dr_%i_HalfWindow.svg',Pinput,ChR2_thr,dr,rmax),'-dsvg');

% %plot profile near probe tip
% rmax =50; % 2*r0; % Transverse plot dimensions(micrometers)
% zmax =100; % 6*r0; % Longitudinal plot dimensions(micrometers)
% dz = 0.05; % z increment (micrometers)
% dr = 0.05; % r increment (micrometers)
% 
% K = 0.000125; % absorption coefficient (um^-1)
% S = 0.00737; % Scattering coefficient (um^-1)
% i=0;
% j=0; %indices
% theta = asin(NA/n0); % Divergence angle
% for z = 0:dz:zmax
%     z0(i+1)=i*dz;
%     i=i+1;
%     for r = -rmax:dr:rmax
%         j=j+1;
%         geometric = r0/(r0 + z*tan(theta)); %Geometric factor
%         gaussian = (1/(2*pi))*exp(-2*(r/(r0 + z*tan(theta)))^2); %Gaussian factor
%         a = 1+(K/S);
%         b = ((a^2)-1)^0.5;
%         abs_scatter = b/((a*sinh(b*S*((r^2+z^2)^0.5)))+(b*cosh(b*S*((r^2+z^2)^0.5)))); %Absorption and scattering
%         Z(j,i)=geometric*gaussian*abs_scatter;
%     end
%     j=0;
% end
% %display image with probe pointing down
% Z=Z';
% Znorm=(Z-min(min(Z)))/(max(max(Z))-min(min(Z)));
% 
% Z_irr=Znorm*Irr0;
% 
% rbias=rmax/dr;
% 
% window_size=100;
% 
% rsmall=((-window_size/(dr*2)):(window_size/(dr*2)))+rbias+1;
% zsmall=((0):(window_size/dz))+1;
% 
% figure;
% imagesc(-window_size/2:window_size/2,0:window_size,Z_irr(zsmall,rsmall));
% %colormap(hot);
% colormap(map);
% ax=ancestor(gca, 'axes');
% xrule=ax.XAxis;
% yrule=ax.YAxis;
% 
% % xL=xlabel('Transverse distance (micons)');
% xL=xlabel(['Axial distance (',char(181),'m)'], 'Interpreter', 'tex');
% 
% yL=ylabel(['Depth (',char(181),'m)'], 'Interpreter', 'tex');
% 
% xrule.FontSize =FontsizeAxis;
% yrule.FontSize =FontsizeAxis;
% xL.FontSize=Fontsize;
% yL.FontSize=Fontsize;
% title('Irradiance profile');
% colorbar;
% % set(gca, 'FontSize', 26);
% height=7;
% width=9;
% set(gcf, 'Position', [100 100 width*100, height*100]); %<- Set size
% set(gcf,'color','w');
% brush('on');
% print('-painters','-bestfit',sprintf('%f_mW_zoomed_profile.pdf',Pinput),'-dpdf');
% print('-painters',sprintf('%f_mW_zoomed_profile.svg',Pinput),'-dsvg');
% figure;
% [m,n]=size(B);
% B0=B(round(m/2),:)/max(B(round(m/2),:));
% % Normalized value of the longitudinal profile
% for r=0;
%     plot(z0,B0);
%     xlabel('Axial distance (microns)');
%     ylabel('Normalized intensity at fiber tip center (N.U.)');
%     title('Longitudinal Irradiance profile');
% end

end
%%
function RGB=colormatch(wavelength)
% Converts Wavelength to RGB based on the color matching functions of
% Stiles and Burch (1959) and Speranskaya (1959)
% Modified from Spectral and XYZ Color Functions by Jeff Mather
% Copyright (c) 2010, Jeff Mather
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are
% met:
% 
%     * Redistributions of source code must retain the above copyright
%       notice, this list of conditions and the following disclaimer.
%     * Redistributions in binary form must reproduce the above copyright
%       notice, this list of conditions and the following disclaimer in
%       the documentation and/or other materials provided with the distribution
% 
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
% POSSIBILITY OF SUCH DAMAGE.

        cmf = [360,  0.000000122200,  0.000000013398,  0.000000535027
               361,  0.000000185138,  0.000000020294,  0.000000810720
               362,  0.000000278830,  0.000000030560,  0.000001221200
               363,  0.000000417470,  0.000000045740,  0.000001828700
               364,  0.000000621330,  0.000000068050,  0.000002722200
               365,  0.000000919270,  0.000000100650,  0.000004028300
               366,  0.000001351980,  0.000000147980,  0.000005925700
               367,  0.000001976540,  0.000000216270,  0.000008665100
               368,  0.000002872500,  0.000000314200,  0.000012596000
               369,  0.000004149500,  0.000000453700,  0.000018201000
               370,  0.000005958600,  0.000000651100,  0.000026143700
               371,  0.000008505600,  0.000000928800,  0.000037330000
               372,  0.000012068600,  0.000001317500,  0.000052987000
               373,  0.000017022600,  0.000001857200,  0.000074764000
               374,  0.000023868000,  0.000002602000,  0.000104870000
               375,  0.000033266000,  0.000003625000,  0.000146220000
               376,  0.000046087000,  0.000005019000,  0.000202660000
               377,  0.000063472000,  0.000006907000,  0.000279230000
               378,  0.000086892000,  0.000009449000,  0.000382450000
               379,  0.000118246000,  0.000012848000,  0.000520720000
               380,  0.000159952000,  0.000017364000,  0.000704776000
               381,  0.000215080000,  0.000023327000,  0.000948230000
               382,  0.000287490000,  0.000031150000,  0.001268200000
               383,  0.000381990000,  0.000041350000,  0.001686100000
               384,  0.000504550000,  0.000054560000,  0.002228500000
               385,  0.000662440000,  0.000071560000,  0.002927800000
               386,  0.000864500000,  0.000093300000,  0.003823700000
               387,  0.001121500000,  0.000120870000,  0.004964200000
               388,  0.001446160000,  0.000155640000,  0.006406700000
               389,  0.001853590000,  0.000199200000,  0.008219300000
               390,  0.002361600000,  0.000253400000,  0.010482200000
               391,  0.002990600000,  0.000320200000,  0.013289000000
               392,  0.003764500000,  0.000402400000,  0.016747000000
               393,  0.004710200000,  0.000502300000,  0.020980000000
               394,  0.005858100000,  0.000623200000,  0.026127000000
               395,  0.007242300000,  0.000768500000,  0.032344000000
               396,  0.008899600000,  0.000941700000,  0.039802000000
               397,  0.010870900000,  0.001147800000,  0.048691000000
               398,  0.013198900000,  0.001390300000,  0.059210000000
               399,  0.015929200000,  0.001674000000,  0.071576000000
               400,  0.019109700000,  0.002004400000,  0.086010900000
               401,  0.022788000000,  0.002386000000,  0.102740000000
               402,  0.027011000000,  0.002822000000,  0.122000000000
               403,  0.031829000000,  0.003319000000,  0.144020000000
               404,  0.037278000000,  0.003880000000,  0.168990000000
               405,  0.043400000000,  0.004509000000,  0.197120000000
               406,  0.050223000000,  0.005209000000,  0.228570000000
               407,  0.057764000000,  0.005985000000,  0.263470000000
               408,  0.066038000000,  0.006833000000,  0.301900000000
               409,  0.075033000000,  0.007757000000,  0.343870000000
               410,  0.084736000000,  0.008756000000,  0.389366000000
               411,  0.095041000000,  0.009816000000,  0.437970000000
               412,  0.105836000000,  0.010918000000,  0.489220000000
               413,  0.117066000000,  0.012058000000,  0.542900000000
               414,  0.128682000000,  0.013237000000,  0.598810000000
               415,  0.140638000000,  0.014456000000,  0.656760000000
               416,  0.152893000000,  0.015717000000,  0.716580000000
               417,  0.165416000000,  0.017025000000,  0.778120000000
               418,  0.178191000000,  0.018399000000,  0.841310000000
               419,  0.191214000000,  0.019848000000,  0.906110000000
               420,  0.204492000000,  0.021391000000,  0.972542000000
               421,  0.217650000000,  0.022992000000,  1.038900000000
               422,  0.230267000000,  0.024598000000,  1.103100000000
               423,  0.242311000000,  0.026213000000,  1.165100000000
               424,  0.253793000000,  0.027841000000,  1.224900000000
               425,  0.264737000000,  0.029497000000,  1.282500000000
               426,  0.275195000000,  0.031195000000,  1.338200000000
               427,  0.285301000000,  0.032927000000,  1.392600000000
               428,  0.295143000000,  0.034738000000,  1.446100000000
               429,  0.304869000000,  0.036654000000,  1.499400000000
               430,  0.314679000000,  0.038676000000,  1.553480000000
               431,  0.324355000000,  0.040792000000,  1.607200000000
               432,  0.333570000000,  0.042946000000,  1.658900000000
               433,  0.342243000000,  0.045114000000,  1.708200000000
               434,  0.350312000000,  0.047333000000,  1.754800000000
               435,  0.357719000000,  0.049602000000,  1.798500000000
               436,  0.364482000000,  0.051934000000,  1.839200000000
               437,  0.370493000000,  0.054337000000,  1.876600000000
               438,  0.375727000000,  0.056822000000,  1.910500000000
               439,  0.380158000000,  0.059399000000,  1.940800000000
               440,  0.383734000000,  0.062077000000,  1.967280000000
               441,  0.386327000000,  0.064737000000,  1.989100000000
               442,  0.387858000000,  0.067285000000,  2.005700000000
               443,  0.388396000000,  0.069764000000,  2.017400000000
               444,  0.387978000000,  0.072218000000,  2.024400000000
               445,  0.386726000000,  0.074704000000,  2.027300000000
               446,  0.384696000000,  0.077272000000,  2.026400000000
               447,  0.382006000000,  0.079979000000,  2.022300000000
               448,  0.378709000000,  0.082874000000,  2.015300000000
               449,  0.374915000000,  0.086000000000,  2.006000000000
               450,  0.370702000000,  0.089456000000,  1.994800000000
               451,  0.366089000000,  0.092947000000,  1.981400000000
               452,  0.361045000000,  0.096275000000,  1.965300000000
               453,  0.355518000000,  0.099535000000,  1.946400000000
               454,  0.349486000000,  0.102829000000,  1.924800000000
               455,  0.342957000000,  0.106256000000,  1.900700000000
               456,  0.335893000000,  0.109901000000,  1.874100000000
               457,  0.328284000000,  0.113835000000,  1.845100000000
               458,  0.320150000000,  0.118167000000,  1.813900000000
               459,  0.311475000000,  0.122932000000,  1.780600000000
               460,  0.302273000000,  0.128201000000,  1.745370000000
               461,  0.292858000000,  0.133457000000,  1.709100000000
               462,  0.283502000000,  0.138323000000,  1.672300000000
               463,  0.274044000000,  0.143042000000,  1.634700000000
               464,  0.264263000000,  0.147787000000,  1.595600000000
               465,  0.254085000000,  0.152761000000,  1.554900000000
               466,  0.243392000000,  0.158102000000,  1.512200000000
               467,  0.232187000000,  0.163941000000,  1.467300000000
               468,  0.220488000000,  0.170362000000,  1.419900000000
               469,  0.208198000000,  0.177425000000,  1.370000000000
               470,  0.195618000000,  0.185190000000,  1.317560000000
               471,  0.183034000000,  0.193025000000,  1.262400000000
               472,  0.170222000000,  0.200313000000,  1.205000000000
               473,  0.157348000000,  0.207156000000,  1.146600000000
               474,  0.144650000000,  0.213644000000,  1.088000000000
               475,  0.132349000000,  0.219940000000,  1.030200000000
               476,  0.120584000000,  0.226170000000,  0.973830000000
               477,  0.109456000000,  0.232467000000,  0.919430000000
               478,  0.099042000000,  0.239025000000,  0.867460000000
               479,  0.089388000000,  0.245997000000,  0.818280000000
               480,  0.080507000000,  0.253589000000,  0.772125000000
               481,  0.072034000000,  0.261876000000,  0.728290000000
               482,  0.063710000000,  0.270643000000,  0.686040000000
               483,  0.055694000000,  0.279645000000,  0.645530000000
               484,  0.048117000000,  0.288694000000,  0.606850000000
               485,  0.041072000000,  0.297665000000,  0.570060000000
               486,  0.034642000000,  0.306469000000,  0.535220000000
               487,  0.028896000000,  0.315035000000,  0.502340000000
               488,  0.023876000000,  0.323335000000,  0.471400000000
               489,  0.019628000000,  0.331366000000,  0.442390000000
               490,  0.016172000000,  0.339133000000,  0.415254000000
               491,  0.013300000000,  0.347860000000,  0.390024000000
               492,  0.010759000000,  0.358326000000,  0.366399000000
               493,  0.008542000000,  0.370001000000,  0.344015000000
               494,  0.006661000000,  0.382464000000,  0.322689000000
               495,  0.005132000000,  0.395379000000,  0.302356000000
               496,  0.003982000000,  0.408482000000,  0.283036000000
               497,  0.003239000000,  0.421588000000,  0.264816000000
               498,  0.002934000000,  0.434619000000,  0.247848000000
               499,  0.003114000000,  0.447601000000,  0.232318000000
               500,  0.003816000000,  0.460777000000,  0.218502000000
               501,  0.005095000000,  0.474340000000,  0.205851000000
               502,  0.006936000000,  0.488200000000,  0.193596000000
               503,  0.009299000000,  0.502340000000,  0.181736000000
               504,  0.012147000000,  0.516740000000,  0.170281000000
               505,  0.015444000000,  0.531360000000,  0.159249000000
               506,  0.019156000000,  0.546190000000,  0.148673000000
               507,  0.023250000000,  0.561180000000,  0.138609000000
               508,  0.027690000000,  0.576290000000,  0.129096000000
               509,  0.032444000000,  0.591500000000,  0.120215000000
               510,  0.037465000000,  0.606741000000,  0.112044000000
               511,  0.042956000000,  0.622150000000,  0.104710000000
               512,  0.049114000000,  0.637830000000,  0.098196000000
               513,  0.055920000000,  0.653710000000,  0.092361000000
               514,  0.063349000000,  0.669680000000,  0.087088000000
               515,  0.071358000000,  0.685660000000,  0.082248000000
               516,  0.079901000000,  0.701550000000,  0.077744000000
               517,  0.088909000000,  0.717230000000,  0.073456000000
               518,  0.098293000000,  0.732570000000,  0.069268000000
               519,  0.107949000000,  0.747460000000,  0.065060000000
               520,  0.117749000000,  0.761757000000,  0.060709000000
               521,  0.127839000000,  0.775340000000,  0.056457000000
               522,  0.138450000000,  0.788220000000,  0.052609000000
               523,  0.149516000000,  0.800460000000,  0.049122000000
               524,  0.161041000000,  0.812140000000,  0.045954000000
               525,  0.172953000000,  0.823330000000,  0.043050000000
               526,  0.185209000000,  0.834120000000,  0.040368000000
               527,  0.197755000000,  0.844600000000,  0.037839000000
               528,  0.210538000000,  0.854870000000,  0.035384000000
               529,  0.223460000000,  0.865040000000,  0.032949000000
               530,  0.236491000000,  0.875211000000,  0.030451000000
               531,  0.249633000000,  0.885370000000,  0.028029000000
               532,  0.262972000000,  0.895370000000,  0.025862000000
               533,  0.276515000000,  0.905150000000,  0.023920000000
               534,  0.290269000000,  0.914650000000,  0.022174000000
               535,  0.304213000000,  0.923810000000,  0.020584000000
               536,  0.318361000000,  0.932550000000,  0.019127000000
               537,  0.332705000000,  0.940810000000,  0.017740000000
               538,  0.347232000000,  0.948520000000,  0.016403000000
               539,  0.361926000000,  0.955600000000,  0.015064000000
               540,  0.376772000000,  0.961988000000,  0.013676000000
               541,  0.391683000000,  0.967540000000,  0.012308000000
               542,  0.406594000000,  0.972230000000,  0.011056000000
               543,  0.421539000000,  0.976170000000,  0.009915000000
               544,  0.436517000000,  0.979460000000,  0.008872000000
               545,  0.451584000000,  0.982200000000,  0.007918000000
               546,  0.466782000000,  0.984520000000,  0.007030000000
               547,  0.482147000000,  0.986520000000,  0.006223000000
               548,  0.497738000000,  0.988320000000,  0.005453000000
               549,  0.513606000000,  0.990020000000,  0.004714000000
               550,  0.529826000000,  0.991761000000,  0.003988000000
               551,  0.546440000000,  0.993530000000,  0.003289000000
               552,  0.563426000000,  0.995230000000,  0.002646000000
               553,  0.580726000000,  0.996770000000,  0.002063000000
               554,  0.598290000000,  0.998090000000,  0.001533000000
               555,  0.616053000000,  0.999110000000,  0.001091000000
               556,  0.633948000000,  0.999770000000,  0.000711000000
               557,  0.651901000000,  1.000000000000,  0.000407000000
               558,  0.669824000000,  0.999710000000,  0.000184000000
               559,  0.687632000000,  0.998850000000,  0.000047000000
               560,  0.705224000000,  0.997340000000,  0.000000000000
               561,  0.722773000000,  0.995260000000,  0.000000000000
               562,  0.740483000000,  0.992740000000,  0.000000000000
               563,  0.758273000000,  0.989750000000,  0.000000000000
               564,  0.776083000000,  0.986300000000,  0.000000000000
               565,  0.793832000000,  0.982380000000,  0.000000000000
               566,  0.811436000000,  0.977980000000,  0.000000000000
               567,  0.828822000000,  0.973110000000,  0.000000000000
               568,  0.845879000000,  0.967740000000,  0.000000000000
               569,  0.862525000000,  0.961890000000,  0.000000000000
               570,  0.878655000000,  0.955552000000,  0.000000000000
               571,  0.894208000000,  0.948601000000,  0.000000000000
               572,  0.909206000000,  0.940981000000,  0.000000000000
               573,  0.923672000000,  0.932798000000,  0.000000000000
               574,  0.937638000000,  0.924158000000,  0.000000000000
               575,  0.951162000000,  0.915175000000,  0.000000000000
               576,  0.964283000000,  0.905954000000,  0.000000000000
               577,  0.977068000000,  0.896608000000,  0.000000000000
               578,  0.989590000000,  0.887249000000,  0.000000000000
               579,  1.001910000000,  0.877986000000,  0.000000000000
               580,  1.014160000000,  0.868934000000,  0.000000000000
               581,  1.026500000000,  0.860164000000,  0.000000000000
               582,  1.038800000000,  0.851519000000,  0.000000000000
               583,  1.051000000000,  0.842963000000,  0.000000000000
               584,  1.062900000000,  0.834393000000,  0.000000000000
               585,  1.074300000000,  0.825623000000,  0.000000000000
               586,  1.085200000000,  0.816764000000,  0.000000000000
               587,  1.095200000000,  0.807544000000,  0.000000000000
               588,  1.104200000000,  0.797947000000,  0.000000000000
               589,  1.112000000000,  0.787893000000,  0.000000000000
               590,  1.118520000000,  0.777405000000,  0.000000000000
               591,  1.123800000000,  0.766490000000,  0.000000000000
               592,  1.128000000000,  0.755309000000,  0.000000000000
               593,  1.131100000000,  0.743845000000,  0.000000000000
               594,  1.133200000000,  0.732190000000,  0.000000000000
               595,  1.134300000000,  0.720353000000,  0.000000000000
               596,  1.134300000000,  0.708281000000,  0.000000000000
               597,  1.133300000000,  0.696055000000,  0.000000000000
               598,  1.131200000000,  0.683621000000,  0.000000000000
               599,  1.128100000000,  0.671048000000,  0.000000000000
               600,  1.123990000000,  0.658341000000,  0.000000000000
               601,  1.118900000000,  0.645545000000,  0.000000000000
               602,  1.112900000000,  0.632718000000,  0.000000000000
               603,  1.105900000000,  0.619815000000,  0.000000000000
               604,  1.098000000000,  0.606887000000,  0.000000000000
               605,  1.089100000000,  0.593878000000,  0.000000000000
               606,  1.079200000000,  0.580781000000,  0.000000000000
               607,  1.068400000000,  0.567653000000,  0.000000000000
               608,  1.056700000000,  0.554490000000,  0.000000000000
               609,  1.044000000000,  0.541228000000,  0.000000000000
               610,  1.030480000000,  0.527963000000,  0.000000000000
               611,  1.016000000000,  0.514634000000,  0.000000000000
               612,  1.000800000000,  0.501363000000,  0.000000000000
               613,  0.984790000000,  0.488124000000,  0.000000000000
               614,  0.968080000000,  0.474935000000,  0.000000000000
               615,  0.950740000000,  0.461834000000,  0.000000000000
               616,  0.932800000000,  0.448823000000,  0.000000000000
               617,  0.914340000000,  0.435917000000,  0.000000000000
               618,  0.895390000000,  0.423153000000,  0.000000000000
               619,  0.876030000000,  0.410526000000,  0.000000000000
               620,  0.856297000000,  0.398057000000,  0.000000000000
               621,  0.836350000000,  0.385835000000,  0.000000000000
               622,  0.816290000000,  0.373951000000,  0.000000000000
               623,  0.796050000000,  0.362311000000,  0.000000000000
               624,  0.775610000000,  0.350863000000,  0.000000000000
               625,  0.754930000000,  0.339554000000,  0.000000000000
               626,  0.733990000000,  0.328309000000,  0.000000000000
               627,  0.712780000000,  0.317118000000,  0.000000000000
               628,  0.691290000000,  0.305936000000,  0.000000000000
               629,  0.669520000000,  0.294737000000,  0.000000000000
               630,  0.647467000000,  0.283493000000,  0.000000000000
               631,  0.625110000000,  0.272222000000,  0.000000000000
               632,  0.602520000000,  0.260990000000,  0.000000000000
               633,  0.579890000000,  0.249877000000,  0.000000000000
               634,  0.557370000000,  0.238946000000,  0.000000000000
               635,  0.535110000000,  0.228254000000,  0.000000000000
               636,  0.513240000000,  0.217853000000,  0.000000000000
               637,  0.491860000000,  0.207780000000,  0.000000000000
               638,  0.471080000000,  0.198072000000,  0.000000000000
               639,  0.450960000000,  0.188748000000,  0.000000000000
               640,  0.431567000000,  0.179828000000,  0.000000000000
               641,  0.412870000000,  0.171285000000,  0.000000000000
               642,  0.394750000000,  0.163059000000,  0.000000000000
               643,  0.377210000000,  0.155151000000,  0.000000000000
               644,  0.360190000000,  0.147535000000,  0.000000000000
               645,  0.343690000000,  0.140211000000,  0.000000000000
               646,  0.327690000000,  0.133170000000,  0.000000000000
               647,  0.312170000000,  0.126400000000,  0.000000000000
               648,  0.297110000000,  0.119892000000,  0.000000000000
               649,  0.282500000000,  0.113640000000,  0.000000000000
               650,  0.268329000000,  0.107633000000,  0.000000000000
               651,  0.254590000000,  0.101870000000,  0.000000000000
               652,  0.241300000000,  0.096347000000,  0.000000000000
               653,  0.228480000000,  0.091063000000,  0.000000000000
               654,  0.216140000000,  0.086010000000,  0.000000000000
               655,  0.204300000000,  0.081187000000,  0.000000000000
               656,  0.192950000000,  0.076583000000,  0.000000000000
               657,  0.182110000000,  0.072198000000,  0.000000000000
               658,  0.171770000000,  0.068024000000,  0.000000000000
               659,  0.161920000000,  0.064052000000,  0.000000000000
               660,  0.152568000000,  0.060281000000,  0.000000000000
               661,  0.143670000000,  0.056697000000,  0.000000000000
               662,  0.135200000000,  0.053292000000,  0.000000000000
               663,  0.127130000000,  0.050059000000,  0.000000000000
               664,  0.119480000000,  0.046998000000,  0.000000000000
               665,  0.112210000000,  0.044096000000,  0.000000000000
               666,  0.105310000000,  0.041345000000,  0.000000000000
               667,  0.098786000000,  0.038750700000,  0.000000000000
               668,  0.092610000000,  0.036297800000,  0.000000000000
               669,  0.086773000000,  0.033983200000,  0.000000000000
               670,  0.081260600000,  0.031800400000,  0.000000000000
               671,  0.076048000000,  0.029739500000,  0.000000000000
               672,  0.071114000000,  0.027791800000,  0.000000000000
               673,  0.066454000000,  0.025955100000,  0.000000000000
               674,  0.062062000000,  0.024226300000,  0.000000000000
               675,  0.057930000000,  0.022601700000,  0.000000000000
               676,  0.054050000000,  0.021077900000,  0.000000000000
               677,  0.050412000000,  0.019650500000,  0.000000000000
               678,  0.047006000000,  0.018315300000,  0.000000000000
               679,  0.043823000000,  0.017068600000,  0.000000000000
               680,  0.040850800000,  0.015905100000,  0.000000000000
               681,  0.038072000000,  0.014818300000,  0.000000000000
               682,  0.035468000000,  0.013800800000,  0.000000000000
               683,  0.033031000000,  0.012849500000,  0.000000000000
               684,  0.030753000000,  0.011960700000,  0.000000000000
               685,  0.028623000000,  0.011130300000,  0.000000000000
               686,  0.026635000000,  0.010355500000,  0.000000000000
               687,  0.024781000000,  0.009633200000,  0.000000000000
               688,  0.023052000000,  0.008959900000,  0.000000000000
               689,  0.021441000000,  0.008332400000,  0.000000000000
               690,  0.019941300000,  0.007748800000,  0.000000000000
               691,  0.018544000000,  0.007204600000,  0.000000000000
               692,  0.017241000000,  0.006697500000,  0.000000000000
               693,  0.016027000000,  0.006225100000,  0.000000000000
               694,  0.014896000000,  0.005785000000,  0.000000000000
               695,  0.013842000000,  0.005375100000,  0.000000000000
               696,  0.012862000000,  0.004994100000,  0.000000000000
               697,  0.011949000000,  0.004639200000,  0.000000000000
               698,  0.011100000000,  0.004309300000,  0.000000000000
               699,  0.010311000000,  0.004002800000,  0.000000000000
               700,  0.009576880000,  0.003717740000,  0.000000000000
               701,  0.008894000000,  0.003452620000,  0.000000000000
               702,  0.008258100000,  0.003205830000,  0.000000000000
               703,  0.007666400000,  0.002976230000,  0.000000000000
               704,  0.007116300000,  0.002762810000,  0.000000000000
               705,  0.006605200000,  0.002564560000,  0.000000000000
               706,  0.006130600000,  0.002380480000,  0.000000000000
               707,  0.005690300000,  0.002209710000,  0.000000000000
               708,  0.005281900000,  0.002051320000,  0.000000000000
               709,  0.004903300000,  0.001904490000,  0.000000000000
               710,  0.004552630000,  0.001768470000,  0.000000000000
               711,  0.004227500000,  0.001642360000,  0.000000000000
               712,  0.003925800000,  0.001525350000,  0.000000000000
               713,  0.003645700000,  0.001416720000,  0.000000000000
               714,  0.003385900000,  0.001315950000,  0.000000000000
               715,  0.003144700000,  0.001222390000,  0.000000000000
               716,  0.002920800000,  0.001135550000,  0.000000000000
               717,  0.002713000000,  0.001054940000,  0.000000000000
               718,  0.002520200000,  0.000980140000,  0.000000000000
               719,  0.002341100000,  0.000910660000,  0.000000000000
               720,  0.002174960000,  0.000846190000,  0.000000000000
               721,  0.002020600000,  0.000786290000,  0.000000000000
               722,  0.001877300000,  0.000730680000,  0.000000000000
               723,  0.001744100000,  0.000678990000,  0.000000000000
               724,  0.001620500000,  0.000631010000,  0.000000000000
               725,  0.001505700000,  0.000586440000,  0.000000000000
               726,  0.001399200000,  0.000545110000,  0.000000000000
               727,  0.001300400000,  0.000506720000,  0.000000000000
               728,  0.001208700000,  0.000471110000,  0.000000000000
               729,  0.001123600000,  0.000438050000,  0.000000000000
               730,  0.001044760000,  0.000407410000,  0.000000000000
               731,  0.000971560000,  0.000378962000,  0.000000000000
               732,  0.000903600000,  0.000352543000,  0.000000000000
               733,  0.000840480000,  0.000328001000,  0.000000000000
               734,  0.000781870000,  0.000305208000,  0.000000000000
               735,  0.000727450000,  0.000284041000,  0.000000000000
               736,  0.000676900000,  0.000264375000,  0.000000000000
               737,  0.000629960000,  0.000246109000,  0.000000000000
               738,  0.000586370000,  0.000229143000,  0.000000000000
               739,  0.000545870000,  0.000213376000,  0.000000000000
               740,  0.000508258000,  0.000198730000,  0.000000000000
               741,  0.000473300000,  0.000185115000,  0.000000000000
               742,  0.000440800000,  0.000172454000,  0.000000000000
               743,  0.000410580000,  0.000160678000,  0.000000000000
               744,  0.000382490000,  0.000149730000,  0.000000000000
               745,  0.000356380000,  0.000139550000,  0.000000000000
               746,  0.000332110000,  0.000130086000,  0.000000000000
               747,  0.000309550000,  0.000121290000,  0.000000000000
               748,  0.000288580000,  0.000113106000,  0.000000000000
               749,  0.000269090000,  0.000105501000,  0.000000000000
               750,  0.000250969000,  0.000098428000,  0.000000000000
               751,  0.000234130000,  0.000091853000,  0.000000000000
               752,  0.000218470000,  0.000085738000,  0.000000000000
               753,  0.000203910000,  0.000080048000,  0.000000000000
               754,  0.000190350000,  0.000074751000,  0.000000000000
               755,  0.000177730000,  0.000069819000,  0.000000000000
               756,  0.000165970000,  0.000065222000,  0.000000000000
               757,  0.000155020000,  0.000060939000,  0.000000000000
               758,  0.000144800000,  0.000056942000,  0.000000000000
               759,  0.000135280000,  0.000053217000,  0.000000000000
               760,  0.000126390000,  0.000049737000,  0.000000000000
               761,  0.000118100000,  0.000046491000,  0.000000000000
               762,  0.000110370000,  0.000043464000,  0.000000000000
               763,  0.000103150000,  0.000040635000,  0.000000000000
               764,  0.000096427000,  0.000038000000,  0.000000000000
               765,  0.000090151000,  0.000035540500,  0.000000000000
               766,  0.000084294000,  0.000033244800,  0.000000000000
               767,  0.000078830000,  0.000031100600,  0.000000000000
               768,  0.000073729000,  0.000029099000,  0.000000000000
               769,  0.000068969000,  0.000027230700,  0.000000000000
               770,  0.000064525800,  0.000025486000,  0.000000000000
               771,  0.000060376000,  0.000023856100,  0.000000000000
               772,  0.000056500000,  0.000022333200,  0.000000000000
               773,  0.000052880000,  0.000020910400,  0.000000000000
               774,  0.000049498000,  0.000019580800,  0.000000000000
               775,  0.000046339000,  0.000018338400,  0.000000000000
               776,  0.000043389000,  0.000017177700,  0.000000000000
               777,  0.000040634000,  0.000016093400,  0.000000000000
               778,  0.000038060000,  0.000015080000,  0.000000000000
               779,  0.000035657000,  0.000014133600,  0.000000000000
               780,  0.000033411700,  0.000013249000,  0.000000000000
               781,  0.000031315000,  0.000012422600,  0.000000000000
               782,  0.000029355000,  0.000011649900,  0.000000000000
               783,  0.000027524000,  0.000010927700,  0.000000000000
               784,  0.000025811000,  0.000010251900,  0.000000000000
               785,  0.000024209000,  0.000009619600,  0.000000000000
               786,  0.000022711000,  0.000009028100,  0.000000000000
               787,  0.000021308000,  0.000008474000,  0.000000000000
               788,  0.000019994000,  0.000007954800,  0.000000000000
               789,  0.000018764000,  0.000007468600,  0.000000000000
               790,  0.000017611500,  0.000007012800,  0.000000000000
               791,  0.000016532000,  0.000006585800,  0.000000000000
               792,  0.000015521000,  0.000006185700,  0.000000000000
               793,  0.000014574000,  0.000005810700,  0.000000000000
               794,  0.000013686000,  0.000005459000,  0.000000000000
               795,  0.000012855000,  0.000005129800,  0.000000000000
               796,  0.000012075000,  0.000004820600,  0.000000000000
               797,  0.000011345000,  0.000004531200,  0.000000000000
               798,  0.000010659000,  0.000004259100,  0.000000000000
               799,  0.000010017000,  0.000004004200,  0.000000000000
               800,  0.000009413630,  0.000003764730,  0.000000000000
               801,  0.000008847900,  0.000003539950,  0.000000000000
               802,  0.000008317100,  0.000003329140,  0.000000000000
               803,  0.000007819000,  0.000003131150,  0.000000000000
               804,  0.000007351600,  0.000002945290,  0.000000000000
               805,  0.000006913000,  0.000002770810,  0.000000000000
               806,  0.000006501500,  0.000002607050,  0.000000000000
               807,  0.000006115300,  0.000002453290,  0.000000000000
               808,  0.000005752900,  0.000002308940,  0.000000000000
               809,  0.000005412700,  0.000002173380,  0.000000000000
               810,  0.000005093470,  0.000002046130,  0.000000000000
               811,  0.000004793800,  0.000001926620,  0.000000000000
               812,  0.000004512500,  0.000001814400,  0.000000000000
               813,  0.000004248300,  0.000001708950,  0.000000000000
               814,  0.000004000200,  0.000001609880,  0.000000000000
               815,  0.000003767100,  0.000001516770,  0.000000000000
               816,  0.000003548000,  0.000001429210,  0.000000000000
               817,  0.000003342100,  0.000001346860,  0.000000000000
               818,  0.000003148500,  0.000001269450,  0.000000000000
               819,  0.000002966500,  0.000001196620,  0.000000000000
               820,  0.000002795310,  0.000001128090,  0.000000000000
               821,  0.000002634500,  0.000001063680,  0.000000000000
               822,  0.000002483400,  0.000001003130,  0.000000000000
               823,  0.000002341400,  0.000000946220,  0.000000000000
               824,  0.000002207800,  0.000000892630,  0.000000000000
               825,  0.000002082000,  0.000000842160,  0.000000000000
               826,  0.000001963600,  0.000000794640,  0.000000000000
               827,  0.000001851900,  0.000000749780,  0.000000000000
               828,  0.000001746500,  0.000000707440,  0.000000000000
               829,  0.000001647100,  0.000000667480,  0.000000000000
               830,  0.000001553140,  0.000000629700,  0.000000000000];
lambda = cmf(:, 1)'; xFcn = cmf(:, 2)'; yFcn = cmf(:, 3)'; zFcn = cmf(:, 4)';
XYZ = interp1(lambda', [xFcn; yFcn; zFcn]', wavelength, 'pchip', 0);
XYZ2sRGB = makecform('xyz2srgb');
RGB = applycform(XYZ, XYZ2sRGB);
end