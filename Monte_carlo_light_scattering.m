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
