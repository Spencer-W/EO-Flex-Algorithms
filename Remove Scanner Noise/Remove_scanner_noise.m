%% Created by Spencer Ward
% 11/14/2020
%Timings are for data from 5/31/19 run for eo-flex
clear; close all;
filename='run13';
load(strcat(filename, '.mat'));

r_tE=mean(diff(tE));
tE_to_tL=(tL(end))/(length(tE)-1);
tL_exp=(0:1:(length(tE)-1))'*tE_to_tL;
laser_exp=interp1(tL, Laser, tL_exp);
laser_thr=laser_exp;
lsr_thr=0.95*max(laser_exp);
laser_thr(laser_thr>lsr_thr)=4.9;
laser_thr(laser_thr<=lsr_thr)=0;
dy=diff(laser_thr);
dy=[0; dy];
rising_ind=find(dy==4.9);
falling_ind=find(dy==-4.9);
pulse_width=(falling_ind(1)-rising_ind(1));
peak_ind=rising_ind;
y=base;
y_factored=factor(length(y))
%% Extract Signal based on finding minimums in electrical data
% num_iter=round(length(y)/613);
num_iter=round(length(y)/y_factored(end));
window=length(y)/num_iter;
window_ind=1:window;
ksr_win=1:round(window);
% ksr_r=tE(window)/(length(ksr_win)-1);
fitsize=num_iter*window_ind(end);
yfit=zeros(num_iter*length(ksr_win),1);
xfit=zeros(num_iter*length(ksr_win),1);

tic
for i=1:num_iter
    r=ksr(tE(window_ind+(i-1)*length(window_ind)),abs(y(window_ind+(i-1)*length(window_ind))),0.005, length(ksr_win));
    yfit(ksr_win+(i-1)*length(ksr_win))=r.f;
    xfit(ksr_win+(i-1)*length(ksr_win))=r.x;
    toc, i
end
y2=y;
y2(yfit>0.021)=0;
%% Plot figures of extracted signal
%% trim data
base_trimmed=base;
tE_trimmed=tE;
Laser_trimmed=laser_exp;
tL_trimmed=tL_exp;

tE_old=tE;
Laser_old=laser_exp;
base_old=base;
tL_old=tL_exp;

r_tL=mean(diff(tL_old));
r_tE=mean(diff(tE_old));

trim_array=yfit>0.021;
% trim_array(expanded_widths)=0;


% base_trimmed(trim_array)=[];
% tE_trimmed(trim_array)=[];
% Laser_trimmed(trim_array)=[];
% tL_trimmed(trim_array)=[];

base_trimmed(trim_array)=NaN;
tE_trimmed(trim_array)=NaN;
Laser_trimmed(trim_array)=NaN;
tL_trimmed(trim_array)=NaN;

tE=tE_trimmed;
tL=tL_trimmed;
base=base_trimmed;
Laser=Laser_trimmed;


tE(trim_array)=[];
tL(trim_array)=[];
tE=(0:(length(tE)-1))'*r_tE;
tL=(0:(length(tL)-1))'*r_tL;
base(trim_array)=[];
Laser(trim_array)=[];

figure;
hold on;
ax1=subplot(2,1,1);
plot(tE_trimmed, base_trimmed);
set(gcf,'color','w');
title('Minimums extracted')
ylabel('Activity (mV)');
xlabel('Time (s)');
ax2=subplot(2,1,2);
plot(tL_trimmed,Laser_trimmed);
title('Stimulus');
ylabel('Stimulus');
xlabel('Time (s)');
set(gca, 'ytick', [0.05,4.9], 'yticklabel', {'Off', 'On'});
linkaxes([ax1,ax2], 'x');

save(strcat(filename,'_trimmed.mat'),'tE','base','tL','Laser', 'base_trimmed','tE_trimmed','Laser_trimmed','tL_trimmed',...
    'tE_old','Laser_old','base_old','tL_old');

%% Generate other plots to show difference between original electrical 
%    signal recorded from EO-Flex vs filtered signal after removing
%    Becquerel effect due to two-photon laser passing over probe tip.
%    This plot was used to create the panel for SI fig 8

figure;
ax1=subplot(3,1,1);
title('original signal');
plot(tE,base);
ax2=subplot(3,1,2);
title('smoothing algorithm applied to base signal');
plot(xfit,yfit);
ax3=subplot(3,1,3);
plot(tE,y2);
linkaxes([ax1,ax2,ax3], 'x');

%% Plot Figures
%use 4-6s as zoomed in region
y2=y;
tE2=tE;
tE2(yfit>0.021)=NaN;
y2(yfit>0.021)=NaN;
figure;
hold on;
ax1=subplot(3,1,1);
plot(tE, y*1000);
set(gcf,'color','w');
% title('Original Signal')
ylabel('Activity (\mum)');
% xlabel('Time (s)');
set(gca,'box','off');
 h = gca; h.XAxis.Visible = 'off';

ax2=subplot(3,1,2);
plot(tE2, y2*1000);
set(gcf,'color','w');
% title('Minimums extracted')
ylabel('Activity (\mum)');
% xlabel('Time (s)');
set(gca,'box','off');
 h = gca; h.XAxis.Visible = 'off';
 
ax3=subplot(3,1,3);
plot(tL,Laser);
% title('Stimulus');
ylabel('Stimulus');
xlabel('Time (s)');
set(gca,'box','off');

set(gca, 'ytick', [0.05,4.9], 'yticklabel', {'Off', 'On'});
linkaxes([ax1,ax2,ax3], 'x');
print('-painters',strcat(filename,'_trimmed_3panel.pdf'),'-dpdf');
%% Plot Figure for 4-6s
%use 4-6s as zoomed in region

region=(round(4/r_tE):1:round(6/r_tE))';

y2=y;
tE2=tE;
tE2(yfit>0.021)=NaN;
y2(yfit>0.021)=NaN;


figure;
hold on;
ax1=subplot(3,1,1);
plot(tE(region), y(region)*1000);
set(gcf,'color','w');
% title('Original Signal')
ylabel('Activity (\mum)');
% xlabel('Time (s)');
set(gca,'box','off');
 h = gca; h.XAxis.Visible = 'off';

ax2=subplot(3,1,2);
plot(tE2(region), y2(region)*1000);
set(gcf,'color','w');
% title('Minimums extracted')
ylabel('Activity (\mum)');
% xlabel('Time (s)');
set(gca,'box','off');
 h = gca; h.XAxis.Visible = 'off';
 
ax3=subplot(3,1,3);
plot(tL_exp(region),laser_thr(region));
% title('Stimulus');
ylabel('Stimulus');
xlabel('Time (s)');
set(gca,'box','off');

set(gca, 'ytick', [0.05,4.9], 'yticklabel', {'Off', 'On'});
linkaxes([ax1,ax2,ax3], 'x');
xlim([4 6])
print('-painters',strcat(filename,'_trimmed_3panel_zoomed.pdf'),'-dpdf');
%% Plot 2 panel Figures
%use 4-6s as zoomed in region

y2=y;
tE2=tE;
tE2(yfit>0.021)=NaN;
y2(yfit>0.021)=NaN;


figure;
hold on;
ax1=subplot(2,1,1);
plot(tE, y*1000);
set(gcf,'color','w');
% title('Original Signal')
ylabel('Activity (\mum)');
% xlabel('Time (s)');
set(gca,'box','off');
 h = gca; h.XAxis.Visible = 'off';

ax2=subplot(2,1,2);
plot(tE2, y2*1000);
set(gcf,'color','w');
% title('Minimums extracted')
ylabel('Activity (\mum)');
% xlabel('Time (s)');
set(gca,'box','off');
xlabel('Time (s)');
linkaxes([ax1,ax2,ax3], 'x');

print('-painters',strcat(filename,'_trimmed_2panel.pdf'),'-dpdf');
%% Plot 2 panel Figure for 4-6s
%use 4-6s as zoomed in region

region=(round(4/r_tE):1:round(6/r_tE))';

y2=y;
tE2=tE;
tE2(yfit>0.021)=NaN;
y2(yfit>0.021)=NaN;


figure;
hold on;
ax1=subplot(2,1,1);
plot(tE(region), y(region)*1000);
set(gcf,'color','w');
% title('Original Signal')
ylabel('Activity (\mum)');
% xlabel('Time (s)');
set(gca,'box','off');
 h = gca; h.XAxis.Visible = 'off';

ax2=subplot(2,1,2);
plot(tE2(region), y2(region)*1000);
set(gcf,'color','w');
% title('Minimums extracted')
ylabel('Activity (\mum)');
% xlabel('Time (s)');
set(gca,'box','off');
xlabel('Time (s)');
linkaxes([ax1,ax2,ax3], 'x');
xlim([4 6])

print('-painters',strcat(filename,'_trimmed_2panel_zoomed.pdf'),'-dpdf');


%%
base=base_old;
tE=tE_old;
Laser=Laser_old;
tL=tL_old;