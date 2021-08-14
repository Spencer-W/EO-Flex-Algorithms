%plotting for suite 2p data
clear; close all;
%for laptop
% pwd1='C:\Users\Spencer\Desktop\Documents\Spencer-salk-comp\2019-5-31-EO-flex_imaging_exps\ProcessedData\Run_13\suite2p\plane0';
%for tower
pwd1='C:\Users\Spencer\Dropbox\UCSD Research\Spencer-salk-comp\2019-5-31-EO-flex_imaging_exps\ProcessedData\Run_13\suite2p\plane0';

filename='run13';

load('run13_trimmed.mat');
load('Fall.mat');

base=base_old;
Laser=Laser_old;
tE=tE_old;
tL=tL_old;
%% find rising edges
tE_to_tL=(tL(end))/(length(tE)-1);
tL_exp=(0:1:(length(tE)-1))'*tE_to_tL;
r_tE=tE(2)-tE(1);

laser_exp=interp1(tL, Laser, tL_exp);
laser_thr=laser_exp;

lsr_thr=0.95*max(laser_exp);
laser_thr(laser_thr>lsr_thr)=4.9; 
laser_thr(laser_thr<=lsr_thr)=0;
% 
% dy=diff(laser_thr);
% dy=[0; dy];
% rising_ind=find(dy==4.9);
% falling_ind=find(dy==-4.9);


%% shorten data set - First 4 frames of fluorescent trace are trimmed off.
fps=8.14; %2-photon scanning rate
trimmed_frames=4; %first 4 frames
time_to_trim=trimmed_frames/fps;
tE_rem=round(20000*time_to_trim);

tE_shortened=tE((tE_rem+1):end,1)-tE(tE_rem+1);
tL_shortened=tL_exp((tE_rem+1):end,1)-tL_exp(tE_rem+1);
Laser_shortened=laser_thr((tE_rem+1):end,1);
base_shortened=base_trimmed((tE_rem+1):end,1);


%% plot old vs trimmed tE/tL
ax1=subplot(4,1,1);
plot(tE, base)
ax2=subplot(4,1,2);
plot(tE_shortened, base_shortened)
ax3=subplot(4,1,3);
plot(tL, Laser)
ax4=subplot(4,1,4);
plot(tL_shortened, Laser_shortened)

linkaxes([ax1,ax2,ax3,ax4],'x');

%%


F_orig=F-0.7*Fneu;
time=(0:1:size(F,2)-1)*1/fps;
plot(time,F_orig);

F_sort=sortrows(F_orig, 'descend');
figure;
plot(time, F_sort((1:5),:));
legend('1', '2', '3', '4', '5');


F0=min(F_orig,[],2);
% F0=mean(F_orig(:,1:ceil((rising_ind(1)*tE_to_tL)*fps)),2);

DelF=(F_orig-F0)./F0;

DelF_norm=(DelF-min(DelF,[],2))./(max(DelF,[],2)-min(DelF,[],2));

figure;
plot(time,DelF);

% F_sort=sortrows(DelF, 'descend');
% figure;
% plot(time, F_sort((1:5),:));
% legend('1', '2', '3', '4', '5');



% 
% figure;
% plot(time, F_sort((6:10),:));
% legend('6', '7', '8', '9', '10');
% 
% figure;
% plot(time, F_sort((11:15),:));
% legend('11', '12', '13', '14', '15');
%%
%plot roi with F0 placed at mean vs median


for i=1:10
    roi=i;
trace=F_orig(roi,:);
figure;
hold on;
plot(time,trace);
plot(time,DelF(roi,:))
plot(time,DelF_norm(roi,:))
threshline = refline([0 mean(trace)]);
threshline.Color = 'r';

threshline = refline([0 median(trace)]);
threshline.Color = 'b';

Outlier= (abs(trace-min(trace))>2*std(trace));
ftrim=trace(~Outlier);

threshline = refline([0 F0(roi,:)]);
threshline.Color = 'g';
hold off;
end
%% remove rois below median of the average deltaF/F
close all;

ROIs=(0:1:size(DelF,1)-1)'; %ROI numbering starts at 0
F_rois=[ROIs DelF];
Frois_norm=DelF_norm;

F_mean=mean(F_orig,2);
% F_max=max(DelF,[],2);
sigma=1;

F_rois(~(F_mean>sigma*median(F_mean)),:)=[];
% F_rois(~(F_mean>1.5*median(F_mean)),:)=[];
Frois_norm(~(F_mean>sigma*median(F_mean)),:)=[];

 
ROI_ind=F_rois(:,1);
F_rois(:,1)=[];
plot(time,F_rois);

DefF_mean=mean(DelF,2);

ytick=DefF_mean(ROI_ind+1); %zeros(1,length(ROI_ind));
incr=4*round(mean(ytick));

F_test=F_rois;
for i=1:length(ROI_ind)
    F_test(i,:)=F_rois(i,:)+incr*i;
    
    ytick(i,:)=ytick(i)+incr*i;
end

plot(time,F_test,'k');
xlabel('Time (s)');
ylabel('ROI Number');
set(gca, 'ytick', ytick, 'yticklabel', ROI_ind);
%set(gca, 'ytick', ROI_ind, 'yticklabel', {'Off', 'On'});

set(gcf,'color','w');
set(gca,'box','off')


figure;
% ytick=((length(ROIs)-1):-1:0);
ytick=(0:(length(ROI_ind)-1));

imagesc(time(1):time(end),ytick,Frois_norm)
colormap(parula);
colorbar;
xlabel('Time (s)');
ylabel('ROI Number');

width=8;
height=width/2;
set(gcf, 'Position', [0 0 width*100, height*100]); %<- Set size
set(gca, 'FontSize', 16);
set(gca, 'ytick', ytick, 'yticklabel', ROI_ind);
set(gcf,'color','w');
set(gca,'box','off')
%% 
F_auto_sum=sum(Frois_norm);

ax1=subplot(4,1,1);
imagesc(time(1):time(end),ytick,Frois_norm)
colormap(parula);
set(gca, 'ytick', ytick, 'yticklabel', ROI_ind);
% c=colorbar;
% c.FontSize=16;
% c.Label.String='\DeltaF/F';
set(gcf,'color','w');
set(gca,'box','off')
h = gca; h.XAxis.Visible = 'off'; 

ax2=subplot(4,1,2);
plot(time,F_auto_sum);
h = gca; h.XAxis.Visible = 'off'; 
h = gca; h.YAxis.Visible = 'off'; 


ax3=subplot(4,1,3);
plot(tE_trimmed,base_trimmed*1000);
set(gcf,'color','w');
set(gca,'box','off')
ylabel('Potential (\muV)');
h = gca; h.XAxis.Visible = 'off'; 


ax4=subplot(4,1,4);
plot(tL_old,Laser_old);
set(gcf,'color','w');
set(gca,'box','off')
ylabel('Stimulus');
xlabel('Time (s)');
set(gca, 'ytick', [0.05,4.9], 'yticklabel', {'Off', 'On'});
linkaxes([ax1,ax2,ax3,ax4],'x');

print('-painters','-bestfit',sprintf('%s_%isigma_ROI_4Panel.pdf',filename,sigma),'-dpdf');
print('-painters',sprintf('%s_%isigma_ROI_4panel.svg',filename,sigma),'-dsvg');

%% sorted by flourescence
F_max=max(max(F_rois));
[max_x, max_y]=find(F_orig==F_max);

[F_sort, F_ind]=sortrows(F_rois, max_y, 'descend');
ytick=F_mean(ROI_ind(F_ind)+1); %zeros(1,length(ROI_ind));
F_test=F_sort;
for i=1:length(F_ind)
    F_test(i,:)=F_sort(i,:)+100*i;
    
    ytick(i,:)=ytick(i)+100*i;
end

plot(time,F_test);
xlabel('Time (s)');
ylabel('ROI Number');
set(gca, 'ytick', ytick, 'yticklabel', ROI_ind(F_ind));
set(gcf,'color','w');
set(gca,'box','off')

%% custom ROI list
ROIs=[30 4 9 0 2 1 11 5];
ROIs=sort(ROIs);

F_selected=DelF(ROIs+1,:);
F_test=F_selected;



ytick=mean(F_selected,2); %zeros(1,length(ROI_ind));
incr=4*round(mean(ytick));

for i=1:length(ROIs)
    F_test(i,:)=F_selected(i,:)+incr*i;
    ytick(i,:)=ytick(i)+incr*i;
end

plot(time,F_test,'k');
xlabel('Time (s)');
ylabel('ROI Number');
width=8;
height=width/2;
pos = get(gcf, 'Position');
set(gcf, 'Position', [pos(1) pos(2) width*100, height*100]); %<- Set size
set(gca, 'FontSize', 16);
set(gca, 'ytick', ytick, 'yticklabel', ROIs);
set(gcf,'color','w');
set(gca,'box','off')
xlim([time(1) time(end)]);
print('-painters','-bestfit',sprintf('%s_selected_ROI_traces.pdf',filename),'-dpdf');
print('-painters',sprintf('%s_selected_ROI_traces.svg',filename),'-dsvg');

F_selected=DelF_norm(ROIs+1,:);

figure;
% ytick=((length(ROIs)-1):-1:0);
ytick=(0:(length(ROIs)-1));

imagesc(time(1):time(end),ytick,F_selected)
colormap(parula);
c=colorbar;
c.FontSize=16;
c.Label.String='\DeltaF/F';

xlabel('Time (s)');
ylabel('ROI Number');

width=8;
height=width/2;
set(gcf, 'Position', [pos(1) pos(2) width*100, height*100]); %<- Set size
set(gca, 'FontSize', 16);
set(gca, 'ytick', ytick, 'yticklabel', ROIs);
set(gcf,'color','w');
set(gca,'box','off')

print('-painters','-bestfit',sprintf('%s_selected_ROI_heatmap.pdf',filename),'-dpdf');
print('-painters',sprintf('%s_selected_ROI_heatmap.svg',filename),'-dsvg');
%% plot heatmap with electrical recording, 4-panel
Fontsize=28;
FontsizeAxis=18;

F_selected_sum=sum(F_selected);

ax1=subplot(4,1,1);
imagesc(time(1):time(end),ytick,F_selected)
colormap(parula);


ax=ancestor(gca, 'axes');
yrule=ax.YAxis;
yrule.FontSize =FontsizeAxis;
yL1=ylabel('ROI Number', 'Interpreter', 'tex');
yL1.FontSize=Fontsize;
yL1.Position(1)= -10

set(gca, 'ytick', ytick, 'yticklabel', ROIs);
% c=colorbar;
% c.FontSize=16;
% c.Label.String='\DeltaF/F';
set(gcf,'color','w');
set(gca,'box','off')
h = gca; h.XAxis.Visible = 'off'; 
width=8;
height=10;
set(gcf, 'Position', [10 10 width*100, height*100]); %<- Set size
ax2=subplot(4,1,2);

plot(time,F_selected_sum);
h = gca; h.XAxis.Visible = 'off'; 
h = gca; h.YAxis.Visible = 'off'; 
set(gca,'YTick',[]);
set(gca,'box','off')

ax=ancestor(gca, 'axes');
yrule=ax.YAxis;
yrule.FontSize =FontsizeAxis;
yL2=ylabel(['(',char(916),'F/F)'], 'Interpreter', 'tex');
yL2.FontSize=Fontsize;
yL2.Position(1)= -10

ax3=subplot(4,1,3);
plot(tE_shortened,base_shortened*1000);
set(gcf,'color','w');
set(gca,'box','off')




ax=ancestor(gca, 'axes');
yrule=ax.YAxis;
yrule.FontSize =FontsizeAxis;
yL3=ylabel(['Potential (',char(181),'V)'], 'Interpreter', 'tex');
yL3.FontSize=Fontsize;
yL3.Position(1)= -10

h = gca; h.XAxis.Visible = 'off'; 
xlim([tE_old(1) tE_old(end)]);


ax4=subplot(4,1,4);
plot(tL_shortened,Laser_shortened);
set(gcf,'color','w');
set(gca,'box','off')




xL=xlabel('Time (s)');
% xticks([PCA_time(1) round(PCA_time(end))]);



ax=ancestor(gca, 'axes');
xrule=ax.XAxis;
yrule=ax.YAxis;
xrule.FontSize =FontsizeAxis;
yrule.FontSize =FontsizeAxis;
yL4=ylabel(['Stimulus'], 'Interpreter', 'tex');
yL4.FontSize=Fontsize;
yL4.Position
yL4.Position(1)= -10
yL4.Position(2)= 2.5

set(gca, 'ytick', [0.05,4.9], 'yticklabel', {'Off', 'On'});
xlim([tL_old(1) tL_old(end)]);
xL.FontSize=Fontsize;


linkaxes([ax1,ax2,ax3,ax4],'x');


print('-painters','-bestfit',sprintf('%s_selected_ROI_4Panel.pdf',filename),'-dpdf');
print('-painters',sprintf('%s_selected_ROI_4panel.svg',filename),'-dsvg');
%% plot heatmap with electrical recording 3 panel, no summed fluorescence
Fontsize=28;
FontsizeAxis=18;
Ylabelpos=-9;

F_selected_sum=sum(F_selected);

ax1=subplot(3,1,1);
imagesc(time(1):time(end),ytick,F_selected)
colormap(parula);


ax=ancestor(gca, 'axes');
yrule=ax.YAxis;
yrule.FontSize =FontsizeAxis;
yL1=ylabel('ROI Number', 'Interpreter', 'tex');
yL1.FontSize=Fontsize;
yL1.Position(1)= Ylabelpos

set(gca, 'ytick', ytick, 'yticklabel', ROIs);
% c=colorbar;
% c.FontSize=16;
% c.Label.String='\DeltaF/F';
set(gcf,'color','w');
set(gca,'box','off')
h = gca; h.XAxis.Visible = 'off'; 
width=8;
height=10;
set(gcf, 'Position', [10 10 width*100, height*100]); %<- Set size

ax2=subplot(3,1,2);
plot(tE_shortened,base_shortened*1000);
set(gcf,'color','w');
set(gca,'box','off')




ax=ancestor(gca, 'axes');
yrule=ax.YAxis;
yrule.FontSize =FontsizeAxis;
yL2=ylabel(['Potential (',char(181),'V)'], 'Interpreter', 'tex');
yL2.FontSize=Fontsize;
yL2.Position(1)= Ylabelpos

h = gca; h.XAxis.Visible = 'off'; 
xlim([tE_old(1) tE_old(end)]);


ax3=subplot(3,1,3);
plot(tL_shortened,Laser_shortened);
set(gcf,'color','w');
set(gca,'box','off')




xL=xlabel('Time (s)');
% xticks([PCA_time(1) round(PCA_time(end))]);



ax=ancestor(gca, 'axes');
xrule=ax.XAxis;
yrule=ax.YAxis;
xrule.FontSize =FontsizeAxis;
yrule.FontSize =FontsizeAxis;
yL3=ylabel(['Stimulus'], 'Interpreter', 'tex');
yL3.FontSize=Fontsize;
yL3.Position
yL3.Position(1)= Ylabelpos
yL3.Position(2)= 2.5

set(gca, 'ytick', [0.05,4.9], 'yticklabel', {'Off', 'On'});
xlim([tL_old(1) tL_old(end)]);
xL.FontSize=Fontsize;


linkaxes([ax1,ax2,ax3],'x');


print('-painters','-bestfit',sprintf('%s_selected_ROI_3Panel.pdf',filename),'-dpdf');
print('-painters',sprintf('%s_selected_ROI_3panel.svg',filename),'-dsvg');

%% plot heatmap with electrical recording 3 panel, no summed fluorescence
Fontsize=28;
FontsizeAxis=18;
Ylabelpos=-9;

F_selected_sum=sum(F_selected);

ax1=subplot(3,1,1);
imagesc(time(1):time(end),ytick,F_selected)
colormap(parula);


ax=ancestor(gca, 'axes');
yrule=ax.YAxis;
yrule.FontSize =FontsizeAxis;
yL1=ylabel('ROI Number', 'Interpreter', 'tex');
yL1.FontSize=Fontsize;
yL1.Position(1)= Ylabelpos

set(gca, 'ytick', ytick, 'yticklabel', ROIs);
% c=colorbar;
% c.FontSize=16;
% c.Label.String='\DeltaF/F';
set(gcf,'color','w');
set(gca,'box','off')
h = gca; h.XAxis.Visible = 'off'; 
width=8;
height=10;
set(gcf, 'Position', [10 10 width*100, height*100]); %<- Set size

ax2=subplot(3,1,2);
plot(tL_shortened,Laser_shortened);
set(gcf,'color','w');
set(gca,'box','off')




xL=xlabel('Time (s)');
% xticks([PCA_time(1) round(PCA_time(end))]);



ax=ancestor(gca, 'axes');
xrule=ax.XAxis;
yrule=ax.YAxis;
xrule.FontSize =FontsizeAxis;
yrule.FontSize =FontsizeAxis;
yL3=ylabel(['Stimulus'], 'Interpreter', 'tex');
yL3.FontSize=Fontsize;
yL3.Position
yL3.Position(1)= Ylabelpos
yL3.Position(2)= 2.5

set(gca, 'ytick', [0.05,4.9], 'yticklabel', {'Off', 'On'});
xlim([tL_old(1) tL_old(end)]);
xL.FontSize=Fontsize;

ax3=subplot(3,1,3);
plot(tE_shortened,base_shortened*1000);
set(gcf,'color','w');
set(gca,'box','off')




ax=ancestor(gca, 'axes');
yrule=ax.YAxis;
yrule.FontSize =FontsizeAxis;
yL2=ylabel(['Potential (',char(181),'V)'], 'Interpreter', 'tex');
yL2.FontSize=Fontsize;
yL2.Position(1)= Ylabelpos

h = gca; h.XAxis.Visible = 'off'; 
xlim([tE_old(1) tE_old(end)]);

linkaxes([ax1,ax2,ax3],'x');


print('-painters','-bestfit',sprintf('%s_selected_ROI_3Panel_diff_order.pdf',filename),'-dpdf');
print('-painters',sprintf('%s_selected_ROI_3panel_diff_order.svg',filename),'-dsvg');
%% individual plots for above 
imagesc(time(1):time(end),ytick,F_selected)
colormap(parula);
set(gca, 'ytick', ytick, 'yticklabel', ROIs);
% c=colorbar;
% c.FontSize=16;
% c.Label.String='\DeltaF/F';
set(gcf,'color','w');
set(gca,'box','off')
xlabel('Time (s)');
width=8;
height=width/2;
set(gcf, 'Position', [pos(1) pos(2) width*100, height*100]); %<- Set size
set(gca, 'FontSize', 16);
% h = gca; h.XAxis.Visible = 'off'; 
print('-painters','-bestfit',sprintf('%s_selected_ROI_heatmap_nobar.pdf',filename),'-dpdf');
print('-painters',sprintf('%s_selected_ROI_heatmap_nobar.svg',filename),'-dsvg');

figure;
plot(tE_shortened,base_shortened*1000);
xlim([tE_old(1) tE_old(end)]);
set(gcf,'color','w');
set(gca,'box','off')
ylabel('Potential (\muV)');
xlabel('Time (s)');
set(gcf,'color','w');
set(gca,'box','off')
width=8;
height=width/2;
set(gcf, 'Position', [pos(1) pos(2) width*100, height*100]); %<- Set size
set(gca, 'FontSize', 16);
% h = gca; h.XAxis.Visible = 'off'; 
print('-painters','-bestfit',sprintf('%s_selected_ROI_Eflex_trimmed.pdf',filename),'-dpdf');
print('-painters',sprintf('%s_selected_ROI_Eflex_trimmed.svg',filename),'-dsvg');
%% old sorting of rois using mean flourescence of F_orig
ROIs=(0:1:size(F_orig,1)-1)'; %ROI numbering starts at 0
F_rois=[ROIs F_orig];

F_mean=mean(F_orig,2);

% F_rois(~(F_mean>1.5*median(F_mean)),:)=[];
F_rois(~(F_mean>1.5*median(F_mean)),:)=[];
 
ROI_ind=F_rois(:,1);
F_rois(:,1)=[];
plot(time,F_rois);

ytick=F_mean(ROI_ind+1); %zeros(1,length(ROI_ind));
F_test=F_rois;
for i=1:length(ROI_ind)
    F_test(i,:)=F_rois(i,:)+100*i;
    
    ytick(i,:)=ytick(i)+100*i;
end

plot(time,F_test);
xlabel('Time (s)');
ylabel('ROI Number');
set(gca, 'ytick', ytick, 'yticklabel', ROI_ind);
%set(gca, 'ytick', ROI_ind, 'yticklabel', {'Off', 'On'});

set(gcf,'color','w');
set(gca,'box','off')

%%  3-Panel with selected ROIs and electrical recording 
Fontsize=28;
FontsizeAxis=18;
Ylabelpos=-9;

F_selected_sum=sum(F_selected);

ylabel_rois=(1:length(ROIs));
incr=1;

for i=1:length(ROIs)
    F_test(i,:)=F_selected(i,:)+incr*i;
%     ylabel_rois(:,1)=ylabel_rois(1,i)+incr*i;
end

ax1=subplot(3,1,1);
% imagesc(time(1):time(end),ytick,F_selected)
plot(time,F_test,'k');
% colormap(parula);


ax=ancestor(gca, 'axes');
yrule=ax.YAxis;
yrule.FontSize =FontsizeAxis;
yL1=ylabel('ROI Number', 'Interpreter', 'tex');
yL1.FontSize=Fontsize;
yL1.Position(1)= Ylabelpos

set(gca, 'ytick', ylabel_rois, 'yticklabel', ROIs);
% c=colorbar;
% c.FontSize=16;
% c.Label.String='\DeltaF/F';
set(gcf,'color','w');
set(gca,'box','off')
h = gca; h.XAxis.Visible = 'off'; 
width=8;
height=10;
set(gcf, 'Position', [10 10 width*100, height*100]); %<- Set size

ax2=subplot(3,1,2);
plot(tE_shortened,base_shortened*1000);
set(gcf,'color','w');
set(gca,'box','off')




ax=ancestor(gca, 'axes');
yrule=ax.YAxis;
yrule.FontSize =FontsizeAxis;
yL2=ylabel(['Potential (',char(181),'V)'], 'Interpreter', 'tex');
yL2.FontSize=Fontsize;
yL2.Position(1)= Ylabelpos

h = gca; h.XAxis.Visible = 'off'; 
xlim([tE_old(1) tE_old(end)]);


ax3=subplot(3,1,3);
plot(tL_shortened,Laser_shortened);
set(gcf,'color','w');
set(gca,'box','off')




xL=xlabel('Time (s)');
% xticks([PCA_time(1) round(PCA_time(end))]);



ax=ancestor(gca, 'axes');
xrule=ax.XAxis;
yrule=ax.YAxis;
xrule.FontSize =FontsizeAxis;
yrule.FontSize =FontsizeAxis;
yL3=ylabel(['Stimulus'], 'Interpreter', 'tex');
yL3.FontSize=Fontsize;
yL3.Position
yL3.Position(1)= Ylabelpos
yL3.Position(2)= 2.5

set(gca, 'ytick', [0.05,4.9], 'yticklabel', {'Off', 'On'});
xlim([tL_old(1) tL_old(end)]);
xL.FontSize=Fontsize;


linkaxes([ax1,ax2,ax3],'x');


print('-painters','-bestfit',sprintf('%s_selected_ROI_3Panel_traces.pdf',filename),'-dpdf');
print('-painters',sprintf('%s_selected_ROI_3panel_traces.svg',filename),'-dsvg');
%% 3-panel plot stimmulus first and then selected ROIs and electrical activity
Fontsize=28;
FontsizeAxis=18;
Ylabelpos=-9;

ax1=subplot(3,1,1);
plot(tL_shortened,Laser_shortened);
set(gcf,'color','w');
set(gca,'box','off')


% xticks([PCA_time(1) round(PCA_time(end))]);

ax=ancestor(gca, 'axes');

yrule=ax.YAxis;

yrule.FontSize =FontsizeAxis;
yL3=ylabel(['Stimulus'], 'Interpreter', 'tex');
yL3.FontSize=Fontsize;
yL3.Position
yL3.Position(1)= Ylabelpos
yL3.Position(2)= 2.5

set(gca, 'ytick', [0.05,4.9], 'yticklabel', {'Off', 'On'});
xlim([tL_shortened(1) tL_shortened(end)]);
h = gca; h.XAxis.Visible = 'off'; 

F_selected_sum=sum(F_selected);

ylabel_rois=(1:length(ROIs));
incr=1;

for i=1:length(ROIs)
    F_test(i,:)=F_selected(i,:)+incr*i;
%     ylabel_rois(:,1)=ylabel_rois(1,i)+incr*i;
end

ax2=subplot(3,1,2);
% imagesc(time(1):time(end),ytick,F_selected)
plot(time,F_test,'k');
% colormap(parula);


ax=ancestor(gca, 'axes');
yrule=ax.YAxis;
yrule.FontSize =FontsizeAxis;
yL1=ylabel('ROI Number', 'Interpreter', 'tex');
yL1.FontSize=Fontsize;
yL1.Position(1)= Ylabelpos

set(gca, 'ytick', ylabel_rois, 'yticklabel', ROIs);
% c=colorbar;
% c.FontSize=16;
% c.Label.String='\DeltaF/F';
set(gcf,'color','w');
set(gca,'box','off')
h = gca; h.XAxis.Visible = 'off'; 
width=8;
height=10;
set(gcf, 'Position', [10 10 width*100, height*100]); %<- Set size

ax3=subplot(3,1,3);
plot(tE_shortened,base_shortened*1000);
set(gcf,'color','w');
set(gca,'box','off')

xL=xlabel('Time (s)');


ax=ancestor(gca, 'axes');
xrule=ax.XAxis;
yrule=ax.YAxis;
xrule.FontSize =FontsizeAxis;
yrule.FontSize =FontsizeAxis;
yL2=ylabel(['Potential (',char(181),'V)'], 'Interpreter', 'tex');
yL2.FontSize=Fontsize;
yL2.Position(1)= Ylabelpos

% h = gca; h.XAxis.Visible = 'off'; 
xlim([tE_shortened(1) tE_shortened(end)]);
xL.FontSize=Fontsize;


linkaxes([ax1,ax2,ax3],'x');


print('-painters','-bestfit',sprintf('%s_selected_ROI_3Panel_traces_stimulus_first.pdf',filename),'-dpdf');
print('-painters',sprintf('%s_selected_ROI_3panel_traces_stimulus_first.svg',filename),'-dsvg');
