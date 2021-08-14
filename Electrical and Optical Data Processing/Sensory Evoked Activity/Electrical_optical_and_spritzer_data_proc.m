%% This Matlab scripts take raw electrical (base) and optical (Laser) data 
%with their respective time points tE and tL to compute all the necessary
%figures for Ward et al, "Electro-optical mechanically flexible coaxial microprobes-
%for minimally invasive interfacing with intrinsic neural circuits"

%Created by Spencer Ward
%12/20/2019

%This script will require the data in the form of a 'runX.mat' file.
%This includes the data for base, tE, tL, and Laser; X is a number between 
%36 and 52. Run 36 represents the baseline parameters from which all other
%variable variations (Duty, Freq, Power, and Probe Retraction) deviation
%from. Run 36 was an optical pulse train of 20 pules over 1 second 
%(20 Hz stimulation), with a 10% duty cycle, and this would repeat for 1s 
%of stimulationon (light on) and 1s recovery (light off). Total laser power
%was 208 uW (DPSS laser reading of 4.03.

%this code uses a shade function created by Javier Montalt on
%mathworks.com:

%Javier Montalt Tordera (2021). Filled area plot 
%(https://www.mathworks.com/matlabcentral/fileexchange/69652-filled-area-plot),
%MATLAB Central File Exchange. Retrieved March 21, 2020.

clear; close all;

prompt1='Would you like to run data in...\nDuty cycle varied runs (1)\nFrequency varied runs (2)\nPower varied runs (3)\nProbe retraction (4)\nSpontaneous recordings (5)?\nPlease enter 1,2,3,4,or 5: ';

str1 = input(prompt1, 's');
switch str1
    case '1'
        pwd_folder='Duty cycle varied runs\';
        %  filename='run42+41+36+43+44'; %for duty
        prompt2='Which trial run would you like to process? Please enter a number (36 or 41-44): ';
        str2 = input(prompt2, 's');
        filename= strcat(pwd_folder, 'run',str2); 
    case '2'
        pwd_folder='Frequncy varied runs\';
        % filename='run40+36+37+38+39'; %for frequency
        prompt2='Which trial run would you like to process? Please enter a number (36-40): ';
        str2 = input(prompt2, 's');
        filename= strcat(pwd_folder, 'run',str2); 
    case '3'
        pwd_folder='Power varied runs\';
        % % % filename='run36+45+46+47+48+49'; %for power
        prompt2='Which trial run would you like to process? Please enter a number (36 or 45-49): ';
        str2 = input(prompt2, 's');
        filename= strcat(pwd_folder, 'run',str2); 
    case '4'
        pwd_folder='Probe retraction\';
        % % filename='run36+51+52'; %for retraction
        prompt2='Which trial run would you like to process? Please enter a number (36 or 51-52): ';
        str2 = input(prompt2, 's');
        filename= strcat(pwd_folder, 'run',str2); 
    case '5'
        pwd_folder='Spontaneous recordings\';
        %run 13 was the data used for Fig 2 in the EO-Flex manuscript
        prompt2='Which trial run would you like to process? Please enter a number (1,5, or 13): ';
        str2 = input(prompt2, 's');
        filename= strcat(pwd_folder, 'run',str2); 
end

% base=y;
%%
clear; close all;
files=dir('*.mat');
filename=files.name;
pwd_folder='Recording 1';

%cd(strcat(pwd,'\', pwd_folder)) 

load(strcat(pwd,'\', filename));
filename=filename(1:(end-4));

base=amplifier;
base=base*1000;
Laser=func_gen;
tL=tE;
%enable following line if you don't want any tL data to be deleted
% if (strcmp(filename,'Recording 4000')==1)
%2/25/2021 exps
% if (strcmp(filename,'Recording 4')==1 || strcmp(filename,'Recording 10')==1|| strcmp(filename,'Recording 14')==1)
%4/12/2021
if (strcmp(filename,'Recording 1')==1 || strcmp(filename,'Recording 2')==1)
%4/11/2021 animal 1
% if (strcmp(filename,'Recording 2')==1|| strcmp(filename,'Recording 3')==1)
%4/11/2021 animal 2
% if (strcmp(filename,'Recording 1')==1|| strcmp(filename,'Recording 2')==1|| strcmp(filename,'Recording 3')==1 ...
%     || strcmp(filename,'Recording 4')==1|| strcmp(filename,'Recording 5')==1|| strcmp(filename,'Recording 6')==1 ...
%     || strcmp(filename,'Recording 9')==1|| strcmp(filename,'Recording 15')==1)

%4/10/2021
% if (strcmp(filename,'Recording 2')==1|| strcmp(filename,'Recording 3')==1|| strcmp(filename,'Recording 4')==1 ...
%    || strcmp(filename,'Recording 5')==1|| strcmp(filename,'Recording 6')==1|| strcmp(filename,'Recording 7')==1 ...
%    || strcmp(filename,'Recording 8')==1|| strcmp(filename,'Recording 9')==1|| strcmp(filename,'Recording 24')==1 ...
%    || strcmp(filename,'Recording 25')==1|| strcmp(filename,'Recording 26')==1|| strcmp(filename,'Recording 27')==1 ...
%    || strcmp(filename,'Recording 32')==1)

%4/16/2021
% if (strcmp(filename,'Recording 1')==1 || strcmp(filename,'Recording 2')==1|| strcmp(filename,'Recording 140')==1)

    clear tL;
end
% tL=picospritz;

% stimulus_rect_facecolor='#4DBEEE'; %cyan
% stimulus_rect_edgecolor='#4DBEEE'; %cyan

stimulus_rect_facecolor=[0.8 0.8 0.2549]; %yellow

r_tE=mean(diff(tE));
% sd=2.8; %standard deviation used for 5/31 exps
%  sd=4; %standard deviation used for 2/20 exps
sd=6; %6 used for chronic 2/27/21 exps

% ymax_volt=round(max(base),-1);
% ymin_volt=round(min(base),-1);

ymax_volt=800;
ymin_volt=-400;


pca_x_lim=[-2 2];
pca_y_lim=[-2 2];

% thresh=sd*(median(abs(base)/0.6745));
thresh=sd*(nanmedian(abs(base)/0.6745));


%Perform bandpass filtering
%first expand laser stimulus to have same length as base and tE
close all;
if exist('tL','var')==1
    
tE_to_tL=(tL(end))/(length(tE)-1);
tL_exp=(0:1:(length(tE)-1))'*tE_to_tL;

laser_exp=interp1(tL, Laser, tL_exp);
laser_thr=laser_exp./2;

lsr_thr=0.9*max(laser_thr);
laser_thr(laser_thr>lsr_thr)=4.9;
% laser_thr(laser_thr>lsr_thr)=4.9;


laser_thr(laser_thr<=lsr_thr)=0;
% plot(tL, laser_thr)

figure;
hold on;
axb=subplot(2,1,1);
plot(tL,picospritz, 'k');
%title('Spontaneous Activity');
xlabel('Time (s)');
ylabR=ylabel('Picospritz', 'Interpreter', 'tex');
set(gca,'box','off');
axs=subplot(2,1,2);
plot(tL_exp,laser_thr, 'Color', stimulus_rect_facecolor);
xlabel('Time (s)');
ylabS=ylabel('Function Gen');
%ylabS.Position= [-4.3410 2.5 -1.000]
%set(gca, 'ytick', [0.05,4.9], 'yticklabel', {'Off', 'On'});
set(gcf,'color','w');
set(gca,'box','off')
%title('Optical Stimulus');
linkaxes([axb,axs], 'x');
hold off;
print('-painters',strcat(filename,'_funcgen_vs_picospritz.pdf'),'-dpdf');

end

[vel, sampleTime, vel_res, frames_res]= locomotion_stab_2ch(ball_enc,20000,2,tE);


% Plot electrical recording, spritzer, and ball encoder
figure;
hold on;
axa=subplot(3,1,1);
plot(tE, base, 'k');
threshline = refline([0 thresh]);
threshline.Color = 'r';
%title('Spontaneous Activity');
xlabel('Time (s)');
ylabR=ylabel('Amplifier', 'Interpreter', 'tex');
set(gca,'box','off');
axb=subplot(3,1,2);
plot(tE,picospritz, 'Color', stimulus_rect_facecolor);
xlabel('Time (s)');
ylabS=ylabel('Picospritzer');
axc=subplot(3,1,3);
plot(tE,vel, 'Color', 'k');
xlabel('Time (s)');
ylabS=ylabel('Ball Encoder (mm/s)');
%ylabS.Position= [-4.3410 2.5 -1.000]
%set(gca, 'ytick', [0.05,4.9], 'yticklabel', {'Off', 'On'});
set(gcf,'color','w');
set(gca,'box','off')
%title('Optical Stimulus');
linkaxes([axa,axb,axc], 'x');
hold off;
print('-painters',strcat(filename,'_amp_pico_ballenc.pdf'),'-dpdf');
 


%%
close all;

spikes_thr_t=(base>thresh);
rising_edge_spikes_t=find(diff(spikes_thr_t)==1);
spike_array=zeros(length(base),1);

seconds_prespk=0.5/1000; %ms converted to s
seconds_postspk=2.5/1000; %ms converted to s, 2ms used for awhile
index_prespk=round(seconds_prespk/r_tE);
index_postspk=seconds_postspk/seconds_prespk*index_prespk; %ratio needs to be whole number

spike_index=(index_prespk+index_postspk+1);
PCA_time=(0:(spike_index-1))*r_tE*10^3; %time scale in ms

spikes_PCA_unnorm=zeros(spike_index,length(rising_edge_spikes_t));
spikes_discarded=cell(1,1);
spike_times_discarded=zeros(1,length(rising_edge_spikes_t));

shifted_max_times=zeros(1,length(rising_edge_spikes_t));
peak_window=round((500/r_tE)*10^-6); %maximum indicies around waveform to find max spike value (based on 500 micro-s)

%for all spikes locate the max and align window around that location, if
%this max location was previously recorded, then discard the spike
for i=1:length(rising_edge_spikes_t)
    if ((rising_edge_spikes_t(i)+index_postspk)<=length(tE) && (rising_edge_spikes_t(i)-index_postspk)>0)
        spike_zone=base((rising_edge_spikes_t(i)-index_prespk):(rising_edge_spikes_t(i)+index_postspk));
        [max_value, max_PCA_loc]=max(spike_zone); %if spike occurs, find maximum peak
        delta_to_max=max_PCA_loc-index_prespk-2;
        shifted_index=rising_edge_spikes_t(i)+delta_to_max;
            if (ismember(shifted_index,shifted_max_times)==0)
                %ensure shifted peak still exists within the data set
                if (shifted_index+index_postspk)<=length(tE)
                    %shift to new spikes location
                    spike_zone=base((shifted_index-index_prespk):(shifted_index+index_postspk));
                    spikes_PCA_unnorm(:,i)=spike_zone(:,1);
                    shifted_max_times(1,i)=shifted_index; %shifted_index;   
                else
                    %spike has a new maximum after shifted, 'discard' current one
                    spike_times_discarded(1,i)=rising_edge_spikes_t(i);
                    %record the spike
                    spikes_discarded{1,1}=[spikes_discarded{1,1} spike_zone(:,1)];
                end
            else
                %spike has a new maximum after shifted, 'discard' current one
               spike_times_discarded(1,i)=rising_edge_spikes_t(i);
                %record the spike
               spikes_discarded{1,1}=[spikes_discarded{1,1} spike_zone(:,1)];
            end
    end
end

%ensure arrays are of same length
spikes_PCA_unnorm=spikes_PCA_unnorm(:,~any(spike_times_discarded,1)); %removes discarded spikes
rising_edge_spikes=rising_edge_spikes_t(~any(spike_times_discarded,1),1);
spike_array(rising_edge_spikes)=1;
spikes_s=tE(rising_edge_spikes)';

%plot all of the discarded spikes
if isempty(spikes_discarded{1,1})~=1
    figure;
    plot(PCA_time,spikes_discarded{1,1})
end

Fontsize=30;
FontsizeAxis=22;
figure;
plot(PCA_time,spikes_PCA_unnorm, 'Color', 'k')
ax=ancestor(gca, 'axes');
xrule=ax.XAxis;
yrule=ax.YAxis;
ylim([ymin_volt ymax_volt]);
xL=xlabel('Time (ms)');
yL=ylabel(['Potential (',char(181),'V)'], 'Interpreter', 'tex');

xrule.FontSize =FontsizeAxis;
yrule.FontSize =FontsizeAxis;
xL.FontSize=Fontsize;
yL.FontSize=Fontsize;
width=4;
height=width;
pos = get(gcf, 'Position');
set(gcf, 'Position', [pos(1) pos(2) width*100, height*100]); %<- Set size
set(gcf,'color','w');
set(gca,'box','off')

ISI_all_spikes=diff(spikes_s);
[coeff,score,latent, tsquared, explained, mu, spikes_PCA]=CalculatePCA(spikes_PCA_unnorm,sd,thresh,filename,pca_x_lim,pca_y_lim);

%%
% find optical trials, use the expanded laser time signal (tL_exp) to find locations of optical stimulus
close all;
if exist('tL','var')==1
    dy=diff(laser_thr);
    dy=[0; dy];
    rising_ind=find(dy==4.9);
    falling_ind=find(dy==-4.9);
    
    %for 2/25/2021 exps
%     if strcmp(filename,'Recording 6')
%         falling_ind(1:2)=[];
%         rising_ind(1)=[];
%     end
    r_tL=tE_to_tL; %sampling time, s same as tE due to resampling
    %calculate epochs
   
 
    pulse_width=(falling_ind(1)-rising_ind(1))*r_tL*10^3; %in ms

    if strcmp(pwd_folder,'Frequncy varied runs\')
    %used for frequency varied runs
        pre_trial_time=5*10^-3; %time before pulse, s
        post_trial_time=15*10^-3; %time after pulse, s 
    else
        pre_trial_time=10*10^-3; %time before pulse, s
        post_trial_time=30*10^-3; %time after pulse, s 
    end

%     pre_trial=round(pre_trial_time/r_tE); %ind width to save before pulse
%     post_trial=round(post_trial_time/r_tE)+(round(pulse_width)*10^-3/r_tE); %ind width to save after pulse
    
%for fixed trial window
    pre_trial=round(pre_trial_time/r_tE); %ind width to save before pulse
    post_trial=round(post_trial_time/r_tE)+(round(100)*10^-3/r_tE); %ind width to save after pulse

    trial_time_ms=((-pre_trial:1:post_trial)*r_tE*10^3)'; %time scale for each trial in ms after pulse
    trial_time_s=((-pre_trial:1:post_trial)*r_tE)'; %time scale for each trial in s after pulse

    trials=zeros(length(trial_time_s),length(rising_ind));
    trial_times=zeros(length(trial_time_s),length(rising_ind));

    for i=1:length(rising_ind)
        if (rising_ind(i))<=(length(tE)-post_trial)
            trials(:,i)=spike_array((rising_ind(i)-pre_trial):(rising_ind(i)+post_trial));
            trial_times(:,i)=tE((rising_ind(i)-pre_trial):(rising_ind(i)+post_trial));
        else
            tmp=spike_array((rising_ind(i)-pre_trial):(length(tE)));
            trials(:,i)=[tmp; zeros(length(trial_time_ms)-length(tmp),1)];
            tmp=tE((rising_ind(i)-pre_trial):(length(tE)));
            trial_times(:,i)=[tmp; zeros(length(trial_time_ms)-length(tmp),1)];
        end
    end

    num_trials=size(trials,2);
    spikes_raster=cell(num_trials,1);
    AllFiringRate = zeros(length(trial_time_s),num_trials);

    spike_times_trials_abs=cell(1,1);
    spike_times_ms_after_stim=cell(1,1);

    for j=1:num_trials
        abs_time=trial_times(:,j);
        spk_array_log=logical(trials(:,j)); %find spikes by looking for rising edges 
        spikes_ms=trial_time_ms(spk_array_log);
        spikes_raster{j}=spikes_ms';    
        spike_times_ms_after_stim{1}=[spike_times_ms_after_stim{1} spikes_ms'];
        spike_times_trials_abs{1}=[spike_times_trials_abs{1} abs_time(spk_array_log)'];
        
        spikes_raster_s=trial_time_s(spk_array_log);

        beta=length(spikes_raster_s)^(4/5);
        %Using Alpha and Beta from: Ahmadi et al, "Esimation of neuronal-
        % firing rate using Bayesian Adaptive Kernel Smoother"
        alpha=4;
        [FiringRate, h] = BAKS(spikes_ms,trial_time_ms,alpha,beta);
        AllFiringRate(:,j)=FiringRate;
    end
    avg_firing=mean(AllFiringRate,2);
end


% generate psth plot. calculate the location of spikes for each trial
Fontsize=30;
FontsizeAxis=22;
if exist('tL','var')==1
 
    figure;
    rectangle('Position',[0 0.5 pulse_width num_trials+10],'FaceColor', stimulus_rect_facecolor, 'EdgeColor',  stimulus_rect_facecolor);
    hold on;

    plotRaster(spikes_raster, trial_time_ms, 'k', 0,1);
    hold off;

    ax=ancestor(gca, 'axes');
    xrule=ax.XAxis;
    yrule=ax.YAxis;

    xL=xlabel('Pulse time (ms)');

    yL=ylabel('Pulse number');

    xrule.FontSize =FontsizeAxis;
    yrule.FontSize =FontsizeAxis;
    xL.FontSize=Fontsize;
    yL.FontSize=Fontsize;

    width=4;
    height=width;
    pos = get(gcf, 'Position');
    set(gcf, 'Position', [pos(1) pos(2) width*100, height*100]); %<- Set size
    % set(gca, 'FontSize', Fontsize);
    set(gcf,'color','w');
    set(gca,'box','off')
%     if strcmp(pwd_folder,'Frequncy varied runs\')
%         xticks([-10 -5 0 5 10 15])
%         xticklabels({'-10','-5','0','5','10','15'})
%     else
%         xticks([-10 0 10 20 30])
%         xticklabels({'-10','0','10','20','30'})
%     end

%     xlim([-pre_trial_time*10^3 post_trial_time*10^3]);
    % saveas(gcf, strcat(filename,sprintf('_%i_uV_%isd_raster.pdf',round(thresh,0),sd)));
    % saveas(gcf, strcat(filename,sprintf('_%i_uV_%isd_raster.eps',round(thresh,0),sd)),'epsc');
    print('-painters',strcat(filename,sprintf('_%i_uV_%isd_raster.pdf',round(thresh,0),sd)),'-dpdf');

    %plot spikes according to time after stimulus
    figure;
    Fontsize=22;
    FontsizeAxis=16;
    blue_rect_height=300;
    rectangle('Position',[0 0 pulse_width blue_rect_height],'FaceColor', stimulus_rect_facecolor, 'EdgeColor',  stimulus_rect_facecolor);
    hold on;
    h=histogram(spike_times_ms_after_stim{1},100);
    h.FaceColor = 'k';
    h.EdgeColor = 'k';
    ax=ancestor(gca, 'axes');
    xrule=ax.XAxis;
    yrule=ax.YAxis;

    xL=xlabel('Pulse Time (ms)');

    yL=ylabel('Spikes');

    xrule.FontSize =FontsizeAxis;
    yrule.FontSize =FontsizeAxis;
    xL.FontSize=Fontsize;
    yL.FontSize=Fontsize;
    %for 4/11/2021 runs
    %ylim([0 300]);
    %otherwise
    ylim([0 100]);
    width=4;
    height=width;
    pos = get(gcf, 'Position');
    set(gcf, 'Position', [pos(1) pos(2) width*100, height*100]); %<- Set size
    set(gcf,'color','w');
    set(gca,'box','off')
    print('-painters',strcat(filename,sprintf('_%i_uV_%isd_totalled_raster.pdf',round(thresh,0),sd)),'-dpdf');
    hold off;
%
    %plot BAKS calculation for firing rate after each pulse
    figure;
    blue_rect_height=max(avg_firing); %800 for power variations
    if blue_rect_height==0
        blue_rect_height=100;
    end
    
    rectangle('Position',[0 0 pulse_width blue_rect_height],'FaceColor', stimulus_rect_facecolor, 'EdgeColor',  stimulus_rect_facecolor);
    hold on;
    %plot(ty_ms,y_ms, 'Color', 'k');
    plot(trial_time_ms,avg_firing, 'Color', 'k', 'LineWidth', 2);

    ax=ancestor(gca, 'axes');
    xrule=ax.XAxis;
    yrule=ax.YAxis;

    xL=xlabel('Pulse time (ms)');
    yL=ylabel('Spike rate (kHz)');

    xrule.FontSize =FontsizeAxis;
    yrule.FontSize =FontsizeAxis;
    xL.FontSize=Fontsize;
    yL.FontSize=Fontsize;

    width=4;
    height=width;
    pos = get(gcf, 'Position');
    set(gcf, 'Position', [pos(1) pos(2) width*100, height*100]); %<- Set size
    % set(gca, 'FontSize', Fontsize);
%     xlim([-pre_trial_time*10^3 post_trial_time*10^3]);
%     if strcmp(pwd_folder,'Frequncy varied runs\')
%         xticks([-10 -5 0 5 10 15])
%         xticklabels({'-10','-5','0','5','10','15'})
%     else
%         xticks([-10 0 10 20 30])
%         xticklabels({'-10','0','10','20','30'})
%     end
    ylim([0 blue_rect_height]);
    %title(sprintf('KBO Spike Rate with threshold at %i std',sd));
    set(gcf,'color','w');
    set(0,'DefaultAxesTitleFontWeight','normal');
    hold off;

    % saveas(gcf, strcat(filename,sprintf('_%i_uV_%isd_BAKS.pdf',round(thresh,0),sd)));
    print('-painters',strcat(filename,sprintf('_%i_uV_%isd_BAKS.pdf',round(thresh,0),sd)),'-dpdf');

    %plot SCALED text size, BAKS calculation for firing rate after each pulse
    Fontsize=22;
    FontsizeAxis=16;
    figure;
    rectangle('Position',[0 0 pulse_width blue_rect_height],'FaceColor', stimulus_rect_facecolor, 'EdgeColor',  stimulus_rect_facecolor);
    hold on;
    %plot(ty_ms,y_ms, 'Color', 'k');
    plot(trial_time_ms,avg_firing, 'Color', 'k', 'LineWidth', 2);

    ax=ancestor(gca, 'axes');
    xrule=ax.XAxis;
    yrule=ax.YAxis;

    xL=xlabel('Pulse time (ms)');
    % xticks([PCA_time(1) round(PCA_time(end))]);
    yL=ylabel('Spike rate (kHz)');

    xrule.FontSize =FontsizeAxis;
    yrule.FontSize =FontsizeAxis;
    xL.FontSize=Fontsize;
    yL.FontSize=Fontsize;

    width=4;
    height=width;
    pos = get(gcf, 'Position');
    set(gcf, 'Position', [pos(1) pos(2) width*100, height*100]); %<- Set size
    % set(gca, 'FontSize', Fontsize);
%     xlim([-pre_trial_time*10^3 post_trial_time*10^3]);
%     if strcmp(pwd_folder,'Frequncy varied runs\')
%         xticks([-10 -5 0 5 10 15])
%         xticklabels({'-10','-5','0','5','10','15'})
%     else
%         xticks([-10 0 10 20 30])
%         xticklabels({'-10','0','10','20','30'})
%     end
    ylim([0 blue_rect_height]);
    %title(sprintf('KBO Spike Rate with threshold at %i std',sd));
    set(gcf,'color','w');
    set(0,'DefaultAxesTitleFontWeight','normal');
    hold off;

    % saveas(gcf, strcat(filename,sprintf('_%i_uV_%isd_BAKS.pdf',round(thresh,0),sd)));
    print('-painters',strcat(filename,sprintf('_%i_uV_%isd_BAKS_scaled.pdf',round(thresh,0),sd)),'-dpdf');


    % plot waveforms of spontaneous activity before stimulation
    figure;
    % tmp=tE<=(tE(rising_ind(1)));
    if round(tE(rising_ind(1))) > (tE(rising_ind(1)))
        tmp=tE<=(round(tE(rising_ind(1)))-1);
    else
    %     tmp=tE<=(round(tE(rising_ind(1))));
        tmp=tE<=(round(tE(rising_ind(1)),1));
    end
    spont_time=tE(tmp);
    spont_base=base(tmp);
    hold on;

    plot(spont_time,spont_base, 'Color', 'k');
    threshline = refline([0 thresh]);
    threshline.Color = 'r';
    ax=ancestor(gca, 'axes');
    xrule=ax.XAxis;
    yrule=ax.YAxis;

    xL=xlabel('Time (s)');
    % xticks([PCA_time(1) round(PCA_time(end))]);
    yL=ylabel(['Potential (',char(181),'V)'], 'Interpreter', 'tex');
    xrule.FontSize =FontsizeAxis;
    yrule.FontSize =FontsizeAxis;
    xL.FontSize=Fontsize;
    yL.FontSize=Fontsize;

    width=4;
    height=width;
    pos = get(gcf, 'Position');
    set(gcf, 'Position', [pos(1) pos(2) width*100, height*100]); %<- Set size
    % set(gca, 'FontSize', Fontsize);
%     xlim([0 8]);
%     xticks([0 2 4 6 8])
%     xticklabels({'0','2','4','6','8'})
%     ylim([ymin_volt/2 ymax_volt/2]);
    set(gcf,'color','w');
    set(0,'DefaultAxesTitleFontWeight','normal');
    hold off;

    % saveas(gcf, strcat(filename,sprintf('_%i_uV_%isd_spont_waveforms.pdf',round(thresh,0),sd)));
    print('-painters',strcat(filename,sprintf('_%i_uV_%isd_spont_waveforms.pdf',round(thresh,0),sd)),'-dpdf');

    hold off;
end

%% Plot epoch/stimulus trains
% Plot example of pulse train only if optical data exists
Fontsize=30;
FontsizeAxis=22;
if exist('tL','var')==1
     tmp=diff(falling_ind);
    
    if (std(tmp) < 0.1*mean(tmp))
        %if std is small relative to the average, then the entire series of
        %pulses is a train.  I.e. there is only 1 train
        train_end=falling_ind(tmp>mean(tmp));
        train_end=train_end(end);
        train_start=rising_ind(1); %add first train index
        test='only one pulse train'
    else
        train_end=falling_ind(tmp>mean(tmp)); %end index of pulse train
        train_end=[train_end; falling_ind(end)]; %add last train index
        tmp=diff(rising_ind);
        train_start=rising_ind(logical([0; tmp>mean(tmp)])); %startindex of pulse train
        train_start=[rising_ind(1); train_start]; %add first train index
    end
    
    if (length(train_end) ~= length(train_start))
       train_width=train_end(1:length(train_start))-train_start;
    else
       train_width=train_end-train_start(1:length(train_end));
    end
    %remove trains that are not full
    bad_trains=train_width < abs(median(train_width)-std(train_width));
    train_start(bad_trains)=[];
    train_end(bad_trains)=[];
    if (length(train_end) ~= length(train_start))
       train_width=train_end(1:length(train_start))-train_start;
    else
       train_width=train_end-train_start(1:length(train_end));
    end
    %using cleaned trains, determine pulses in each train
    tmp=rising_ind(logical((rising_ind<=train_end(1)).*(rising_ind>=train_start(1))));
    pulses_in_train=([tmp falling_ind(logical((falling_ind<=train_end(1)).*(falling_ind>=train_start(1))))]-tmp(1))*r_tL;

    
    prepost_train=round(0.2/r_tE); %seconds recorded before and after train
    num_trains=length(train_start);

    total_epoch_length=round(median(train_width))+2*prepost_train;
    epochs=zeros(total_epoch_length+1,num_trains);
    
    epoch_waveforms=zeros(total_epoch_length+1,num_trains);
    
    epoch_spike_times=cell(num_trains,1);
    epoch_spike_times_all=cell(1,1);
    epoch_spike_times_abs=cell(num_trains,1);
    epoch_spike_times_abs_all=cell(1,1);
    
    epoch_time=((-prepost_train:1:(total_epoch_length-prepost_train))*r_tL)';
    for j=1:num_trains
        if (train_start(j)-prepost_train+total_epoch_length)<length(tE)
            epoch_waveforms(:,j)=base((train_start(j)-prepost_train):(train_start(j)-prepost_train+total_epoch_length));
            tmp=tE((train_start(j)-prepost_train):(train_start(j)-prepost_train+total_epoch_length));
            epochs(:,j)=spike_array((train_start(j)-prepost_train):(train_start(j)-prepost_train+total_epoch_length));
            epoch_spike_times{j,1}=epoch_time(logical(epochs(:,j)),1)';
            epoch_spike_times_abs{j,1}=tmp(logical(epochs(:,j)),1)';
        end
    end
    
    %any unused epochs are removed due to not having enough time after
    %recording ends
    %commented out for processing ball vs no ball rolling whisker stim runs
%     epochs=epochs(:,any(epochs));
%     epoch_spike_times=epoch_spike_times(any(epochs));
%     epoch_spike_times_abs=epoch_spike_times_abs(any(epochs));
    
    epoch_spike_times_abs_all{1,1}=[epoch_spike_times_abs{:}]';
    epoch_spike_times_all{1,1}=[epoch_spike_times{:}]';
    
    %    num_trains=size(epochs, 2);

    
    if tE(end)<=60
            endpoint=60;
    else
        endpoint=tE(end);
    end
    first_train=(train_start(1)-prepost_train):1:(train_end(1)+prepost_train);
    first_train_time=(0:1:(length(first_train)-1))'.*r_tL;

    figure;
    hold on;
    blue_rect_height=1;
    for i=1:length(pulses_in_train)
        rectangle('Position',[pulses_in_train(i,1) 0 (pulses_in_train(i,2)-pulses_in_train(i,1)) blue_rect_height],'FaceColor', stimulus_rect_facecolor, 'EdgeColor', stimulus_rect_facecolor);
    end
   
    ax=ancestor(gca, 'axes');
    xrule=ax.XAxis;
    yrule=ax.YAxis;
%     xlim([epoch_time(1) 1.2]);
%     xticks([-0.4 0 0.4 0.8 1.2])
%     xticklabels({'-0.4','0','0.4','0.8','1.2'})
    xL=xlabel('Pulse train time (s)');
    yL=ylabel('Stimulus');
    set(gca, 'ytick', [0 1], 'yticklabel', {'Off', 'On'});
    xrule.FontSize =FontsizeAxis;
    yrule.FontSize =FontsizeAxis;
    xL.FontSize=Fontsize;
    yL.FontSize=Fontsize;

    width=4;
    height=width;
    pos = get(gcf, 'Position');
    set(gcf, 'Position', [pos(1) pos(2) width*100, height*100]); %<- Set size
    set(gcf,'color','w');
    set(gca,'box','off')
    hold off;
    print('-painters',strcat(filename,'_single_pulse_train.pdf'),'-dpdf');

    figure;
    hold on;
    blue_rect_height=1;
    for i=1:length(falling_ind)
        rectangle('Position',[tL_exp(rising_ind(i)) 0 tL_exp((falling_ind(i)-rising_ind(i))) blue_rect_height],'FaceColor', stimulus_rect_facecolor, 'EdgeColor',  stimulus_rect_facecolor);
    end

    width=8;
    height=width/6;
    set(gcf, 'Position', [pos(1) pos(2) width*100, height*100]); %<- Set size
    if endpoint<60
        xticks([0 20 40 floor(tL_exp(end))])
        xticklabels({'0','20','40',sprintf('%i',floor(tL_exp(end)))})%
        xlim([0 tL(end)]);
    end
    xlim([0 endpoint])
    set(gca, 'ytick', [0 1], 'yticklabel', {'Off', 'On'});
    set(gca, 'FontSize', 16);
    set(gca, 'FontSize', 16);
    set(gcf,'color','w');
    xlabel('Time (s)');
    ylabel('Stimulus');

    set(gca,'box','off')
    hold off;
    print('-painters',strcat(filename,'_full_pulse_train.pdf'),'-dpdf');
    
    %plot waveforms of a single pulse train
    figure;

    hold on;
    % plot(epoch_time,epoch_waveforms(:,24), 'Color', 'k');
%     plot(epoch_time,epoch_waveforms(:,22), 'Color', 'k');
     plot(epoch_time,epoch_waveforms(:,1), 'Color', 'k');

    threshline = refline([0 thresh]);
    threshline.Color = 'r';
    ax=ancestor(gca, 'axes');
    xrule=ax.XAxis;
    yrule=ax.YAxis;
    xL=xlabel('Pulse train time (s)');
    % xticks([PCA_time(1) round(PCA_time(end))]);
    yL=ylabel(['Potential (',char(181),'V)'], 'Interpreter', 'tex');
    xrule.FontSize =FontsizeAxis;
    yrule.FontSize =FontsizeAxis;
    xL.FontSize=Fontsize;
    yL.FontSize=Fontsize;

    width=4;
    height=width;
    pos = get(gcf, 'Position');
    set(gcf, 'Position', [pos(1) pos(2) width*100, height*100]); %<- Set size
    % set(gca, 'FontSize', Fontsize);
%     xlim([epoch_time(1) 1.2]);
%     xticks([-0.4 0 0.4 0.8 1.2])
%     xticklabels({'-0.4','0','0.4','0.8','1.2'})
%     ylim([ymin_volt/2 ymax_volt/2]);
    set(gcf,'color','w');
    set(0,'DefaultAxesTitleFontWeight','normal');
    hold off;

    % saveas(gcf, strcat(filename,sprintf('_%i_uV_%isd_train_waveforms.pdf',round(thresh,0),sd)));
    print('-painters',strcat(filename,sprintf('_%i_uV_%isd_train_waveforms.pdf',round(thresh,0),sd)),'-dpdf');
    
    Fontsize=30;
    FontsizeAxis=22;
    %plot psth for each pulse train
    epoch_AllFiringRate = zeros(length(epoch_time),num_trains);
    train_spikes=cell(num_trains,1);

    if train_width*r_tE<10
        
        for j=1:num_trains
            spk_array=logical(epochs(:,j));

            tmp=epoch_time(spk_array);
            train_spikes{j}=tmp; 

            beta=length(tmp)^(4/5);
            %from Ahmadi et al
            alpha=4;
            [FiringRate, h] = BAKS(tmp,epoch_time,alpha,beta);
            epoch_AllFiringRate(:,j)=FiringRate;
        end
    else
        FiringRate=[0];
        h=[0];
        epoch_AllFiringRate(:,j)=FiringRate;
    end
    epoch_avg_firing=mean(epoch_AllFiringRate,2);

    figure;
    blue_rect_height=max(epoch_avg_firing);
    
    if blue_rect_height==0
        blue_rect_height=100;
    end
    hold on;
    for i=1:length(pulses_in_train)
        rectangle('Position',[pulses_in_train(i,1) 0 (pulses_in_train(i,2)-pulses_in_train(i,1)) blue_rect_height],'FaceColor', stimulus_rect_facecolor, 'EdgeColor',  stimulus_rect_facecolor);
    end

    plot(epoch_time,epoch_avg_firing, 'Color', 'k', 'LineWidth', 2);

    ax=ancestor(gca, 'axes');
    xrule=ax.XAxis;
    yrule=ax.YAxis;

    xL=xlabel('Pulse train time (s)');
    yL=ylabel('Spike rate (Hz)');

    xrule.FontSize =FontsizeAxis;
    yrule.FontSize =FontsizeAxis;
    xL.FontSize=Fontsize;
    yL.FontSize=Fontsize;

    width=4;
    height=width;
    pos = get(gcf, 'Position');
    set(gcf, 'Position', [pos(1) pos(2) width*100, height*100]); %<- Set size
%     xlim([epoch_time(1) 1.2]);
%     xticks([-0.4 0 0.4 0.8 1.2])
%     xticklabels({'-0.4','0','0.4','0.8','1.2'})
    ylim([0 blue_rect_height]);
    set(gcf,'color','w');
    set(0,'DefaultAxesTitleFontWeight','normal');
    hold off;

    print('-painters',strcat(filename,sprintf('_%i_uV_%isd_epochs_BAKS.pdf',round(thresh,0),sd)),'-dpdf');

%==================================================
    % Scaled text on axis of BAKS epochs plot for fig3
    Fontsize=22;
    FontsizeAxis=16;
    blue_rect_height=20;

    figure;
    hold on;
    for i=1:length(pulses_in_train)
        rectangle('Position',[pulses_in_train(i,1) 0 (pulses_in_train(i,2)-pulses_in_train(i,1)) blue_rect_height],'FaceColor', stimulus_rect_facecolor, 'EdgeColor',  stimulus_rect_facecolor);
    end

    plot(epoch_time,epoch_avg_firing, 'Color', 'k', 'LineWidth', 2);
    ax=ancestor(gca, 'axes');
    xrule=ax.XAxis;
    yrule=ax.YAxis;

    xL=xlabel('Pulse train time (s)');
    % xticks([PCA_time(1) round(PCA_time(end))]);
    yL=ylabel('Spike rate (Hz)');


    xrule.FontSize =FontsizeAxis;
    yrule.FontSize =FontsizeAxis;
    xL.FontSize=Fontsize;
    yL.FontSize=Fontsize;

    width=4;
    height=width;
    pos = get(gcf, 'Position');
    set(gcf, 'Position', [pos(1) pos(2) width*100, height*100]); %<- Set size
    % set(gca, 'FontSize', Fontsize);
%     xlim([epoch_time(1) 1.2]);
%     xticks([-0.4 0 0.4 0.8 1.2])
%     xticklabels({'-0.4','0','0.4','0.8','1.2'})
    ylim([0 blue_rect_height]);
    set(gcf,'color','w');
    set(0,'DefaultAxesTitleFontWeight','normal');
    hold off;

    % saveas(gcf, strcat(filename,sprintf('_%i_uV_%isd_epochs_BAKS.pdf',round(thresh,0),sd)));
    print('-painters',strcat(filename,sprintf('_%i_uV_%isd_epochs_BAKS_scaled.pdf',round(thresh,0),sd)),'-dpdf');

    % plot raster of pulse train
    figure;
    hold on;
    for i=1:length(pulses_in_train)
        rectangle('Position',[pulses_in_train(i,1) 0 (pulses_in_train(i,2)-pulses_in_train(i,1)) size(epochs,2)+1],'FaceColor', stimulus_rect_facecolor, 'EdgeColor',  stimulus_rect_facecolor);
        %rectangle('Position',[(rising_ind(i)-rising_ind(1))*r_tL 0 (falling_ind(i)-rising_ind(i))*r_tL blue_rect_height],'FaceColor','#4DBEEE', 'EdgeColor', '#4DBEEE');
    end
    % plotSpikeRaster(train_spikes,'PlotType','vertline','XLimForCell',[epoch_time(1) epoch_time(end)]);
    plotRaster ( train_spikes , epoch_time, 'k', 0,1 )
    ax=ancestor(gca, 'axes');
    xrule=ax.XAxis;
    yrule=ax.YAxis;

    xL=xlabel('Pulse train time (s)');
    % xticks([PCA_time(1) round(PCA_time(end))]);
    yL=ylabel('Pulse train');

    xrule.FontSize =FontsizeAxis;
    yrule.FontSize =FontsizeAxis;
    xL.FontSize=Fontsize;
    yL.FontSize=Fontsize;

    width=4;
    height=width;
    pos = get(gcf, 'Position');
    set(gcf, 'Position', [pos(1) pos(2) width*100, height*100]); %<- Set size
    % set(gca, 'FontSize', Fontsize);
    set(gcf,'color','w');
    set(gca,'box','off')
%     xlim([-0.2 1.2]);
%     xticks([0 0.4 0.8 1.2])
%     xticklabels({'0','0.4','0.8','1.2'})
    % saveas(gcf, strcat(filename,sprintf('_%i_uV_%isd_epochs_raster.pdf',round(thresh,0),sd)));
    hold off;
    print('-painters',strcat(filename,sprintf('_%i_uV_%isd_epochs_raster.pdf',round(thresh,0),sd)),'-dpdf');

end
 
    
%%
% Single plots for stimulus, recording
% close all;
if tE(end)<=60
        endpoint=60;
else
    endpoint=tE(end);
end
Fontsize=20;
FontsizeAxis=16;
tmp={spikes_s};
plotRaster(tmp,tE, 'k', 1,1);
h = gca; h.YAxis.Visible = 'off'; 
xlabel('Time (s)');
width=8;
height=width/8;
set(gca, 'FontSize', 16);
set(gcf, 'Position', [100 100 width*100, height*100]); %<- Set size
xlim([0 endpoint]);
print('-painters',strcat(filename,sprintf('_%i_uV_%isd_raster_full.pdf',round(thresh,0),sd)),'-dpdf');

%plot ISI for full raster
ISI=diff(spikes_s)*10^3;
figure;
hs=histogram(ISI,[0:0.5:10]);
hs.FaceColor = 'k';
hs.EdgeColor = 'k';
ax=ancestor(gca, 'axes');
xrule=ax.XAxis;
yrule=ax.YAxis;

xL=xlabel('Spike Interval (ms)');

yL=ylabel('Count');

xrule.FontSize =FontsizeAxis;
yrule.FontSize =FontsizeAxis;
xL.FontSize=Fontsize;
yL.FontSize=Fontsize;

width=4;
height=width;
pos = get(gcf, 'Position');
set(gcf, 'Position', [pos(1) pos(2) width*100, height*100]); %<- Set size
set(gcf,'color','w');
set(gca,'box','off')
print('-painters',strcat(filename,sprintf('_%i_uV_%isd_ISI_hist.pdf',round(thresh,0),sd)),'-dpdf');
%
figure;
plot(tE,base, 'k');
threshline = refline([0 thresh]);
threshline.Color = 'r';

ax=ancestor(gca, 'axes');
xrule=ax.XAxis;
yrule=ax.YAxis;

xL=xlabel('Time (s)');
yL=ylabel(['Potential (',char(181),'V)'], 'Interpreter', 'tex');
xrule.FontSize =FontsizeAxis;
yrule.FontSize =FontsizeAxis;
xL.FontSize=Fontsize;
yL.FontSize=Fontsize;

width=8;
height=width/3;
%pos = get(gcf, 'Position');
% set(gca, 'FontSize', 16);
ylim([ymin_volt ymax_volt]);
set(gcf, 'Position', [100 100 width*100, height*100]); %<- Set size
%title('10 Hz Stimulation');
% xlabel('Time (s)');
% ylabel(['Potential (',char(181),'V)'], 'Interpreter', 'tex');
set(gca,'box','off');
set(gcf,'color','w')
xlim([0 endpoint]);

%filename='run36+50+51+52'; %for retraction
%filename='run42+41+36+43+44'; %for duty
% filename='run36+45+46+47+48+49'; %for power
%filename='run40+36+37+38+39'; %for frequency
if strcmp(filename,'run36+50+51+52') || strcmp(filename,'run42+41+36+43+44')...
        || strcmp(filename,'run36+45+46+47+48+49') || strcmp(filename,'run40+36+37+38+39') || strcmp(filename,'run36+51+52')
    h = gca; h.XAxis.Visible = 'off';
    test=1
end
    
% saveas(gcf, strcat(filename,sprintf('_%i_uV_%isd_full_w_thr.pdf',round(thresh,0),sd)));
print('-painters',strcat(filename,sprintf('_%i_uV_%isd_full_w_thr.pdf',round(thresh,0),sd)),'-dpdf');

if exist('tL_exp','var')==1
    figure;
%plot(tL_exp(100000:155000),laser_exp(100000:155000), 'Color', '[0.4 0.4 0.4]');
    plot(tL_exp,laser_exp, 'Color', '[0.4 0.4 0.4]');
    width=8;
    height=width/5;
    pos = get(gcf, 'Position');
    set(gcf, 'Position', [pos(1) pos(2) width*100, height*100]); %<- Set size
    xlabel('Time (s)');
    set(gca, 'FontSize', 16);
    ylabel('Stimulus');
    set(gca, 'ytick', [0.05,4.9], 'yticklabel', {'Off', 'On'});
    set(gcf,'color','w');
    set(gca,'box','off')
    %xlim([5 tL_exp(155000)]);
    xlim([0 endpoint]);
    ylim([0 4.9]);
%     saveas(gcf, strcat(filename,'_stimulus.pdf'));
    print('-painters',strcat(filename,'_stimulus.pdf'),'-dpdf');

end

% plot moving window of waveforms 
figure;
tmp=(12.5<=tE & tE<=(13.5));

spont_time=tE(tmp);
spont_base=base(tmp);
hold on;

plot(spont_time,spont_base, 'Color', 'k');
threshline = refline([0 thresh]);
threshline.Color = 'r';

ax=ancestor(gca, 'axes');
xrule=ax.XAxis;
yrule=ax.YAxis;

xL=xlabel('Time (s)');
% xticks([PCA_time(1) round(PCA_time(end))]);
yL=ylabel(['Potential (',char(181),'V)'], 'Interpreter', 'tex');
xrule.FontSize =FontsizeAxis;
yrule.FontSize =FontsizeAxis;
xL.FontSize=Fontsize;
yL.FontSize=Fontsize;

width=8;
height=width/2;
pos = get(gcf, 'Position');
set(gcf, 'Position', [pos(1) pos(2) width*100, height*100]); %<- Set size
% set(gca, 'FontSize', Fontsize);
% xlim([0 spont_time(end)]);
% xticks([10 11 12 13 14 15])
% xticklabels({'0','1','2','3','4','5'})
xticks([12.5 13 13.5])
xticklabels({'0','0.5','1'})

ylim([ymin_volt ymax_volt]);
set(gcf,'color','w');
set(0,'DefaultAxesTitleFontWeight','normal');
hold off;

% saveas(gcf, strcat(filename,sprintf('_%i_uV_%isd_spont_waveforms.pdf',round(thresh,0),sd)));
print('-painters',strcat(filename,sprintf('_%i_uV_%isd_window_waveforms.pdf',round(thresh,0),sd)),'-dpdf');

hold off;

%% fit to gaussian clutering

ColOrd = [0,    0.7,    0; 0.3,    0.2,    0.3;0.2314,    0.4784,    0.0157;       0.7686,0.6784,0.5;      0.635,0.078,0.184; 0.494,0.184,0.556; 0.466,0.674,0.188;    0.307,0.47,0.36;0.85,0.325,0.098;         0.929,0.694,0.125;0.301,0.745,0.733;...
    1,0,0;0,1,0;...
     0,0,0;1,0,1;...
     1,1,0; 0,0,1;...
    0,0.447,0.741;0.85,0.325,0.098;0.929,0.694,0.125;...
    0.494,0.184,0.556;0.466,0.674,0.188;0.635,0.078,0.184;...
    0.307,0.47,0.36; 0,0,0;1,0,0; 1,0,1;1,1,0;...
    0.301,0.745,0.733; 0,1,0; 0,0,1;...
    ];

close all;
PCx=1;
PCy=2;

PCA_time=(0:(spike_index-1))*r_tE*10^3; %time scale in ms

PCA_scores=[score(:,PCx),score(:,PCy)];
%eva_gmdist=evalclusters(PCA_scores,'gmdistribution','CalinskiHarabasz','KList',(1:15));
%original value used for optical spike sorting: 15
eva_gmdist=evalclusters(PCA_scores,'gmdistribution','CalinskiHarabasz','KList',(1:6));

i=eva_gmdist.OptimalK

%use this option to force cluster to a certain amount if there are errors
%with the estimation above.
 i=3;

gm=fitgmdist(PCA_scores,i,'Options', statset('Display','final','MaxIter',100,'TolFun',1e-5));
% gm=fitgmdist(PCA_scores,i,'Options', statset('Display','final','MaxIter',1000,'TolFun',1e-5));


figure;
scatter(score(:,PCx),score(:,PCy),10,'.') % Scatter plot with points of size 10
hold on;
gmPDF = @(x,y)reshape(pdf(gm,[x(:) y(:)]),size(x));
fcontour(gmPDF,[-1 1 -1 1])

[idx,nlogL,P,logpdf] = cluster(gm,PCA_scores);

mean(PCA_scores(idx==1),1)
mean(PCA_scores(idx==2),1)
if mean(PCA_scores(idx==1),1) > mean(PCA_scores(idx==2),1)
    ColOrd = [0,    0.7,    0; 0.3,    0.2,    0.3;0.2314,    0.4784,    0.0157;       0.7686,0.6784,0.5;      0.635,0.078,0.184; 0.494,0.184,0.556; 0.466,0.674,0.188;    0.307,0.47,0.36;0.85,0.325,0.098;         0.929,0.694,0.125;0.301,0.745,0.733;...
    1,0,0;0,1,0;...
     0,0,0;1,0,1;...
     1,1,0; 0,0,1;...
    0,0.447,0.741;0.85,0.325,0.098;0.929,0.694,0.125;...
    0.494,0.184,0.556;0.466,0.674,0.188;0.635,0.078,0.184;...
    0.307,0.47,0.36; 0,0,0;1,0,0; 1,0,1;1,1,0;...
    0.301,0.745,0.733; 0,1,0; 0,0,1;...
    ];
else
    ColOrd = [0.3,    0.2,    0.3;0,    0.7,    0;0.2314,    0.4784,    0.0157;       0.7686,0.6784,0.5;      0.635,0.078,0.184; 0.494,0.184,0.556; 0.466,0.674,0.188;    0.307,0.47,0.36;0.85,0.325,0.098;         0.929,0.694,0.125;0.301,0.745,0.733;...
    1,0,0;0,1,0;...
     0,0,0;1,0,1;...
     1,1,0; 0,0,1;...
    0,0.447,0.741;0.85,0.325,0.098;0.929,0.694,0.125;...
    0.494,0.184,0.556;0.466,0.674,0.188;0.635,0.078,0.184;...
    0.307,0.47,0.36; 0,0,0;1,0,0; 1,0,1;1,1,0;...
    0.301,0.745,0.733; 0,1,0; 0,0,1;...
    ];
end
hold off;

figure;
gscatter(score(:,PCx),score(:,PCy),idx,ColOrd);


% use optimal cluster size to sort PCA
%close all;



[PCA_avg, PCA_std, PCA_spike_locations]=CalcandPlotClusWavef(PCA_time,spikes_PCA_unnorm,idx,gm,...
    i,spikes_s,ColOrd,PCA_scores,pca_x_lim,pca_y_lim,ymin_volt, ymax_volt,PCx,PCy,filename,sd,thresh);

peak_amps=[PCA_avg(12,:); PCA_std(12,:)]';

if exist('tL','var')==1
    [PCA_trial_spike_locations, PCA_epoch_spike_locations,PCA_trial_spikes_colored, PCA_epoch_spikes_colored]=ColorSpikes(spikes_raster,...
        epoch_spike_times,spike_times_trials_abs, epoch_spike_times_abs_all, spikes_s, i,idx);
end
% Plot each component in the PCA+Waveforms plot by itself
%close all;

rem_cluster=0; %cluster number to remove. 0=none
if exist('tL','var')==1
    PlotClusterWvefms(tE,PCA_time,PCA_avg,PCA_std,PCA_spike_locations,score,idx,gm,ColOrd,pwd_folder,...
        filename, thresh, sd,ymin_volt, ymax_volt,PCx,PCy,pca_x_lim,pca_y_lim,rem_cluster,tL,trial_time_ms,pulse_width,...
        pulses_in_train,epoch_time,epochs,PCA_trial_spike_locations,PCA_epoch_spike_locations)
else
    PlotClusterWvefms(tE,PCA_time,PCA_avg,PCA_std,PCA_spike_locations,score,idx,gm,ColOrd,pwd_folder,...
    filename, thresh, sd,ymin_volt, ymax_volt,PCx,PCy,pca_x_lim,pca_y_lim,rem_cluster)

end
%plot colored raster
if exist('tL','var')==1
    plotColoredRaster(PCA_epoch_spikes_colored,epoch_time,PCA_trial_spikes_colored, trial_time_ms,...
            thresh,sd, pulse_width, pre_trial_time,post_trial_time, pulses_in_train, ColOrd,stimulus_rect_facecolor, pwd_folder,filename)
end

%% BAKS for full run figure
close all;
if tE(end)<=60
        endpoint=60;
else
    endpoint=tE(end);
end
Fontsize=20;
FontsizeAxis=16;
% if strcmp(filename,'run36+50+51+52') || strcmp(filename,'run42+41+36+43+44')...
%         || strcmp(filename,'run36+45+46+47+48+49') || strcmp(filename,'run40+36+37+38+39')
%     
% else
        %calc overall spike rate
    
    spk_array=logical(spike_array);%find spikes by looking for rising edges 
    spikes_s=tE(spk_array);
    beta=length(spikes_s)^(4/5);
    %from Ahmadi et al
    alpha=4;
    [FiringRate, h] = BAKS(spikes_s,tE,alpha,beta);

    %avg_firing=mean(FiringRate,2);
%%
    plot(tE,FiringRate, 'Color', 'k', 'LineWidth', 1);

    ax=ancestor(gca, 'axes');
    xrule=ax.XAxis;
    yrule=ax.YAxis;

    xL=xlabel('Time (s)');
    % xticks([PCA_time(1) round(PCA_time(end))]);
    yL=ylabel('Spike rate (Hz)');

    xrule.FontSize =FontsizeAxis;
    yrule.FontSize =FontsizeAxis;
    xL.FontSize=Fontsize;
    yL.FontSize=Fontsize;
    
    ylim([0 max(FiringRate)]);
    width=8;
    height=width/3;
    pos = get(gcf, 'Position');
    set(gcf, 'Position', [pos(1) pos(2) width*100, height*100]); %<- Set size
   
    set(gcf,'color','w');
    set(gca,'box','off')
    xlim([0 endpoint]);
%     xlim([tE(1) tE(end)]);
%     saveas(gcf, strcat(filename,sprintf('_%i_uV_%isd_BAKS_full.pdf',round(thresh,0),sd)));
    print('-painters',strcat(filename,sprintf('_%i_uV_%isd_BAKS_full.pdf',round(thresh,0),sd)),'-dpdf');

% end 
figure;
plot(tE,FiringRate, 'Color', 'k', 'LineWidth', 1);

xlabel('Time (s)');
ylabel('Spike Rate (Hz)');
ylim([0 max(FiringRate)]);
width=4;
height=2*width/3;
pos = get(gcf, 'Position');
set(gcf, 'Position', [pos(1) pos(2) width*100, height*100]); %<- Set size
set(gca, 'FontSize', 16);
set(gcf,'color','w');
set(gca,'box','off')
xlim([0 endpoint]);
%     xlim([tE(1) tE(end)]);
%     saveas(gcf, strcat(filename,sprintf('_%i_uV_%isd_BAKS_full.pdf',round(thresh,0),sd)));
print('-painters',strcat(filename,sprintf('_%i_uV_%isd_BAKS_full_plot_half_width.pdf',round(thresh,0),sd)),'-dpdf');

%% plot Full recording+raster+laser stimulus
figure;
% rectangle('Position',[0 0 12.2 size(trials,1)],'FaceColor','#4DBEEE', 'EdgeColor', '#4DBEEE');
hold on;
pos1=[0.1 0.65 0.7 0.25];
%axb=subplot(3,1,1);
axb=subplot('Position', pos1);
plot(tE,base, 'k');
threshline = refline([0 thresh]);
threshline.Color = 'r';
%title('Spontaneous Activity');
xlabel('Time (s)');
ylabel(['Potential (',char(181),'V)'], 'Interpreter', 'tex');
ylim([ymin_volt ymax_volt]);
set(gca,'box','off');


pos2=[0.1 0.425 0.7 0.1];
axr=subplot('Position', pos2);

%plotSpikeRaster(spike_times_trials_abs,'PlotType','vertline','XLimForCell',[tE(1) tE(end)]);
plotRaster({spikes_s},tE, 'k', 1,1);
h = gca; h.YAxis.Visible = 'off'; 

%title(sprintf('Raster plot for Run 14, Threshold at %i std',sd));

xlabel('Time (s)');
ylabel('Raster', 'Color', 'k');
set(gcf,'color','w');
set(gca,'box','off')
hold on;

pos3=[0.1 0.1 0.7 0.2];
axs=subplot('Position', pos3);
plot(tL_exp,laser_thr, 'Color', 'k');
xlabel('Time (s)');
ylabel('Stimulus');
set(gca, 'ytick', [0.05,4.9], 'yticklabel', {'Off', 'On'});
set(gcf,'color','w');
set(gca,'box','off')
%title('Optical Stimulus');
linkaxes([axb, axr,axs], 'x');

hold off;

% p=get(gca, 'Position');
% p(2)=p(2)*2;
% set(gca, 'Position',p);

% saveas(gcf, strcat(filename,sprintf('_%i_uV_%isd_record+stim+raster.pdf',round(thresh,0),sd)));
print('-painters',strcat(filename,sprintf('_%i_uV_%isd_record+stim+raster.pdf',round(thresh,0),sd)),'-dpdf');

%Only stimulation and activity

figure;
% rectangle('Position',[0 0 12.2 size(trials,1)],'FaceColor','#4DBEEE', 'EdgeColor', '#4DBEEE');
hold on;
% pos1=[0.1 0.65 0.7 0.25];
%axb=subplot(3,1,1);
% axb=subplot('Position', pos1);
axb=subplot(2,1,1);
plot(tE,base, 'k');
threshline = refline([0 thresh]);
threshline.Color = 'r';
%title('Spontaneous Activity');
xlabel('Time (s)');
ylabR=ylabel(['Potential (',char(181),'V)'], 'Interpreter', 'tex');

% ylim([ymin_volt/1.5 ymax_volt/1.5]);
set(gca,'box','off');


% pos2=[0.1 0.425 0.7 0.1];
% axs=subplot('Position', pos2);
axs=subplot(2,1,2);
plot(tL_exp,laser_thr, 'Color', stimulus_rect_facecolor);
xlabel('Time (s)');
ylabS=ylabel('Whisker Stimulus');
ylabS.Position= [-0.05*tE(end) 2.5 -1.000]
set(gca, 'ytick', [0.05,4.9], 'yticklabel', {'Off', 'On'});
set(gcf,'color','w');
set(gca,'box','off')
%title('Optical Stimulus');
linkaxes([axb,axs], 'x');

hold off;

% p=get(gca, 'Position');
% p(2)=p(2)*2;
% set(gca, 'Position',p);

% saveas(gcf, strcat(filename,sprintf('_%i_uV_%isd_record+stim+raster.pdf',round(thresh,0),sd)));
print('-painters',strcat(filename,sprintf('_%i_uV_%isd_record+stim.pdf',round(thresh,0),sd)),'-dpdf');

%%
Fontsize=20;
FontsizeAxis=16;
% small_selection of stim+recording
% start=150000;
% stop=230000; 
% for 4/11/2021 run, animal 2
start=150000;
stop=400000; 

figure;
% rectangle('Position',[0 0 12.2 size(trials,1)],'FaceColor','#4DBEEE', 'EdgeColor', '#4DBEEE');
hold on;
% pos1=[0.1 0.65 0.7 0.25];
%axb=subplot(3,1,1);
% axb=subplot('Position', pos1);
axb=subplot(2,1,1);
plot(tE(start:stop),base(start:stop), 'k');
    ax=ancestor(gca, 'axes');
    xrule=ax.XAxis;
    yrule=ax.YAxis;
    
width=5;
    height=8;
    pos = get(gcf, 'Position');
    set(gcf, 'Position', [0 0 width*100, height*100]); %<- Set size
threshline = refline([0 thresh]);
threshline.Color = 'r';
%title('Spontaneous Activity');
xL=xlabel('Time (s)');
ylabR1=ylabel(['Potential (',char(181),'V)'], 'Interpreter', 'tex');

    xrule.FontSize =FontsizeAxis;
    yrule.FontSize =FontsizeAxis;
    xL.FontSize=Fontsize;
    ylabR1.FontSize=Fontsize;
    
ylabR1.Position
ylim([ymin_volt/1.5 ymax_volt/1.5]);
set(gca,'box','off');


% pos2=[0.1 0.425 0.7 0.1];
% axs=subplot('Position', pos2);
axs=subplot(2,1,2);
plot(tL_exp(start:stop),laser_thr(start:stop), 'Color', stimulus_rect_facecolor);
    ax=ancestor(gca, 'axes');
    xrule=ax.XAxis;
    yrule=ax.YAxis;

xL=xlabel('Time (s)');
ylabS1=ylabel('Whisker Stimulus');
    xrule.FontSize =FontsizeAxis;
    yrule.FontSize =FontsizeAxis;
    xL.FontSize=Fontsize;
    ylabS1.FontSize=Fontsize;

set(gca, 'ytick', [0.05,4.9], 'yticklabel', {'Off', 'On'});
set(gcf,'color','w');
set(gca,'box','off')
%title('Optical Stimulus');
linkaxes([axb,axs], 'x');
xlim([7.5 20])
%ylabS1.Position=ylabR1.Position;
%ylabS1.Position= [6.6382 -24.9999 -1.000]
hold off;

% p=get(gca, 'Position');
% p(2)=p(2)*2;
% set(gca, 'Position',p);

% saveas(gcf, strcat(filename,sprintf('_%i_uV_%isd_record+stim+raster.pdf',round(thresh,0),sd)));
print('-painters',strcat(filename,sprintf('_%i_uV_%isd_record+stim_selection.pdf',round(thresh,0),sd)),'-dpdf');

%%=========================================================================================
%%===================================================================================================================
%%===================================================================================================================
%%===================================================================================================================
%% ========================== function definitions begin====================================================================================
%%===================================================================================================================
%%===================================================================================================================

function [coeff,score,latent, tsquared, explained, mu, spikes_PCA]=CalculatePCA(spikes_PCA_unnorm,sd,thresh,filename,pca_x_lim,pca_y_lim)

spikes_PCA=(spikes_PCA_unnorm-min(min(spikes_PCA_unnorm)))/(max(max(spikes_PCA_unnorm))-min(min(spikes_PCA_unnorm)));

[coeff,score,latent, tsquared, explained, mu]=pca(double(spikes_PCA'));

%generate scree plot
figure;
explained_sum=cumsum(explained);
hold on;
plot(1:length(explained), explained./100, 'Color', 'k')
plot(1:length(explained), explained_sum./100,'k--')
hold off;
xlabel('Principle component');
ylabel('Proportion of variance');
set(gca,'box','off');
set(gcf,'color','w')
width=4;
height=width;
pos = get(gcf, 'Position');
set(gcf, 'Position', [pos(1) pos(2) width*100, height*100]); %<- Set size
set(gca, 'FontSize', 16);
% saveas(gcf, strcat(filename,sprintf('_%i_uV_%isd_scree_plot.pdf',round(thresh,0),sd)));
print('-painters',strcat(filename,sprintf('_%i_uV_%isd_scree_plot.pdf',round(thresh,0),sd)),'-dpdf');

%plot principle componenets
figure;

PCx=1;
PCy=2;

Fontsize=30;
FontsizeAxis=22;
MarkerS=4; %marker size
plot(score(:,PCx),score(:,PCy),'.','MarkerSize',MarkerS, 'Color', 'k'); %blue
ax=ancestor(gca, 'axes');
xrule=ax.XAxis;
yrule=ax.YAxis;

xL=xlabel(sprintf('PC %i',PCx));
yL=ylabel(sprintf('PC %i',PCy));

xrule.FontSize =FontsizeAxis;
yrule.FontSize =FontsizeAxis;
xL.FontSize=Fontsize;
yL.FontSize=Fontsize;
width=4;
height=width;
pos = get(gcf, 'Position');
set(gcf, 'Position', [pos(1) pos(2) width*100, height*100]); %<- Set size
set(gca,'box','off');
set(gcf,'color','w')
xlim(pca_x_lim);
ylim(pca_y_lim);
end
%%
function [PCA_avg, PCA_std, PCA_spike_locations]=CalcandPlotClusWavef(PCA_time,spikes_PCA,idx,gm,num_clust,spikes_s,ColOrd,PCA_scores,pca_x_lim,pca_y_lim,ymin_volt, ymax_volt,PCx,PCy,filename,sd,thresh)

PCA_avg=zeros(size(spikes_PCA,1),num_clust);
PCA_max=zeros(size(spikes_PCA,1),num_clust); %std(spikes_PCA(:,PCA_thr),0,2);
PCA_min=zeros(size(spikes_PCA,1),num_clust);
PCA_std=zeros(size(spikes_PCA,1),num_clust);
PCA_spike_locations=cell(num_clust,1);

if num_clust>1
    for j=1:num_clust
        PCA_avg(:,j)=mean(spikes_PCA(:,(idx==j)),2);
        PCA_max(:,j)=max(spikes_PCA(:,(idx==j)),[],2); %std(spikes_PCA(:,PCA_thr),0,2);
        PCA_min(:,j)=min(spikes_PCA(:,(idx==j)),[],2); 
        PCA_std(:,j)=std(spikes_PCA(:,(idx==j)),[],2);
        PCA_spike_locations{j,1}=spikes_s((idx==j));
    end
else
    PCA_avg(:,1)=mean(spikes_PCA(:,:),2);
    PCA_max(:,1)=max(spikes_PCA(:,:),[],2); %std(spikes_PCA(:,PCA_thr),0,2);
    PCA_min(:,1)=min(spikes_PCA(:,:),[],2); 
    PCA_std(:,1)=std(spikes_PCA(:,:),[],2);
    PCA_spike_locations{1,1}=spikes_s(:);
end

figure;
subplot(2,1,1);
gscatter(PCA_scores(:,PCx),PCA_scores(:,PCy), idx,ColOrd, '.', 6,'off');
% width=4;
% height=1*width;
% pos = get(gcf, 'Position');
% set(gcf, 'Position', [pos(1) pos(2) width*10, height*10]); %<- Set size
hold on;
if isprop(gm, 'AIC')

    gmPDF = @(x,y)reshape(pdf(gm,[x(:) y(:)]),size(x));
    fcontour(gmPDF,[pca_x_lim pca_y_lim]) %pca_x_lim,pca_y_lim
else
    plot(gm(:,1),gm(:,2),'kx', 'MarkerSize',10,'LineWidth',3)
end
%plot(C(:,1),C(:,2),'kx', 'MarkerSize',10,'LineWidth',3)
gmPDF = @(x,y)reshape(pdf(gm,[x(:) y(:)]),size(x));
fcontour(gmPDF,[pca_x_lim pca_y_lim]) %pca_x_lim,pca_y_lim


set(gca, 'FontSize', 16);
xlabel(sprintf('PC %i',PCx));
ylabel(sprintf('PC %i',PCy));
xlim(pca_x_lim);
ylim(pca_y_lim);


set(gca,'box','off')

hold off;
subplot(2,1,2);
width=3;
height=width*2.5;

set(gcf, 'Position', [0 0 width*100, height*100]); %<- Set size
hold on;
for i=1:size(PCA_avg,2)
    plot(PCA_time,PCA_avg(:,i), 'LineWidth', 2, 'Color', ColOrd(i,:));
    shade(PCA_time,(PCA_avg(:,i)-PCA_std(:,i)),PCA_time,(PCA_avg(:,i)+PCA_std(:,i)),'FillType',[1 2;2 1], 'Color', ColOrd(i,:), 'LineStyle', 'none');
end

ylim([ymin_volt ymax_volt]);
set(gca, 'FontSize', 16);
xlabel('Time (ms)');
ylabel(['Potential (',char(181),'V)'], 'Interpreter', 'tex');
set(gca,'box','off')
set(gcf,'color','w')
hold off;
% saveas(gcf, strcat(filename,sprintf('_%i_uV_%isd_PCA+waveforms.pdf',round(thresh,0),sd)));
print('-painters',strcat(filename,sprintf('_%i_uV_%isd_PCA+waveforms.pdf',round(thresh,0),sd)),'-dpdf');

end
%% Calculate kmeans
function [kmean_mins,idx]=KmeansPlot(tot_clusters,PCA_scores,filename,thresh, sd)
avg_cluster_dist=zeros(tot_clusters, 1);

for i=1:tot_clusters
    %[idx, C, sumd,D]=kmeans(PCA_scores,i, 'MaxIter', 1000);
    [idx, C, sumd,D]=kmeans(PCA_scores,i, 'Replicates', 200);
    single_avg_dist=zeros(i,1);
    for j=1:size(C,1) 
        tmp=((PCA_scores((idx==j),:)-C(j,:)));
        single_avg_dist(j)=mean((tmp(:,1).^2+tmp(:,2).^2).^0.5);
    end
    avg_cluster_dist(i)=mean(single_avg_dist);
end
kmean_mins=islocalmin(avg_cluster_dist);
figure;
plot((1:tot_clusters),avg_cluster_dist);
width=4;
height=width;
%pos = get(gcf, 'Position');
set(gcf, 'Position', [100 100 width*100, height*100]); %<- Set size
xlabel('Number of Clusters');
ylabel('Avg Inter Cluster Dist.');
set(gca, 'FontSize', 16);
set(gca,'box','off')
set(gcf,'color','w')

% saveas(gcf, strcat(filename,sprintf('_%i_uV_%isd_kmeans_opt.pdf',round(thresh,0),sd)));
print('-painters',strcat(filename,sprintf('_%i_uV_%isd_kmeans_opt.pdf',round(thresh,0),sd)),'-dpdf');

end
%%
function [PCA_trial_spike_locations, PCA_epoch_spike_locations,...
        PCA_trial_spikes_colored, PCA_epoch_spikes_colored,idx_times_trials,idx_times_epochs]=ColorSpikes(spikes_raster, epoch_spike_times,...
        spike_times_trials_abs, epoch_spike_times_abs_all, spikes_s, i,idx)

    PCA_trial_spike_locations=cell(i,1);
    PCA_epoch_spike_locations=cell(i,1);

    PCA_trial_spikes_colored=cell(length(spikes_raster),i);
    PCA_epoch_spikes_colored=cell(length(epoch_spike_times),i);
    spike_times_during_trials=[spikes_raster{:}];
    [non_trial_spikes, non_trial_locs] = setdiff(spikes_s, spike_times_trials_abs{1,1} );
%    [non_trial_spikes, non_trial_locs] = setdiff(spikes_s, spike_times_trials_abs{1,1},'legacy' );
   
    idx_times_trials=idx;
    idx_times_trials(non_trial_locs)=[];
    for j=1:i
        PCA_trial_spike_locations{j,1}=spike_times_during_trials(1,(idx_times_trials==j));
    end
    for j=1:i
        for k=1:length(spikes_raster)
            %determine the number of spikes found in each trial
            tmp=size([spikes_raster{1:k-1,1}],2);
            %knowing how many spikes, shift the index by tmp for indicating which
            %cluster they belong to
            idx_tmp=idx_times_trials((1:size(spikes_raster{k,1},2))+tmp,1);
            tmp=spikes_raster{k,1};
            PCA_trial_spikes_colored{k,j}=tmp(idx_tmp==j);
        end
    end
    
%    spike_times_during_epochs=epoch_spike_times_all{1};
    spike_times_during_epochs=[epoch_spike_times{:}];
   [non_epoch_spikes, non_epoch_locs] = setdiff(spikes_s, epoch_spike_times_abs_all{1,1} );

    idx_times_epochs=idx;
    idx_times_epochs(non_epoch_locs)=[];
    for j=1:i
        PCA_epoch_spike_locations{j,1}=spike_times_during_epochs(1,(idx_times_epochs==j));
    end
    for j=1:i
        for k=1:length(epoch_spike_times)
            tmp=size([epoch_spike_times{1:k-1,1}],2);
            idx_tmp=idx_times_epochs((1:size(epoch_spike_times{k,1},2))+tmp,1);
            tmp=epoch_spike_times{k,1};
            PCA_epoch_spikes_colored{k,j}=tmp(idx_tmp==j);
        end
    end

end
%%
function []=PlotClusterWvefms(tE,PCA_time,PCA_avg,PCA_std,PCA_spike_locations,score,idx,C,ColOrd,pwd_folder,filename, thresh, sd,ymin_volt, ymax_volt,PCx,PCy,pca_x_lim,pca_y_lim, rem_cluster,tL,trial_time_ms,pulse_width,pulses_in_train,epoch_time,epochs,PCA_trial_spike_locations,PCA_epoch_spike_locations)
%scaled text for figure 3
ylabpos=-4;
Fontsize=27;
FontsizeAxis=20;
figure;
blue_rect_height=1;
height=4;
width=height*1;

%added for spritzer processing
clust_num=size(PCA_spike_locations,1);
ymin_volt=ymin_volt/1.6;
ymax_volt=-ymin_volt;

%pos = get(gcf, 'Position');
set(gcf, 'Position', [0 0 width*100, height*100]); %<- Set size
hold on;
for i=1:size(PCA_avg,2)
    if i~=rem_cluster
        h=plot(PCA_time,PCA_avg(:,i), 'LineWidth', 2, 'Color', ColOrd(i,:));
        shade(PCA_time,(PCA_avg(:,i)-PCA_std(:,i)),PCA_time,(PCA_avg(:,i)+PCA_std(:,i)),'FillType',[1 2;2 1], 'Color', ColOrd(i,:), 'LineStyle', 'none');

    end
end
ylim([ymin_volt ymax_volt]);
 
ax=ancestor(h, 'axes');
xrule=ax.XAxis;
yrule=ax.YAxis;

xL=xlabel('Time (ms)');
% xticks([PCA_time(1) round(PCA_time(end))]);
yL=ylabel(['Potential (',char(181),'V)'], 'Interpreter', 'tex');

xrule.FontSize =FontsizeAxis;
yrule.FontSize =FontsizeAxis;
xL.FontSize=Fontsize;
yL.FontSize=Fontsize;

% set(gca, 'FontSize', 20);
set(gca,'box','off')
set(gcf,'color','w')
% saveas(gcf, strcat(filename,sprintf('_%i_uV_%isd_waveforms.pdf',round(thresh,0),sd)));
print('-painters',strcat(filename,sprintf('_%i_uV_%isd_%i_clust_waveforms_scaled.pdf',round(thresh,0),sd,clust_num)),'-dpdf');

hold off;

Fontsize=30;
FontsizeAxis=22;

figure;
blue_rect_height=1;
height=4;
width=height*1;
%pos = get(gcf, 'Position');
set(gcf, 'Position', [0 0 width*100, height*100]); %<- Set size
hold on;
for i=1:size(PCA_avg,2)
    if i~=rem_cluster
        h=plot(PCA_time,PCA_avg(:,i), 'LineWidth', 2, 'Color', ColOrd(i,:));
        shade(PCA_time,(PCA_avg(:,i)-PCA_std(:,i)),PCA_time,(PCA_avg(:,i)+PCA_std(:,i)),'FillType',[1 2;2 1], 'Color', ColOrd(i,:), 'LineStyle', 'none');
    end
end
ylim([ymin_volt ymax_volt]);

ax=ancestor(h, 'axes');
xrule=ax.XAxis;
yrule=ax.YAxis;

xL=xlabel('Time (ms)');
% xticks([PCA_time(1) round(PCA_time(end))]);
yL=ylabel(['Potential (',char(181),'V)'], 'Interpreter', 'tex');

xrule.FontSize =FontsizeAxis;
yrule.FontSize =FontsizeAxis;
xL.FontSize=Fontsize;
yL.FontSize=Fontsize;

% set(gca, 'FontSize', 20);
set(gca,'box','off')
set(gcf,'color','w')
% saveas(gcf, strcat(filename,sprintf('_%i_uV_%isd_waveforms.pdf',round(thresh,0),sd)));
print('-painters',strcat(filename,sprintf('_%i_uV_%isd_%i_clust_waveforms.pdf',round(thresh,0),sd,clust_num)),'-dpdf');

hold off;



%plot waveforms without shaded std
figure;
blue_rect_height=1;
width=4;
height=width;
%pos = get(gcf, 'Position');
set(gcf, 'Position', [0 0 width*100, height*100]); %<- Set size
hold on;
for i=1:size(PCA_avg,2)
    if i~=rem_cluster
    h=plot(PCA_time,PCA_avg(:,i), 'LineWidth', 2, 'Color', ColOrd(i,:));
%     shade(PCA_time,(PCA_avg(:,i)-PCA_std(:,i)),PCA_time,(PCA_avg(:,i)+PCA_std(:,i)),'FillType',[1 2;2 1], 'Color', ColOrd(i,:), 'LineStyle', 'none');
    end
end
ax=ancestor(h, 'axes');
xrule=ax.XAxis;
yrule=ax.YAxis;
xL=xlabel('Time (ms)');
yL=ylabel(['Potential (',char(181),'V)'], 'Interpreter', 'tex');

xrule.FontSize =FontsizeAxis;
yrule.FontSize =FontsizeAxis;
xL.FontSize=Fontsize;
yL.FontSize=Fontsize;

ylim([ymin_volt ymax_volt]);
% set(gca, 'FontSize', Fontsize);

set(gca,'box','off')
set(gcf,'color','w')
% saveas(gcf, strcat(filename,sprintf('_%i_uV_%isd_nonshaded_waveforms.pdf',round(thresh,0),sd)));
print('-painters',strcat(filename,sprintf('_%i_uV_%isd_%i_clust_nonshaded_waveforms.pdf',round(thresh,0),sd,clust_num)),'-dpdf');

hold off;
%plot each waveform by itself

for i=1:size(PCA_avg,2)
    figure;
    h=plot(PCA_time,PCA_avg(:,i), 'LineWidth', 2, 'Color', ColOrd(i,:));
    hold on;
    shade(PCA_time,(PCA_avg(:,i)-PCA_std(:,i)),PCA_time,(PCA_avg(:,i)+PCA_std(:,i)),'FillType',[1 2;2 1], 'Color', ColOrd(i,:), 'LineStyle', 'none');
    hold off;
    ylim([ymin_volt ymax_volt]);
    width=5;
    height=width;
    %pos = get(gcf, 'Position');
    set(gcf, 'Position', [50 50 width*50, height*100]); %<- Set size
    set(gca, 'FontSize', Fontsize);
    ax=ancestor(h, 'axes');
    xrule=ax.XAxis;
    yrule=ax.YAxis;
    xL=xlabel('Time (ms)');
    yL=ylabel(['Potential (',char(181),'V)'], 'Interpreter', 'tex');

    xrule.FontSize =FontsizeAxis;
    yrule.FontSize =FontsizeAxis;
    xL.FontSize=Fontsize;
    yL.FontSize=Fontsize;

    set(gca,'box','off')
    set(gcf,'color','w')
%     saveas(gcf, strcat(filename,sprintf('_%i_uV_%isd_cluster_%i_waveform.pdf',round(thresh,0),sd,i)));
    print('-painters',strcat(filename,sprintf('_%i_uV_%isd_cluster_%i_waveform.pdf',round(thresh,0),sd,i)),'-dpdf');
end

% plot PCA plot by itself
figure;
gscatter(score(:,PCx),score(:,PCy), idx,ColOrd, '.', 6,'off');
hold on;
if isprop(C, 'AIC')

    gmPDF = @(x,y)reshape(pdf(C,[x(:) y(:)]),size(x));
    fcontour(gmPDF,[pca_x_lim pca_y_lim]) %pca_x_lim,pca_y_lim
else
    plot(C(:,1),C(:,2),'kx', 'MarkerSize',10,'LineWidth',3);
end
set(gca, 'FontSize', Fontsize);

ax=ancestor(gca, 'axes');
xrule=ax.XAxis;
yrule=ax.YAxis;
xL=xlabel(sprintf('PC %i',PCx));
yL=ylabel(sprintf('PC %i',PCy));

xrule.FontSize =FontsizeAxis;
yrule.FontSize =FontsizeAxis;
xL.FontSize=Fontsize;
yL.FontSize=Fontsize;

xlim(pca_x_lim);
ylim(pca_y_lim);
set(gca,'box','off')
set(gcf,'color','w')
width=4;
height=width;
%pos = get(gcf, 'Position');
set(gcf, 'Position', [50 50 width*100, height*100]);
hold off;

% saveas(gcf, strcat(filename,sprintf('_%i_uV_%isd_PCA.pdf',round(thresh,0),sd)));
print('-painters',strcat(filename,sprintf('_%i_uV_%isd_%i_clust_PCA.pdf',round(thresh,0),sd,clust_num)),'-dpdf');

%plot PCA with scaled text for figure 3
Fontsize=27;
FontsizeAxis=20;
figure;
gscatter(score(:,PCx),score(:,PCy), idx,ColOrd, '.', 6,'off');
hold on;
if isprop(C, 'AIC')

    gmPDF = @(x,y)reshape(pdf(C,[x(:) y(:)]),size(x));
    fcontour(gmPDF,[pca_x_lim pca_y_lim]) %pca_x_lim,pca_y_lim
else
    plot(C(:,1),C(:,2),'kx', 'MarkerSize',10,'LineWidth',3);
end
set(gca, 'FontSize', Fontsize);

ax=ancestor(gca, 'axes');
xrule=ax.XAxis;
yrule=ax.YAxis;
xL=xlabel(sprintf('PC %i',PCx));
yL=ylabel(sprintf('PC %i',PCy));

xrule.FontSize =FontsizeAxis;
yrule.FontSize =FontsizeAxis;
xL.FontSize=Fontsize;
yL.FontSize=Fontsize;

xlim(pca_x_lim);
ylim(pca_y_lim);
set(gca,'box','off')
set(gcf,'color','w')
width=4;
height=width;
%pos = get(gcf, 'Position');
set(gcf, 'Position', [50 50 width*100, height*100]);
hold off;

% saveas(gcf, strcat(filename,sprintf('_%i_uV_%isd_PCA.pdf',round(thresh,0),sd)));
print('-painters',strcat(filename,sprintf('_%i_uV_%isd_%i_clust_PCA_scaled.pdf',round(thresh,0),sd,clust_num)),'-dpdf');
    
end
%%
% Bayesian Adaptive Kernel Smoother (BAKS)
% BAKS is a method for estimating firing rate from spike train data that uses kernel smoothing technique 
% with adaptive bandwidth determined using a Bayesian approach
% ---------------INPUT---------------
% - SpikeTimes : spike event times [nSpikes x 1]
% - Time : time at which the firing rate is estimated [nTime x 1]
% - a : shape parameter (alpha) 
% - b : scale paramter (beta)
% ---------------INPUT---------------
% - h : adaptive bandwidth [nTime x 1]
% - FiringRate : estimated firing rate [nTime x 1]
% More information, please refer to "Estimation of neuronal firing rate using Bayesian adaptive kernel smoother (BAKS)"

function [FiringRate, h] = BAKS(SpikeTimes,Time,a,b)
N = length(SpikeTimes);
sumnum = 0; sumdenum = 0;
for i=1:N
    numerator = (((Time-SpikeTimes(i)).^2)./2 + 1./b).^(-a);
    denumerator = (((Time-SpikeTimes(i)).^2)./2 + 1./b).^(-a-0.5);
    sumnum = sumnum + numerator;
    sumdenum = sumdenum + denumerator;
end
h = (gamma(a)/gamma(a+0.5)).*(sumnum./sumdenum);

FiringRate = zeros(length(Time),1);
for j=1:N
    K = (1./(sqrt(2.*pi).*h)).*exp(-((Time-SpikeTimes(j)).^2)./(2.*h.^2));
    FiringRate = FiringRate + K;
end
end
%% plot colored rasters
function []=plotColoredRaster(PCA_epoch_spikes_colored,epoch_time,PCA_trial_spikes_colored, trial_time_ms,...
        thresh,sd, pulse_width, pre_trial_time,post_trial_time, pulses_in_train, ColOrd,stimulus_rect_facecolor, pwd_folder,filename)
    
    Fontsize=30;
    FontsizeAxis=22;
    ylabpos=-1;
    figure;
    rectangle('Position',[0 0.5 pulse_width size(PCA_trial_spikes_colored,1)+10],'FaceColor', stimulus_rect_facecolor, 'EdgeColor',  stimulus_rect_facecolor);
    % pgon=polyshape([0 0 pulse_width pulse_width],[0.5 0.5+num_trials+10 0.5+num_trials+10 0],'FaceColor','#4DBEEE','EdgeColor', '#4DBEEE')
    hold on;
    for j=1:size(PCA_trial_spikes_colored,2)
        plotRaster(PCA_trial_spikes_colored(:,j), trial_time_ms, ColOrd(j,:), 1,1);
    end
    
    ax=ancestor(gca, 'axes');
    xrule=ax.XAxis;
    yrule=ax.YAxis;
%     xL=xlabel('Time after pulse (ms)');
    xL=xlabel('Pulse time (ms)');
    yL=ylabel('Pulse number');

    xrule.FontSize =FontsizeAxis;
    yrule.FontSize =FontsizeAxis;
    xL.FontSize=Fontsize;
    yL.FontSize=Fontsize;

    width=4;
    height=width;
    pos = get(gcf, 'Position');
    set(gcf, 'Position', [pos(1) pos(2) width*100, height*100]); %<- Set size
%     set(gca, 'FontSize', Fontsize);
    set(gcf,'color','w');
    set(gca,'box','off')
    
%      if strcmp(pwd_folder,'Frequncy varied runs\')
%     xticks([-10 -5 0 5 10 15])
%     xticklabels({'-10','-5','0','5','10','15'})
%     else
%     xticks([-10 0 10 20 30])
%     xticklabels({'-10','0','10','20','30'})
%     end
    hold off;
%     xlim([-pre_trial_time*10^3 post_trial_time*10^3]);
%     ylim([0 500]);
%     yticks([ 0 100 200 300 400 500])
%     yticklabels({'0','100','200','300','400','500'})
%         ylim([0 size(PCA_trial_spikes_colored,1)]);

%     saveas(gcf, strcat(filename,sprintf('_%i_uV_%isd_raster_colored.pdf',round(thresh,0),sd)));
    print('-painters',strcat(filename,sprintf('_%i_uV_%isd_raster_colored.pdf',round(thresh,0),sd)),'-dpdf');

    % plot colored raster of pulse train
    figure;
    hold on;
    for i=1:length(pulses_in_train)
        rectangle('Position',[pulses_in_train(i,1) 0 (pulses_in_train(i,2)-pulses_in_train(i,1)) size(epoch_time,1)+1],'FaceColor', stimulus_rect_facecolor, 'EdgeColor',  stimulus_rect_facecolor);
        %rectangle('Position',[(rising_ind(i)-rising_ind(1))*r_tL 0 (falling_ind(i)-rising_ind(i))*r_tL blue_rect_height],'FaceColor','#4DBEEE', 'EdgeColor', '#4DBEEE');
    end
    for j=1:size(PCA_epoch_spikes_colored,2)
        plotRaster(PCA_epoch_spikes_colored(:,j), epoch_time, ColOrd(j,:), 0,1);
    end
 
    ax=ancestor(gca, 'axes');
    xrule=ax.XAxis;
    yrule=ax.YAxis;
    xL=xlabel('Pulse train time (s)');
    yL=ylabel('Pulse train');

    xrule.FontSize =FontsizeAxis;
    yrule.FontSize =FontsizeAxis;
    xL.FontSize=Fontsize;
    yL.FontSize=Fontsize;

    
    width=4;
    height=width;
    pos = get(gcf, 'Position');
    set(gcf, 'Position', [pos(1) pos(2) width*100, height*100]); %<- Set size
%     set(gca, 'FontSize', Fontsize);
    set(gcf,'color','w');
    set(gca,'box','off')
%     ylim([0 25]);
%     xlim([-0.2 1.2]);
%     xticks([0 0.4 0.8 1.2])
%     xticklabels({'0','0.4','0.8','1.2'})
%     saveas(gcf, strcat(filename,sprintf('_%i_uV_%isd_epochs_raster_colored.pdf',round(thresh,0),sd)));
    print('-painters',strcat(filename,sprintf('_%i_uV_%isd_epochs_raster_colored.pdf',round(thresh,0),sd)),'-dpdf');

% make scaled versions for fig 3 of paper    
     Fontsize=22;
    FontsizeAxis=16;
    
    figure;
    rectangle('Position',[0 0.5 pulse_width size(PCA_trial_spikes_colored,1)+10],'FaceColor', stimulus_rect_facecolor, 'EdgeColor',  stimulus_rect_facecolor);
    % pgon=polyshape([0 0 pulse_width pulse_width],[0.5 0.5+num_trials+10 0.5+num_trials+10 0],'FaceColor','#4DBEEE','EdgeColor', '#4DBEEE')
    hold on;
    for j=1:size(PCA_trial_spikes_colored,2)
        plotRaster(PCA_trial_spikes_colored(:,j), trial_time_ms, ColOrd(j,:), 1,1);
    end
    
    ax=ancestor(gca, 'axes');
    xrule=ax.XAxis;
    yrule=ax.YAxis;
%     xL=xlabel('Time after pulse (ms)');
    xL=xlabel('Pulse time (ms)');
    yL=ylabel('Pulse number');

    xrule.FontSize =FontsizeAxis;
    yrule.FontSize =FontsizeAxis;
    xL.FontSize=Fontsize;
    yL.FontSize=Fontsize;
    
    width=4;
    height=width;
    pos = get(gcf, 'Position');
    set(gcf, 'Position', [pos(1) pos(2) width*100, height*100]); %<- Set size
%     set(gca, 'FontSize', Fontsize);
    set(gcf,'color','w');
    set(gca,'box','off')
    
%     if strcmp(pwd_folder,'Frequncy varied runs\')
%     xticks([-10 -5 0 5 10 15])
%     xticklabels({'-10','-5','0','5','10','15'})
%     else
%     xticks([-10 0 10 20 30])
%     xticklabels({'-10','0','10','20','30'})
%     end
%     hold off;
%     xlim([-pre_trial_time*10^3 post_trial_time*10^3]);
%     ylim([0 500]);
%     yticks([ 0 100 200 300 400 500])
%     yticklabels({'0','100','200','300','400','500'})
%         ylim([0 size(PCA_trial_spikes_colored,1)]);

%     saveas(gcf, strcat(filename,sprintf('_%i_uV_%isd_raster_colored.pdf',round(thresh,0),sd)));
    print('-painters',strcat(filename,sprintf('_%i_uV_%isd_raster_colored_scaled.pdf',round(thresh,0),sd)),'-dpdf');

    % plot colored raster of pulse train
    figure;
    hold on;
    for i=1:length(pulses_in_train)
        rectangle('Position',[pulses_in_train(i,1) 0 (pulses_in_train(i,2)-pulses_in_train(i,1)) size(epoch_time,1)+1],'FaceColor', stimulus_rect_facecolor, 'EdgeColor',  stimulus_rect_facecolor);
        %rectangle('Position',[(rising_ind(i)-rising_ind(1))*r_tL 0 (falling_ind(i)-rising_ind(i))*r_tL blue_rect_height],'FaceColor','#4DBEEE', 'EdgeColor', '#4DBEEE');
    end
    for j=1:size(PCA_epoch_spikes_colored,2)
        plotRaster(PCA_epoch_spikes_colored(:,j), epoch_time, ColOrd(j,:), 0,1);
    end
   
    ax=ancestor(gca, 'axes');
    xrule=ax.XAxis;
    yrule=ax.YAxis;
    xL=xlabel('Pulse train time (s)');
    yL=ylabel('Pulse train');

    xrule.FontSize =FontsizeAxis;
    yrule.FontSize =FontsizeAxis;
    xL.FontSize=Fontsize;
    yL.FontSize=Fontsize;
   
    width=4;
    height=width;
    pos = get(gcf, 'Position');
    set(gcf, 'Position', [pos(1) pos(2) width*100, height*100]); %<- Set size
%     set(gca, 'FontSize', Fontsize);
    set(gcf,'color','w');
    set(gca,'box','off')
%     ylim([0 25]);
%     xlim([-0.2 1.2]);
%     xticks([0 0.4 0.8 1.2])
%     xticklabels({'0','0.4','0.8','1.2'})
%     saveas(gcf, strcat(filename,sprintf('_%i_uV_%isd_epochs_raster_colored.pdf',round(thresh,0),sd)));
    print('-painters',strcat(filename,sprintf('_%i_uV_%isd_epochs_raster_colored_scaled.pdf',round(thresh,0),sd)),'-dpdf');   
end
%% plot raster
function [] = plotRaster( spikeMat , tVec, color, offset,linewidth )
% Visualize raster plot
hold all;
if offset~=1
    offset=0;
end

for trialCount = 1:size(spikeMat,1)
    spikePos = spikeMat{trialCount};
    for spikeCount = 1:length(spikePos)
        plot ([ spikePos(spikeCount) spikePos(spikeCount)] , [trialCount-offset-0.4 trialCount+offset+0.4], 'Color', color, 'LineWidth', linewidth) ;
%         plot ([ spikePos(spikeCount) spikePos(spikeCount)] , [trialCount-0.4 trialCount+0.4], 'Color', color) ;
    end
end

set(gca,'box','off')
set(gcf,'color','w')
ylim([0 size(spikeMat,1)+1]);
xlim([tVec(1) round(tVec(end),1)]);

end
%% calc ball enc
function [vel, sampleTime, vel_res, frames_res]= locomotion_stab_2ch(a_data,encoding_freq,filtwin,a_frames,fps,nFrames)

% analyze locomotor activity
% translate V changes into speed; downsampling, filtering and resampling of
% motion data; the resulting vector has the same length as number of recorded frames
%
% INPUT
% input a_data =  1 channel or 2 channels recorded (2nd channel used to indentify motion direction)
% encoding_freq = 10000; % 10000 recommended to avoid undersampling
% fps =  fps of imaging data, if locomotion should be synched, otherwise fps=[]
% a_frames
%
% OUTPUT
% f_loc: figure handle for running speed(mm/s) vs. sampling time (s)
% f_loc_res: resampled locomotion data running speed(mm/s) vs. imaging frames
% vel: velocity, downsampled and filtered 
% sampleTime: corresponding experimental time (s) for velocity
% vel_res: downsampled and filtered, resampled velocity; corresponding to imaging frames
% frames_res: corresponding frames for resampled velocity
 sampleTime=0;
 vel_res=0;
 frames_res=0;

% Parameter
circumference = 373;
cycles_rev = 720;


% remove more artifacts
digit_thr=0; %0.5*(max(a_data(:))+min(a_data(:)));
a_data_dig=a_data;
a_data_dig(a_data_dig>digit_thr)=1;
a_data_dig(a_data_dig<digit_thr)=0;
for a_ch=1:size(a_data,2)
    ind=[strfind(a_data_dig(:,a_ch)', [1 0 1]) strfind(a_data_dig(:,a_ch)', [0 1 0])];
    for ind1=1:length(ind)
        a_data(ind(ind1)+1,a_ch)=a_data(ind(ind1),a_ch);
    end
end


% Find indices the voltage jumps
% Jump criteria = 0.5 of the max/min difference of the voltages
dV = 0.5*(max(a_data(:))-min(a_data(:)));

% one channel, no direction
if size(a_data,2) == 1 
    jump_indices = find( abs(diff(a_data(:))) > dV);
    % Vchanges = 0 except at indices where there is a voltage change
    Vchanges = zeros(size(a_data(:)));
    Vchanges(jump_indices) = 1;
    posvel=find(Vchanges==1)';
    negvel=find(Vchanges==-1)';
% two channels, speed and direction (QEI - quadratur encoder interface)
else 
    jump_indices_rise = find(diff(a_data(:,1)) > dV);
    jump_indices_fall = find(diff(a_data(:,1)) < -dV);
    
    ch2_1 = find (a_data(jump_indices_fall,2)< 0);
    ch2_2 = find (a_data(jump_indices_rise,2)>= 0);
    ch2_3 = find (a_data(jump_indices_fall,2)>= 0);
    ch2_4 = find (a_data(jump_indices_rise,2)< 0);
    
    Vchanges=zeros(1, size(a_data,1));
    Vchanges(jump_indices_fall(ch2_1))=-1; 
    Vchanges(jump_indices_rise(ch2_2))=-1; 
    Vchanges(jump_indices_fall(ch2_3))=1; 
    Vchanges(jump_indices_rise(ch2_4))=1;
    posvel=find(Vchanges==1);
    negvel=find(Vchanges==-1);
end

% calculate speed
if ~isempty(posvel) || ~isempty(negvel)
    posvel(2,:)=ones(1, length( posvel));
    negvel(2,:)=ones(1, length( negvel)) * (-1);
    allvel=sortrows([posvel negvel]')';
    allvel(2,:)=[(allvel(2,2:end)+allvel(2,1:end-1))/2 0]; % set 1 to -1 and -1 to 1 changes to 0
    int0=[diff(allvel(1,:)) 0];
    intall=[zeros(1,allvel(1,1)) repelem( allvel(2,:).*int0 ,  int0) zeros(1,length(Vchanges)-allvel(1,end)) ];
    vel =  (circumference / (2*cycles_rev)) ./ (intall / 1000); % [mm/sec]
    vel(isinf(vel))=0;
else
    vel=zeros(size(Vchanges));
end

% filter speed
% weights = [1/(2*filtwin);(1/filtwin*ones(filtwin-1,1));1/(2*filtwin)];
% % Moving average is given by the convolution of the weights with the data
% vel = conv(vel,weights,'same'); % same

% % test of functioning figure
% figure('Position', [187 409 1391 499]); 
% hold on;
% plot(Vchanges*10,'b'); 
% plot(vel,'r'); 
% plot(vel_filt,'k');
% plot(vel,'g');  
% hold off


% % Convert the times (s) to frame numbers
% if isempty(a_frames)
%     a_frames=transpose(0:1/encoding_freq:(length(a_data)-1)/encoding_freq);
% end
% if ~isempty(fps)
%     frames = fps*a_frames(:,1);
%     sampleTime=[1:size(a_data,1)] / encoding_freq;
% else
%     sampleTime=[];
% end
% 
% % if locomotion should be synchronized to imaging data and recording done in M1
% if ~isempty(fps) && ~isempty(a_frames)
%     % Resample the locomotion
%     % frames_l and vel_l do not have the size of the number of actual frames in tiff file currently
% 
%     frames_res = [ceil(frames(1))]; vel_res = [vel(1)];
%     for i=[2:size(frames,1)]
%         if ceil(frames(i)) ~= ceil(frames(i-1)) % close to a unit frame number
%             frames_res = [frames_res, ceil(frames(i))];
%             vel_res = [vel_res, vel(i)];
%         end
%     end
%     % Pad the missing information with zeros
%     frames_res = [(1:floor(frames(1))),frames_res];
%     mm = frames_res(end);
%     frames_res = [frames_res, (mm+1:nFrames)];
%     vel_res = [zeros(1,floor(frames(1))),vel_res];
%     vel_res = [vel_res, zeros(1,nFrames-(mm))];
%     
%     frdiff=size(vel_res,2) - nFrames;
%     if frdiff>0
%         % remove one/two elements, maybe because some frames were skipped during recording ??
%         frames_res = frames_res(1:(end-frdiff));
%         vel_res = vel_res(1:(end-frdiff));
%         disp(strcat('Frame difference=',num2str(frdiff)));
%     end
%     
%     % plot
% %     f_loc=figure ('Position',[440   568   703   230]);
% %     plot(frames_res,vel_res);
% %     xlabel('Frames');
% %     ylabel('Running speed (mm/s)');
% %     xlim([frames_res(1),frames_res(end)]);
% %     ylim([min(vel)-2,max(vel)+2]);
% %     %ylim([0,120]);
% %     pbaspect([3,1,1]);
%     
% else
%     frames_res=[];
%     vel_res=[];
%     
% %     % plot
% %     sampleTime=sampleTime/60; % in min
% %     f_loc=figure ('Position',[440   568   703   230]);
% %     plot(sampleTime,vel);
% %     xlabel('Sample Time (min)');
% %     ylabel('Running speed (mm/s)');
% %     xlim([sampleTime(1),sampleTime(end)]);
% %     ylim([min(vel)-2,max(vel)+2]);
% %     %ylim([0,120]);
% %     pbaspect([3,1,1]);
% end

end
