%power coupling calculation
% data imported from .mat is coupling average and std at each wavelength in form of [probe1 probe2 probe3]
clear;
load('optial_coupling.mat')
load('optical_simulation_data.mat')

smf_smf_disp=smf_smf_misalign(:,1); %displacement in microns
smf_lfb_disp=smf_lfb_misalign(:,1);%displacement in microns
smf_lfb_wavelengths=smf_lfb_coupling_spectrum(:,1);

perfect_alignment=10*log10(smf_lfb_coupling_spectrum(:,3));
smf_lfb_misalign(:,2:end);

% %plot loss with fixed smf_lfb_misalign(6,:) (500nm misalignment b/w SMF and microfiber)
% loss_misalign=10*log10(smf_lfb_misalign(6,2:end).*smf_smf_misalign([6, 11, 13, 15, 16, 17],2:end));

%plot loss with fixed smf_lfb_misalign(1,:) (0nm misalignment b/w SMF and microfiber)
loss_misalign=10*log10(smf_lfb_misalign(1,2:end).*smf_smf_misalign([6, 11, 13, 15, 16, 17],2:end));

coupling_473_dB=10*log10(coupling_473_max);
coupling_544_dB=10*log10(coupling_544_max);
coupling_600_dB=10*log10(coupling_600_max);

coupling_473_std_of_max=std(coupling_473_dB);
coupling_544_std_of_max=std(coupling_544_dB);
coupling_600_std_of_max=std(coupling_600_dB);

wavelengths_measured=[473 544 600];
wavelengths_simulated=[400 450 473 500 543 550 600 650 700];
all_coupling=[coupling_473_dB' coupling_544_dB' coupling_600_dB'];
all_errors=[coupling_473_std_of_max' coupling_544_std_of_max' coupling_600_std_of_max'];

hold on;

errorbar(wavelengths_measured, all_coupling(1,:), all_errors(1,:), 'o', 'MarkerSize', 10, 'LineWidth',1.5);
plot(smf_lfb_wavelengths, perfect_alignment, '-k', 'Linewidth',1.5);
plot(wavelengths_simulated,loss_misalign, '-k', 'Linewidth',1.5);

hold off;
xlabel('Wavelength (nm)');
ylabel('Coupling (dB)');
width=5;
height=width;
pos = get(gcf, 'Position');
set(gcf, 'Position', [50 50 width*100, height*100]); %<- Set size

xlim([400 700]);
%xticks([-0.4 0 0.4 0.8 1.2])
%xticklabels({'-0.4','0','0.4','0.8','1.2'})
ylim([-10 0]);
set(gca, 'FontSize', 24);
set(gcf,'color','w');
set(gca,'box','off')
set(0,'DefaultAxesTitleFontWeight','normal');
hold off;

print('-painters','Coupling_efficiency_no_smf-lfb_misalign','-dpdf');

%% power coupling calculation, 500nm misalignment between 
% data imported from .mat is coupling average and std at each wavelength in form of [probe1 probe2 probe3]
clear;
load('optial_coupling.mat')
load('optical_simulation_data.mat')

smf_smf_disp=smf_smf_misalign(:,1); %displacement in microns
smf_lfb_disp=smf_lfb_misalign(:,1);%displacement in microns
smf_lfb_wavelengths=smf_lfb_coupling_spectrum(:,1);

perfect_alignment=10*log10(smf_lfb_coupling_spectrum(:,3));
smf_lfb_misalign(:,2:end);

% %plot loss with fixed smf_lfb_misalign(6,:) (500nm misalignment b/w SMF and microfiber)
% loss_misalign=10*log10(smf_lfb_misalign(6,2:end).*smf_smf_misalign([6, 11, 13, 15, 16, 17],2:end));

%plot loss with fixed smf_lfb_misalign(1,:) (0nm misalignment b/w SMF and microfiber)
loss_misalign=10*log10(smf_lfb_misalign(6,2:end).*smf_smf_misalign([6, 11, 13, 15, 16, 17],2:end));

coupling_473_dB=10*log10(coupling_473_max);
coupling_544_dB=10*log10(coupling_544_max);
coupling_600_dB=10*log10(coupling_600_max);

coupling_473_std_of_max=std(coupling_473_dB);
coupling_544_std_of_max=std(coupling_544_dB);
coupling_600_std_of_max=std(coupling_600_dB);

wavelengths_measured=[473 544 600];
wavelengths_simulated=[400 450 473 500 543 550 600 650 700];
all_coupling=[coupling_473_dB' coupling_544_dB' coupling_600_dB'];
all_errors=[coupling_473_std_of_max' coupling_544_std_of_max' coupling_600_std_of_max'];

hold on;

errorbar(wavelengths_measured, all_coupling(1,:), all_errors(1,:), 'o', 'MarkerSize', 10, 'LineWidth',1.5);
plot(smf_lfb_wavelengths, perfect_alignment, '-k', 'Linewidth',1.5);
plot(wavelengths_simulated,loss_misalign, '-k', 'Linewidth',1.5);

hold off;
xlabel('Wavelength (nm)');
ylabel('Coupling (dB)');
width=5;
height=width;
pos = get(gcf, 'Position');
set(gcf, 'Position', [50 50 width*100, height*100]); %<- Set size

xlim([400 700]);
%xticks([-0.4 0 0.4 0.8 1.2])
%xticklabels({'-0.4','0','0.4','0.8','1.2'})
ylim([-10 0]);
set(gca, 'FontSize', 24);
set(gcf,'color','w');
set(gca,'box','off')
set(0,'DefaultAxesTitleFontWeight','normal');
hold off;

print('-painters','Coupling_efficiency_500nm_smf-lfb_misalign','-dpdf');

%% Plot EIS data
clear; close all;
load('EIS_PBS_no_PEDOT.mat');
freq_before_pedot=all_data(:,9);
imp_before_pedot=all_data(:,10)*10^-6;
FontSize=20;
FontSize_zoomed=26;

load('eis_data_IrOx_probes.mat');

avg_imp= flip(mean(impedances,2)*10^-6);

std_imp=flip(std(impedances,0 ,2)*10^-6);
high_pass=1:length(freq);
log_freq=log10(flip(freq));
figure;
hold on;
plot(log_freq(high_pass),avg_imp(high_pass), 'LineWidth', 2, 'Color', 'k');
width=5;
height=width;
pos = get(gcf, 'Position');
set(gcf, 'Position', [0 0 width*100, height*100]); %<- Set size
shade(log_freq(high_pass)', (avg_imp(high_pass)+std_imp(high_pass))',log_freq(high_pass)', (avg_imp(high_pass)-std_imp(high_pass))','FillType',[1 2;2 1], 'Color', [0.5 0.5 0.5], 'LineStyle', 'none');%


xticks([1 2 3 4 5 6])
xticklabels({'10^{1}','10^{2}','10^{3}','10^{4}','10^{5}','10^{6}'});


xL=xlabel('Frequency (Hz)');
yL=ylabel('Impedance (M\Omega)');
set(gca, 'FontSize',FontSize_zoomed);
xL.FontSize=20;
yL.FontSize=20;
set(gca,'box','off')
set(gcf,'color','w')
hold off;

saveas(gcf, 'impedances_after_PEDOT_IrOx.svg');

figure;
hold on;
plot(log_freq(high_pass),avg_imp(high_pass), 'LineWidth', 2, 'Color', 'k');

plot(log10(freq_before_pedot),imp_before_pedot, 'LineWidth', 2, 'Color', 'g');

width=5;
height=width;
pos = get(gcf, 'Position');
set(gcf, 'Position', [0 0 width*100, height*100]); %<- Set size
shade(log_freq(high_pass)', (avg_imp(high_pass)+std_imp(high_pass))',log_freq(high_pass)', (avg_imp(high_pass)-std_imp(high_pass))','FillType',[1 2;2 1], 'Color', [0.5 0.5 0.5], 'LineStyle', 'none');%


xticks([1 2 3 4 5 6])
xticklabels({'10^{1}','10^{2}','10^{3}','10^{4}','10^{5}','10^{6}'})

xL=xlabel('Frequency (Hz)');
yL=ylabel('Impedance (M\Omega)');
xL.FontSize=20;
yL.FontSize=20;

legend('After PEDOT:PSS','Before PEDOT:PSS');
set(gca,'box','off')
set(gcf,'color','w')
hold off;
set(gca, 'FontSize', FontSize);
saveas(gcf, 'impedances_without_and_with_PEDOT_IrOx_legend.svg');

figure;
hold on;
plot(log_freq(high_pass),avg_imp(high_pass), 'LineWidth', 2, 'Color', 'k');

plot(log10(freq_before_pedot),imp_before_pedot, 'LineWidth', 2, 'Color', 'g');

width=5;
height=width;
pos = get(gcf, 'Position');
set(gcf, 'Position', [0 0 width*100, height*100]); %<- Set size
shade(log_freq(high_pass)', (avg_imp(high_pass)+std_imp(high_pass))',log_freq(high_pass)', (avg_imp(high_pass)-std_imp(high_pass))','FillType',[1 2;2 1], 'Color', [0.5 0.5 0.5], 'LineStyle', 'none');%

xlabel('Frequency (Hz)');
xticks([1 2 3 4 5 6])
xticklabels({'10^{1}','10^{2}','10^{3}','10^{4}','10^{5}','10^{6}'})
ylabel('Impedance (M\Omega)');
set(gca,'box','off')
set(gcf,'color','w')
hold off;
set(gca, 'FontSize', FontSize);
saveas(gcf, 'impedances_without_and_with_PEDOT_IrOx_no_legend.svg');

%% log plots
figure;
hold on;
plot(log_freq(high_pass),log10(avg_imp(high_pass)), 'LineWidth', 2, 'Color', 'k');
width=5;
height=width;
pos = get(gcf, 'Position');
set(gcf, 'Position', [0 0 width*100, height*100]); %<- Set size
shade(log_freq(high_pass)', log10((avg_imp(high_pass)+std_imp(high_pass))'),log_freq(high_pass)', log10((avg_imp(high_pass)-std_imp(high_pass))'),'FillType',[1 2;2 1], 'Color', [0.5 0.5 0.5], 'LineStyle', 'none');%


xticks([1 2 3 4 5 6])
xticklabels({'10^{1}','10^{2}','10^{3}','10^{4}','10^{5}','10^{6}'});

yticks([-3 -2 -1 0 1])
yticklabels({'0.001','0.01','0.1','1','10'});

xL=xlabel('Frequency (Hz)');
yL=ylabel('Impedance (M\Omega)');
set(gca, 'FontSize',FontSize_zoomed);
xL.FontSize=20;
yL.FontSize=20;
set(gca,'box','off')
set(gcf,'color','w')
hold off;

saveas(gcf, 'impedances_after_PEDOT_IrOx_log.svg');

figure;
hold on;

plot(log10(freq_before_pedot),log10(imp_before_pedot), 'LineWidth', 2, 'Color', 'g');
plot(log_freq(high_pass),log10(avg_imp(high_pass)), 'LineWidth', 2, 'Color', 'k');



width=5;
height=width;
pos = get(gcf, 'Position');
set(gcf, 'Position', [0 0 width*100, height*100]); %<- Set size
shade(log_freq(high_pass)', log10((avg_imp(high_pass)+std_imp(high_pass))'),log_freq(high_pass)', log10((avg_imp(high_pass)-std_imp(high_pass))'),'FillType',[1 2;2 1], 'Color', [0.5 0.5 0.5], 'LineStyle', 'none');%


xticks([1 2 3 4 5 6])
xticklabels({'10^{1}','10^{2}','10^{3}','10^{4}','10^{5}','10^{6}'})
yticks([-3 -2 -1 0 1 2])
yticklabels({'0.001','0.01','0.1','1','10', '100'});

xL=xlabel('Frequency (Hz)');
yL=ylabel('Impedance (M\Omega)');
xL.FontSize=20;
yL.FontSize=20;

legend('Before','After', 'Location', 'SouthWest');
set(gca,'box','off')
set(gcf,'color','w')
hold off;
set(gca, 'FontSize', FontSize);
saveas(gcf, 'impedances_without_and_with_PEDOT_IrOx_legend_log.svg');

figure;
hold on;
plot(log_freq(high_pass),log10(avg_imp(high_pass)), 'LineWidth', 2, 'Color', 'k');

plot(log10(freq_before_pedot),log10(imp_before_pedot), 'LineWidth', 2, 'Color', 'g');

width=5;
height=width;
pos = get(gcf, 'Position');
set(gcf, 'Position', [0 0 width*100, height*100]); %<- Set size
shade(log_freq(high_pass)', log10((avg_imp(high_pass)+std_imp(high_pass))'),log_freq(high_pass)', log10((avg_imp(high_pass)-std_imp(high_pass))'),'FillType',[1 2;2 1], 'Color', [0.5 0.5 0.5], 'LineStyle', 'none');%

xlabel('Frequency (Hz)');
xticks([1 2 3 4 5 6])
xticklabels({'10^{1}','10^{2}','10^{3}','10^{4}','10^{5}','10^{6}'})
ylabel('Impedance (M\Omega)');
set(gca,'box','off')
set(gcf,'color','w')
hold off;
set(gca, 'FontSize', FontSize);
saveas(gcf, 'impedances_without_and_with_PEDOT_IrOx_no_legend_log.svg');
%% Nyquist plots

avg_zre=flip(mean(Zre,2))*10^-6;
std_zre=flip(std(Zre,0,2))*10^-6;

avg_zim=flip(mean(Zim,2))*10^-6;
std_zim=flip(std(Zim,0,2))*10^-6;

plot(avg_zre, -avg_zim, 'Color', 'k')

xlim([-0.1 1])
set(gca, 'FontSize', 14);
ylabel('-Im(Z) (M\Omega)');
xlabel('Re(Z) (M\Omega)');
width=5;
height=width;
pos = get(gcf, 'Position');
set(gcf, 'Position', [0 0 width*100, height*100]); %<- Set size

set(gca,'box','off')
set(gcf,'color','w')
hold off;
saveas(gcf, 'nyquists_IrOx.svg');