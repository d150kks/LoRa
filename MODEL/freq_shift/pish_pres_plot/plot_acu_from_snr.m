clc
clear all
close all

grid on

%% ================================= Comparison of the methods for the F shift comp
ferr_our = load('fest_our_snr.mat').est_err_from_pre;
ferr_ghan = load('fest_ghanaatian_snr.mat').est_err_from_pre;
ferr_gss = load('fest_gss_snr.mat').est_err_from_pre;
snr = [-12: 0];


% Plots
% xtext = 'Number of the preamble symbols';
% ytext = 'RMSE of the frequency error';
% legendtext = ['old','new'];
xtext = 'Количество передаваемых символов преамбулы';
ytext = 'RMSE ошибки оценки частоты';
legendtext = {'Предложенный метод','метод Ghanaatian', 'Метод Golden Section Search'};

figure(1); hold on
plot(snr, ferr_our, 'bo-')
plot(snr, ferr_ghan, 'ro-')
plot(snr, ferr_gss, 'ko-')
xlabel(xtext);
ylabel(ytext);
legend(legendtext)
% title('F offset only')

% figure(2); hold on
% plot(num_pre_list_old, fterr_old, 'bo-')
% plot(num_pre_list_new, fterr_new, 'ro-')
% xlabel('Number of the preamble symbols');
% ylabel('RMSE of the frequency error');
% legend('old','new')
% title('F and T offsets')

return
%% ================================= Comparison of BER
ferr_BER_old = load('fest_ber_old_std.mat').BER;
ferr_BER_new = load('fest_ber_new_std.mat').BER;

fterr_BER_old = load('ftest_ber_old_std.mat').BER;
fterr_BER_new = load('ftest_ber_new_std.mat').BER;

snr = -16:0;



figure(3); 
semilogy(snr, ferr_BER_old, 'bo-')
hold on
semilogy(snr, ferr_BER_new, 'ro-')
xlabel('Number of the preamble symbols');
ylabel('RMSE of the frequency error');
legend('old','new')
title('F offset only')

figure(4); 
semilogy(snr, fterr_BER_old, 'bo-')
hold on
semilogy(snr, fterr_BER_new, 'ro-')
xlabel('Number of the preamble symbols');
ylabel('RMSE of the frequency error');
legend('old','new')
title('F T offset')



