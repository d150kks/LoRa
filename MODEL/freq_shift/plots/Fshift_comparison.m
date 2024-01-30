clc
clear all
close all

grid on

%% ================================= Comparison of the methods for the F shift comp
ferr_old = load('fest_old_std.mat').est_err_from_pre;
ferr_new = load('fest_ghanaatian_std.mat').est_err_from_pre;
% ferr_new = load('fest_new_std.mat').est_err_from_pre;

% fterr_old = load('ftest_old_std.mat').est_err_from_pre;
% fterr_new = load('ftest_new_std.mat').est_err_from_pre;

num_pre_list_old = 2:1:8;
num_pre_list_new = 2:2:8;


% Plots
% xtext = 'Number of the preamble symbols';
% ytext = 'RMSE of the frequency error';
% legendtext = ['old','new'];
xtext = 'Количество передаваемых символов преамбулы';
ytext = 'RMSE ошибки оценки частоты';
legendtext = {'Предложенный метод','метод Ghanaatian'};

figure(1); hold on
plot(num_pre_list_old, ferr_old, 'bo-')
plot(num_pre_list_old, ferr_new, 'ro-')
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



