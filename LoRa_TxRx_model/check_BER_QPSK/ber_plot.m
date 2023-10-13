clc
clear all
close all

%% ================================= Переменные

% qpsk_ber = load("qpsk_ber.mat").BER;
% qpsk_ldpc_ber = load("qpsk_ldpc_ber").BER;
% snr_no_code = load("snr.mat").snr;
% 
% figure(3)
% semilogy(snr_no_code,qpsk_ber,'color', "#FF0000", 'LineWidth',2);
% hold on
% semilogy(snr_no_code,qpsk_ldpc_ber,'color', "k", 'LineWidth',2);
% grid on
% legend('no code', 'ldpc')
% xlabel('SNR (dB)')
% ylabel('Bit Error Rate')
% 
% hold off
% title('SNR');


%% 
lora_ber = load("lora_ber.mat").BER;
lora_llr_ldpc_ber = load("lora_llr_ldpc_ber").BER;
snr_lora = load("snr_lora.mat").snr;

figure(3)
semilogy(snr_lora,lora_ber,'color', "#FF0000", 'LineWidth',2);
hold on
semilogy(snr_lora,lora_llr_ldpc_ber,'color', "k", 'LineWidth',2);
grid on
legend('no code', 'ldpc')
xlabel('SNR (dB)')
ylabel('Bit Error Rate')

hold off
title('SNR');


% %% ================================= LDPC Config
% qpsk_config_ldpc_ber = load("qpsk_config_ldpc_ber.mat").BER;
% snr_config_ldpc = load("snr_config_ldpc.mat").snr;
% 
% 
% rate = [1/2 2/3 3/4 5/6];
% cdwlen = [648, 1296, 1944]; % Codeword length
% 
% figure(3)
% 
% xlabel('SNR (dB)')
% ylabel('Bit Error Rate')
% 
% hold off
% title('SNR');
% 
% legend_list = [];
% mlist=["o", "+", "*", ".", "x", "_", "|", "square", "diamond", "^", "v", ">"];
% clist=["#0072BD", "#D95319", "#EDB120", "#7E2F8E", "#77AC30", "#4DBEEE", "#A2142F", "#000000", "#FFFF00", "#FF00FF",  "#00FFFF", "#0000FF"];
% 
% 
% counter = 1;
% for nLen = 1:3
%     for nRate = 1:4
%         semilogy(snr_config_ldpc,qpsk_config_ldpc_ber(:,nLen,nRate),strcat('-',mlist(counter)),'color',clist(counter), 'LineWidth',2);
%         legend_list = [legend_list, {strcat(num2str(cdwlen(nLen)), '----', num2str(rate(nRate)))} ];
%         xlim([-20 0])
%         hold on
%         counter = counter+1;
%     end
% end
% semilogy(snr_no_code,qpsk_ber,'color', "#FF0000", 'LineWidth',2);
% legend_list = [legend_list, {'No Code'} ];
% grid on
% legend(legend_list)

