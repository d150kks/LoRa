clc
clear all
close all


%% ================================= CRC Type
lora_crc_type2_aos_ber = load("lora_crc_type2_aos_ber_v2.mat").BER;
snr = load("snr.mat").snr;

% BER = [lora_crc_type2_aos_ber; lora_crc_type2_aos_ber_v2(4:end,:)];
% save('lora_crc_type2_aos_ber_v2.mat','BER')
% return

figure(1)

xlabel('SNR (dB)')
ylabel('Bit Error Rate')

hold off
title('SNR');

legend_list = [];
mlist=["o", "+", "*", ".", "x", "_", "|", "square", "diamond", "^", "v", ">"];
clist=["#0072BD", "#D95319", "#EDB120", "#7E2F8E", "#77AC30", "#4DBEEE", "#A2142F", "#000000", "#FFFF00", "#FF00FF",  "#00FFFF", "#0000FF"];


counter = 1;
for aos = 1:6
    semilogy(snr,lora_crc_type2_aos_ber(aos,:),strcat('-',mlist(counter)),'color',clist(counter), 'LineWidth',2);
    legend_list = [legend_list, {strcat('aos', '----', num2str(aos))} ];
    xlim([-20 0])
    hold on
    counter = counter+1;
end
% semilogy(snr,lora_ber,'color', "#FF0000", 'LineWidth',2);
legend_list = [legend_list, {'No Code'} ];
grid on
legend(legend_list)



% figure(1)
% semilogy(snr, lora_crc_ber,'-*','color','k');
% hold on
% semilogy(snr,lora_crc_type2_ber,'-s','color','r');
% legend('Sort by amp', 'Sort by amp and position')
% grid
% xlabel('SNR (dB)')
% ylabel('Bit Error Rate')
% hold off
% title('SNR');



return
%% ================================= LDPC Config
lora_config_ldpc_ber = load("lora_config_ldpc_ber.mat").BER;
lora_ber = load("lora_ber.mat").BER;
snr = load("snr.mat").snr;


rate = [1/2 2/3 3/4 5/6];
cdwlen = [648, 1296, 1944]; % Codeword length




%% ================================= CRC+LDPC Config
% lora_crc_config_ldpc_ber = load("lora_crc_config_ldpc_ber.mat").BER;
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
%         semilogy(snr_config_ldpc,lora_crc_config_ldpc_ber(:,nLen,nRate),strcat('-',mlist(counter)),'color',clist(counter), 'LineWidth',2);
%         legend_list = [legend_list, {strcat(num2str(cdwlen(nLen)), '----', num2str(rate(nRate)))} ];
%         xlim([-20 0])
%         hold on
%         counter = counter+1;
%     end
% end
% semilogy(snr_crc,lora_crc_ber,'color', "#FF0000", 'LineWidth',2);
% legend_list = [legend_list, {'No Code'} ];
% grid on
% legend(legend_list)

