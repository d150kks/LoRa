clc
clear all
close all


%% ================================= Переменные
% коэффициенты
SF = 7;        % коэффициент расширения спектра (от 7 до 12)
rc_size = 4;
rc = (SF-rc_size);
BW = 2e6;
snr = [-16:1:0];
nIter = 10;

LORA = myLoRaClass_true(SF,BW);
Base = LORA.Base;
downch = LORA.downch;
chirp = LORA.chirp;

num_pre = 4;

num_sym = 1000;
numinfobits = num_sym*rc; 
% numinfobits = num_sym*SF; 
data = randi([0 1],1, numinfobits); 





% m0 = LORA.M0
% 
% return
%% ================================= Mодуляция
[mod_chirp, check_data, check_no_gray] = LORA.lorax_modified_crcrs(data, num_sym);
% [mod_chirp, check_data, check_no_gray] = LORA.lorax_modified(data, num_sym, 1);
tx_preamble = repmat(chirp, 1, num_pre);

%% ================================= BER
% return
tic
% snr = 10;
for n = 1:length(snr)
    fprintf('Iter: %d\n', n) 

    [numErr, NumData] = deal(0);

    for iter = 1:nIter
        
        % АБГШ 
        rxSig = awgn(mod_chirp, snr(n), 'measured');
        rx_preamble = awgn(tx_preamble,snr(n),'measured');

        % демодуляция
%         [soft_bits, hard_bits, sv_decode, sv, fourier] = LORA.delorax_modified( rxSig, num_sym, tx_preamble, rx_preamble);
        [soft_bits, hard_bits, sv_rs, sv, fourier, fourier_rs] = LORA.delorax_crcrs( rxSig, num_sym, tx_preamble, rx_preamble);


%         figure(1)
%         plot(normalize(hard_bits))
%         hold on
%         plot(-normalize(soft_bits))
%         xlim([1 100])
%         return

        % подсчет БЕР с учетом задержки
        err = sum(hard_bits~=data);

        % Increment the error and bit counters
        numErr = numErr + err;        
        NumData = NumData + numinfobits;
    end

    % Estimate the BER for both methods
    BER(n) = numErr/NumData;
end
toc

%%
figure(1)
semilogy(snr,BER,'-*','color','k');
grid
xlabel('SNR (dB)')
ylabel('Bit Error Rate')
hold off
title('SNR');


% save('lora_ber.mat','BER')
% save('lora_rsg_ber.mat','BER')