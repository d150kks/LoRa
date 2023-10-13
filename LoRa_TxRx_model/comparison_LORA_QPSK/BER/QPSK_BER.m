clear all
close all
clc

%% ================================= Variables
fc = 1000e6;
fs = 30e6;

M = 4;   
Pream_len = 4096;
Chip_len = 340;
Data_len = 100000;
snr = [-36:2:-16];

% Формируем преамбулу
x = randi([0, M-1],1,Pream_len); % Случайное сообщение
% Data
data = randi([0, 1],1,Data_len); % Случайное сообщение


%% ================================= Modulation
% Используем 16-точечную КАМ
y = qammod(x, M);

pream = [y, y];

% формируем сигналы для 0 и 1
chip_zero = (2*randi([0,1],1,Chip_len)-1) + 1i*(2*randi([0,1],1,Chip_len)-1);
chip_one = (2*randi([0,1],1,Chip_len)-1) + 1i*(2*randi([0,1],1,Chip_len)-1);


sig_mod = zeros(1,Data_len*Chip_len);
for i = 1:Data_len
    if(data(i) == 1)
        sig_mod(i*Chip_len-Chip_len+1:i*Chip_len) = chip_one;
    end
    
    if(data(i) == 0)
        sig_mod(i*Chip_len-Chip_len+1:i*Chip_len) = chip_zero;
    end
end


% frame_tx = [pream, sig_mod];



%% ================================= BER
tic
for n = 1:length(snr)

    [numErr, NumData] = deal(0);

    while  NumData < 1.1e5
        % АБГШ демодуляция и декодирование
        rxSig = awgn(sig_mod,snr(n),'measured');
        demod_data = QPSK_DEMOD(rxSig, chip_one, chip_zero, Chip_len, Data_len);
%         %% ================================= Demodulation
%         win_one = zeros(1,Chip_len);
%         win_zero = zeros(1,Chip_len);
%         demod_data = zeros(1,Data_len);
%         for i=1:Data_len
%             win_one = abs(sum(rxSig(i*Chip_len-Chip_len+1:i*Chip_len).*conj(chip_one)));
%             win_zero = abs(sum(rxSig(i*Chip_len-Chip_len+1:i*Chip_len).*conj(chip_zero)));
%             if(win_one>win_zero)
%                 demod_data(i) = 1;
%             else
%                 demod_data(i) = 0;
%             end
%         end
        

        
        % подсчет БЕР с учетом задержки
        err = sum(demod_data~=data);

        % Increment the error and bit counters
        numErr = numErr + err;        
        NumData = NumData + Data_len;
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

save('qpsk_ber.mat','BER')
save('qpsk_snr.mat','snr')


%% ================================= Demodulation
function demod_data = QPSK_DEMOD(rxSig, chip_one, chip_zero, Chip_len, Data_len)

    
    win_one = zeros(1,Chip_len);
    win_zero = zeros(1,Chip_len);
    demod_data = zeros(1,Data_len);
    for i=1:Data_len
        win_one = abs(sum(rxSig(i*Chip_len-Chip_len+1:i*Chip_len).*conj(chip_one)));
        win_zero = abs(sum(rxSig(i*Chip_len-Chip_len+1:i*Chip_len).*conj(chip_zero)));
        if(win_one>win_zero)
            demod_data(i) = 1;
        else
            demod_data(i) = 0;
        end
    end

end
