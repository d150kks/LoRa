clc
clear all
close all


tic
rng default
%% ================================= Переменные

% коэффициенты
SF = 9;        % коэффициент расширения спектра (от 7 до 12)
bits2sym = SF-4;
BW = 125e3;
fc = 2200e6;

LORA = myLoRaClass(SF,BW);
Base = LORA.Base;
Ts = LORA.Ts;
chirp = LORA.chirp;
downch = LORA.downch;

num_pre = 8;
numcodebits = 648;
data = randi([0 1],1, numcodebits); 
data_ldpc_code = data;


%% ================================= Mодуляция
snr = 0;
num_sym = numcodebits/SF;
[mod_chirp, check_data, check_no_gray] = LORA.lorax_modified( data_ldpc_code, num_sym, 1);
tx_preamble = repmat(chirp,1,num_pre);
% tx_sig = [downch, tx_preamble, mod_chirp];
tx_sig = [tx_preamble, mod_chirp];

%% ================================= АБГШ 
rx_sig = awgn(tx_sig, snr,'measured');
rx_preamble = rx_sig(1:length(tx_preamble));
rxSig = rx_sig(length(tx_preamble)+1:end);


%% ================================= демодуляция
% aos = 5;
% [sv, product] = delorax_cyclic( num_sym, SF, chirp, rxSig, LORA.grayCode);
% dbits = de2bi(sv,SF).';
% rxbits = dbits(:).';
N=2^SF;

product = zeros(1,N);
% sv = zeros(1)
for i = 1:num_sym
    
    extract_chirp = mod_chirp(i*N-N+1:i*N);
    for j=0:N-1
        j = 20;
        cyclic_chirp = circshift(chirp,-j);
%         cyclic_chirp2 = circshift(chirp,-j);
        product(j+1) = abs(sum(extract_chirp.*conj(cyclic_chirp)));

        
    end
%         [~, indexMax] = max( product ); % находим щелчок  частоты в чирпе
        [~, indexMax] = max( product(LORA.grayCode+1) ); % находим щелчок  частоты в чирпе
        sv(i) = indexMax-1;
    figure(1)
%     stem(product)
%     plot(real(extract_chirp))
    hold on
    plot(real(cyclic_chirp))
    plot(real(cyclic_chirp2))
    return
end

dbits = de2bi(sv,SF).';
rxbits = dbits(:).';
% end


    figure(1)
    stem(product)
%     plot(real(extract_chirp))
%     hold on
%     plot(real(cyclic_chirp))


%% ================================= Debugging
errStats = sum(data~=rxbits);
fprintf('Number of errors = %d\n', errStats)



% fprintf('Bits per Sec = %.2f B/s \n', sum(sv_cor-check_data))
% fprintf('Spectral efficiency = %d\n', errStats)



% figure(1)
% plot(normalize(double(rxbits)))
% hold on
% plot(normalize(data))
% 
% return
function [sv, product] = delorax_cyclic( num_sym, SF, chirp, mod_chirp, G)
    N=2^SF;
    
    product = zeros(1,N);
    for i = 1:num_sym
        
        extract_chirp = mod_chirp(i*N-N+1:i*N);
        for j=0:N-1
            cyclic_chirp = circshift(chirp,-j);
            product(j+1) = abs(sum(extract_chirp.*conj(cyclic_chirp)));
        end 
        [~, indexMax] = max( product(G+1) ); % находим щелчок  частоты в чирпе
        sv(i) = indexMax-1;

        % вычисляем значение кодового слова исходя из базы сигнала
%         [~, indexMax] = max( product ); % находим щелчок  частоты в чирпе
%         sv(i) = indexMax;

        
    end
end