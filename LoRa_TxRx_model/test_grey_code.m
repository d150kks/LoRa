clc
clear all
close all


%% ================================= Переменные

% коэффициенты
SF = 8;        % коэффициент расширения спектра (от 7 до 12)
BW = 125e3;

LORA = myLoRaClass(SF,BW);
Base = LORA.Base;
downch = LORA.downch;

% массивы данных
numinfobits = 648;
numcodebits = 648;
data = randi([0, 1], 1, numinfobits); % формирование массива бит
num_sym = numcodebits/SF;
num_pre = 8;




% %% ================================= Grey coding
% datade = 1:Base;
% databi = de2bi(datade);
% generateGrayCode(SF)


%% ================================= модуляция
[mod_chirp, check_data] = LORA.lorax_modified( data, num_sym, 1);

%% ================================= АБГШ 
snr = 100;
[rxSig, nvar] = awgn( mod_chirp, snr, 'measured');
% preamble = repmat(check_chirp,1,num_pre);
% [rxSig, nvar] = awgn([preamble, mod_chirp], snr, 'measured');
% rx_pre = rxSig(1:length(preamble));
% rxSig = rxSig(length(preamble)+1:end);


%% ================================= Demodulation
[sv_decode, sv, fourier] = LORA.delorax_modified(rxSig, num_sym);
Gdbits = de2bi(sv_decode,SF).';
Gdbits = Gdbits(:).';

dbits = de2bi(sv,SF).';
dbits = dbits(:).';


%% ================================= Debugging

% evm = sqrt( (sum(abs(real(preamble)-real(rx_pre))).^2 + ...
%              sum(abs(imag(preamble)-imag(rx_pre))).^2)/length(preamble));
% evm = std( abs(rx_pre-preamble).^2);
GerrStats = sum(data~=Gdbits);
errStats = sum(data~=dbits);

fprintf('Number of errors Gray = %d\n', GerrStats)
fprintf('Number of errors      = %d\n', errStats)

% figure(1)
% plot(real(rx_pre))
% hold on
% plot(real(preamble),'LineWidth',2)

% figure(1)
% plot(normalize(L))
% hold on 
% plot(normalize(data))

%%%%%%
% scatterplot(mod_sig)
% figure(1)
% plot(normalize(L), 'LineWidth',2)
% hold on
% plot(normalize(demod_sig))
% legend('my llr', 'llr')
% grid on
% title('QPSK LLR')

function grayCode = generateGrayCode(n)
    % Generate Gray code sequence of length 2^n
    grayCode = zeros(1, 2^n);
    
    for i = 0:2^n-1
        grayCode(i+1) = bitxor(i, bitshift(i, -1));
    end
end
