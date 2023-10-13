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
numinfobits = 486;
numcodebits = 648;
data = randi([0, 1], 1, numinfobits); % формирование массива бит
num_sym = numcodebits/SF;
num_pre = 8;





%% ================================= LDPC coding
P = [16 17 22 24  9  3 14 -1  4  2  7 -1 26 -1  2 -1 21 -1  1  0 -1 -1 -1 -1
     25 12 12  3  3 26  6 21 -1 15 22 -1 15 -1  4 -1 -1 16 -1  0  0 -1 -1 -1
     25 18 26 16 22 23  9 -1  0 -1  4 -1  4 -1  8 23 11 -1 -1 -1  0  0 -1 -1
      9  7  0  1 17 -1 -1  7  3 -1  3 23 -1 16 -1 -1 21 -1  0 -1 -1  0  0 -1
     24  5 26  7  1 -1 -1 15 24 15 -1  8 -1 13 -1 13 -1 11 -1 -1 -1 -1  0  0
      2  2 19 14 24  1 15 19 -1 21 -1  2 -1 24 -1  3 -1  2  1 -1 -1 -1 -1  0
    ];

blockSize = 27;
pcmatrix = ldpcQuasiCyclicMatrix(blockSize,P);
cfgLDPCEnc = ldpcEncoderConfig(pcmatrix);
codeword = ldpcEncode(data.', cfgLDPCEnc).';

%% ================================= модуляция
% data = 0:Base-1;
% data = reshape(de2bi(data).',1,[]);
G = LORA.grayCode;

% num_sym = Base;
[mod_chirp, check_data] = LORA.lorax_modified( codeword, num_sym, 1);

% datanogray = de2bi(check_data_no_gray);
% datagray = de2bi(check_data);
% k = find(data==16)
% 

% return

%% ================================= АБГШ 
snr = -0;
preamble = repmat(LORA.chirp,1,num_pre);
[rxSig, nvar] = awgn([preamble, mod_chirp], snr, 'measured');
rx_pre = rxSig(1:length(preamble));
rxSig = rxSig(length(preamble)+1:end);

%% ================================= Mapper
bitmap = de2bi(0:2^SF-1);
M0 = zeros(Base,SF);
M1 = zeros(Base,SF);
A=[];
for nBit=1:SF
    for nSym = 1:Base
        if( bitmap(nSym,nBit) == 0)
            M0(nSym,nBit)= nSym;
        else
            M1(nSym,nBit)= nSym;
        end
    end
end
% LORA.grayCode(

%% ================================= Demodulation
% [sv, fourier] = delorax_modified( length(data), SF, downch, rxSig);
% dbits = de2bi(sv,SF).';
% rxbits = dbits(:).';

num = length(downch);
B=2^SF;

sv = zeros(1,num_sym);
L = [];
for i = 1:num_sym

    d = rxSig(num*i-num+1:num*i).*downch;   % перемножаем входной и опорный ОБРАТНый чирп
    
    fourier = abs(fft(d));            % переводим результат в область частот
    [peakMak, indexMax] = max( fourier ); % находим щелчок  частоты в чирпе

    % вычисляем значение кодового слова исходя из базы сигнала
    sv(i) = indexMax-1;

    % LLR

        for nBit=1:SF
          m0 = M0(:,nBit);
          m0(m0==0)=[];
          m1 = M1(:,nBit);
          m1(m1==0)=[];
%           LLR = -(1/nvar)*(min( (peakMak-fourier(m0)).^2 ) - min( (peakMak-fourier(m1)).^2 ));
            LLR = -(1/nvar)*(min( (peakMak-fourier( LORA.grayCode(m0)+1 )).^2 ) - min( (peakMak-fourier( LORA.grayCode(m1)+1 )).^2 ));
          L = [L LLR];
        end

end

%% ================================= LDPC decoding
maxnumiter = 100;
cfgLDPCDec = ldpcDecoderConfig(pcmatrix);
rxbits = ldpcDecode(L.', cfgLDPCDec, maxnumiter).';

errStats = sum(data~=rxbits);



%% ================================= Debugging

% evm = sqrt( (sum(abs(real(preamble)-real(rx_pre))).^2 + ...
%              sum(abs(imag(preamble)-imag(rx_pre))).^2)/length(preamble));
% evm = std( abs(rx_pre-preamble).^2);
evm = std( abs(rx_pre-preamble).^2);
fprintf('Number of errors = %d\n', errStats)
fprintf('nvar = %f\n', nvar)
fprintf('evm = %f\n', evm)

% figure(1)
% plot(real(rx_pre))
% hold on
% plot(real(preamble),'LineWidth',2)

figure(1)
plot(normalize(L))
hold on 
plot(normalize(data))

%%%%%%
% scatterplot(mod_sig)
% figure(1)
% plot(normalize(L), 'LineWidth',2)
% hold on
% plot(normalize(demod_sig))
% legend('my llr', 'llr')
% grid on
% title('QPSK LLR')


% function [sv, fourier] = delorax_llr( length_data, SF, downch, chirp)
%     num = length(downch);
%     B=2^SF;
% 
%     num_sym = length_data/SF;
%     sv = zeros(1,num_sym);
%     for i = 1:num_sym
% 
%         d = chirp(num*i-num+1:num*i).*downch;   % перемножаем входной и опорный ОБРАТНый чирп
%         
%         fourier = abs(fft(d));            % переводим результат в область частот
%         [~, indexMax] = max( abs(fourier+fourier2+fourier3+fourier4) ); % находим щелчок  частоты в чирпе
%         % css = fourier
%         %%%%%%%
%         % вычисляем значение кодового слова исходя из базы сигнала
%         
%         sv(i) = indexMax-1;
%         % indexMax-1;
%         % if indexMax>B
%         %     sv(i) = B - (num - indexMax) - 1;
%         % else
%         %     sv(i) = indexMax - 1;
%         % end
%     
%     end
% 
% end

function grayCode = generateGrayCode(n)
    % Generate Gray code sequence of length 2^n
    grayCode = zeros(1, 2^n);
    
    for i = 0:2^n-1
        grayCode(i+1) = bitxor(i, bitshift(i, -1));
    end
end