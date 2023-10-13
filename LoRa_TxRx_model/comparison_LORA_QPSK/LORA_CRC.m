clc
clear all
close all

%% ================================= Переменные

% коэффициенты
SF = 12;        % коэффициент расширения спектра (от 7 до 12)
Base = 2^SF;   % База сигнала
num_sym = 1;% Число передаваемых символов
num_pre = 4;   % Число символов Преамбулы
num_check = 2; % Число символов Синхронизации
snr = [-36:2:-16];

% массивы данных
bits_per_sym = (SF-4);
numbits = bits_per_sym*num_sym;          % число бит
data = randi([0, 1], 1, numbits); % формирование массива бит

% частота, фаза, время
Fmin = 0;              % минимальное значение частоты сигнала
Fmax = 30e6+Fmin;     % максимальное значение частоты сигнала
BW = Fmax-Fmin;        % ширина полосы частот
fc = 1000e6;            % центральная частота
band = 30e6;
Ts = (2^SF)/BW;        % длительность сигнала


%% ================================= CRC coding
data_win = zeros(1, bits_per_sym);
for i=1:num_sym
    data_win = data(i*bits_per_sym-bits_per_sym+1:i*bits_per_sym);
    CRC_Bits = CRC4(data_win.').';
    data_code(i*SF-SF+1:i*SF) = [data_win, CRC_Bits].';
end

%% ================================= модуляция
os = 1;
fs = BW*os;
ts=1/fs;

[mod_chirp, check_chirp, downch, check_data] = lorax_modified( data_code, fs, Base, Ts, SF, BW);
% [mod_chirp, check_chirp, downch, check_data] = LORAX_CRC( data_code, fs, Base, Ts, SF, BW);
N = length(downch);

% Forming Signal
tx_chirp = mod_chirp;
tx_length = length(tx_chirp);


%% ================================= Demodulation
rxSig = awgn(mod_chirp, -26,'measured');
% [demod_bits, sv, fourier] = DELORAX_CRC( length(data_code), SF, downch, rxSig);
num = length(downch);
B=2^SF;
demod_bits = [];

for i = 1:length(data_code)/SF

    % Fourier
    d = rxSig(num*i-num+1:num*i).*downch;   % перемножаем входной и опорный ОБРАТНый чирп
    fourier = abs(fft(d));            % переводим результат в область частот
    % fourier = abs(fourier(end/2:end));
    [~, indexMax] = max( fourier ); % находим щелчок  частоты в чирпе

    % вычисляем значение кодового слова исходя из базы сигнала
    if indexMax>B
        sv(i) = B - (num - indexMax) - 1;
    else
        sv(i) = indexMax - 1;
    end
    
    % CRC
    for n=-1:1
        sv_cor(i) = sv(i)+n;
        dbits = de2bi(abs(sv(i)),SF)';
        dbits = dbits(:)';
        check_crc = sum(CRC4( dbits.').');
        if(check_crc==0)
            break
        end
    end
    demod_bits = [demod_bits, dbits];
end

fprintf('sv_tx = %d\n', check_data)
fprintf('sv_rx = %d\n', sv)
fprintf('data_code  = %s\n', num2str(data_code))
fprintf('demod_bits = %s\n', num2str(demod_bits))

figure(1)
plot(abs(fourier))
return

%% ================================= CRC Decoding
data_decode = zeros(1, numbits);
for i=1:num_sym
    data_win = demod_bits(i*SF-SF+1:i*SF);
    data_decode(i*bits_per_sym-bits_per_sym+1:i*bits_per_sym) = data_win(bits_per_sym-bits_per_sym+1:bits_per_sym);
end
% [sv, fourier] = delorax_modified( length(data_code), SF, downch, rxSig);
% dbits = de2bi(sv,SF)';
% dbits = dbits(:)';



% подсчет БЕР с учетом задержки
err = sum(data_decode~=data)





return
%% ================================= BER
tic
for n = 1:length(snr)

    [numErr, NumData] = deal(0);

    while  NumData < 1e5
        % АБГШ демодуляция и декодирование
        rxSig = awgn(mod_chirp,snr(n),'measured');
        [sv, fourier] = delorax_modified( length(data), SF, downch, rxSig);
        dbits = de2bi(sv,SF)';
        dbits = dbits(:)';
        

        
        % подсчет БЕР с учетом задержки
        err = sum(dbits~=data);

        % Increment the error and bit counters
        numErr = numErr + err;        
        NumData = NumData + numbits;
    end

    % Estimate the BER for both methods
    BER(n) = numErr/NumData;
end
toc

%%
% figure(1)
% semilogy(snr,BER,'-*','color','k');
% grid
% xlabel('SNR (dB)')
% ylabel('Bit Error Rate')
% hold off
% title('SNR');

% save('lora_ber.mat','BER')
% save('lora_snr.mat','snr')






% function [upchirp, check_chirp, downchirp, check_data] = LORAX_CRC( data, fs, B, Ts, SF, BW)
%     % размерности
%     N = length(data);       % Размер последовательности данных
%     
%     % время и частота
%     ts = 1/fs;              % время дискретизации
%     m = BW/Ts;              % коэффициент измнения частоты
%     t = 0:ts:Ts-ts;         % полоса времени ОДНОГО чирпа
% 
%     % формирование пустых массивов
%     chirp_len = length(t);
%     chirp_num = N/SF;
%     code_word = zeros(1,SF); % текущее кодовое слово
%     check_data = zeros(1,chirp_num);
%     upchirp = zeros(1,chirp_len*chirp_num);           % чирпы
% %     upchirp = [];           % чирпы
%     
%     %% Создание чирпов
%     check_chirp = exp( 1i * (2*pi*(0*t + (m*t.^2)/2)) );
%     chirp = exp( 1i * (2*pi*(0*t + (m*t.^2)/2)) );
%     downchirp = exp( -1i * (2*pi*(0*t + (m*t.^2)/2)) );
%     
%     %% Модуляция
%     for i = 1:chirp_num
%         code_word = bi2de(data(SF*i-SF+1:SF*i));      % значение кодового слова
%         check_data(i) = code_word; % значение закодированного символа
%         cs = single((code_word/B)*Ts/ts);             % место сдвига
%         
%         chirp1 = chirp(cs+1:end); % первый чирп
%         chirp2 = chirp(1:cs); % второй чирп
%         
%         %%%%%%%
%         upchirp(i*chirp_len-chirp_len+1:i*chirp_len) = cat(2,chirp1,chirp2);
% %         upchirp = [upchirp, chirp1, chirp2];          % выходной чирп
%     end
% 
% end
function [demod_bits, sv, fourier] = DELORAX_CRC( length_data, SF, downch, chirp)
    num = length(downch);
    B=2^SF;
    demod_bits = [];

    for i = 1:length_data/SF

        % Fourier
        d = chirp(num*i-num+1:num*i).*downch;   % перемножаем входной и опорный ОБРАТНый чирп
        fourier = abs(fft(d));            % переводим результат в область частот
        % fourier = abs(fourier(end/2:end));
        [~, indexMax] = max( fourier ); % находим щелчок  частоты в чирпе

        % вычисляем значение кодового слова исходя из базы сигнала
        if indexMax>B
            sv(i) = B - (num - indexMax) - 1;
        else
            sv(i) = indexMax - 1;
        end
        
        % CRC
        for n=-1:1
            sv(i) = sv(i)+n;
            dbits = de2bi(abs(sv(i)),SF)';
            dbits = dbits(:)';
            check_crc = sum(CRC4( dbits.').');
            if(check_crc==0)
                break
            end
        end
        demod_bits = [demod_bits, dbits];
    end

end