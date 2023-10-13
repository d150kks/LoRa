clc
clear all
close all


tic
% rng(1345)

%% ================================= Переменные

% коэффициенты
SF = 8;        % коэффициент расширения спектра (от 7 до 12)
bits2sym = SF-4;
BW = 30e6;

LORA = myLoRaClass(SF,BW);
Base = LORA.Base;
downch = LORA.downch;

num_pre = 8;

% num_sym = 10;
nbits = bits2sym*23; 
data = randi([0 1],1, nbits); 

%% ================================= Rate matching
[dataRM, numcodebitsRM, num_sym, zeros2end, flagRM] = LORA.RM(data);


%% ================================= CRC coding
[data_code] = LORA.codeCRC(dataRM, num_sym);

% return
%% ================================= Mодуляция
snr = -0;
[mod_chirp, check_data, check_no_gray] = LORA.lorax_modified( data_code, num_sym, 1);
tx_preamble = repmat(LORA.chirp,1,num_pre);
tx_sig = [tx_preamble, mod_chirp];
% tx_sig = conv(tx_sig, [0.3, -0.2, 0.6, -0.5], 'same');


%% ================================= АБГШ 
fps = BW/Base;
freq_shift = -fps*1.5;
% freq_shift = 0;
dphi=freq_shift*2*pi*(1/BW);% сдвиг

% вводим частотный сдвиг
for j=1:length(tx_sig)
    shift_sig(j)=tx_sig(j)*exp(1i*dphi*j);
end

rx_sig = awgn(shift_sig, snr,'measured');
rx_preamble = rx_sig(1:length(tx_preamble));
rxSig = rx_sig(length(tx_preamble)+1:end);


fcheck = abs(fft(rx_preamble(1:Base).*downch));
figure(1)
stem(fcheck);
%% ================================= демодуляция
aos = 2;
% [soft_bits, hard_bits, sv, sv_cor, fourier] = LORA.DELORAX_CRC( rxSig, num_sym, tx_preamble, rx_preamble, aos);

% Param
Gcode = LORA.grayCode;
hard_bits = [];
aos_win = -aos:aos;
nvar = std( abs(tx_preamble-rx_preamble).^2);
sv = zeros(1,num_sym);
sv_cor = zeros(1,num_sym);
soft_bits = [];

% Demodulation
for i = 1:num_sym

    d = rxSig(Base*i-Base+1:Base*i).*downch;   % перемножаем входной и опорный ОБРАТНый чирп
                
    fourier1 = abs(fft(d));            % переводим результат в область частот
    fourier = fourier1;
%     fourier = fourier1(Gcode+1);
    [~, indexMax] = max( fourier ); % находим щелчок  частоты в чирпе

    % CRC
% % % % % % % % % % % % % % % % % %     sv(i) = indexMax-1;
    sv(i) = indexMax;
% sv(i) = Gcode(indexMax);
    peak_win = indexMax+aos_win;
%     peak_win
    for pk=1:length(peak_win)
        if(peak_win(pk)<=0)
            peak_win(pk)=peak_win(pk)+Base;
        end
        if(peak_win(pk)>Base)
            peak_win(pk)=peak_win(pk)-Base;
        end
    end
%     peak_win(peak_win<=0)=[];
%     peak_win(peak_win>Base)=[];
% peak_win
    peaks_amp = fourier(peak_win);
    [sort_amp, sort_idx] = sort(peaks_amp,'descend');
    peak_sort = peak_win(sort_idx);
%     peak_sort = Gcode(peak_sort)+1;
% find(Gcode==63)
% return
    for n=1:length(peak_sort)
        sv_cor(i) = find(Gcode==peak_sort(n)-1)-1;
%         sv_cor(i) = peak_sort(n)-1;
        dbits = de2bi(sv_cor(i),SF)';
        dbits = dbits(:)';
        check_crc = sum(LORA.CRC4(dbits.').');
        peakMakcor = sort_amp(n);
        if(check_crc==0)
            break
        end
    end
    hard_bits = [hard_bits, dbits];

%     fourier = abs(fft(d));            % переводим результат в область частот
%     [~, indexMax] = max( fourier ); % находим щелчок  частоты в чирпе
%     [~, indexMaxGray] = max( fourier(Gcode+1) ); % находим щелчок  частоты в чирпе
% 
%     % вычисляем значение кодового слова исходя из базы сигнала
%     sv(i) = indexMax-1;
%     sv_decode(i) = indexMaxGray-1;
%     dbits = de2bi(sv_decode(i),SF)';
%     dbits = dbits(:)';
%     hard_bits = [hard_bits, dbits];

end
peak_win
[data_decode] = LORA.decodeCRC(hard_bits, num_sym, zeros2end, 0);

%% ================================= Debugging
errStats = sum(data~=data_decode);
fprintf('Number of errors = %d\n', errStats)


fprintf('check G = %s\n', num2str(check_data))
fprintf('sv G    = %s\n\n', num2str(sv))

fprintf('check  = %s\n', num2str(check_no_gray))
fprintf('sv_cor = %s\n', num2str(sv_cor))

% fprintf('Spectral efficiency = %d\n', errStats)

% figure(1)
% stem(fourier)
% hold on
% stem(fourier1)
% legend('no gray', 'gray')

% figure(1)
% plot(normalize(double(rxbits)))
% hold on
% plot(normalize(data))
% 
% return
