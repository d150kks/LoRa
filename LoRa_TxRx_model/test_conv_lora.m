clc
clear all
close all


tic

%% ================================= Переменные

% коэффициенты
SF = 9;        % коэффициент расширения спектра (от 7 до 12)
Base = 2^SF;   % База сигнала
num_sym = 10;% Число передаваемых символов
num_pre = 8;   % Число символов Преамбулы
num_check = 2; % Число символов Синхронизации
snr = 10;      % ОСШ

% массивы данных
numbits = 486;
data = randi([0, 1], 1, numbits);

% частота, фаза, время
Fmin = 0;              % минимальное значение частоты сигнала
Fmax = 125e3+Fmin;     % максимальное значение частоты сигнала
BW = Fmax-Fmin;        % ширина полосы частот
fc = 2.4e9;            % центральная частота
Ts = (Base)/BW;        % длительность сигнала


%% ================================= Mодуляция
os = 1;
fs = BW*os;
ts=1/fs;

[mod_chirp, check_chirp, downch, check_data] = lorax_modified( data, fs, Base, Ts, SF, BW);
N = length(downch);


% %% ================================= Demodulation
% figure(1)
% for i=1:Base
%     [sv, fourier] = delorax_modified2( SF, SF, downch, mod_chirp(i*Base-Base+1:Base*i));
%     
% %     subplot(211)
%     plot(real(fourier(i)), imag(fourier(i)),'o')
%     hold on
%     drawnow
%     
% %     subplot(212)
% %     stem(angle(fourier))
% %     return
% %     fourier(i)
% end
% 
% return
%% ================================= Convolutional coding
trellis = poly2trellis(7,[171 133]);
puncpat = [1;1;0];
K = log2(trellis.numInputSymbols); % Number of input streams
N = log2(trellis.numOutputSymbols); % Number of output streams
unpunc_coderate = K/N; % Unpunctured code rate
punc_coderate = (K/N)*length(puncpat)/sum(puncpat); % Punctured code rate
fprintf('K is %d and N is %d. The unpunctured code rate is %3.2f and the punctured code rate is %3.2f.\n',K,N,unpunc_coderate,punc_coderate)


ConstraintLength = (log2(trellis.numStates) + 1);
codedData = convenc(data,trellis,puncpat);
tbdepth = (ConstraintLength-1)/(1-punc_coderate); % Traceback depth for Viterbi decoder

codedData(10) = ~codedData(10);
codedData(131) = ~codedData(131);
codedData(256) = ~codedData(256);
codedData(366) = ~codedData(366);
codedData(406) = ~codedData(406);
codedData(455) = ~codedData(455);
codedData(544) = ~codedData(544);




%% ================================= Convolutional decoding
% rxbits = vitdec(codedData,trellis,tbdepth,'trunc','hard');
rxbits = vitdec(codedData,trellis,tbdepth,'trunc','hard',puncpat);

errStats = sum(data~=rxbits);
fprintf('Number of errors = %d\n', errStats)
