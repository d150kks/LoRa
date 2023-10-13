function [upchirp, check_chirp, downchirp, check_data] = lorax_modified( data, fs, B, Ts, SF, BW)
% размерности
N = length(data);       % Размер последовательности данных

% время и частота
% fs = 1e6;              % частота дискретизации
ts = 1/fs;              % время дискретизации
m = BW/Ts;              % коэффициент измнения частоты
t = 0:ts:Ts-ts;         % полоса времени ОДНОГО чирпа
% t = 0:127;

% формирование пустых массивов
code_word = zeros(1,SF); % текущее кодовое слово
upchirp = [];           % чирпы

%% Создание чирпов
% upchirp = zeros(1,samples);
% chirp = A*cos(2*pi*(0*t + (m*t.^2)/2));
% downchirp = flip(chirp);
% check_chirp = (1/sqrt(B))*exp( 1i * (2*pi*(0*t + (m*t.^2)/2)) );
% chirp = (1/sqrt(B))*exp( 1i * (2*pi*(0*t + (m*t.^2)/2)) );
% downchirp = (1/sqrt(B))*exp( -1i * (2*pi*(0*t + (m*t.^2)/2)) );

check_chirp = exp( 1i * (2*pi*(0*t + (m*t.^2)/2)) );
chirp = exp( 1i * (2*pi*(0*t + (m*t.^2)/2)) );
downchirp = exp( -1i * (2*pi*(0*t + (m*t.^2)/2)) );

%% Модуляция
for i = 1:N/SF
code_word = bi2de(data(SF*i-SF+1:SF*i));      % значение кодового слова
check_data(i) = code_word; % значение закодированного символа
cs = single((code_word/B)*Ts/ts);             % место сдвига

chirp1 = chirp(cs+1:end); % первый чирп
chirp2 = chirp(1:cs); % второй чирп

%%%%%%%
upchirp = [upchirp, chirp1, chirp2];          % выходной чирп
end

end