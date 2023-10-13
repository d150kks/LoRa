clear all
close all
clc

fc = 1000e6;
fs = 30e6;
Buf_len_koef = 5;
dec_bit = 255;


f0 = 15e6;
f2 = f0+15e6;
f1 = f0-15e6;
% f3 = f1+5;
BW = f2-f1;
sf = 11;
N = fs*2^sf/BW;
T = N/fs;
t = (0:N-1)/fs;

% sf=log2(BW*N/fs);


mu = 2*pi*BW/T;
chirp = sin(2*pi*f1*t+mu*t.^2/2);
n_bit = sf;


% формируем массив временных отрезков кратных sf

n_start = (N/2^sf)*dec_bit;
chirp_shift = circshift(chirp,[0,-n_start+N/2^sf]);
pream = chirp(end:-1:1);
pream = pream/max(pream);




n_start1 = (N/2^sf)*25;
chirp_shift1 = circshift(chirp,[0,-n_start1+N/2^sf]);
pream = chirp(end:-1:1);
pream = pream/max(pream);




frame_tx = [pream,pream, chirp_shift, chirp_shift1];
t_sdvig = 1e-99;

%  frame_tx = time_sdvig(frame_tx, t_sdvig, fs);
frame_tx = complex(frame_tx);

figure
plot(real(frame_tx));
Tx_Data = frame_tx.';



%% Передатчик

return
tx = sdrtx('Pluto','RadioID','usb:0'); % Определяем устройство передачи 
tx.ShowAdvancedProperties = true;
tx.UseCustomFilter = false;
tx.CenterFrequency = fc; % Задаем несущую частоту
tx.BasebandSampleRate = fs; % Задаем полосу частот
tx.Gain =-45; % -89.75 ... 0
transmitRepeat(tx, Tx_Data); % Передача данных Tx_Data
% tx.ShowAdvancedProperties = true;
% 

%% Приемник

rx = sdrrx('Pluto','RadioID','usb:1'); % Определяем устройство передачи
rx.CenterFrequency = fc; % Выбираем несущую частоту
rx.BasebandSampleRate = fs; % Выбираем полосу частот
rx.Gain = 10; % -4 ... 71

rx.ShowAdvancedProperties = true;
rx.EnableQuadratureCorrection = true;
rx.EnableRFDCCorrection = true;
rx.EnableBasebandDCCorrection = true;

rx.OutputDataType = 'double';
rx.SamplesPerFrame = length(Tx_Data)*Buf_len_koef; % Размер буффера



rx_data_in_time=[];
rx_data_in_time(:,1) = rx(); % Принятые данные
rx_data_in_time = rx_data_in_time.';

% figure(1)
% plot(real(rx_data_in_time))
% return
korr = xcorr(rx_data_in_time, [pream,pream]);

[A,B] = max(korr(1:length(rx_data_in_time)+length(rx_data_in_time)/2));

start = B-length(rx_data_in_time) +1;

Pream_len = N;
 sym1 = rx_data_in_time(start:start+Pream_len-1);
 sym2 = rx_data_in_time(start+Pream_len:start+2*Pream_len-1);
figure
plot(abs(korr));

corr=xcorr(sym2,sym1);

[a,b]=max(corr);
figure
plot(abs(corr));

% Оценка фазового набега
Dphi=angle(a);

dt=1/fs;
tau=(1000)*dt;

ocen_freq=Dphi/(2*pi*tau)


% расчитываем фазовый набег за интервал дискретизации по полученной оценке
d_phi_ocen=(ocen_freq*2*pi)/fs;


% компенсируем частотный вводя фазовый набег но с отрицательным знаком
for j=1:length(rx_data_in_time)
rx_data_in_time_comp(j)=rx_data_in_time(j)*exp(1i*j*(-d_phi_ocen));
end
figure
plot(abs(korr));

chirp_shift_channel = rx_data_in_time_comp(start+N*2:start+3*N-1);
chirp_inv=chirp(end:-1:1);
 chirp_pr  =(chirp_shift_channel.*chirp_inv);

  
%   dw = 
%   chirp_pr_1 = cos(
  
  
 chirp_pr_sp  =fft(chirp_pr);


 
 figure
plot(((chirp_pr)));
set(gca,'FontName','Times New Roman','fontsize',12);
grid on

ylabel('Амплитуда');
xlabel('Время, с');
% 
%  
  figure
 plot((abs(chirp_pr_sp)));
 grid on
set(gca,'FontName','Times New Roman','fontsize',12);
grid on

ylabel('Амплитуда');
xlabel('Частота, Гц');
% 


[AA,BB] = max(abs(chirp_pr_sp(1:N/2)));
if BB > 2^sf
BB = BB-2^sf;
end
BB = BB+1;

% 
 demod_bi = de2bi(BB);
 demod_bit = zeros(1,sf);
 demod_bit(1:length(demod_bi)) = demod_bi;

 spec_eff = N/sf 

 chirp_shift_channel = rx_data_in_time_comp(start+N*3:start+4*N-1);
chirp_inv=chirp(end:-1:1);
 chirp_pr  =(chirp_shift_channel.*chirp_inv);

  
%   dw = 
%   chirp_pr_1 = cos(
  
  
 chirp_pr_sp  =fft(chirp_pr);


 
 figure
plot(((chirp_pr)));
set(gca,'FontName','Times New Roman','fontsize',12);
grid on

ylabel('Амплитуда');
xlabel('Время, с');
% 
%  
  figure
 plot((abs(chirp_pr_sp)));
 grid on
set(gca,'FontName','Times New Roman','fontsize',12);
grid on

ylabel('Амплитуда');
xlabel('Частота, Гц');
% 


[AA,CC] = max(abs(chirp_pr_sp(1:N/2)));
if CC > 2^sf
CC = CC-2^sf;
end
CC = CC+1;

% 
% 
% 





% sig_inf_reciv = rx_data_in_time(start+Pream_len*2:start+Pream_len*2+Chip_len*Data_len-1);
% 
% 
% corr=xcorr(sig_inf_reciv,chip_one);
% figure
% plot(abs(corr));
% 
% chip_matr = reshape(sig_inf_reciv, Chip_len, Data_len);
% 
% chip_matr = chip_matr.';
% 
% for i = 1: Data_len
% korr_koeff_one(i) = abs(sum(chip_matr (i,:).*conj(chip_one)));
% korr_koeff_zero(i) = abs(sum(chip_matr (i,:).*conj(chip_zero)));
% if (korr_koeff_one (i)>korr_koeff_zero(i))
% bits_demod (i) = 1;
% else
%     bits_demod (i) = 0;
% end
% end
% 
% [number,ratio] = biterr(data,bits_demod)
