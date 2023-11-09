clc
clear all
close all


%% ================================= Переменные
% коэффициенты
M = 4;
snr = 5;

numinfobits = 100;          % число бит
data = randi([0, 1], 1, numinfobits,'int8'); % формирование массива бит


%% ================================= модуляция
mod_sig = pskmod(data.', M,InputType='bit').';

%% ================================= АБГШ 
[rxSig, nvar] = awgn(mod_sig,snr,'measured');
% nvar
%% ================================= демодуляция
demod_sig = pskdemod(rxSig.', M, OutputType='approxllr',NoiseVariance=nvar).';

S = [0,0;...
     0,1;...
     1,0;...
     1,1];

S0 = [1+1i*0, 0+1i, 0-1i, -1+1i*0];

for i=1:2
    idx_s0 = find(S(:,i)==0);
    idx_s1 = find(S(:,i)==1);
    sx0(i,:) = real(S0(idx_s0));
    sy0(i,:) = imag(S0(idx_s0));
    sx1(i,:) = real(S0(idx_s1));
    sy1(i,:) = imag(S0(idx_s1));
end

L = [];
for b=1:numinfobits/2
    x = real(rxSig(b));
    y = imag(rxSig(b));

    for i=1:2
%     LLR1 = log( exp(-( (x-sx01).^2+(y-sy01).^2 ))/exp(-( (x-sx11).^2+(y-sy11).^2 )) );
%     LLR2 = log( exp(-( (x-sx02).^2+(y-sy02).^2 ))/exp(-( (x-sx12).^2+(y-sy12).^2 )) );
%     LLR1 =  -(min( (x-sx01).^2+(y-sy01).^2 ) - min( (x-sx11).^2+(y-sy11).^2 ));
%     LLR2 =  -(min( (x-sx02).^2+(y-sy02).^2 ) - min( (x-sx12).^2+(y-sy12).^2 ));
      LLR = -(min( (x-sx0(i,:)).^2+(y-sy0(i,:)).^2 ) - min( (x-sx1(i,:)).^2+(y-sy1(i,:)).^2 ));
%       return
      L = [L LLR];
    end
    
%     L = [L LLR1 LLR2];
end  
%%%%%%
% scatterplot(mod_sig)
figure(1)
plot(normalize(L), 'LineWidth',2)
hold on
plot(normalize(demod_sig))
legend('my llr', 'llr')
grid on
title('QPSK LLR')