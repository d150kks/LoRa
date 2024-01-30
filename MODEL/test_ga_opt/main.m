clc
clear all
close all


fun = @(x) MAIN_GA(x(1));
% fun = ga( @(x) MAIN_GA(x(1), x(2), x(3)));
% fun = @(x) x.^2;
% fun = @ps_example;

% PARAM_ALPHA, PARAM_WINSIZE, PARAM_NVAR_LIMIT
A = [];
b = [];
Aeq = [];
beq = [];
lb = [0.1];
ub = [10];
nonlcon = [];
% IntCon = 2;

options = optimoptions('ga','ConstraintTolerance',1e-6, ...
                    'CrossoverFraction', 0.8, ...
                    'MaxGenerations', 300, ...
                    'MigrationInterval', 20, ...
                    'PlotFcn', @gaplotbestf, ...
                    'PopulationSize', 50);
% options.CrossoverFraction = 0.8; % {0.8}
% options.MaxGenerations = 300; % {100*numberOfVariables}
% options.MigrationInterval = 20; % default 20
% options.PlotFcn = @gaplotbestf; % plots the best score value and mean score versus generation.
% options.PopulationSize = 50; %min(max(10*nvars,40),100)



% [x,fval,exitflag,output,population,scores] = ga(fun, 1, A, b, Aeq, beq, lb, ub, nonlcon, IntCon, options);
[x,fval,exitflag,output,population,scores] = ga(fun, 1, A, b, Aeq, beq, lb, ub, nonlcon, options);
% 2.7346    7.0000   10.7295 % BAD
% 2.8334    6.0000   14.6313
% 2.5
best_pars = x;


function BER = MAIN_GA(PARAM_ALPHA)
% fprintf('%.2f\n', PARAM_ALPHA)
% fprintf('%.2f\n', PARAM_WINSIZE)
% fprintf('%.2f\n', PARAM_NVAR_LIMIT)
%% ================================= Переменные
% коэффициенты
SF = 7;        % коэффициент расширения спектра (от 7 до 12)
rc_size = 4;
rc = (SF-rc_size);
BW = 2e6;
nIter = 50;

LORA = myLoRaClass_RSG_GA(SF, BW, PARAM_ALPHA); 
Base = LORA.Base;
downch = LORA.downch;
chirp = LORA.chirp;
% Ts = LORA.Ts;
num_pre = 4;


%% ================================= Data bit
% Number of message bits
numinfobits = 1200; 

% Message/Iformation bits
data = randi([0 1],1, numinfobits); 
num_sym = length(data)/rc;

%% ================================= Mодуляция
% RS
[mod_chirp, ~, ~] = LORA.lorax_modified_crcrs(data, num_sym);
tx_preamble = repmat(chirp,1,num_pre);
tx_downch = repmat(downch,1,num_pre);

tx_chirp = [ tx_downch, tx_preamble, mod_chirp];
tx_length = length(tx_chirp);

%% ================================= CHANNEL

% Channel
h11 = load('h.mat').h;
h11 = [h11 zeros(1, tx_length-length(h11))];
H11 = fft(h11);
tx_chirp_h = ifft( fft(tx_chirp).*H11 );


%% ================================= BER
snr = -8;
[numErr, NumData] = deal(0);
for iter = 1:nIter
    
    % АБГШ 
    rxSig = awgn(tx_chirp_h, snr, 'measured');

    % Freq Sync
    rx_preamble = rxSig(num_pre*Base+1:num_pre*2*Base);
    corrected_signal = rxSig(num_pre*2*Base+1:end);

    % Demodulation
    % RS
    [~, hard_bits, ~, ~, ~, ~] = LORA.delorax_crcrs( corrected_signal, num_sym, tx_preamble, rx_preamble);

    % подсчет БЕР с учетом задержки
    err = sum(hard_bits~=data);

    % Increment the error and bit counters
    numErr = numErr + err;        
    NumData = NumData + numinfobits;
end

% Estimate the BER for both methods
BER = numErr/NumData;

end