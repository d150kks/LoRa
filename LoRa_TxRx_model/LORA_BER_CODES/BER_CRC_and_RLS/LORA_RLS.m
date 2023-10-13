function [FilterOutput, w] = LORA_RLS(N, DesiredSignal, TrainSignal, InputSignal);

%% =============================================== Входные аргументы функции
%N              - % Порядок фильтра
%lambda         - % коэффициент забывания 0< lambda <1
%DesiredSignal  - % Опорный сигнал
%InputSignal    - % Входной искаженный сигнал
%NumIter        - % Количество итераций
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
%% =============================================== Выходные аргументы функции
%FilterOutput   - % Выходной сигнал фильтра
%ErrorSignal    - % Сигнал ошибки фильтра
%w              - % Вектор коэффициентов
%Error_Plane    - % Сигнал ошибки фильтра
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %





%% =============================================== Adaptive Filter
% Параметры фильтра
lambda = 1;
laminv  = 1/lambda;                     % forgetting factor
delta   = 1000.0/(std(TrainSignal).^2); % initialization parameter
K       = zeros(N,1);                   % vector
I       = eye(N);
P       = delta*I;                 % inverse correlation matrix
w       = zeros(N,1);                   % filter coefficients
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

% Входные сигналы                               
delay_start = ceil((N-1)/2);  % Задержка сдвигового регистра
delay_end   = floor((N-1)/2); % Задержка сдвигового регистра

d = DesiredSignal.';         % Задержанный опорный сигнал
x = [zeros(1,delay_start), TrainSignal, zeros(1,delay_end)].';       % Входной искаженный тренировочный сигнал
nIter = length(d);
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

% Регистры фильтра
% e = 0;                          % Регистр вектора ошибок
FilterOutput = zeros(1, nIter); % Регистр вектора выходного сигнала фильтра
% ErrorSignal  = zeros(1, nIter); % Регистр вектора выходного сигнала фильтра
% ErrorSquare  = zeros(1, nIter);   % Регистр квадрата вектора ошибок

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

% Вычисление коэффициентов фильтра
for i = 1:nIter
    u = x(i+N-1:-1:i);
    e = (d(i) - w'*u);

%     % Update inverse correlation matrix
%     P = (1/lambda) * (P - (P * u * u' * P) / (lambda + u' * P * u));
%     
%     % Update equalizer weights
%     w = w + P * u * e;

    K = (P*u)/(lambda + u'*P*u); 
    P = laminv * (P - K*u'*P);

    w = w + K*conj(e);
end

%% Filtration
xx = [zeros(1,delay_start), InputSignal, zeros(1,delay_end)].';       % Входной искаженный тренировочный сигнал
for i = 1:length(InputSignal)
    u = xx(i+N-1:-1:i);
    FilterOutput(i) = w'*u;
end
% for i = 1:nIter
%     up = x(i+N-1:-1:i);
%     FilterPreamb(i) = w'*up;
% end
% FilterOutput = filter(w, 1, InputSignal); %%%%%%%%%%%%%%%%%
% FilterPreamb = conv(TrainSignal, conj(w), 'same'); %%%%%%%%%%%%%%%%%


% Output Signals
% ErrorSignal = DesiredSignal-FilterOutput;
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
end