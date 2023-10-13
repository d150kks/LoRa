clc
clear all
close all

N = 100;    
User_id = 1;
N_users = 2;
k=2;
bits1 = randi([0, 1], 1, N);
bits2 = randi([0, 1], 1, N);

mod_bits1 = 2*bits1-1;
mod_bits2 = 2*bits2-1;

mtx=hadamard(N_users);
% ##Произведите кодирование сформированных данных, используя функции Уолша.
j=1;
for i=1:N
   if mod_bits1 (i) == 1
       mod_bits_1_code (j:j+k-1)= mtx(1,:);
   else
       mod_bits_1_code (j:j+k-1)= mtx(1,:)*(-1);
   end

   if mod_bits2 (i) == 1
       mod_bits_2_code (j:j+k-1)= mtx(2,:);
   else
       mod_bits_2_code (j:j+k-1)= mtx(2,:)*(-1);
   end

  j=j+k;  
end

figure(1)
plot(mod_bits_1_code)
hold on
plot(mod_bits_2_code)
return
sum_of_signals = mod_bits_1_code+mod_bits_2_code+mod_bits_3_code+mod_bits_4_code;
sum_of_signals_noise = awgn(sum_of_signals,50,'measured');

B_long = zeros(1,length(mod_bits_1_code));
B_long(1:k)=mtx(User_id,:); 
korr = ifft(fft(sum_of_signals_noise).*conj(fft(B_long)));

j=1;
for i=1:k:N*k
   korr_decim(j) =  korr(i)/k;
   j=j+1;
end   

for i=1:N
   if real(korr_decim(i)) < 0;
     demod(i) = 0;
   else
     demod(i) = 1;
   end
end


figure
subplot(3,2,1)
plot(mod_bits1)
% ##title('1 user');

subplot(3,2,2)
plot(mod_bits2)
% ##title('2 user');

subplot(3,2,3)
plot(mod_bits3)
% ##title('3 user');

subplot(3,2,4)
plot(mod_bits4)
% ##title('4 user');

subplot(3,2,5)
plot(sum_of_signals_noise)
% ##title('sum signal');

subplot(3,2,6)
plot(real(korr_decim))

num_err = biterr(demod,bits3);