
function [RESULT,ERROR]=CRC_CODEC(DATA_CODEC,CRC_Type,ENCODER)


%% Тип CRC
if(CRC_Type == 8)
    crc_poly = [1,1,0,0,1,1,0,1,1];                                 crc_len  = 8;
elseif(CRC_Type == 16)
    crc_poly = [1,0,0,0,1,0,0,0,0,0,0,1,0,0,0,0,1];                 crc_len  = 16;
elseif(CRC_Type == '24A')
    crc_poly = [1,1,0,0,0,0,1,1,0,0,1,0,0,1,1,0,0,1,1,1,1,1,0,1,1]; crc_len  = 24;
elseif(CRC_Type == '24B')
    crc_poly = [1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,1,1]; crc_len  = 24; 
end;


crc_rem = zeros(1,crc_len+1);
if (ENCODER==1)  tmp_array = [DATA_CODEC, zeros(1,crc_len)]; else tmp_array=DATA_CODEC; end;

%% Вычисление CRC
for n=1:length(tmp_array)
    for m=1:crc_len crc_rem(m) = crc_rem(m+1); end;
    crc_rem(crc_len+1) = tmp_array(n);

    if(crc_rem(1) ~= 0)
        for(m=1:crc_len+1) crc_rem(m) = mod(crc_rem(m)+crc_poly(m), 2); end; 
    end; 
end;
       
   CRC_Bits=crc_rem(2:end);
   
   if (ENCODER==1) RESULT=[DATA_CODEC, CRC_Bits]; ERROR=0; 
   else
       if ((ENCODER==0)&&(bi2de(CRC_Bits))~=0)  RESULT=DATA_CODEC(1:(length(tmp_array)-crc_len));  ERROR=1; 
         
   else
       if ((ENCODER==0)&&(bi2de(CRC_Bits))==0) RESULT=DATA_CODEC(1:(length(tmp_array)-crc_len));  ERROR=0;
       
       end;
   end;
   end;
end
