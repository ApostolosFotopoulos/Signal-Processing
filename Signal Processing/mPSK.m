#Apostolos Fotopoulos 2202
#Eb=1

clear;
pkg load signal;

%Inputs
N = 10000;                %Number  of bits.
Rb = 1000;                %Bit rate (bits/sec).
Oversampling = 5;
SNR = [0 2 4 6 8 10];     %Signal to noise ratio (db).
BitsPerSym = 2;           %Bits transmitted per symbol. Possible values: 1 2 3 or 4.

encodedData = randi([0,1],1,N);
%encodedData = [0 0 0 1 1 0 0 0 1 1 0 1 0 1 0 1 0 1 0 1];

constellation = [];
equal = [];
for k = 1:(2^BitsPerSym)
  if (BitsPerSym==1)
    constellation(k) = exp((k-1)*i*pi/2^(BitsPerSym-1));
  else
    constellation(k) = sqrt((BitsPerSym*1)/2)*exp((k-1)*i*pi/2^(BitsPerSym-1));
  endif
  
  if (abs(real(constellation(k)))<1e-10)
    constellation(k) = 1i*imag(constellation(k));
  endif
  if (abs(imag(constellation(k)))<1e-10)
    constellation(k) = real(constellation(k));
  endif
  
  equal = [equal,dec2bin(k-1,BitsPerSym)-'0'];
endfor

tst = constellation
Is = [];
Qs = [];
for k = 1:(size(encodedData,2)/BitsPerSym)
  for j = 1:(size(equal,2)/BitsPerSym)
    if (encodedData(BitsPerSym*k-(BitsPerSym-1):BitsPerSym*k) == equal(BitsPerSym*j-(BitsPerSym-1):BitsPerSym*j))
      Is(k) = real(constellation(j));
      Qs(k) = imag(constellation(j));
    endif
  endfor
endfor

t = linspace(0,(N-1)*1/Rb,N/BitsPerSym);

Isfilt = [];
Qsfilt = [];
for i = 1:columns(Is)
  Isfilt((i-1)*Oversampling+1:Oversampling*i) = Is(i)/sqrt(Oversampling);
  Qsfilt((i-1)*Oversampling+1:Oversampling*i) = Qs(i)/sqrt(Oversampling);
endfor
constellation=constellation/sqrt(Oversampling);

complexOs = complex(Isfilt,Qsfilt);
tfilt = linspace(0,(N-1)*1/Rb,size(Isfilt,2));

SNRLin = [];
ber = [];
errors = zeros(1,size(SNR,2));
for k = 1:columns(SNR)
  SNRLin(k) = 10^(SNR(k)/10)
  Variance = 1/(2*(SNRLin(k)));
  Inoise = [];
  Qnoise = [];
  Inoise = sqrt(Variance)*randn(1,size(Isfilt,2));
  Qnoise = sqrt(Variance)*randn(1,size(Qsfilt,2));

  tst= Isfilt;
  tst = Qsfilt;
  noise = complex(Inoise,Qnoise);
  noisyY = complexOs + noise;
  
  Iy = [];
  Qy = [];
  Iy = real(noisyY);
  Qy = imag(noisyY);

  Iyfilt = [];
  Qyfilt = [];
  IyConv = [];
  QyConv = [];

  for j = 1:(size(encodedData,2)/BitsPerSym)
    IyConv = zeros(1,Oversampling) + abs(Is(j))/sqrt(Oversampling);
    QyConv = zeros(1,Oversampling) + abs(Qs(j))/sqrt(Oversampling);
    IconvBit = Iy((j-1)*Oversampling+1:Oversampling*j);
    IconvRes = conv(IconvBit,IyConv);
    QconvBit = Qy((j-1)*Oversampling+1:Oversampling*j);
    QconvRes = conv(QconvBit,QyConv);
    
    if (j==1)
      Iyfilt = IconvRes;
      Qyfilt = QconvRes;
    else
      Iyfilt((Oversampling*(j-1))+1:(Oversampling*(j-1)+Oversampling-1)) = Iyfilt((Oversampling*(j-1))+1:(Oversampling*(j-1)+Oversampling-1)) +  IconvRes(1:Oversampling-1);
      Iyfilt = [Iyfilt,IconvRes(Oversampling:(Oversampling*2-1))];
      Qyfilt((Oversampling*(j-1))+1:(Oversampling*(j-1)+Oversampling-1)) = Qyfilt((Oversampling*(j-1))+1:(Oversampling*(j-1)+Oversampling-1)) +  QconvRes(1:Oversampling-1);
      Qyfilt = [Qyfilt,QconvRes(Oversampling:(Oversampling*2-1))];
    endif
  endfor
    
  tConv = linspace(0,(N-1)*1/Rb,size(Iyfilt,2));
  
  IySamp = [];
  QySamp = [];
  for j = 1:(size(encodedData,2)/BitsPerSym)
    IySamp(j) = Iyfilt(Oversampling*j);
    QySamp(j) = Qyfilt(Oversampling*j);
  endfor

  complexNum = complex(IySamp,QySamp);
  
  ipHat = [];
  for n = 1:(size(encodedData,2)/BitsPerSym)
    des = 0;
    prev = 100000;
    for j = 1:columns(constellation)
      d = norm(complexNum(n) - constellation(j));
      if (d<prev)
        des = j;
        prev = d;
      endif
    endfor
    if (BitsPerSym == 1)
      ipHat(n) = equal(des); 
    else
      ipHat = [ipHat,equal(des*BitsPerSym-(BitsPerSym-1):des*BitsPerSym)];
    endif
  endfor
  
  for j = 1:size(ipHat,2)
    if (encodedData(j) != ipHat(j))
      errors(k) =  errors(k) + 1;
    endif  
  endfor
  if (errors(k)==0)
    errors(k) =  errors(k) + 1e-10;
  endif
  ber(k) = errors(k)/size(encodedData,2);
endfor

figure(1);
subplot(3,4,1);
stem(encodedData(BitsPerSym+1:BitsPerSym+12),"r");
title("encodedData[t]");

subplot(3,4,2);
stem(t(1:12/BitsPerSym),Is(2:12/BitsPerSym + 1),"r");
title("I - s[t]");

subplot(3,4,3);
stem(t(1:12/BitsPerSym),Qs(2:12/BitsPerSym + 1),"r");
title("Q - s[t]");

subplot(3,4,4);
stem(tfilt(1:12*Oversampling/BitsPerSym),Isfilt(Oversampling+1:12*Oversampling/BitsPerSym + Oversampling),"r");
title("I - sfilt[t]");

subplot(3,4,5);
stem(tfilt(1:12*Oversampling/BitsPerSym),Qsfilt(Oversampling+1:12*Oversampling/BitsPerSym + Oversampling),"r");
title("Q - sfilt[t]");

subplot(3,4,6);
stem(tfilt(1:12*Oversampling/BitsPerSym),Iy(Oversampling+1:12*Oversampling/BitsPerSym + Oversampling),"r");
title("I - y[t]");

subplot(3,4,7);
stem(tfilt(1:12*Oversampling/BitsPerSym),Qy(Oversampling+1:12*Oversampling/BitsPerSym + Oversampling),"r");
title("Q - y[t]");

subplot(3,4,8);
stem(tfilt(1:12*Oversampling/BitsPerSym),Iyfilt(Oversampling+fix(Oversampling/2)+1:12*Oversampling/BitsPerSym + Oversampling + fix(Oversampling/2)),"r");
title("I - yFilt[t]");

subplot(3,4,9);
stem(IySamp(2:12/BitsPerSym + 1),"r");
title("I - ySamp[t]");

subplot(3,4,10);
stem(QySamp(2:12/BitsPerSym + 1),"r");
title("Q - ySamp[t]"); 

subplot(3,4,11);
plot(real(constellation),imag(constellation),"*r");
title("constellation[t]");

subplot(3,4,12);
stem(ipHat(BitsPerSym+1:BitsPerSym+12),"r");
title("ipHat[t]");

%{
figure(2);
stem(tConv,Qyfilt,"r");
title("Q - yFilt[t]");
%}

%{
figure(2);
stem(errors,"r");
title("Errors");
%}

errTheo = 0.5*erfc(sqrt(SNRLin));
figure(2);
if (BitsPerSym<3)
  semilogy(SNR, ber, 'bx', SNR, errTheo, 'r', 'MarkerSize', 10, 'LineWidth', 2);
  xlabel("SNR");
  ylabel("BER");  
  title("Ber vs SNR");
  legend("Simulation","Theory");
  grid on
else
  semilogy(SNR, ber,"bx");
  xlabel('Eb/No (dB)');
  ylabel('Bit error rate');
  grid on
  title("Ber vs SNR");  
endif

%{
figure(5);
stem(ber,"r");
title("BER");
figure(6);
stem(errTheo,"r");
title("BER THEO");
%}




