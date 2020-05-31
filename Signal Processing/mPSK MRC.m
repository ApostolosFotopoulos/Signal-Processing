#Apostolos Fotopoulos 2202
#Eb=1

clear;
pkg load signal;

%Inputs
N = 10000;                %Number  of bits.
Rb = 1000;                %Bit rate (bits/sec).
Oversampling = 5;         %Number of samples after Tx Filter.
SNR = [0 2 4 6 8 10 12];  %Signal to noise ratio (db).
BitsPerSym = 1;           %Bits transmitted per symbol. Possible values: 1 2 3 or 4.
T=1;                      %How often h changes value (secs).
c=1;                      %A real number that gets added to h.
enableC=0;                %Choose whether to enable c or not.
antennas=1;               %Ammount of antennas that the receiver has.

encodedData = randi([0,1],1,N);
%encodedData = [0 0 0 1 1 0 0 0 1 1 0 1 0 1 0 1 0 1 0 1];

constellation = [];
equal = [];
%%%%%%%%%%Constellation%%%%%%%%%%%%
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

Is = [];
Qs = [];
%%%%%%%%%%Modulation%%%%%%%%%%%%
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
%%%%%%%%%%Oversampling%%%%%%%%%%%%
for i = 1:columns(Is)
  Isfilt((i-1)*Oversampling+1:Oversampling*i) = Is(i)/sqrt(Oversampling);
  Qsfilt((i-1)*Oversampling+1:Oversampling*i) = Qs(i)/sqrt(Oversampling);
endfor
constellation=constellation/sqrt(Oversampling);

complexOs = complex(Isfilt,Qsfilt);
tfilt = linspace(0,(N-1)*1/Rb,size(Isfilt,2));

Fs = Rb*Oversampling;
Rxx = xcorr(complexOs,'biased');
Rxxdft = abs(fft(Rxx));
freq = -Fs/2:Fs/length(Rxx):Fs/2-(Fs/length(Rxx));
figure(3);
plot(freq,Rxxdft);
title("x - Transmit spectrum");

SNRLin = [];
ber = [];
errors = zeros(1,size(SNR,2));
hValues=N/(T*Rb);
h=[];
Iy = [];
if (rem(size(Isfilt,2),hValues)!=0)
  Iy = zeros(columns(SNR),size(Isfilt,2)-(size(Isfilt,2)-fix(size(Isfilt,2)/hValues)*hValues),antennas);
  Qy = zeros(columns(SNR),size(Isfilt,2)-(size(Isfilt,2)-fix(size(Isfilt,2)/hValues)*hValues),antennas);
else
  Iy = zeros(columns(SNR),size(Isfilt,2),antennas);
  Qy = zeros(columns(SNR),size(Isfilt,2),antennas);
endif

%%%%%%%%%%Transmission%%%%%%%%%%%%
for m = 1:antennas
  IyPart2 = [];
  QyPart2 = [];
  n=1;
  r=1;
  for n = 1:hValues
    h(m,n) = complex(sqrt(1/2)*randn(1,1),sqrt(1/2)*randn(1,1));
    %h(m,n) = 2*m;
    %h(m,n) = complex(sqrt(2),sqrt(2));
    IyPart = [];
    QyPart = [];
    
    for k = 1:columns(SNR)
      SNRLin(k) = 10^(SNR(k)/10);
      Variance = 1/(2*(SNRLin(k)));
      Inoise = [];
      Qnoise = [];
      Inoise = sqrt(Variance)*randn(1,fix(size(Isfilt,2)/hValues));
      Qnoise = sqrt(Variance)*randn(1,fix(size(Qsfilt,2)/hValues));
      
      noise = complex(Inoise,Qnoise);
      noisyY = complexOs(1+(n-1)*(fix(size(Isfilt,2)/hValues)):n*(fix(size(Isfilt,2)/hValues)))*h(m,n) + noise;
      
      IyPart(k,:)=real(noisyY);
      QyPart(k,:)=imag(noisyY);
    endfor
    IyPart2 = [IyPart2,IyPart];
    QyPart2 = [QyPart2,QyPart];  
  endfor
  Iy(:,:,m)=IyPart2;
  Qy(:,:,m)=QyPart2;
endfor
complexY = complex(Iy,Qy);

%{
figure(5);
plot(complexY(1,:));
title("y");
grid on
%}

Fs = Rb*Oversampling;
Rxx = xcorr(complexY(1,:),'biased');
Rxxdft = abs(fft(Rxx));
freq = -Fs/2:Fs/length(Rxx):Fs/2-(Fs/length(Rxx));
figure(4);
plot(freq,Rxxdft);
title("y - Transmit spectrum");

%%%%%%%%%%Approximation of c%%%%%%%%%%%%
if (enableC==1)
  hNew = h+c;
  cApprox = real(sum(hNew,2)/size(h,2));
  hApprox = hNew - cApprox;
  %hApprox = hNew;   
endif

IyRecv = [];
QyRecv = [];

%%%%%%%%%%Channel Inversion%%%%%%%%%%%%
for n = 1:hValues
  IyRecvPart = [];
  QyRecvPart = [];
  for k = 1:columns(SNR)
    complexYpart = zeros(antennas,fix(size(complexY,2)/(hValues)));
    if (enableC==1)
      for m = 1:antennas
        complexYpart(m,:) = complexY(k,1+(n-1)*(fix(size(complexY,2)/(hValues))):n*(fix(size(complexY,2)/(hValues))),m);
      endfor
      
      complexYrecv = (ctranspose(hApprox(:,n))*complexYpart)/(ctranspose(hApprox(:,n))*hApprox(:,n));
    else
      for m = 1:antennas
        complexYpart(m,:) = complexY(k,1+(n-1)*(fix(size(complexY,2)/(hValues))):n*(fix(size(complexY,2)/(hValues))),m);
      endfor
      
      complexYrecv = (ctranspose(h(:,n))*complexYpart)/(ctranspose(h(:,n))*h(:,n));
    endif
    
    IyRecvPart(k,:) = real(complexYrecv);
    QyRecvPart(k,:) = imag(complexYrecv); 
  endfor
  
  IyRecv = [IyRecv,IyRecvPart];
  QyRecv = [QyRecv,QyRecvPart];
endfor

%{
for n = 1:(hValues)
  IyRecvPart = [];
  QyRecvPart = [];
  for k = 1:columns(SNR)
    if (enableC==1)
      complexYrecv = (complexY(k,1+(n-1)*(size(complexY,2)/(hValues)):n*(size(complexY,2)/(hValues)))*conj(hApprox(n)))/(hApprox(n)*conj(hApprox(n)));
    else
      complexYrecv = (complexY(k,1+(n-1)*(size(complexY,2)/(hValues)):n*(size(complexY,2)/(hValues)))*conj(h(m,n)))/(h(m,n)*conj(h(m,n)));
    endif
    
    IyRecvPart(k,:) = real(complexYrecv);
    QyRecvPart(k,:) = imag(complexYrecv); 
  endfor
  
  IyRecv = [IyRecv,IyRecvPart];
  QyRecv = [QyRecv,QyRecvPart];
endfor
%}

%%%%%%%%%%Demodulation%%%%%%%%%%%%
for k = 1:columns(SNR)
  Iyfilt = [];
  Qyfilt = [];
  IyConv = [];
  QyConv = [];

  for j = 1:(fix(size(IyRecv,2)/Oversampling))
    IyConv = zeros(1,Oversampling) + abs(Is(j))/sqrt(Oversampling);
    QyConv = zeros(1,Oversampling) + abs(Qs(j))/sqrt(Oversampling);
    IconvBit = IyRecv(k,(j-1)*Oversampling+1:Oversampling*j);
    IconvRes = conv(IconvBit,IyConv);
    QconvBit = QyRecv(k,(j-1)*Oversampling+1:Oversampling*j);
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
  for j = 1:(fix(size(IyRecv,2)/Oversampling))
    IySamp(j) = Iyfilt(Oversampling*j);
    QySamp(j) = Qyfilt(Oversampling*j);
  endfor

  complexNum = complex(IySamp,QySamp);
  
  ipHat = [];
  for n = 1:(fix(size(IyRecv,2)/Oversampling))
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
subplot(4,4,1);
stem(encodedData(BitsPerSym+1:BitsPerSym+12),"r");
title("encodedData[t]");

subplot(4,4,2);
stem(t(1:12/BitsPerSym),Is(2:12/BitsPerSym + 1),"r");
title("I - s[t]");

subplot(4,4,3);
stem(t(1:12/BitsPerSym),Qs(2:12/BitsPerSym + 1),"r");
title("Q - s[t]");

subplot(4,4,4);
stem(tfilt(1:12*Oversampling/BitsPerSym),Isfilt(Oversampling+1:12*Oversampling/BitsPerSym + Oversampling),"r");
title("I - sfilt[t]");

subplot(4,4,5);
stem(tfilt(1:12*Oversampling/BitsPerSym),Qsfilt(Oversampling+1:12*Oversampling/BitsPerSym + Oversampling),"r");
title("Q - sfilt[t]");

subplot(4,4,6);
stem(tfilt(1:12*Oversampling/BitsPerSym),Iy(k,Oversampling+1:12*Oversampling/BitsPerSym + Oversampling),"r");
title("I - y[t]");

subplot(4,4,7);
stem(tfilt(1:12*Oversampling/BitsPerSym),Qy(k,Oversampling+1:12*Oversampling/BitsPerSym + Oversampling),"r");
title("Q - y[t]");

subplot(4,4,8);
stem(tfilt(1:12*Oversampling/BitsPerSym),IyRecv(k,Oversampling+1:12*Oversampling/BitsPerSym + Oversampling),"r");
title("I - yRecv[t]");

subplot(4,4,9);
stem(tfilt(1:12*Oversampling/BitsPerSym),QyRecv(k,Oversampling+1:12*Oversampling/BitsPerSym + Oversampling),"r");
title("Q - yRecv[t]");

subplot(4,4,10);
stem(tfilt(1:12*Oversampling/BitsPerSym),Iyfilt(Oversampling+fix(Oversampling/2)+1:12*Oversampling/BitsPerSym + Oversampling + fix(Oversampling/2)),"r");
title("I - yFilt[t]");

subplot(4,4,11);
stem(tfilt(1:12*Oversampling/BitsPerSym),Qyfilt(Oversampling+fix(Oversampling/2)+1:12*Oversampling/BitsPerSym + Oversampling + fix(Oversampling/2)),"r");
title("Q - yFilt[t]");

subplot(4,4,12);
stem(IySamp(2:12/BitsPerSym + 1),"r");
title("I - ySamp[t]");

subplot(4,4,13);
stem(QySamp(2:12/BitsPerSym + 1),"r");
title("Q - ySamp[t]"); 

subplot(4,4,14);
plot(real(constellation),imag(constellation),"*r");
title("constellation[t]");

subplot(4,4,15);
stem(ipHat(BitsPerSym+1:BitsPerSym+12),"r");
title("ipHat[t]");

subplot(4,4,16);
stem(SNR,errors,"r");
title("Errors");

errTheo = 0.5*erfc(sqrt(SNRLin));
figure(2);
if (BitsPerSym<3)
  semilogy(SNR, ber, 'bx', SNR, errTheo, 'r', 'MarkerSize', 10, 'LineWidth', 2);
  xlabel("SNR");
  ylabel("BER");  
  title("Ber vs SNR");
  legend("Simulation","AWGN Theory");
  grid on
else
  semilogy(SNR, ber,"bx");
  xlabel('Eb/No (dB)');
  ylabel('Bit error rate');
  grid on
  title("Ber vs SNR");  
endif


