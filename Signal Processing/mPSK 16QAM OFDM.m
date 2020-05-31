#Apostolos Fotopoulos 2202
#Eb=1

clear;
close all;
pkg load signal;

%Inputs
N = 10000;                %Number  of bits.
Rb = 1000;                %Bit rate (bits/sec).
Oversampling = 1;         %Number of samples after Tx Filter.
SNR = [0 2 4 6 8 10 12];  %Signal to noise ratio (db).
BitsPerSym = 1;           %Bits transmitted per symbol. Possible values: 1 2 3 or 4.
enableQAM = 1;            %Choose whether to enable the 16-QAM Modulation when transmitting 4 bits per symbol.
Subcarriers = 64;         %Number of subcarriers.
CyclicPrefix = 16;        %Length of cyclic prefix.
enableLTI = 0;            %Choose whether to enable the LTI chanel.
LTIvals = [complex(0.9,0.9),complex(0.6,0.6),complex(0.3,0.3)];

encodedData = randi([0,1],1,N);
%encodedData = [0 0 0 1 1 0 1 1 1 1 0 1 0 1 0 1 0 1 1 1];

if (enableLTI==1)
  CyclicPrefix = 4; 
endif

constellation = [];
equal = [];
%%%%%%%%%%Constellation%%%%%%%%%%%%
if (enableQAM==1 && BitsPerSym==4)
  a=sqrt(16/40);   %Opou a to platos tou 1ou shmeiou tou eswterikou kyklou tou 16QAM
  
  k=1;
  m=-3*a;
  for mm = 1:4
    n=-3*a;
    for nn = 1:4
      constellation(k) = complex(m,n);
      equal = [equal,dec2bin(k-1,BitsPerSym)-'0'];
      n=n+2*a;
      k=k+1;
    endfor
    m=m+2*a;
  endfor
else
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
endif

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

SNRLin = [];
ber = [];
packets = fix(size(Isfilt,2)/Subcarriers);
errors = zeros(1,size(SNR,2));
Iy = [];
Qy = [];
if (rem(size(Isfilt,2),packets)!=0)
  Iy = zeros(columns(SNR),packets,Subcarriers);
  Qy = zeros(columns(SNR),packets,Subcarriers);
else
  Iy = zeros(columns(SNR),packets,Subcarriers);
  Qy = zeros(columns(SNR),packets,Subcarriers);
endif

%%%%%%%%%%Subcarriers%%%%%%%%%%%%
for n = 1:packets
  IyPart = [];
  QyPart = [];
  
  for k = 1:columns(SNR)
    part=[];
    part = complexOs(1+(n-1)*(fix(size(Isfilt,2)/packets)):n*(fix(size(Isfilt,2)/packets)));
    
    IyPart(k,:)=real(part);
    QyPart(k,:)=imag(part);
  endfor
  Iy(:,n,:)=IyPart;
  Qy(:,n,:)=QyPart;
endfor
complexYsub = complex(Iy,Qy);

%%%%%%%%%%IDFT%%%%%%%%%%%%
complexYifft = zeros(columns(SNR),packets,Subcarriers);
for k = 1:columns(SNR)
  for n = 1:packets
    if (enableQAM==1 && BitsPerSym==4)
      complexYifft(k,n,1:Subcarriers) = ifft(fftshift(complexYsub(k,n,:)));
    else
      complexYifft(k,n,1:Subcarriers) = (Subcarriers/(sqrt(Subcarriers)))*ifft(fftshift(complexYsub(k,n,:)));
    endif
  endfor
endfor


%%%%%%%%%%Add Cyclic Prefix%%%%%%%%%%%%
complexYcp = complexYifft;
for k = 1:CyclicPrefix
  complexYcp(:,:,Subcarriers+k) = complexYcp(:,:,k);
endfor

complexYplot = [];
for n = 1:packets
  complexYpart = [];
  for k = 1:columns(SNR)
    part = [];
    part = complexYcp(k,n,:);
    
    complexYpart(k,:)=part;
  endfor
  complexYplot = [complexYplot,complexYpart];
endfor

Fs=Rb;
Ts=1/Fs;
t=0:Ts:(N-1)*Ts;
Rxx = xcorr(complexYplot(7,:),'biased');
Rxxdft = abs(fft(Rxx));
%Rxxdft = abs(fftshift(fft(Rxx)));
freq = -Fs/2:Fs/length(Rxx):Fs/2-(Fs/length(Rxx));
figure(1)
plot(freq,Rxxdft);

complexYlti = [];
%%%%%%%%%%LTI Channel%%%%%%%%%%%%
if (enableLTI==1)
  for n = 1:packets
    for k = 1:columns(SNR)
      complexYlti(k,n,:) = conv(squeeze(complexYcp(k,n,:)),LTIvals);
    endfor
  endfor
endif

%%%%%%%%%%Transmission%%%%%%%%%%%%
for k = 1:columns(SNR)
  SNRLin(k) = 10^(SNR(k)/10);
  Variance = 1/(2*(SNRLin(k)));
  Inoise = [];
  Qnoise = [];  
  noisyY=[];
  
  if (enableLTI==1)
    Inoise = sqrt(Variance)*randn(1,packets,fix(size(complexYcp,3)));
    Qnoise = sqrt(Variance)*randn(1,packets,fix(size(complexYcp,3))); 
    noise = complex(Inoise,Qnoise);
    
    complexY(k,:,:) = complexYlti(k,:,1:Subcarriers+CyclicPrefix) + noise;
  else
    Inoise = sqrt(Variance)*randn(1,packets,fix(size(complexYcp,3)));
    Qnoise = sqrt(Variance)*randn(1,packets,fix(size(complexYcp,3))); 
    noise = complex(Inoise,Qnoise);
    
    complexY(k,:,:) = complexYcp(k,:,:) + noise;
  endif 
endfor

%%%%%%%%%%Remove Cyclic Prefix%%%%%%%%%%%%
complexYrecv = complexY;
for k = 1:CyclicPrefix
  complexYrecv(:,:,Subcarriers+CyclicPrefix-k+1) = [];
endfor

%%%%%%%%%%DFT%%%%%%%%%%%%
complexYfft = zeros(columns(SNR),packets,Subcarriers);
for k = 1:columns(SNR)
  for n = 1:packets
    complexYfft(k,n,1:Subcarriers) = fftshift(fft(complexYrecv(k,n,:)));
  endfor
endfor

complexYps = [];
%%%%%%%%%%Subcarriers%%%%%%%%%%%%
for n = 1:packets
  complexYpart = [];
  for k = 1:columns(SNR)
    part = [];
    part = complexYfft(k,n,:);
    
    complexYpart(k,:)=part;
  endfor
  complexYps = [complexYps,complexYpart];
endfor
IyRecv = real(complexYps);
QyRecv = imag(complexYps);

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

figure(2);
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
figure(3);
if (BitsPerSym<3)
  semilogy(SNR, ber, 'bx', SNR, errTheo, 'r', 'MarkerSize', 10, 'LineWidth', 2);
  xlabel("SNR");
  ylabel("BER");  
  title("Ber vs SNR");
  legend("Simulation","AWGN Theory");
  grid on
else
  semilogy(SNR, ber,"bx");
  xlabel("SNR");
  ylabel("BER");
  grid on
  title("Ber vs SNR");   
endif
%}


