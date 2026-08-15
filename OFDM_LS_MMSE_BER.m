% OFDM_LS_MMSE_BER
% OFDM system with LS channel estimation.
% Compare ZF and MMSE equalization under different SNR conditions.
%
% Processing chain:
% Random Bits
% -> 16-QAM
% -> Pilot Insertion
% -> IFFT
% -> Cyclic Prefix
% -> Multipath Channel
% -> AWGN
% -> CP Removal
% -> FFT
% -> LS Channel Estimation
% -> ZF / MMSE Equalization
% -> QAM Demodulation
% -> BER
SNR=0:2:30;
BER_ZF=zeros(length(SNR),1);
BER_MMSE=zeros(length(SNR),1);
MSE_H = zeros(length(SNR),1);
numBits=2^20;
srcBits=randi([0,1],numBits,1);
modOrder=16;
numCarr=64;
cycPrefLen=5;
qamModOut=qammod(srcBits,modOrder,InputType="bit");
qamPower = mean(abs(qamModOut(:)).^2);

pilot=sqrt(10)*ones(numCarr,1);%H^ls=H+​W/Xp,Xp[k]=1功率过小，信噪比过大，导致噪声项过大，信道估计有误差
%sqrt(10)使|Xk|^2=10,为16QAM符号平均功率

qamModOut = reshape(qamModOut,numCarr,[]);
qamModOut=[pilot qamModOut];
ifftModOut=ifft(qamModOut,numCarr);
cp = ifftModOut(end-cycPrefLen+1:end,:);
ofdmModOut = [cp;ifftModOut];%注意
ofdmModOut = ofdmModOut(:);%实际 OFDM 信号应该是串行的
mpChan=[0,0.1,0.6,0.5,0.8]';
mpChanOut=filter(mpChan,1,ofdmModOut);
Hhat=zeros(numCarr,length(SNR));
for i=1:length(SNR)
snr=SNR(i);
% Calculate noise power
snrLinear = 10^(snr/10);
noiseVarFreq = qamPower / snrLinear;
noiseVarTime = noiseVarFreq / numCarr;
%Generate complex AWGN
noise = sqrt(noiseVarTime/2) * (randn(size(mpChanOut)) + 1j*randn(size(mpChanOut)));
chanOut=mpChanOut+noise;

rx = reshape(chanOut,numCarr+cycPrefLen,[]);
rxNoCP = rx(cycPrefLen+1:end,:);
fftModOut=fft(rxNoCP,numCarr,1);
pilotRx =fftModOut(:,1);
fftModOut=fftModOut(:,2:end);
Hhat(:,i) = pilotRx ./ pilot;%存在信道估计误差
MSE_H(i) = mean(abs(Hhat(:,i)-Htrue).^2);

eqOut_ZF=fftModOut./Hhat(:,i);
qamdeModOut_ZF=qamdemod(eqOut_ZF,modOrder,OutputType="bit");
BER_ZF(i)=nnz(qamdeModOut_ZF(:)~=srcBits)/numBits;

noiseVarRatio=1/snrLinear ;
mmseWeight = conj(Hhat(:,i)) ./ (abs(Hhat(:,i)).^2 + noiseVarRatio);
eqOut_MMSE=fftModOut.*mmseWeight;
qamdeModOut_MMSE=qamdemod(eqOut_MMSE,modOrder,OutputType="bit");
BER_MMSE(i)=nnz(qamdeModOut_MMSE(:)~=srcBits)/numBits;
end
hold on
semilogy(SNR,BER_ZF,"o-r")%在 x 轴上使用线性刻度、在 y 轴上使用以 10 为底的对数刻度来绘制 x 和 y 坐标
semilogy(SNR,BER_MMSE,"*--b");
legend("ZF","MMSE");
xlabel("SNR(dB)");
ylabel("BER");
title("OFDM BER Performance: ZF vs MMSE");
hold off

Htrue=fft(mpChan,numCarr);
figure;
plot(abs(Htrue),"o-r");
hold on;
plot(abs(Hhat(:,10/2+1)),"x-b");
plot(abs(Hhat(:,30/2+1)),"+-g");
grid on;
xlabel("Subcarrier");
ylabel("|H[k]|");
legend("True Channel","LS Estimate 10dB","LS Estimate 30dB");
title("True Channel vs LS Channel Estimate");
figure;
semilogy(SNR,MSE_H,"o-");
grid on;
xlabel("SNR (dB)");
ylabel("Channel Estimation MSE");
title("LS Channel Estimation MSE");