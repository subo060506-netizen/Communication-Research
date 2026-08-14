%最简单的方法是:
%发送一个接收机事先知道的 OFDM 符号
%这个特殊符号叫：Pilot（导频）
%假设某个子载波上发送的导频是：Xp=[k]
%Yp​[k]=H[k]Xp​[k]+W[k]
%Y[k]≈X[k]H^[k]
%误差J=∣Y[k]−Xp[k]H^[k]∣^2
%求极值，H^[k]=Yp[k]/Xp[k]
%取Xp(k)=1,H^(k)=Yp(k)
SNR=0:2:50;
BER=zeros(length(SNR),1);
numBits=2^20;
srcBits=randi([0,1],numBits,1);
modOrder=16;
numCarr=64;
cycPrefLen=5;
qamModOut=qammod(srcBits,modOrder,InputType="bit");
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
for i=1:length(SNR)
snr=SNR(i);
chanOut=awgn(mpChanOut,snr,"measured");
rx = reshape(chanOut,numCarr+cycPrefLen,[]);
rxNoCP = rx(cycPrefLen+1:end,:);
fftModOut=fft(rxNoCP,numCarr,1);
pilotRx =fftModOut(:,1);
Hhat = pilotRx ./ pilot;%存在信道估计误差
fftModOut=fftModOut(:,2:end);
eqOut=fftModOut./Hhat;
qamdeModOut=qamdemod(eqOut,modOrder,OutputType="bit");
BER(i)=nnz(qamdeModOut(:)~=srcBits)/numBits;
end
semilogy(SNR,BER,"o-")%在 x 轴上使用线性刻度、在 y 轴上使用以 10 为底的对数刻度来绘制 x 和 y 坐标
xlabel("SNR(dB)")
ylabel("BER")
title("OFDM BER Performance with LS Channel Estimation and ZF Equalization");
scatterplot(eqOut(:));
%LS 信道估计的误差与 Pilot 的信噪比直接相关。LS 在低 SNR 下不够好
%对比曲线：∣H[k]∣和|H^LS[k]|
% Htrue=fft(mpChan,numCarr);
% figure;
% plot(abs(Htrue),"o-");
% hold on;
% plot(abs(Hhat),"x-");
% grid on;
% xlabel("Subcarrier");
% ylabel("|H[k]|");
% legend("True Channel","LS Estimate");
% title("True Channel vs LS Channel Estimate");