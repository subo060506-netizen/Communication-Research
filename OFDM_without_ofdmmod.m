% OFDM system with ZF equalization
% This script implements a basic OFDM communication system
% with 16-QAM, multipath channel, AWGN and ZF equalization.
SNR=0:2:50;
BER=zeros(length(SNR),1);
numBits=2^20;
srcBits=randi([0,1],numBits,1);
modOrder=16;
numCarr=64;
cycPrefLen=5;
qamModOut=qammod(srcBits,modOrder,InputType="bit");
qamModOut = reshape(qamModOut,numCarr,[]);%可以指定 [] 的单个维度大小，以便自动计算维度大小，以使 B 中的元素数与 A 中的元素数相匹配。例如，如果 A 是一个 10×10 矩阵，则 reshape(A,2,2,[]) 将 A 的 100 个元素重构为一个 2×2×25 数组。注意：与位置有关
ifftModOut=ifft(qamModOut,numCarr);
cp = ifftModOut(end-cycPrefLen+1:end,:);
ofdmModOut = [cp;ifftModOut];%注意
ofdmModOut = ofdmModOut(:);%实际 OFDM 信号应该是串行的
% ofdmModOut=zeros(floor(length(ifftModOut))/numCarr*cycPrefLen+length(ifftModOut),1);length() 对矩阵不是返回元素总数，而是返回最大维度。
% for idx=1:length(ifftModOut)
%     n=floor((idx-1)/numCarr)+1;
%     y=mod(idx-1,numCarr);
%     b=numCarr-cycPrefLen;
%     if y>b
%        ofdmModOut(n*numCarr+n*cycPrefLen+y-b)=ifftModOut(idx);
%        ofdmModOut(n*numCarr+n*cycPrefLen+cycPrefLen+y)=ifftModOut(idx);
%     else
%        ofdmModOut(n*numCarr+n*cycPrefLen+cycPrefLen+y)=ifftModOut(idx); 
%     end
% end
mpChan=[0,0.1,0.6,0.5,0.8]';
mpChanOut=filter(mpChan,1,ofdmModOut);
for i=1:length(SNR)
snr=SNR(i);
chanOut=awgn(mpChanOut,snr,"measured");
rx = reshape(chanOut,numCarr+cycPrefLen,[]);
rxNoCP = rx(cycPrefLen+1:end,:);
fftModOut=fft(rxNoCP,numCarr,1);
mpChanFreq=fft(mpChan,numCarr);
eqOut=fftModOut./mpChanFreq;
qamdeModOut=qamdemod(eqOut,modOrder,OutputType="bit");
BER(i)=nnz(qamdeModOut(:)~=srcBits)/numBits;
end
semilogy(SNR,BER,"o-")%在 x 轴上使用线性刻度、在 y 轴上使用以 10 为底的对数刻度来绘制 x 和 y 坐标
xlabel("SNR(dB)")
ylabel("BER")
title("OFDM BER Performance with ZF Equalization");
figure;
plot(0:numCarr-1,abs(mpChanFreq),'o-');
grid on;
xlabel("Subcarrier");
ylabel("|H[k]|");
title("Channel Frequency Response");