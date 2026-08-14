%Y[k]=H[k]X[k]+W[k]
%X^ZF​[k]=H[k]Y[k]​​
%X^ZF​[k]=X[k]+H[k]W[k]
% 当∣H[k]∣→0，噪声就会被巨大放大（低信噪比影响更明显）
% MMSE在“消除信道影响”和“抑制噪声增强”之间做一个折中。
% MMSE 均衡器：G^MMSE​[k]=H*[k] / (|H[k]|^2+σw^2/σx^2)​
% X^[k]=GMMSE​[k]Y[k]
% H[k]：信道频率响应
%H*[k]：信道频率响应的共轭
%σw^2：噪声功率
%σx^2：发送信号功率​
% ZF：G^ZF[K]=1/H[k]    MMSE：G^MMSE[K]
%信噪比高时H*[k] / |H[k]|^2=1/H[k],MMSE->LS,噪声影响忽略，只消除信道影响
%H[k]->0,G^MMSE=σx^2/σw^2=SNR,信噪比低时，抑制噪声增强，但信道影响无法消除
SNR=0:2:30;
BER_ZF=zeros(length(SNR),1);
BER_MMSE=zeros(length(SNR),1);
numBits=2^20;
srcBits=randi([0,1],numBits,1);
modOrder=16;
numCarr=64;
cycPrefLen=5;
qamModOut=qammod(srcBits,modOrder,InputType="bit");
qamModOut = reshape(qamModOut,numCarr,[]);
ifftModOut=ifft(qamModOut,numCarr);
cp = ifftModOut(end-cycPrefLen+1:end,:);
ofdmModOut = [cp;ifftModOut];%注意
ofdmModOut = ofdmModOut(:);%实际 OFDM 信号应该是串行的
mpChan=[0,0.1,0.6,0.5,0.8]';
mpChanOut=filter(mpChan,1,ofdmModOut);

qamPower = mean(abs(qamModOut(:)).^2);
MSE_ZF = zeros(length(SNR),1);
MSE_MMSE = zeros(length(SNR),1);
for i=1:length(SNR)
snr=SNR(i);
%chanOut=awgn(mpChanOut,snr,"measured");
%awgn中SNR信号功率是mpChanOut信号，与MMSE中不同
% Calculate noise power
snrLinear = 10^(snr/10);
noiseVarFreq = qamPower / snrLinear;
noiseVarTime = noiseVarFreq / numCarr;
%Generate complex AWGN
noise = sqrt(noiseVarTime/2) * (randn(size(mpChanOut)) + 1j*randn(size(mpChanOut)));%randn(size(mpChanOut))：生成与mpChanOut维度相同的实高斯噪声（均值=0，方差=1）
%σ² = noiseVarTime/2 才能使总噪声功率等于 noiseVarTime
chanOut = mpChanOut + noise;
rx = reshape(chanOut,numCarr+cycPrefLen,[]);
rxNoCP = rx(cycPrefLen+1:end,:);
fftModOut=fft(rxNoCP,numCarr,1);
mpChanFreq=fft(mpChan,numCarr);
noiseVarRatio = 1/snrLinear;
mmseWeight = conj(mpChanFreq) ./ (abs(mpChanFreq).^2 + noiseVarRatio);
eqOut_ZF=fftModOut ./mpChanFreq;
eqOut_MMSE= fftModOut .* mmseWeight;
%MMSE估计的幅度收缩α[k]=∣H[k]∣^2/(|H[k]|^2​+λ)→QAM判决偏差
%理解：MMSE允许一定的信号衰减,换取减少噪声增强
qamdeModOut_ZF=qamdemod(eqOut_ZF,modOrder,OutputType="bit");
qamdeModOut_MMSE=qamdemod(eqOut_MMSE,modOrder,OutputType="bit");
BER_ZF(i)=nnz(qamdeModOut_ZF(:)~=srcBits)/numBits;
BER_MMSE(i)=nnz(qamdeModOut_MMSE(:)~=srcBits)/numBits;
MSE_ZF(i) = mean(abs(eqOut_ZF(:) - qamModOut(:)).^2);
MSE_MMSE(i) = mean(abs(eqOut_MMSE(:) - qamModOut(:)).^2);
end
hold on
semilogy(SNR,BER_ZF,"o-r")%在 x 轴上使用线性刻度、在 y 轴上使用以 10 为底的对数刻度来绘制 x 和 y 坐标
semilogy(SNR,BER_MMSE,"*--b");
legend("ZF","MMSE");
xlabel("SNR(dB)");
ylabel("BER");
title("OFDM BER Performance: ZF vs MMSE");
hold off
figure;
semilogy(SNR,MSE_ZF,"o-r");
hold on;
semilogy(SNR,MSE_MMSE,"*--b");
grid on;
legend("ZF","MMSE");
xlabel("SNR (dB)");
ylabel("MSE");
title("ZF vs MMSE: MSE Performance");
%MMSE最小化的是均方误差 MSE，而不是直接最小化BER,MSE更小不代表BER一定更低​
%MSE问的是：“你估计得有多接近？”BER问的是：“你最终判错了多少个比特？”