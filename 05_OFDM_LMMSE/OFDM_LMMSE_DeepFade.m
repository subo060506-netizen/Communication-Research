SNR=0:2:30;
BER_Htrue=zeros(length(SNR),1);
BER_LS=zeros(length(SNR),1);
BER_LMMSE=zeros(length(SNR),1);
numBits=2^20;
srcBits=randi([0,1],numBits,1);
modOrder=16;
numCarr=64;
cycPrefLen=5;

%% 信道构建
alphaList=[0.50 0.80 0.95 0.99];
mpChan = zeros(2,length(alphaList));
%功率归一
for i = 1:length(alphaList)

    alpha = alphaList(i);

    mpChan(:,i) = ...
        [1; -alpha] / sqrt(1 + alpha^2);

end
Htrue=fft(mpChan,numCarr,1);
depth=min(abs(Htrue),[],1);
hold on
for i = 1:length(alphaList)  
    
    plot(abs(Htrue(:,i)),"o-");
end
hold off
xlabel("Subcarrier");
ylabel("|H[K]|");
legend( ...
    "\alpha = 0.50", ...
    "\alpha = 0.80", ...
    "\alpha = 0.95", ...
    "\alpha = 0.99");
title("Frequency Response Under Different Deep-Fading Levels");

%% 解释深衰落影响
%深衰落主要是在后面的均衡阶段把 LS 的误差放大

%% 求 R_HH
% Estimate R_HH using Monte Carlo
numRealizations = 10000;
R_HH= zeros(numCarr,numCarr,length(alphaList));
for i=1:length(alphaList)
    for m = 1:numRealizations

        % 生成一次随机信道冲激响应
        alpha = alphaList(i);
        g = sqrt(1/2) * (randn + 1j*randn);
        h = zeros(numCarr,1);
        h(1) = g / sqrt(1 + alpha^2);
        h(2) = -alpha * g / sqrt(1 + alpha^2);
        %第二径相对于第一径的幅度比例，并具有固定的 π 相对相位。
        %构造的是一个固定深衰落结构 + 随机整体复幅度的统计模型
        H = fft(h,numCarr);

        % 累加 H*H^H
        R_HH(:,:,i)=  R_HH(:,:,i)+ H*H';%H'对于复数向量表示共轭转置。
    end
    R_HH(:,:,i)= R_HH(:,:,i)/ numRealizations;
end
%%
qamModOut=qammod(srcBits,modOrder,InputType="bit");
pilot=sqrt(10)*ones(numCarr,1);
qamModOut=reshape(qamModOut,numCarr,[]);
qamModOut=[pilot qamModOut];
ifftmodOut=ifft(qamModOut,numCarr);
cp= ifftmodOut(numCarr-cycPrefLen+1:end,:);
ofdmmodOut=[cp;ifftmodOut];
ofdmmodOut = ofdmmodOut(:);%实际 OFDM 信号应该是串行的
qamPower = mean(abs(qamModOut(:)).^2);
pilotPower = mean(abs(pilot).^2);
for a=1:length(alphaList)
    mpChanOut=filter(mpChan(:,a),1,ofdmmodOut);
    for i=1:length(SNR)
        snr=SNR(i);
        snrLinear = 10^(snr/10);
        noiseVarFreq = qamPower / snrLinear;
        noiseVarTime = noiseVarFreq / numCarr;
        noiseVarRatio=1/snrLinear ;
        %Generate complex AWGN
        noise = sqrt(noiseVarTime/2) * (randn(size(mpChanOut)) + 1j*randn(size(mpChanOut)));
        %复高斯噪声由实部和虚部两个独立的高斯分量组成
        chanOut = mpChanOut + noise;
        rx = reshape(chanOut,numCarr+cycPrefLen,[]);
        rxNoCP = rx(cycPrefLen+1:end,:);
        fftModOut=fft(rxNoCP,numCarr,1);
        pilotRx =fftModOut(:,1);
        fftModOut=fftModOut(:,2:end);
        Hhat = pilotRx ./ pilot;%存在信道估计误差
        H_LMMSE = R_HH(:,:,a) * ((R_HH (:,:,a)+noiseVarFreq / pilotPower *eye(numCarr)) \ Hhat);
        %LMMSE中这里的 σx2应该对应导频功率。
        %MMSE中信噪比信号代表的是QAM
        mmseWeight_LS=conj(Hhat)./(abs(Hhat).^2+noiseVarRatio);
        eqOut_LS=fftModOut.*mmseWeight_LS;
        qamdeModOut_LS=qamdemod(eqOut_LS,modOrder,"OutputType","bit");
        BER_LS(i)=nnz(qamdeModOut_LS(:)~=srcBits)/numBits;

        mmseWeight_LMMSE=conj(H_LMMSE)./(abs(H_LMMSE).^2+noiseVarRatio);
        eqOut_LMMSE=fftModOut.*mmseWeight_LMMSE;
        qamdeModOut_LMMSE=qamdemod(eqOut_LMMSE,modOrder,"OutputType","bit");
        BER_LMMSE(i)=nnz(qamdeModOut_LMMSE(:)~=srcBits)/numBits;

        mmseWeight_Htrue=conj(Htrue(:,a))./(abs(Htrue(:,a)).^2+noiseVarRatio);
        eqOut_Htrue=fftModOut.*mmseWeight_Htrue;
        qamdeModOut_Htrue=qamdemod(eqOut_Htrue,modOrder,"OutputType","bit");
        BER_Htrue(i)=nnz(qamdeModOut_Htrue(:)~=srcBits)/numBits;
    end
    figure
    hold on
    semilogy(SNR,BER_LMMSE,"ro-");
    semilogy(SNR,BER_LS,"b*-");
    semilogy(SNR,BER_Htrue,"g+-");
    set(gca,'YScale','log');
    legend("LMMSE","LS","Perfect_CSI");
    xlabel("SNR(dB)");
    ylabel("BER");
    %title(['OFDM BER Performance(alpha=',alphaList(a),'): LMMSE vs LS vs Perfect_CSI']);
    %方法一（推荐，使用 + 拼接）：
    title('OFDM BER Performance (alpha = ' + string(alphaList(a)) + '): LMMSE vs LS vs Perfect CSI');
    %方法二（使用 sprintf 格式化）：title(sprintf('OFDM BER Performance (alpha = %.2f): LMMSE vs LS', alphaList(a)));
    hold off
    
end
%在所构建的两径深衰落信道中，随着第二径幅度比例 α 增大，频域深衰落逐渐加剧，系统 BER 性能明显下降。
% LMMSE 信道估计在不同深衰落条件下均表现出优于 LS 的性能，并且其 BER 曲线与 Perfect CSI 条件下的结果高度接近，表明 LMMSE 已能够有效降低信道估计误差。
% 然而，在极端深衰落条件下，即使采用 Perfect CSI，系统仍存在较高 BER，
% 说明此时系统性能损失主要来源于信道本身的深衰落，而非信道估计误差。