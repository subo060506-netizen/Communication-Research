% LMMSE 信道估计的意义，就是：已知 H^LS很 noisy，能不能利用信道的统计特性，对这个 H^LS再进行一次最优线性修正？
%LMMSE 的全称：Linear Minimum Mean Square Error即：线性最小均方误差估计。最小化：E[||H-H^||^2] 
%|| ||叫做范数（norm），可以先简单理解成：用来表示一个向量“有多大”。
% e=[e1 e2 e3]它的二范数（Euclidean norm）是：||e||=sqrt(|e1|^2+|e2|^2+|e3|^2)
%既然我们要求“Linear”，那么可以假设：H^=AH^LS,线性是相对H^LS这个观测值
% 估计误差：e=H−AH^LS 均方误差：J(A)=E[||e||^2]=E[(H−A*H^LS​)^H*(H−A*H^LS​)] H共轭转置：
%推导最优 A:z=H^LS​,H^=Az,e=H−Az
%根据正交性原理，最优线性 MMSE 估计满足：E[(H−Az)z^H]=0
%如果误差还和观测量 z 有相关性，就说明：z 中还有一些可以被利用的信息。那就还可以继续调整 A，使误差进一步降低。
%因此最优情况下必须满足：E[Hz^H]−AE[zz^H]=0 
% ​A=E[Hz^H](E[zz^H])^(−1)  H^LMMSE​=E[H*H^LS^H]*E[H^LS​*H^LS^H]^(−1)*H^LS
% H^LS​=H+N假设：E[HN^H]=0也就是：信道和噪声互不相关。E[WW^H]=σw^2​*I   X^(−1)*X^(−1)H=σx^2​*I I为单位矩阵
% 最终得到 LMMSE H^LMMSE​=R^HH​*(R^HH​+​σw^2/σx^2​​*I)^(−1)*H^LS​ R^HH​=E[H*H^H]
% R^HH信道自相关矩阵,统计先验信息​,R^HH​ 描述了不同子载波信道之间的统计相关性,对角线→各子载波自身功率,非对角线→不同子载波之间的相关性 
%注意化简得R^HH是真实信道自相关矩阵 这个矩阵实际仿真里我们可以：方法 1：理论计算 方法 2：蒙特卡洛估计
% ​噪声趋近于 0​,H^LMMSE​→H^LS 噪声非常大 H^LMMSE​→0 SNR 高 → 相信实际观测；SNR 低 → 相信统计先验。​
SNR=0:2:30;
BER_LS=zeros(length(SNR),1);
BER_LMMSE=zeros(length(SNR),1);
numBits=2^20;
srcBits=randi([0,1],numBits,1);
modOrder=16;
numCarr=64;
cycPrefLen=5;
qamModOut=qammod(srcBits,modOrder,InputType="bit");
pilot=sqrt(10)*ones(numCarr,1);
qamModOut=reshape(qamModOut,numCarr,[]);
qamModOut=[pilot qamModOut];
ifftmodOut=ifft(qamModOut,numCarr);
cp= ifftmodOut(numCarr-cycPrefLen+1:end,:);
ofdmmodOut=[cp;ifftmodOut];
ofdmmodOut = ofdmmodOut(:);%实际 OFDM 信号应该是串行的
mpChan=[0,0.1,0.6,0.5,0.8]';
mpChanOut=filter(mpChan,1,ofdmmodOut);
qamPower = mean(abs(qamModOut(:)).^2);
% Estimate R_HH using Monte Carlo
numRealizations = 10000;
R_HH = zeros(numCarr);
for m = 1:numRealizations

    % 生成一次随机信道冲激响应
    h = zeros(numCarr,1);
    %假设E[|h(k)|^2]=[0 0.01 0.36 0.25 0.64 ]
    h(2) = sqrt(0.01/2) * (randn + 1j*randn);
    h(3) = sqrt(0.36/2) * (randn + 1j*randn);
    h(4) = sqrt(0.25/2) * (randn + 1j*randn);
    h(5) = sqrt(0.64/2) * (randn + 1j*randn);
    % 蒙特卡洛方法本质只是：按照一个已经确定的概率模型，其包含信道特征，反复随机采样（即随机生成h）然后统计平均，以此作为优化H_LS依据
    % Monte Carlo 本身并不产生信道模型。正确逻辑应该是：物理环境→信道模型→p(H)→Monte Carlo采样→RHH
	% 0.8/0.4/0.2/0.1 的代码只是教学用的统计信道模型，不应该直接把它当成真实工程信道。一般使用经典模型/标准模型/实测数据模型​​
    % 转换到频域
    H = fft(h,numCarr);

    % 累加 H*H^H
    R_HH = R_HH + H*H';
end
R_HH = R_HH / numRealizations;

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
    H_LMMSE = R_HH * ((R_HH + noiseVarRatio*eye(numCarr)) \ Hhat);

    mmseWeight_LS=conj(Hhat)./(abs(Hhat).^2+noiseVarRatio);
    eqOut_LS=fftModOut.*mmseWeight_LS;
    qamdeModOut_LS=qamdemod(eqOut_LS,modOrder,"OutputType","bit");
    BER_LS(i)=nnz(qamdeModOut_LS(:)~=srcBits)/numBits;
    mmseWeight_LMMSE=conj(H_LMMSE)./(abs(H_LMMSE).^2+noiseVarRatio);
    eqOut_LMMSE=fftModOut.*mmseWeight_LMMSE;
    qamdeModOut_LMMSE=qamdemod(eqOut_LMMSE,modOrder,"OutputType","bit");
    BER_LMMSE(i)=nnz(qamdeModOut_LMMSE(:)~=srcBits)/numBits;
end
hold on
semilogy(SNR,BER_LMMSE,"ro-");
semilogy(SNR,BER_LS,"b*-");
legend("LMMSE","LS");
xlabel("SNR(dB)");
ylabel("BER");
title("OFDM BER Performance: LMMSE vs LS");
hold off

%% Visualize R_HH

figure;

imagesc(real(R_HH));
colorbar;
xlabel('Subcarrier index');
ylabel('Subcarrier index');
title('Real(R_{HH})');


figure;

imagesc(imag(R_HH));
colorbar;
xlabel('Subcarrier index');
ylabel('Subcarrier index');
title('Imag(R_{HH})');


figure;

imagesc(abs(R_HH));
colorbar;
xlabel('Subcarrier index');
ylabel('Subcarrier index');
title("|R_{HH}|");
%% Check Hermitian property检查R^HH是否为自共轭矩阵

hermitianError = norm(R_HH - R_HH', 'fro');

fprintf('Hermitian error = %.6e\n', hermitianError);
%Hermitian error足够小时，认为满足自共轭，矩阵第i行第j列的元素都与第j行第i列的元素的共轭相等
%% Check positive semidefinite property检查 R^HH是否半正定

eig_RHH = eig(R_HH);

fprintf('Minimum eigenvalue = %.6e\n', min(real(eig_RHH)));
figure;

plot(real(eig_RHH), 'o-');
xlabel('Eigenvalue index');
ylabel('Eigenvalue');
title('Eigenvalues of R_{HH}');
grid on;
%% Diagonal of R_HH
%每个 OFDM 子载波平均经历多大的信道功率。
figure;

plot(0:numCarr-1, real(diag(R_HH)), 'o-');
xlabel('Subcarrier index');
ylabel('E[|H[k]|^2]');
title('Diagonal of R_{HH}');
grid on;
%% Theoretical R_HH

pathDelay = [1 2 3 4];
pathPower = [0.01 0.36 0.25 0.64];

R_HH_theory = zeros(numCarr,numCarr);

for k = 0:numCarr-1

    for m = 0:numCarr-1

        value = 0;

        for l = 1:length(pathDelay)

            value = value + ...
                pathPower(l) * ...
                exp(-1j*2*pi*(k-m)*pathDelay(l)/numCarr);

        end

        R_HH_theory(k+1,m+1) = value;

    end

end
error_RHH = norm(R_HH - R_HH_theory, 'fro') ...
            / norm(R_HH_theory, 'fro');
% 并不是“理论模型误差”，而是Monte Carlo 有限采样误差​
fprintf('Relative error = %.6e\n', error_RHH);
figure;

imagesc(abs(R_HH - R_HH_theory));
colorbar;

xlabel('Subcarrier index');
ylabel('Subcarrier index');

title('|R_{HH}^{MC}-R_{HH}^{Theory}|');
% LS：H^LS->min|Y−X*H^LS|^2它关心:“我怎样最贴合这一次收到的 Y？”
%LMMSE：H^LMMSE->minE[||H-H^||^2]​它关心：“我怎样利用统计规律，使长期平均的真实信道估计误差最小？”
% LMMSE 利用 E[ ] 把问题从“一次观测的最优拟合”提升到了“随机信道和噪声条件下的平均最优估计”​​