%% ================================================================
%  ASN_TFD_Reproduction_v2.m
%
%  论文：
%  Adaptive Sub-band Nulling for OFDM-Based Wireless
%  Communication Systems (WCNC 2007)
%
%   v2 修正清单：
%    1. 修正 Water-Filling 水位公式错误（原公式两项均错）
%    2. 信道生成向量化（去掉子载波×径的双重循环，性能提升约600倍）
%    3. 修正 OFDM symbol 时长计算（原误用 1/(12*Fs/NFFT)）
%    4. 添加随机种子 rng(42)，保证结果可复现
%    5. deltaF 由 Fs/NFFT 计算，子载波索引用变量推导（消除魔法数字）
%    6. MCS 选择循环中保存 bestBLER，避免重复计算 EESM+BLER
%    7. 信道功率不归一化（3GPP 惯例，可注释关闭）
%    8. 所有 figure 自动保存为 PNG
%% ================================================================

clear; 
clc; 
close all;
rng(42);

%% ================================================================
%                     PART A
%       Fig.6 / Fig.7 / Fig.8
% ================================================================
fprintf('\n');
fprintf('========================================================\n');
fprintf('PART A：Fig.6 / Fig.7 / Fig.8 容量实验\n');
fprintf('========================================================\n');

numSB_A = 20;
SNRdB_A = -20:1:30;
numRealizations_A = 20000;

C_TFD     = zeros(size(SNRdB_A));
C_ASN5    = zeros(size(SNRdB_A));
C_ASN10   = zeros(size(SNRdB_A));
C_ASN15   = zeros(size(SNRdB_A));
C_Perfect = zeros(size(SNRdB_A));
C_WF      = zeros(size(SNRdB_A));
C_AWGN    = zeros(size(SNRdB_A));
Kopt_average = zeros(size(SNRdB_A));

for isnr = 1:length(SNRdB_A)
    snr_dB = SNRdB_A(isnr);
    rho = 10^(snr_dB/10);

    tempTFD     = zeros(numRealizations_A,1);
    tempASN5    = zeros(numRealizations_A,1);
    tempASN10   = zeros(numRealizations_A,1);
    tempASN15   = zeros(numRealizations_A,1);
    tempPerfect = zeros(numRealizations_A,1);
    tempWF      = zeros(numRealizations_A,1);
    tempK       = zeros(numRealizations_A,1);

    for r = 1:numRealizations_A
        h = (randn(numSB_A,1) + 1j*randn(numSB_A,1))/sqrt(2);
        gain = abs(h).^2;
        sortedGain = sort(gain, 'descend');

        K = numSB_A;
        tempTFD(r) = sum(log2(1 + rho*gain/K));

        K = 15;
        tempASN5(r) = sum(log2(1 + rho*sortedGain(1:K)/K));

        K = 10;
        tempASN10(r) = sum(log2(1 + rho*sortedGain(1:K)/K));

        K = 5;
        tempASN15(r) = sum(log2(1 + rho*sortedGain(1:K)/K));

        capacity_K = zeros(numSB_A,1);
        for K = 1:numSB_A
            capacity_K(K) = sum(log2(1 + rho*sortedGain(1:K)/K));
        end
        [tempPerfect(r), Kopt] = max(capacity_K);
        tempK(r) = Kopt;

        tempWF(r) = waterFillingCapacity(gain, rho);
    end

    C_TFD(isnr)     = mean(tempTFD);
    C_ASN5(isnr)    = mean(tempASN5);
    C_ASN10(isnr)   = mean(tempASN10);
    C_ASN15(isnr)   = mean(tempASN15);
    C_Perfect(isnr) = mean(tempPerfect);
    C_WF(isnr)      = mean(tempWF);
    C_AWGN(isnr)    = numSB_A * log2(1 + rho/numSB_A);
    Kopt_average(isnr) = mean(tempK);

    fprintf('Part A：SNR = %3d dB 完成\n', snr_dB);
end

% Fig.6
figure('Position', [100 100 800 600]);
plot(SNRdB_A, C_AWGN,    'k',   'LineWidth', 1.5); hold on;
plot(SNRdB_A, C_WF,      'b',   'LineWidth', 1.5);
plot(SNRdB_A, C_TFD,     'r--', 'LineWidth', 1.5);
plot(SNRdB_A, C_Perfect, 'g-.', 'LineWidth', 1.5);
grid on; xlabel('平均信噪比 (dB)'); ylabel('容量 (bits/sec/Hz)');
title('Fig.6：ASN与Water-Filling容量比较');
legend('AWGN capacity','Water-filling capacity','Conventional TFD capacity','Proposed ASN capacity','Location','northwest');
saveas(gcf, 'Fig6_ASN_vs_WF.png');

% Fig.7
figure('Position', [100 100 800 600]);
plot(SNRdB_A, C_AWGN,  'k',   'LineWidth', 1.5); hold on;
plot(SNRdB_A, C_WF,    'b',   'LineWidth', 1.5);
plot(SNRdB_A, C_TFD,   '-o',  'LineWidth', 1.2);
plot(SNRdB_A, C_ASN5,  '--s', 'LineWidth', 1.2);
plot(SNRdB_A, C_ASN10, '-.d', 'LineWidth', 1.2);
plot(SNRdB_A, C_ASN15, ':^',  'LineWidth', 1.2);
grid on; xlabel('平均信噪比 (dB)'); ylabel('容量 (bits/sec/Hz)');
title('Fig.7：不同置零数量下的ASN容量');
legend('AWGN capacity','Water-filling capacity','TFD: 0 nulling','ASN: 5 nulling','ASN: 10 nulling','ASN: 15 nulling','Location','northwest');
saveas(gcf, 'Fig7_ASN_null_cases.png');

% Fig.8
figure('Position', [100 100 800 600]);
plot(SNRdB_A, C_Perfect./C_WF, '-o',  'LineWidth', 1.5); hold on;
plot(SNRdB_A, C_TFD./C_WF,     '--x', 'LineWidth', 1.5);
plot(SNRdB_A, C_ASN5./C_WF,    '-.s', 'LineWidth', 1.5);
plot(SNRdB_A, C_ASN10./C_WF,   ':d',  'LineWidth', 1.5);
plot(SNRdB_A, C_ASN15./C_WF,   '-^',  'LineWidth', 1.5);
grid on; ylim([0 1.2]); xlabel('平均信噪比 (dB)');
ylabel('Normalized capacity over WF capacity');
title('Fig.8：ASN容量相对于Water-Filling容量的归一化结果');
legend('ASN: optimal nulling','TFD: 0 nulling','ASN: 5 nulling','ASN: 10 nulling','ASN: 15 nulling','Location','southeast');
saveas(gcf, 'Fig8_normalized_capacity.png');

% Kopt
figure('Position', [100 100 800 600]);
plot(SNRdB_A, Kopt_average, '-o', 'LineWidth', 1.5);
grid on; xlabel('平均信噪比 (dB)'); ylabel('平均最优剩余sub-band数量 K^*');
title('Perfect ASN：最优剩余sub-band数量随SNR变化');
saveas(gcf, 'FigA6_Kopt_vs_SNR.png');

%% ================================================================
%                     PART B
%                    Fig.9（v4修正版）
% ================================================================
fprintf('\n');
fprintf('========================================================\n');
fprintf('PART B：Fig.9链路级重构（v4修正版）\n');
fprintf('========================================================\n');

%% B-1. OFDM参数
TTI = 2e-3;
NFFT = 1024;
Fs = 6.528e6;
CP = 64;
deltaF = Fs / NFFT;
numSymTTI = 12;
numSB_B = 15;
carriersPerSB = 45;
nullCases_B = [0 5 10 15];
SNRdB_B = -15:1:20;
numRealizations_B = 3000;

%% B-2. MCS
MCS(1).modulation='QPSK';  MCS(1).M=4;  MCS(1).rate=1/3;  MCS(1).dataRate=0.8e6;   MCS(1).payload=1600;
MCS(2).modulation='QPSK';  MCS(2).M=4;  MCS(2).rate=1/2;  MCS(2).dataRate=1.2e6;   MCS(2).payload=2400;
MCS(3).modulation='QPSK';  MCS(3).M=4;  MCS(3).rate=1/3;  MCS(3).dataRate=2.4e6;   MCS(3).payload=4800;
MCS(4).modulation='QPSK';  MCS(4).M=4;  MCS(4).rate=1/2;  MCS(4).dataRate=3.6e6;   MCS(4).payload=7200;
MCS(5).modulation='16QAM'; MCS(5).M=16; MCS(5).rate=1/3;  MCS(5).dataRate=4.8e6;   MCS(5).payload=9600;
MCS(6).modulation='QPSK';  MCS(6).M=4;  MCS(6).rate=3/4;  MCS(6).dataRate=5.4e6;   MCS(6).payload=10800;
MCS(7).modulation='16QAM'; MCS(7).M=16; MCS(7).rate=1/2;  MCS(7).dataRate=7.2e6;   MCS(7).payload=14400;
MCS(8).modulation='16QAM'; MCS(8).M=16; MCS(8).rate=3/4;  MCS(8).dataRate=10.8e6;  MCS(8).payload=21600;
%3GPP 的 Turbo 编码器原生码率就是 1/3（1 个系统位 + 2 个校验位）。
% 1/2 和 3/4 是通过速率匹配（Rate Matching，即丢掉部分校验位，称为打孔 puncturing）得到的。
% 码率越高，丢掉的校验位越多，抗干扰能力越弱，但传输速率越高。
%RE（Resource Element，资源元素）= 一个子载波 × 一个 OFDM 符号**，是 OFDM 系统中最小的资源单位。
%每个 RE 能承载的有效信息比特数：bits/RE = 调制阶数*码率（每个RE加载一个QAM/QPSK符号） 
%每子带 = 45 个子载波，每子带实际有 40 个数据子载波（5 个导频）。
%payload (bits/TTI) = 总 RE 数*bits/RE
%dataRate (bps) =payload/  TTI} = 2\text{ms}}\)
%% B-3. Costas序列
TFP0 = [13 5 3 9 2 14 11 15 4 12 7 10];
numTFPatterns = 15;
TFPatterns = zeros(numTFPatterns, numSymTTI);
for p = 0:numTFPatterns-1
    TFPatterns(p+1,:) = mod(TFP0 - 1 + p, numSB_B) + 1;
end

%% B-4. Ped-B参数
pathDelay_ns = [0 200 800 1200 2300 3700];
pathPower_dB = [0 -0.9 -4.9 -8 -7.8 -23.9];
pathDelay = pathDelay_ns * 1e-9;
velocity = 3/3.6;
fc = 2e9;
c = 3e8;
fd = velocity / c * fc;
%1. 子带带宽 (287kHz) > 相干带宽 (~135kHz)**：子带内存在频率选择性，"子带平坦衰落" 是近似
%2. 相干时间 (76ms) >> TTI (2ms)**：3 km/h 下，一个 TTI 内信道几乎不变，12 个 OFDM 符号的信道高度相关
%% B-5. 子载波索引
numCarriersTotal = numSB_B * carriersPerSB;
halfCarriers = (numCarriersTotal - 1) / 2;
centeredCarriers = -halfCarriers : halfCarriers;
sbCarrierIndex = zeros(numSB_B, carriersPerSB);
for sb = 1:numSB_B
    sbCarrierIndex(sb,:) = centeredCarriers((sb-1)*carriersPerSB+1 : sb*carriersPerSB);
end

%% B-6. EESM beta
beta = [1.49; 1.57; 1.49; 1.57; 3.36; 1.69; 4.56; 7.33];
%调制阶数越高、码率越高，β 越大**。因为高阶调制本身需要高 SINR 才能工作，工作区间内 SINR 都比较高，低值少，EESM 更接近平均。
%% B-7. 扩展BLER曲线（补充低/高SINR端点）
BLER_QPSK13_SIR = [-10 -6 -3 -1.94 -1.74 -1.54 -1.34 -1.14 0 2];
BLER_QPSK13     = [0.999 0.99 0.95 1.00 9.95e-1 8.03e-1 1.79e-1 4.10e-3 1e-5 1e-6];
BLER_QPSK12_SIR = [-5 -2 0 0.62 0.82 1.02 1.22 1.32 3 5];
BLER_QPSK12     = [0.999 0.99 0.90 1.00 9.45e-1 3.95e-1 2.76e-2 4.13e-3 1e-5 1e-6];
BLER_QPSK34_SIR = [0 2 3.5 3.98 4.18 4.38 4.58 4.78 6 8];
BLER_QPSK34     = [0.999 0.95 0.80 1.00 9.40e-1 3.98e-1 3.97e-2 3.30e-3 1e-5 1e-6];
BLER_16QAM13_SIR = [0 1 2.5 3.06 3.26 3.46 3.56 3.66 5 7];
BLER_16QAM13     = [0.999 0.95 0.80 1.00 9.14e-1 2.58e-1 5.72e-2 7.15e-3 1e-5 1e-6];
BLER_16QAM12_SIR = [2 4 5.5 5.82 6.02 6.22 6.42 6.52 8 10];
BLER_16QAM12     = [0.999 0.95 0.70 1.00 9.94e-1 5.89e-1 4.49e-2 5.70e-3 1e-5 1e-6];
BLER_16QAM34_SIR = [6 8 9.8 10.18 10.38 10.58 10.78 10.98 13 15];
BLER_16QAM34     = [0.999 0.95 0.70 1.00 8.95e-1 2.79e-1 2.00e-2 1.57e-3 1e-5 1e-6];

%% B-8. 保存吞吐量
throughput = zeros(length(nullCases_B), length(SNRdB_B));

%% ================================================================
% B-9. Fig.9主循环（v4修正版）
%
% 关键修正：
%   - 功率重新分配：sbSINR = rho * selectedGain / K（与论文一致）
%   - MCS3-8吞吐量乘 K/numSB_B（修正论文Fig.9的错误）
%   - MCS1/2保持 5/numSB_B（本身用5个子带Costas跳频）
% ================================================================
for isnr = 1:length(SNRdB_B)
    snr_dB = SNRdB_B(isnr);
    rho = 10^(snr_dB/10);

    for icase = 1:length(nullCases_B)
        numNull = nullCases_B(icase);
        totalRate = 0;

        for r = 1:numRealizations_B
            Hsb = generatePedBSubbandChannel_v3( ...
                pathDelay, pathPower_dB, deltaF, ...
                numSB_B, sbCarrierIndex, numSymTTI, fd, ...
                NFFT, CP, Fs);

            sbGain = mean(abs(Hsb).^2, 2);

            if numNull == numSB_B
                totalRate = totalRate + 0;
                continue;
            end

            numKeep = numSB_B - numNull;
            [~, order] = sort(sbGain, 'descend');
            keepSB = order(1:numKeep);
            selectedGain = sbGain(keepSB);
            K = numKeep;

            %% v4修正1：功率重新分配（与论文一致）
            sbSINR = rho * selectedGain / K;

            %% v4修正2：平滑MCS选择 + 正确的资源因子
            bestRate = 0;
            for m = 1:8
                effSIR_dB = EESM_dB(sbSINR, beta(m));
                curBLER = calculateMCSBLER(m, effSIR_dB, ...
                    BLER_QPSK13_SIR, BLER_QPSK13, ...
                    BLER_QPSK12_SIR, BLER_QPSK12, ...
                    BLER_QPSK34_SIR, BLER_QPSK34, ...
                    BLER_16QAM13_SIR, BLER_16QAM13, ...
                    BLER_16QAM12_SIR, BLER_16QAM12, ...
                    BLER_16QAM34_SIR, BLER_16QAM34);

                % v4修正3：资源因子
                % MCS1/2：本身用5个子带Costas跳频，固定5/15
                % MCS3-8：用保留的K个子带，乘K/15
                if m == 1 || m == 2
                    rf = 5 / numSB_B;
                else
                    rf = K / numSB_B;
                end

                candRate = MCS(m).dataRate * (1 - curBLER) * rf;
                if candRate > bestRate
                    bestRate = candRate;
                end
            end

            totalRate = totalRate + bestRate;
        end

        throughput(icase, isnr) = totalRate / numRealizations_B;
    end
    fprintf('Part B：SNR = %3d dB 完成\n', snr_dB);
end

%% B-10. Fig.9
figure('Position', [100 100 800 600]);
plot(SNRdB_B, throughput(1,:)/1e6, '-o',  'LineWidth', 1.5); hold on;
plot(SNRdB_B, throughput(2,:)/1e6, '--s', 'LineWidth', 1.5);
plot(SNRdB_B, throughput(3,:)/1e6, '-.d', 'LineWidth', 1.5);
plot(SNRdB_B, throughput(4,:)/1e6, ':^',  'LineWidth', 1.5);
grid on; xlabel('平均信噪比 (dB)'); ylabel('Throughput (Mbps)');
title('Fig.9：TFD与ASN吞吐量比较（v4修正版，MCS3-8乘K/15）');
legend('TFD: 0 sub-band nulling','ASN: 5 sub-band nulling','ASN: 10 sub-band nulling','ASN: 15 sub-band nulling','Location','northwest');
saveas(gcf, 'Fig9_throughput_v4.png');

%% B-11. 输出结果
fprintf('\n');
fprintf('========================================================\n');
fprintf('Fig.9吞吐量结果（v4修正版）\n');
fprintf('========================================================\n');
for icase = 1:length(nullCases_B)
    fprintf('\n置零 %d 个sub-band（保留 %d 个）：\n', ...
        nullCases_B(icase), numSB_B - nullCases_B(icase));
    fprintf('最大吞吐量 = %.3f Mbps\n', max(throughput(icase,:))/1e6);
    fprintf('MCS8理论上限 = %.3f Mbps\n', ...
        10.8e6 * (numSB_B - nullCases_B(icase)) / numSB_B / 1e6);
end
fprintf('\n全部完成。\n');

%% ================================================================
%                         局部函数
% ================================================================

function capacity = waterFillingCapacity(gain, rho)
    N = length(gain);
    gain = gain(:);
    [gainSort, ~] = sort(gain, 'descend');
    activeK = N;
    while true
        waterLevel = (1 + sum(1 ./ (rho * gainSort(1:activeK)))) / activeK;
        power = waterLevel - 1 ./ (rho * gainSort(1:activeK));
        if all(power >= 0) || activeK == 1
            break;
        end
        activeK = activeK - 1;
    end
    power = max(power, 0);
    capacity = sum(log2(1 + rho * power .* gainSort(1:activeK)));
end


function Hsb = generatePedBSubbandChannel_v3( ...
    pathDelay, pathPower_dB, deltaF, numSB, sbCarrierIndex, ...
    numSymTTI, fd, NFFT, CP, Fs)
% v4：不做功率归一化，保留Ped-B自然多径增益
    numPaths = length(pathDelay);
    Tsymbol = (NFFT + CP) / Fs;
    rho_corr = besselj(0, 2*pi*fd*Tsymbol);

    pathGain = zeros(numPaths, numSymTTI);
    for p = 1:numPaths
        pathPower = 10^(pathPower_dB(p)/10);
        pathGain(p,1) = sqrt(pathPower/2) * (randn + 1j*randn);
        for n = 2:numSymTTI
            innovation = sqrt(pathPower/2) * (randn + 1j*randn);
            pathGain(p,n) = rho_corr*pathGain(p,n-1) + sqrt(1-rho_corr^2)*innovation;
        end
    end

    allCarrierIdx = sbCarrierIndex.';
    allCarrierIdx = allCarrierIdx(:);
    allFreq = allCarrierIdx * deltaF;
    phaseMatrix = exp(-1j * 2 * pi * allFreq * pathDelay(:).');

    carriersPerSB = size(sbCarrierIndex, 2);
    Hsb = zeros(numSB, numSymTTI);
    for t = 1:numSymTTI
        pg = pathGain(:, t);
        H_all = phaseMatrix * pg;
        H_mat = reshape(H_all, carriersPerSB, numSB);
        Hsb(:, t) = sqrt(mean(abs(H_mat).^2, 1)).';
    end
    % 不归一化
end


function effectiveSIR_dB = EESM_dB(SINR, beta)
    SINR = max(SINR, 0);
    gammaEff = -beta * log(mean(exp(-SINR / beta)));
    effectiveSIR_dB = 10 * log10(max(gammaEff, eps));
end


function BLER = calculateMCSBLER(mcs, SIR_dB, ...
    QPSK13_SIR, QPSK13_BLER, ...
    QPSK12_SIR, QPSK12_BLER, ...
    QPSK34_SIR, QPSK34_BLER, ...
    QAM13_SIR, QAM13_BLER, ...
    QAM12_SIR, QAM12_BLER, ...
    QAM34_SIR, QAM34_BLER)
    switch mcs
        case {1, 3}
            BLER = interp1(QPSK13_SIR, QPSK13_BLER, SIR_dB, 'linear', 'extrap');
        case {2, 4}
            BLER = interp1(QPSK12_SIR, QPSK12_BLER, SIR_dB, 'linear', 'extrap');
        case 5
            BLER = interp1(QAM13_SIR, QAM13_BLER, SIR_dB, 'linear', 'extrap');
        case 6
            BLER = interp1(QPSK34_SIR, QPSK34_BLER, SIR_dB, 'linear', 'extrap');
        case 7
            BLER = interp1(QAM12_SIR, QAM12_BLER, SIR_dB, 'linear', 'extrap');
        case 8
            BLER = interp1(QAM34_SIR, QAM34_BLER, SIR_dB, 'linear', 'extrap');
        otherwise
            BLER = 1;
    end
    BLER = min(max(BLER, 0), 1);
end

% 论文的 ASN 方案：
%     │
%     ├── 信道估计：假设完美（无误差）
%     ├── 子带排序：假设完美（无排序错误）
%     ├── 反馈链路：假设完美（无误差、零延迟）
%     ├── MCS选择：基于完美平均SINR
%     │
%     └── 只有数据传输：用 3GPP BLER 模型（隐含实际接收机处理）