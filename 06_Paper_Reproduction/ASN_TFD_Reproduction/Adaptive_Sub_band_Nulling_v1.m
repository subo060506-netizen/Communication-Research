%% ================================================================
%  ASN_TFD_Reproduction_v1.m
%
%  论文：
%  Adaptive Sub-band Nulling for OFDM-Based Wireless
%  Communication Systems
%
%  本程序分成两个完全独立的实验：
%
%  Part A：
%       复现论文 Fig.6 / Fig.7 / Fig.8
%       使用：
%           20个sub-band
%           独立Rayleigh fading
%
%  Part B：
%       重构论文 Fig.9
%       使用：
%           15个sub-band
%           每个sub-band 45个相邻子载波
%           1024 FFT
%           12 OFDM symbols / TTI
%           Pedestrian-B
%           MCS 1~8
%           Costas T-F mapping
%           EESM + 3GPP参考BLER模型
%
%  重要说明：
%
%  Fig.6~8：
%       可以严格按照论文Section III的容量公式复现。
%
%  Fig.9：
%       论文引用3GPP TR 25.892作为link-level simulator依据。
%       由于论文正文没有给出作者完整的原始链路级程序，
%       因此这里采用3GPP公开的EESM + AWGN BLER参考曲线
%       进行系统级重构。
%
% ================================================================

clear;
clc;
close all;

rng(2026);

%% ================================================================
%                     PART A
%       Fig.6 / Fig.7 / Fig.8
% ================================================================

fprintf('\n');
fprintf('========================================================\n');
fprintf('PART A：Fig.6 / Fig.7 / Fig.8 容量实验\n');
fprintf('========================================================\n');

%% ------------------------------------------------
% A-1. Fig.6~8参数
% -------------------------------------------------

numSB_A = 20;

% SNR范围
SNRdB_A = -20:1:30;

% Monte-Carlo次数
numRealizations_A = 20000;

% 固定置零数量
nullCases_A = [0 5 10 15];

% 保存容量
C_TFD = zeros(size(SNRdB_A));
C_ASN5 = zeros(size(SNRdB_A));
C_ASN10 = zeros(size(SNRdB_A));
C_ASN15 = zeros(size(SNRdB_A));

% Perfect ASN
C_Perfect = zeros(size(SNRdB_A));

% Water-Filling
C_WF = zeros(size(SNRdB_A));

% AWGN
C_AWGN = zeros(size(SNRdB_A));

% 保存平均最优剩余sub-band数量
Kopt_average = zeros(size(SNRdB_A));

%% ------------------------------------------------
% A-2. Monte-Carlo容量计算
% -------------------------------------------------

for isnr = 1:length(SNRdB_A)

    snr_dB = SNRdB_A(isnr);

    % 平均SNR
    rho = 10^(snr_dB/10);

    tempTFD = zeros(numRealizations_A,1);
    tempASN5 = zeros(numRealizations_A,1);
    tempASN10 = zeros(numRealizations_A,1);
    tempASN15 = zeros(numRealizations_A,1);

    tempPerfect = zeros(numRealizations_A,1);
    tempWF = zeros(numRealizations_A,1);

    tempK = zeros(numRealizations_A,1);

    for r = 1:numRealizations_A

        %% --------------------------------------------
        % 生成20个独立Rayleigh子带信道
        %
        % h_n ~ CN(0,1)
        % |h_n|^2服从指数分布
        % ---------------------------------------------

        h = (randn(numSB_A,1) + ...
             1j*randn(numSB_A,1))/sqrt(2);

        gain = abs(h).^2;

        %% --------------------------------------------
        % TFD
        %
        % 所有20个sub-band都使用
        % 总功率平均分配
        % ---------------------------------------------

        K = numSB_A;

        tempTFD(r) = ...
            sum(log2(1 + rho*gain/K));

        %% --------------------------------------------
        % ASN-5
        %
        % 保留最好的15个
        % ---------------------------------------------

        sortedGain = sort(gain,'descend');

        K = 15;

        selectedGain = sortedGain(1:K);

        tempASN5(r) = ...
            sum(log2(1 + rho*selectedGain/K));

        %% --------------------------------------------
        % ASN-10
        %
        % 保留最好的10个
        % ---------------------------------------------

        K = 10;

        selectedGain = sortedGain(1:K);

        tempASN10(r) = ...
            sum(log2(1 + rho*selectedGain/K));

        %% --------------------------------------------
        % ASN-15
        %
        % 保留最好的5个
        % ---------------------------------------------

        K = 5;

        selectedGain = sortedGain(1:K);

        tempASN15(r) = ...
            sum(log2(1 + rho*selectedGain/K));

        %% --------------------------------------------
        % Perfect ASN
        %
        % K = 1,...,20
        %
        % 对每一个K：
        %   选择最好的K个sub-band
        %   将总功率平均分给K个
        %   计算容量
        %
        % 最后：
        %
        % K* = argmax(C_K)
        % ---------------------------------------------

        capacity_K = zeros(numSB_A,1);

        for K = 1:numSB_A

            selectedGain = sortedGain(1:K);

            capacity_K(K) = ...
                sum(log2(1 + rho*selectedGain/K));

        end

        [tempPerfect(r),Kopt] = max(capacity_K);

        tempK(r) = Kopt;

        %% --------------------------------------------
        % Water-Filling
        %
        % 对20个信道进行最优功率分配
        %
        % C_WF = sum log2(1 + p_n*g_n*rho)
        %
        % 其中：
        %
        % sum(p_n) = 1
        % ---------------------------------------------

        tempWF(r) = ...
            waterFillingCapacity(gain,rho);

    end

    %% --------------------------------------------
    % Monte-Carlo平均
    % ---------------------------------------------

    C_TFD(isnr) = mean(tempTFD);
    C_ASN5(isnr) = mean(tempASN5);
    C_ASN10(isnr) = mean(tempASN10);
    C_ASN15(isnr) = mean(tempASN15);

    C_Perfect(isnr) = mean(tempPerfect);

    C_WF(isnr) = mean(tempWF);

    % AWGN容量
    %
    % AWGN中所有20个子带信道增益均为1
    %
    C_AWGN(isnr) = ...
        numSB_A * log2(1 + rho/numSB_A);

    Kopt_average(isnr) = mean(tempK);

    fprintf('Part A：SNR = %3d dB 完成\n',snr_dB);

end

%% ================================================================
% A-3. Fig.6
% ================================================================

figure;

plot(SNRdB_A,C_AWGN,'k','LineWidth',1.5);
hold on;

plot(SNRdB_A,C_WF,'b','LineWidth',1.5);
plot(SNRdB_A,C_TFD,'r--','LineWidth',1.5);
plot(SNRdB_A,C_Perfect,'g-.','LineWidth',1.5);

grid on;

xlabel('平均信噪比 (dB)');
ylabel('容量 (bits/sec/Hz)');

title('Fig.6：ASN与Water-Filling容量比较');

legend( ...
    'AWGN capacity', ...
    'Water-filling capacity', ...
    'Conventional TFD capacity', ...
    'Proposed ASN capacity', ...
    'Location','northwest');

%% ================================================================
% A-4. Fig.7
% ================================================================

figure;

plot(SNRdB_A,C_AWGN,'k','LineWidth',1.5);
hold on;

plot(SNRdB_A,C_WF,'b','LineWidth',1.5);

plot(SNRdB_A,C_TFD,'-o','LineWidth',1.2);
plot(SNRdB_A,C_ASN5,'--s','LineWidth',1.2);
plot(SNRdB_A,C_ASN10,'-.d','LineWidth',1.2);
plot(SNRdB_A,C_ASN15,':^','LineWidth',1.2);

grid on;

xlabel('平均信噪比 (dB)');
ylabel('容量 (bits/sec/Hz)');

title('Fig.7：不同置零数量下的ASN容量');

legend( ...
    'AWGN capacity', ...
    'Water-filling capacity', ...
    'TFD: 0 sub-band nulling', ...
    'ASN: 5 sub-band nulling', ...
    'ASN: 10 sub-band nulling', ...
    'ASN: 15 sub-band nulling', ...
    'Location','northwest');

%% ================================================================
% A-5. Fig.8
%
% 注意：
% 这里必须除以Water-Filling容量。
%
% 这是上一版最严重的错误之一。
% ================================================================

figure;

plot(SNRdB_A,C_Perfect./C_WF,'-o','LineWidth',1.5);
hold on;

plot(SNRdB_A,C_TFD./C_WF,'--x','LineWidth',1.5);
plot(SNRdB_A,C_ASN5./C_WF,'-.s','LineWidth',1.5);
plot(SNRdB_A,C_ASN10./C_WF,':d','LineWidth',1.5);
plot(SNRdB_A,C_ASN15./C_WF,'-^','LineWidth',1.5);

grid on;

ylim([0 1.2]);

xlabel('平均信噪比 (dB)');
ylabel('Normalized capacity over WF capacity');

title('Fig.8：ASN容量相对于Water-Filling容量的归一化结果');

legend( ...
    'ASN: optimal number of nulled sub-bands', ...
    'TFD: 0 sub-band nulling', ...
    'ASN: 5 sub-band nulling', ...
    'ASN: 10 sub-band nulling', ...
    'ASN: 15 sub-band nulling', ...
    'Location','southeast');

%% ================================================================
% A-6. 最优剩余sub-band数量
% ================================================================

figure;

plot(SNRdB_A,Kopt_average,'-o','LineWidth',1.5);

grid on;

xlabel('平均信噪比 (dB)');
ylabel('平均最优剩余sub-band数量 K^*');

title('Perfect ASN：最优剩余sub-band数量随SNR变化');

%% ================================================================
%                     PART B
%                    Fig.9
% ================================================================

fprintf('\n');
fprintf('========================================================\n');
fprintf('PART B：Fig.9链路级重构\n');
fprintf('========================================================\n');

%% ------------------------------------------------
% B-1. OFDM参数
% -------------------------------------------------

TTI = 2e-3;

NFFT = 1024;

Fs = 6.528e6;

CP = 64;

deltaF = 6.375e3;

numSymTTI = 12;

numSB_B = 15;

carriersPerSB = 45;

% 每个TTI固定置零数量
nullCases_B = [0 5 10 15];

% SNR范围
SNRdB_B = -15:1:20;

% Monte-Carlo次数
numRealizations_B = 3000;

%% ------------------------------------------------
% B-2. MCS
% -------------------------------------------------

MCS(1).modulation = 'QPSK';
MCS(1).M = 4;
MCS(1).rate = 1/3;
MCS(1).dataRate = 0.8e6;
MCS(1).payload = 1600;

MCS(2).modulation = 'QPSK';
MCS(2).M = 4;
MCS(2).rate = 1/2;
MCS(2).dataRate = 1.2e6;
MCS(2).payload = 2400;

MCS(3).modulation = 'QPSK';
MCS(3).M = 4;
MCS(3).rate = 1/3;
MCS(3).dataRate = 2.4e6;
MCS(3).payload = 4800;

MCS(4).modulation = 'QPSK';
MCS(4).M = 4;
MCS(4).rate = 1/2;
MCS(4).dataRate = 3.6e6;
MCS(4).payload = 7200;

MCS(5).modulation = '16QAM';
MCS(5).M = 16;
MCS(5).rate = 1/3;
MCS(5).dataRate = 4.8e6;
MCS(5).payload = 9600;

MCS(6).modulation = 'QPSK';
MCS(6).M = 4;
MCS(6).rate = 3/4;
MCS(6).dataRate = 5.4e6;
MCS(6).payload = 10800;

MCS(7).modulation = '16QAM';
MCS(7).M = 16;
MCS(7).rate = 1/2;
MCS(7).dataRate = 7.2e6;
MCS(7).payload = 14400;

MCS(8).modulation = '16QAM';
MCS(8).M = 16;
MCS(8).rate = 3/4;
MCS(8).dataRate = 10.8e6;
MCS(8).payload = 21600;

%% ------------------------------------------------
% B-3. 3GPP Costas序列
%
% TR 25.892参数集2：
%
% TFP0 =
% [13 5 3 9 2 14 11 15 4 12 7 10]
%
% 其他TFP：
% 由TFP0在频域进行不同的循环移位得到。
% -------------------------------------------------

TFP0 = [13 5 3 9 2 14 11 15 4 12 7 10];

numTFPatterns = 15;

TFPatterns = zeros(numTFPatterns,numSymTTI);

for p = 0:numTFPatterns-1

    TFPatterns(p+1,:) = ...
        mod(TFP0 - 1 + p,numSB_B) + 1;

end

%% ------------------------------------------------
% B-4. 显示TFP0和TFP1
% -------------------------------------------------

fprintf('\nTFP0：\n');
disp(TFPatterns(1,:));

fprintf('TFP1：\n');
disp(TFPatterns(2,:));

%% ------------------------------------------------
% B-5. 画15个正交TF pattern
% -------------------------------------------------

figure;

TFGrid = zeros(numSB_B,numSymTTI);

for p = 1:numTFPatterns

    for t = 1:numSymTTI

        sb = TFPatterns(p,t);

        TFGrid(sb,t) = p;

    end

end

imagesc(TFGrid);

xlabel('OFDM symbol');
ylabel('sub-band');

title('3GPP Costas T-F pattern集合');

colorbar;

%% ------------------------------------------------
% B-6. Pedestrian-B参数
% -------------------------------------------------

pathDelay_ns = [0 200 800 1200 2300 3700];

pathPower_dB = [0 -0.9 -4.9 -8 -7.8 -23.9];

pathDelay = pathDelay_ns * 1e-9;

% 3 km/h
velocity = 3/3.6;

% 假设载频2 GHz
fc = 2e9;

c = 3e8;

fd = velocity/c*fc;

%% ------------------------------------------------
% B-7. 建立15个sub-band的45个子载波
% -------------------------------------------------

centeredCarriers = -337:337;

sbCarrierIndex = zeros(numSB_B,carriersPerSB);

for sb = 1:numSB_B

    startIndex = ...
        (sb-1)*carriersPerSB + 1;

    endIndex = ...
        sb*carriersPerSB;

    sbCarrierIndex(sb,:) = ...
        centeredCarriers(startIndex:endIndex);

end

%% ------------------------------------------------
% B-8. EESM beta参数
%
% 3GPP TR 25.892给出的参考beta：
%
% QPSK 1/3 = 1.49
% QPSK 1/2 = 1.57
% QPSK 3/4 = 1.69
%
% 16QAM 1/3 = 3.36
% 16QAM 1/2 = 4.56
% 16QAM 3/4 = 7.33
%
% -------------------------------------------------

beta = zeros(8,1);

beta(1) = 1.49;
beta(2) = 1.57;
beta(3) = 1.49;
beta(4) = 1.57;
beta(5) = 3.36;
beta(6) = 1.69;
beta(7) = 4.56;
beta(8) = 7.33;

%% ------------------------------------------------
% B-9. 3GPP AWGN参考BLER曲线
%
% 这里只列出论文MCS真正需要的：
%
% QPSK 1/3
% QPSK 1/2
% QPSK 3/4
% 16QAM 1/3
% 16QAM 1/2
% 16QAM 3/4
%
% 这些曲线来自TR 25.892的参考AWGN TTI BLER曲线。
% -------------------------------------------------

% QPSK 1/3
BLER_QPSK13_SIR = ...
    [-1.94 -1.74 -1.54 -1.34 -1.14];

BLER_QPSK13 = ...
    [1.00 9.95e-1 8.03e-1 1.79e-1 4.10e-3];

% QPSK 1/2
BLER_QPSK12_SIR = ...
    [0.62 0.82 1.02 1.22 1.32];

BLER_QPSK12 = ...
    [1.00 9.45e-1 3.95e-1 2.76e-2 4.13e-3];

% QPSK 3/4
BLER_QPSK34_SIR = ...
    [3.98 4.18 4.38 4.58 4.78];

BLER_QPSK34 = ...
    [1.00 9.40e-1 3.98e-1 3.97e-2 3.30e-3];

% 16QAM 1/3
BLER_16QAM13_SIR = ...
    [3.06 3.26 3.46 3.56 3.66];

BLER_16QAM13 = ...
    [1.00 9.14e-1 2.58e-1 5.72e-2 7.15e-3];

% 16QAM 1/2
BLER_16QAM12_SIR = ...
    [5.82 6.02 6.22 6.42 6.52];

BLER_16QAM12 = ...
    [1.00 9.94e-1 5.89e-1 4.49e-2 5.70e-3];

% 16QAM 3/4
BLER_16QAM34_SIR = ...
    [10.18 10.38 10.58 10.78 10.98];

BLER_16QAM34 = ...
    [1.00 8.95e-1 2.79e-1 2.00e-2 1.57e-3];

%% ------------------------------------------------
% B-10. 保存吞吐量
% -------------------------------------------------

throughput = zeros( ...
    length(nullCases_B), ...
    length(SNRdB_B));

%% ================================================================
% B-11. Fig.9主循环
% ================================================================

for isnr = 1:length(SNRdB_B)

    snr_dB = SNRdB_B(isnr);

    rho = 10^(snr_dB/10);

    for icase = 1:length(nullCases_B)

        numNull = nullCases_B(icase);

        totalRate = 0;

        for r = 1:numRealizations_B

            %% ----------------------------------------------------
            % 生成Pedestrian-B频域信道
            %
            % 输出：
            %
            % Hcarrier：
            % 1024个FFT位置中的实际频域响应
            %
            % 这里最终关心15个sub-band。
            % -----------------------------------------------------

            Hsb = generatePedBSubbandChannel_v3( ...
                pathDelay, ...
                pathPower_dB, ...
                deltaF, ...
                numSB_B, ...
                sbCarrierIndex, ...
                numSymTTI, ...
                fd);

            %% ----------------------------------------------------
            % 计算每个sub-band在12个OFDM symbol中的平均增益
            % -----------------------------------------------------

            sbGain = mean(abs(Hsb).^2,2);

            %% ----------------------------------------------------
            % 固定ASN：
            %
            % 0个null：
            %   所有15个SB
            %
            % 5个null：
            %   保留最好的10个
            %
            % 10个null：
            %   保留最好的5个
            %
            % 15个null：
            %   没有数据SB
            % -----------------------------------------------------

            if numNull == 15

                totalRate = totalRate + 0;

                continue;

            end

            numKeep = numSB_B - numNull;

            [~,order] = sort(sbGain,'descend');

            keepSB = order(1:numKeep);

            %% ----------------------------------------------------
            % 计算剩余sub-band上的平均接收SNR
            % -----------------------------------------------------

            selectedGain = sbGain(keepSB);

            % ASN把原来被置零sub-band的功率
            % 重新分配给剩余sub-band。
            %
            % 因此：
            %
            % 每个剩余sub-band功率 ∝ 1/K
            %
            % 平均SNR：
            % rho * gain / K
            %
            K = numKeep;

            sbSINR = rho * selectedGain / K;

            %% ----------------------------------------------------
            % MCS选择
            %
            % 论文：
            %
            % BS根据MS反馈的average received SINR
            % 选择合适的MCS。
            %
            % 这里通过所有候选MCS的BLER模型寻找：
            %
            % “满足BLER <= 10%的最高速率MCS”
            %
            % 这比上一版人工写固定SINR门限更接近
            % 3GPP的系统级仿真思想。
            % -----------------------------------------------------

            bestMCS = 0;

            bestRate = 0;

            for m = 1:8

                effectiveSIR_dB = ...
                    EESM_dB( ...
                    sbSINR, ...
                    beta(m));

                currentBLER = ...
                    calculateMCSBLER( ...
                    m, ...
                    effectiveSIR_dB, ...
                    BLER_QPSK13_SIR, ...
                    BLER_QPSK13, ...
                    BLER_QPSK12_SIR, ...
                    BLER_QPSK12, ...
                    BLER_QPSK34_SIR, ...
                    BLER_QPSK34, ...
                    BLER_16QAM13_SIR, ...
                    BLER_16QAM13, ...
                    BLER_16QAM12_SIR, ...
                    BLER_16QAM12, ...
                    BLER_16QAM34_SIR, ...
                    BLER_16QAM34);

                % BLER <= 10%认为当前MCS可以使用
                if currentBLER <= 0.10

                    if MCS(m).dataRate > bestRate

                        bestRate = MCS(m).dataRate;

                        bestMCS = m;

                    end

                end

            end

            %% ----------------------------------------------------
            % MCS1/2特殊资源结构
            %
            % 论文明确：
            %
            % MCS1/2：
            %   5个OFDM data units
            %   15个sub-band中只有5个用于一个用户
            %   每个OFDM symbol的5个sub-band
            %   按Costas sequence变化
            %
            % MCS3~8：
            %   全部15个sub-band
            %
            % 因此这里对MCS1/2增加5/15的资源占用因子。
            %
            % 注意：
            % 这是系统级资源占用模型，而不是逐bit Turbo
            % 编码器的完整实现。
            % -----------------------------------------------------

            if bestMCS == 1 || bestMCS == 2

                resourceFactor = 5/15;

            else

                resourceFactor = 1;

            end

            %% ----------------------------------------------------
            % 计算有效吞吐量
            %
            % 采用：
            %
            % throughput =
            %       MCS data rate
            %       × (1-BLER)
            %       × resourceFactor
            %
            % 作为系统级重构指标。
            % -----------------------------------------------------

            if bestMCS == 0

                currentRate = 0;

            else

                effectiveSIR_dB = ...
                    EESM_dB(sbSINR,beta(bestMCS));

                currentBLER = ...
                    calculateMCSBLER( ...
                    bestMCS, ...
                    effectiveSIR_dB, ...
                    BLER_QPSK13_SIR, ...
                    BLER_QPSK13, ...
                    BLER_QPSK12_SIR, ...
                    BLER_QPSK12, ...
                    BLER_QPSK34_SIR, ...
                    BLER_QPSK34, ...
                    BLER_16QAM13_SIR, ...
                    BLER_16QAM13, ...
                    BLER_16QAM12_SIR, ...
                    BLER_16QAM12, ...
                    BLER_16QAM34_SIR, ...
                    BLER_16QAM34);

                currentRate = ...
                    MCS(bestMCS).dataRate ...
                    * (1-currentBLER) ...
                    * resourceFactor;

            end

            totalRate = totalRate + currentRate;

        end

        %% --------------------------------------------------------
        % Monte-Carlo平均
        % --------------------------------------------------------

        throughput(icase,isnr) = ...
            totalRate / numRealizations_B;

    end

    fprintf('Part B：SNR = %3d dB 完成\n',snr_dB);

end

%% ================================================================
% B-12. Fig.9
% ================================================================

figure;

plot( ...
    SNRdB_B, ...
    throughput(1,:)/1e6, ...
    '-o', ...
    'LineWidth',1.5);

hold on;

plot( ...
    SNRdB_B, ...
    throughput(2,:)/1e6, ...
    '--s', ...
    'LineWidth',1.5);

plot( ...
    SNRdB_B, ...
    throughput(3,:)/1e6, ...
    '-.d', ...
    'LineWidth',1.5);

plot( ...
    SNRdB_B, ...
    throughput(4,:)/1e6, ...
    ':^', ...
    'LineWidth',1.5);

grid on;

xlabel('平均信噪比 (dB)');
ylabel('Throughput (Mbps)');

title('Fig.9：TFD与ASN吞吐量比较');

legend( ...
    'TFD: 0 sub-band nulling', ...
    'ASN: 5 sub-band nulling', ...
    'ASN: 10 sub-band nulling', ...
    'ASN: 15 sub-band nulling', ...
    'Location','northwest');

%% ================================================================
% B-13. 输出结果
% ================================================================

fprintf('\n');
fprintf('========================================================\n');
fprintf('Fig.9吞吐量结果\n');
fprintf('========================================================\n');

for icase = 1:length(nullCases_B)

    fprintf('\n置零 %d 个sub-band：\n', ...
        nullCases_B(icase));

    fprintf( ...
        '最大吞吐量 = %.3f Mbps\n', ...
        max(throughput(icase,:))/1e6);

end


%% ================================================================
%                         局部函数
% ================================================================


function capacity = waterFillingCapacity(gain,rho)

% ================================================================
% Water-Filling容量
%
% 给定：
%
%   gain(n) = |h_n|^2
%
%   rho = 平均SNR
%
% 求：
%
%   max sum log2(1 + rho*p_n*gain_n)
%
% s.t.
%
%   sum(p_n)=1
%   p_n>=0
%
% ================================================================

    N = length(gain);

    gain = gain(:);

    % 将信道增益从大到小排列
    [gainSort,~] = sort(gain,'descend');

    activeK = N;

    while true

        invGain = 1 ./ gainSort(1:activeK);

        waterLevel = ...
            (N + rho*sum(invGain)) / ...
            (rho*activeK);

        power = ...
            waterLevel - ...
            1./(rho*gainSort(1:activeK));

        if all(power >= 0) || activeK == 1

            break;

        end

        activeK = activeK - 1;

    end

    power = max(power,0);

    capacity = ...
        sum(log2( ...
        1 + rho*power.*gainSort(1:activeK)));

end


function Hsb = generatePedBSubbandChannel_v3( ...
    pathDelay, ...
    pathPower_dB, ...
    deltaF, ...
    numSB, ...
    sbCarrierIndex, ...
    numSymTTI, ...
    fd)

% ================================================================
% 生成Pedestrian-B频率响应
%
% 输出：
%
%   Hsb：
%       numSB × numSymTTI
%
% 每一个元素代表：
%
%       一个sub-band
%       在一个OFDM symbol上的等效信道增益
%
% ================================================================

    numPaths = length(pathDelay);

    %% ------------------------------------------------------------
    % OFDM symbol之间的时间间隔
    % -------------------------------------------------------------

    Tsymbol = 1/(12*6.528e6/1024);

    rho = besselj(0,2*pi*fd*Tsymbol);

    %% ------------------------------------------------------------
    % 生成时变路径系数
    % -------------------------------------------------------------

    pathGain = zeros(numPaths,numSymTTI);

    for p = 1:numPaths

        pathPower = 10^(pathPower_dB(p)/10);

        pathGain(p,1) = ...
            sqrt(pathPower/2) * ...
            (randn + 1j*randn);

        for n = 2:numSymTTI

            innovation = ...
                sqrt(pathPower/2) * ...
                (randn + 1j*randn);

            pathGain(p,n) = ...
                rho*pathGain(p,n-1) + ...
                sqrt(1-rho^2)*innovation;

        end

    end

    %% ------------------------------------------------------------
    % 计算频域响应
    % -------------------------------------------------------------

    Hsb = zeros(numSB,numSymTTI);

    for sb = 1:numSB

        carrierIndex = sbCarrierIndex(sb,:);

        frequency = carrierIndex * deltaF;

        for t = 1:numSymTTI

            Hcarrier = zeros(1,length(frequency));

            for k = 1:length(frequency)

                Htemp = 0;

                for p = 1:numPaths

                    Htemp = Htemp + ...
                        pathGain(p,t) * ...
                        exp(-1j*2*pi* ...
                        frequency(k)*pathDelay(p));

                end

                Hcarrier(k) = Htemp;

            end

            % sub-band等效信道增益
            %
            % 论文Section III的核心假设：
            % 一个sub-band内信道增益近似相同。
            Hsb(sb,t) = ...
                sqrt(mean(abs(Hcarrier).^2));

        end

    end

end


function effectiveSIR_dB = EESM_dB(SINR,beta)

% ================================================================
% EESM：
%
%   gamma_eff =
%       -beta * ln(
%           mean(exp(-gamma/beta))
%       )
%
% gamma使用线性值。
%
% 最后转换为dB。
% ================================================================

    SINR = max(SINR,0);

    gammaEff = ...
        -beta * log(mean( ...
        exp(-SINR/beta)));

    effectiveSIR_dB = ...
        10*log10(max(gammaEff,eps));

end


function BLER = calculateMCSBLER( ...
    mcs, ...
    SIR_dB, ...
    QPSK13_SIR, ...
    QPSK13_BLER, ...
    QPSK12_SIR, ...
    QPSK12_BLER, ...
    QPSK34_SIR, ...
    QPSK34_BLER, ...
    QAM13_SIR, ...
    QAM13_BLER, ...
    QAM12_SIR, ...
    QAM12_BLER, ...
    QAM34_SIR, ...
    QAM34_BLER)

% ================================================================
% 根据MCS调用对应的3GPP参考AWGN BLER曲线
%
% 论文使用的8种MCS中：
%
% MCS1 = QPSK 1/3
% MCS2 = QPSK 1/2
% MCS3 = QPSK 1/3
% MCS4 = QPSK 1/2
% MCS5 = 16QAM 1/3
% MCS6 = QPSK 3/4
% MCS7 = 16QAM 1/2
% MCS8 = 16QAM 3/4
%
% 对于表中没有直接提供的MCS，只能按照相同
% modulation/code-rate link mode进行映射。
% ================================================================

    switch mcs

        case {1,3}

            BLER = interp1( ...
                QPSK13_SIR, ...
                QPSK13_BLER, ...
                SIR_dB, ...
                'linear','extrap');

        case {2,4}

            BLER = interp1( ...
                QPSK12_SIR, ...
                QPSK12_BLER, ...
                SIR_dB, ...
                'linear','extrap');

        case 5

            BLER = interp1( ...
                QAM13_SIR, ...
                QAM13_BLER, ...
                SIR_dB, ...
                'linear','extrap');

        case 6

            BLER = interp1( ...
                QPSK34_SIR, ...
                QPSK34_BLER, ...
                SIR_dB, ...
                'linear','extrap');

        case 7

            BLER = interp1( ...
                QAM12_SIR, ...
                QAM12_BLER, ...
                SIR_dB, ...
                'linear','extrap');

        case 8

            BLER = interp1( ...
                QAM34_SIR, ...
                QAM34_BLER, ...
                SIR_dB, ...
                'linear','extrap');

        otherwise

            BLER = 1;

    end

    BLER = min(max(BLER,0),1);

end