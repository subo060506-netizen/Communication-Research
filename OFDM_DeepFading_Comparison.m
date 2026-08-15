%% OFDM_DeepFading_Scan.m
% 深衰落程度扫描实验
% 比较不同两径信道下 ZF 与 MMSE 的 BER 性能

clear;
clc;
close all;

%% ==================== 参数设置 ====================

numBits = 2^20;          % 比特数
modOrder = 16;           % 16-QAM
numCarr = 64;            % OFDM 子载波数
cycPrefLen = 5;          % 循环前缀长度

SNR = 0:2:30;             % SNR 扫描范围

% 不同程度的深衰落信道
alphaList = [0.50, 0.80, 0.95, 0.99];

% 保存 BER
BER_ZF = zeros(length(alphaList), length(SNR));
BER_MMSE = zeros(length(alphaList), length(SNR));

% 保存每个信道的最小频响
minH = zeros(length(alphaList),1);

%% ==================== 产生发送信号 ====================

% 随机比特
srcBits = randi([0 1], numBits, 1);

% 16-QAM 调制
qamModOut = qammod(srcBits, modOrder, InputType="bit");

% 整理成 OFDM 符号
qamModOut = reshape(qamModOut, numCarr, []);

% IFFT
ifftModOut = ifft(qamModOut, numCarr, 1);

% 添加 CP
cp = ifftModOut(end-cycPrefLen+1:end, :);

ofdmModOut = [cp; ifftModOut];

% 串行化
ofdmModOut = ofdmModOut(:);

%% ==================== 开始扫描不同信道 ====================

for a = 1:length(alphaList)

    alpha = alphaList(a);

    % 当前信道
    mpChan = [1; -alpha];

    % 计算真实信道频率响应
    Htrue = fft(mpChan, numCarr);

    % 记录最小信道增益
    minH(a) = min(abs(Htrue));

    fprintf("\nalpha = %.2f\n", alpha);
    fprintf("最小 |H[k]| = %.4f\n", minH(a));

    %% ==================== SNR 扫描 ====================

    for i = 1:length(SNR)

        snr = SNR(i);

        %% ---------- 多径信道 ----------

        chanOut = filter(mpChan, 1, ofdmModOut);

        %% ---------- 产生 AWGN ----------

        % 原始 QAM 符号平均功率
        qamPower = mean(abs(qamModOut(:)).^2);

        % SNR 转换为线性值
        snrLinear = 10^(snr/10);

        % 频域噪声方差
        noiseVarFreq = qamPower / snrLinear;

        % IFFT 后时域信号功率缩小 numCarr 倍
        noiseVarTime = noiseVarFreq / numCarr;

        % 复高斯白噪声
        noise = sqrt(noiseVarTime/2) * ...
            (randn(size(chanOut)) + 1j*randn(size(chanOut)));

        % 加噪
        rx = chanOut + noise;

        %% ---------- OFDM 解调 ----------

        % 恢复 OFDM 符号矩阵
        rx = reshape(rx, numCarr + cycPrefLen, []);

        % 去 CP
        rxNoCP = rx(cycPrefLen+1:end, :);

        % FFT
        fftModOut = fft(rxNoCP, numCarr, 1);

        %% ==================== ZF ====================

        eqOut_ZF = fftModOut ./ Htrue;

        % QAM 解调
        qamDeModOut_ZF = qamdemod( ...
            eqOut_ZF, ...
            modOrder, ...
            OutputType="bit");

        % 转成列向量
        rxBits_ZF = qamDeModOut_ZF(:);

        % BER
        BER_ZF(a,i) = ...
            nnz(rxBits_ZF ~= srcBits) / numBits;

        %% ==================== MMSE ====================

        % MMSE 权重
        noiseVarRatio = 1 / snrLinear;

        mmseWeight = conj(Htrue) ./ ...
            (abs(Htrue).^2 + noiseVarRatio);

        % MMSE 均衡
        eqOut_MMSE = fftModOut .* mmseWeight;

        % QAM 解调
        qamDeModOut_MMSE = qamdemod( ...
            eqOut_MMSE, ...
            modOrder, ...
            OutputType="bit");

        % 转成列向量
        rxBits_MMSE = qamDeModOut_MMSE(:);

        % BER
        BER_MMSE(a,i) = ...
            nnz(rxBits_MMSE ~= srcBits) / numBits;

    end
end

%% =========================================================
%  绘制不同深衰落信道的频率响应
% ==========================================================

figure;

for a = 1:length(alphaList)

    alpha = alphaList(a);

    mpChan = [1; -alpha];

    Htrue = fft(mpChan, numCarr);

    plot(abs(Htrue), "o-");
    hold on;

end

grid on;

xlabel("Subcarrier");
ylabel("|H[k]|");

legend( ...
    "\alpha = 0.50", ...
    "\alpha = 0.80", ...
    "\alpha = 0.95", ...
    "\alpha = 0.99");

title("Frequency Response Under Different Deep-Fading Levels");

%% =========================================================
%  绘制 ZF BER
% ==========================================================

figure;

for a = 1:length(alphaList)

    semilogy( ...
        SNR, ...
        BER_ZF(a,:), ...
        "o-");

    hold on;

end

grid on;

xlabel("SNR (dB)");
ylabel("BER");

legend( ...
    "\alpha = 0.50", ...
    "\alpha = 0.80", ...
    "\alpha = 0.95", ...
    "\alpha = 0.99");

title("ZF BER Under Different Deep-Fading Levels");

%% =========================================================
%  绘制 MMSE BER
% ==========================================================

figure;

for a = 1:length(alphaList)

    semilogy( ...
        SNR, ...
        BER_MMSE(a,:), ...
        "o-");

    hold on;

end

grid on;

xlabel("SNR (dB)");
ylabel("BER");

legend( ...
    "\alpha = 0.50", ...
    "\alpha = 0.80", ...
    "\alpha = 0.95", ...
    "\alpha = 0.99");

title("MMSE BER Under Different Deep-Fading Levels");

%% =========================================================
%  ZF 与 MMSE 对比
% ==========================================================

figure;

for a = 1:length(alphaList)

    semilogy( ...
        SNR, ...
        BER_ZF(a,:), ...
        "o-");

    hold on;

    semilogy( ...
        SNR, ...
        BER_MMSE(a,:), ...
        "x--");

end

grid on;

xlabel("SNR (dB)");
ylabel("BER");

legend( ...
    "ZF, \alpha=0.50", ...
    "MMSE, \alpha=0.50", ...
    "ZF, \alpha=0.80", ...
    "MMSE, \alpha=0.80", ...
    "ZF, \alpha=0.95", ...
    "MMSE, \alpha=0.95", ...
    "ZF, \alpha=0.99", ...
    "MMSE, \alpha=0.99");

title("ZF vs MMSE Under Different Deep-Fading Levels");

%% =========================================================
%  输出最小信道增益
% ==========================================================

fprintf("\n====================================\n");
fprintf("不同信道的最小频率响应\n");
fprintf("====================================\n");

for a = 1:length(alphaList)

    fprintf( ...
        "alpha = %.2f   min|H| = %.4f\n", ...
        alphaList(a), ...
        minH(a));

end