%%SINR vs SNR
%SNR=Ps/Pn描述信号与噪声的功率关系
%SINR=Ps/（Pn+Pi) S:Signal  I:Interference N:Noise
%在这篇 ASN 论文里,SINR 的重要性在于：同一个 OFDM symbol 的不同频率资源,信道不同,SINR不同

%MCS:Modulation and Coding Scheme调制方法+信道编码方式
%SINR->MCS,即AMC,信道好,16/64QAMA,高码率,每个资源传更多bit;信道差,QPSK,低码率
%AMC:Adaptive Modulation and Coding

%BLER:Block Error Rate 数据包传输错误率
%在译码后进行CRC(循环冗余校验),看是否传错,再计算误码率

%EESM：Effective Exponential SNR Mapping
%把不同sub_band的SINR=[y1 y2 y3 …yn]压缩成y^eff,用对数以此来确定MCS

%TFD=Time_Frequency Domain
%(t,f)二维资源映射

%Costas sequence:一种资源映射规则

%Turbo Coding一种信道编码

%HARQ重传机制

%Throughput吞吐量实际系统有效传输速率
%除信道容量外,还考虑MCS、信道编码、重传等影响,一般Troughput<capacity
%Throughput与R^MCS*(1-BLER) R^MCS为MCS确定的数据传输效率