# Communication-Research
This repository records my undergraduate research training in communication signal processing, with OFDM as the main research platform.

The project focuses on gradually developing the ability to independently implement, understand, analyze, and improve communication algorithms using MATLAB.

---
## 1. Research Roadmap

The current research path is:

OFDM Fundamentals  
→ LS Channel Estimation  
→ MMSE Equalization  
→ LS / MMSE Performance Comparison  
→ LMMSE Channel Estimation  
→ Deep-Fading Analysis  
→ Perfect CSI Benchmark  
→ Literature Review  
→ Paper Reproduction  
→ Algorithm Improvement  
→ Research Paper

---
## 2. Repository Structure

Communication-Research/
│
├── README.md
│
├── 01_OFDM_Basics/
│   └── OFDM_without_ofdmmod.m
│
├── 02_OFDM_LS/
│   └── OFDM_LS_BER.m
│
├── 03_OFDM_MMSE/
│   └── OFDM_MMSE_BER.m
│
├── 04_OFDM_LS_MMSE/
│   ├── OFDM_4Methods_Comparison.m
│   ├── OFDM_DeepFading_Comparison.m
│   └── OFDM_LS_MMSE_BER.m
│
└── 05_OFDM_LMMSE/
    ├── OFDM_LMMSE_DeepFade.m
    └── OFDM_LMMSE_MMSE_BER.m

---
## 3. Research Progress

01 OFDM Basics

Folder: 01_OFDM_Basics/

File: OFDM_without_ofdmmod.m

Main contents:

QAM modulation and demodulation
OFDM symbol generation
IFFT and FFT
Cyclic prefix
Multipath channel
AWGN
Basic OFDM BER simulation

The purpose of this stage was to understand the complete OFDM transmission and reception process and gradually move away from relying entirely on high-level communication toolbox functions.

02 LS Channel Estimation

Folder: 02_OFDM_LS/

File: OFDM_LS_BER.m

Main contents:

Pilot-based channel estimation
LS channel estimation
Channel frequency response
Channel equalization
BER performance analysis

This stage established the basic pilot-based channel estimation framework for later LMMSE research.

03 MMSE Equalization

Folder: 03_OFDM_MMSE/

File: OFDM_MMSE_BER.m

Main contents:

ZF equalization
MMSE equalization
Noise amplification
BER performance comparison

The purpose of this stage was to understand how channel fading and noise affect OFDM equalization and BER performance.

04 LS / MMSE Comparison

Folder: 04_OFDM_LS_MMSE/

Files:

OFDM_4Methods_Comparison.m
OFDM_DeepFading_Comparison.m
OFDM_LS_MMSE_BER.m

Main contents:

LS channel estimation
ZF equalization
MMSE equalization
Different channel estimation and equalization combinations
BER performance comparison
Initial deep-fading experiments

This stage was used to compare different OFDM receiver methods and investigate the relationship between channel fading, channel estimation, equalization, and BER.

05 LMMSE Channel Estimation

Folder: 05_OFDM_LMMSE/

Files:

OFDM_LMMSE_MMSE_BER.m
OFDM_LMMSE_DeepFade.m

This is the current research stage.

Main contents:

LMMSE channel estimation
Channel correlation matrix estimation
Monte Carlo estimation of channel statistics
LS and LMMSE BER comparison
Deep-fading channel construction
Deep-fading frequency response analysis
Perfect CSI benchmark

The main purpose of this stage is to investigate whether LMMSE can effectively reduce the performance loss caused by channel estimation errors under deep-fading conditions.

---

## 4. Deep-Fading Experiment

The current deep-fading experiment uses a normalized two-path channel.

The relative amplitude of the second path is controlled by the parameter alpha.

Four different fading conditions are investigated:

alpha = 0.50
alpha = 0.80
alpha = 0.95
alpha = 0.99

The experiment contains three main parts:

Analyze the frequency response under different fading conditions.
Compare LS and LMMSE channel estimation through BER performance.
Introduce Perfect CSI as an ideal benchmark.
Main observations

As alpha increases, the frequency response develops deeper frequency-selective fading.

The BER results show that LMMSE performs better than LS under the tested conditions.

The LMMSE BER performance is also very close to the Perfect CSI benchmark.

This indicates that LMMSE can effectively reduce the channel estimation error in the current simulation system.

However, under severe deep fading, even Perfect CSI cannot completely eliminate the BER degradation.

Therefore, the remaining performance limitation is mainly related to the deep fading of the channel itself rather than the channel estimation error.

---

## 5. Current Research Problem

The current experiments suggest that further improving the LMMSE channel estimator may not be the most effective way to solve the remaining performance degradation.

The research problem is therefore being extended from:

"How to improve channel estimation?"

to:

"How to mitigate the performance degradation caused by deep-faded OFDM subcarriers?"

Possible research directions include:

Deep-fade subcarrier detection
Adaptive subcarrier selection
Adaptive modulation
Adaptive power allocation
Joint subcarrier and power adaptation
Deep-fade-aware equalization

These directions will first be investigated through literature review and paper reproduction before developing a new algorithm.

---

## 6. Current Status

Completed
 OFDM basic transmission and reception
 QAM modulation and demodulation
 IFFT and FFT
 Cyclic prefix
 Multipath channel simulation
 AWGN simulation
 LS channel estimation
 ZF equalization
 MMSE equalization
 LS and MMSE BER experiments
 LMMSE channel estimation
 Channel correlation matrix analysis
 Monte Carlo estimation of channel statistics
 LS and LMMSE BER comparison
 Deep-fading frequency response analysis
 Deep-fading BER experiment
 Perfect CSI benchmark

Next Steps
 Literature review
 Select representative papers
 Reproduce existing algorithms
 Analyze limitations of existing methods
 Identify a research gap
 Design an improved algorithm
 Conduct comparative simulations
 Analyze experimental results
 Prepare a research paper

 ---

## 7. Tools

The main tools used in this project are:

MATLAB
Git
GitHub

MATLAB is primarily used for communication system simulation, algorithm implementation, and performance evaluation.

Git and GitHub are used to manage the research code, experiments, and development history.

---

## 8. Research Goal

The long-term goal of this project is to develop the ability to independently conduct communication signal processing research.

The research process will gradually transition from:

Learning existing methods
→ Implementing existing methods
→ Understanding existing methods
→ Reproducing research results
→ Identifying limitations
→ Designing improved methods
→ Evaluating proposed methods
→ Writing a research paper

The ultimate goal is to obtain complete undergraduate research experience and develop a research topic and paper based on OFDM, LMMSE, and deep-fading-related problems.

---