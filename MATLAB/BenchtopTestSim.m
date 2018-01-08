% July 2017
% Brown University Larson Lab
% Chester Kilfoyle
% Read in benchtop data, filter, apply PLL and bit slicing


%% Requires:
% trial and run
    tr = 21;
    rn = 4;
% IF guess, expects  ~0.5 - 1MHz accuracy
    IF = 26e6;
% Simulation duration (full dataset or est. # of bits)
    runFull = 0; %bool
    runBits = 2000;

%% ARW 9/8/17 rev: add path management, assumes working above three folders
% unzipped from GitHub
clc
addpath('Simulink');
addpath('Matlab');
addpath('Datasets');

%% ____________________________________________
% Open simulink system
open_system('BenchtopSkinpatch_SingleTone');

%% File Prep
% ARW 9/8/17: generate filenames and paths from trial and run numbers
    % set file paths
        % fp is metadata file
        % fn is data file
        % pathStr is the file path
    fp = ['Datasets' filesep sprintf('trial%d',tr) filesep sprintf('T%d.%d', tr, rn) ... 
            filesep [sprintf('T%d.%d', tr, rn) '.parameters.mat'] ];
    pathStr = fileparts(fp);
    fn = [pathStr filesep 'IQ.dat'];

    load(fp)% gets metadata

% ARW 9/8/17: process data
    % map interleaved binary file 
    % extract I and Q 
    % convert binary int16 to double 
    % scale to actual capture voltage levels
m = memmapfile(fn, 'Format', 'int16');
Iin = m.Data(1:2:end); Iin = double(Iin)*2^-15*(s.rangeADC);
Qin = m.Data(2:2:end); Qin = double(Qin)*2^-15*(s.rangeADC);

%% Simulation Parameters
% Sampling Freq and time vector
Fs = s.fs*1e6;
Ts = 1/Fs;
time = (0:Ts:Ts*(length(Iin) - 1))';

%FFT length
CFC_FFTlength = 512;

%bit/symbol rate estimate
Fbit = IF/3;
Tsym = 1/(Fbit);

%LP filter cutoff
Fcutoff = Fbit;

%Sim duration either full duration or fixed number of bits (approx.)
if(runFull)
    sim_time = time(end);
else
    sim_time = runBits/Fbit; % ARW 9/8/17: increment is samples for slx run, not bits
end

%create timeseries
cc = 1;%cos(-2*pi*time*IF);
ss = 0;%sin(-2*pi*time*IF);
I = timeseries(cc.*Iin - ss.*Qin,time);
Q = timeseries(cc.*Qin + ss.*Iin,time);

%% -- UNCOMMENT this to view data on spec. analyzer
h1 = dsp.SpectrumAnalyzer(1);
    h1.SampleRate=Fs;
    h1.SpectralAverages=50;
    h1.ReferenceLoad = 50;
    h1.FrequencySpan = 'Span and center frequency';
    h1.Span = 125e6;
    h1.CenterFrequency = 0;
    h1.YLimits = [-100 -60];
    h1.ChannelNames = {sprintf('IQ DDC at %1.0f MHz IF',1e-6*IF)};
    h1.ShowLegend = true;
    h1.Title = sprintf('IQ DDC spectrum at 50 Ohm for Trial %d, Run %d, LO %1.0f MHz',tr,rn,s.LO);

%%
iq = I.data + 1i*Q.data;
step(h1, iq)

%% filter
  I = idealfilter(I-mean(I),[15e6 45e6],'pass');
  Q = idealfilter(Q-mean(Q),[15e6 45e6],'pass');

%% Run the simulation
sim('BenchtopSkinpatch_SingleTone');

%% process output bits
bits = outBits.data(bitStrobe.data > 0);
bits = bits(:);

%31 bit PN sequence
bitseq = [1 1 0 0 0 1 1 1 1 1 0 0 1 1 0 1 0 0 1 0 0 0 0 1 0 1 0 1 1 1 0];

%compute error for all 31 phases of bit sequence
for i = 1:31
    seq = repmat([bitseq(i:end) bitseq(1:i-1)]',floor(length(bits)/31),1); 
    e(i) = sum(bits(1:length(seq)) == seq); 
end

%account for incorrect polarity
[vmin, imin] = min(e);
[vmax, imax] = max(e);

%choose best overall match
if(vmin < (length(seq) - vmax))
    phase = imin;
    errors = biterrors(bits, bitseq', phase);
else
    phase = imax;
    errors = biterrors(bits, bitseq', phase);
    errors = ~errors;
end

%bit error rate
BER = sum(errors)/length(errors)

%bit error rate ignoring first 10% of simulation to allow lock-on
os = round(length(errors)/10);
BER_offset = sum(errors(os:end))/length(errors(os:end))