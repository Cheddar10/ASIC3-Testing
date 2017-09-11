% July 2017
% Brown University Larson Lab
% Chester Kilfoyle
% Read in benchtop data, filter, apply PLL and bit slicing


%% Requires:
% trial and run
    tr = 9;
    rn = 1;
% IF guess
    IF = 31.6e6;

%% ARW 9/8/17 rev: add path management, assumes working above three folders
% unzipped from GitHub
clc
addpath('Simulink');
addpath('Matlab');
addpath('Datasets');

%% ____________________________________________
% Open simulink system
open_system('BenchtopSkinpatch_SingleTone');

%% -- EDIT THIS SECTION --
% Read in data and produce two data arrays Iin, Qin
% These arrays should represent voltage levels directly, eg. not int16
                % tr = 'T9.1';
                % fp = strcat('Datasets\trial9\',tr,'\');
                % fn = 'IQ.dat';
                % load(strcat(fp,tr,'.parameters.mat'));

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
% -----------------------

%% -- EDIT THIS SECTION --
% Sampling Freq, Tx Freq, and IF estimate (requires ~.5-1MHz accuracy)
Fs = s.fs*1e6;
RF = s.LO*1e6;

% -----------------------

% a few more sim parameters
Ts = 1/Fs;
Fbit = IF/3;
Tsym = 1/(Fbit);
Fcutoff = Fbit;
time = (0:Ts:Ts*(length(Iin) - 1))';

% -- EDIT THIS SECTION --
%choose to run for ~nbits duration or full capture length
nbits = 1e4;
sim_time = nbits/Fbit; % ARW 9/8/17: increment is samples for slx run, not bits
%sim_time = time(end);
% -----------------------

%% downconvert from IF
cc = cos(2*pi*IF*time);
ss = sin(2*pi*IF*time);
I = Iin.*cc - Qin.*ss;
Q = Iin.*ss + Qin.*cc;

%create timeseries
I = timeseries(I,time);
Q = timeseries(Q,time);

%% -- UNCOMMENT this to view data on spec. analyzer
h1 = dsp.SpectrumAnalyzer(1);
    h1.SampleRate=Fs;
    h1.SpectralAverages=50;
    h1.ReferenceLoad = 50;
    h1.FrequencySpan = 'Span and center frequency';
    h1.Span = 100e6;
    h1.CenterFrequency = 0;
    h1.YLimits = [-90 -60];
    h1.ChannelNames = {sprintf('IQ DDC at %1.0f MHz IF',1e-6*IF)};
    h1.ShowLegend = true;
    h1.Title = sprintf('IQ DDC spectrum at 50 Ohm for Trial %d, Run %d, LO %1.0f MHz',tr,rn,s.LO);

%%
iq = I.data+j.*Q.data;
step(h1, iq)

%% filter
uI = mean(I);
uQ = mean(Q);
interval = [0 Fcutoff];
I = idealfilter(I-uI,interval,'pass');
Q = idealfilter(Q-uQ,interval,'pass');
iqFilt = I.data + j.*Q.data;

%%
h2 = dsp.SpectrumAnalyzer(1);
    h2.SampleRate=Fs;
    h2.SpectralAverages=50;
    h2.ReferenceLoad = 50;
    h2.FrequencySpan = 'Span and center frequency';
    h2.Span = 30e6;
    h2.CenterFrequency = 0;
    h2.YLimits = [-90 -60];
    h2.ChannelNames = {sprintf('IQ DDC at %1.0f MHz IF',1e-6*IF)};
    h2.ShowLegend = true;
    h2.Title = sprintf('Filtered IQ DDC spectrum at 50 Ohm for Trial %d, Run %d, LO %1.0f MHz',tr,rn,s.LO);

%%
step(h2,iqFilt)

%% Run the simulation
sim('BenchtopSkinpatch_SingleTone');

%% process output bits
bits = outBits.data(bitStrobe.data > 0);
bits = bits(:);

%31 bit PN sequence
bitseq = [1 1 0 0 0 1 1 1 1 1 0 0 1 1 0 1 0 0 1 0 0 0 0 1 0 1 0 1 1 1 0];

%compute error for all phases of bit sequence
for i = 1:31
    seq = repmat([bitseq(i:end) bitseq(1:i-1)]',floor(length(bits)/31),1); 
    e(i) = sum(bits(1:length(seq)) == seq); 
end

%account for incorrect polarity
[vmin, imin] = min(e);
[vmax, imax] = max(e);

%find best match
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