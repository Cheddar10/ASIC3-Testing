% July 2017
% Brown University Larson Lab
% Chester Kilfoyle
% Read in benchtop data, filter, apply PLL and bit slicing

% Open simulink system
open_system('BenchtopSkinpatch_SingleTone');

% -- EDIT THIS SECTION --
% Read in data and produce two data arrays Iin, Qin
% These arrays should represent voltage levels directly, eg. not int16
tr = 'T9.1';
fp = strcat('Datasets\trial9\',tr,'\');
fn = 'IQ.dat';
load(strcat(fp,tr,'.parameters.mat'));
m = memmapfile(strcat(fp,fn), 'Format', 'int16');
Iin = m.Data(1:2:end); Iin = double(Iin)*2^-15*.2;
Qin = m.Data(2:2:end); Qin = double(Qin)*2^-15*.2;
% -----------------------

% -- EDIT THIS SECTION --
% Sampling Freq, Tx Freq, and IF estimate (requires ~.5-1MHz accuracy)
Fs = s.fs*1e6;
RF = s.LO*1e6;
IF = 31.6e6;
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
sim_time = nbits/Fbit;
%sim_time = time(end);
% -----------------------

% downconvert from IF
cc = cos(2*pi*IF*time);
ss = sin(2*pi*IF*time);
I = Iin.*cc - Qin.*ss;
Q = Iin.*ss + Qin.*cc;

%create timeseries
I = timeseries(I,time);
Q = timeseries(Q,time);

% -- UNCOMMENT this to view data on spec. analyzer
h1 = dsp.SpectrumAnalyzer;
h1.SampleRate=Fs;
h1.SpectralAverages=25;
h1.ReferenceLoad = 50;
iq = I.data+1i*Q.data;
step(h1,iq)

% filter
uI = mean(I);
uQ = mean(Q);
interval = [0 Fcutoff];
I = idealfilter(I-uI,interval,'pass');
Q = idealfilter(Q-uQ,interval,'pass');

% Run the simulation
sim('BenchtopSkinpatch_SingleTone');

%process output bits
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