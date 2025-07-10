%%
addpath('D:/Dropbox/projects/matlab2/pulseq150/matlab/');
addpath('D:\Dropbox\projects\matlab2\safe_pns_prediction\');


%%
clearvars
clc


%%
res = load('diff_safe.mat');

%%

% Set system limits
sys = mr.opts('MaxGrad',800,'GradUnit','mT/m',...
    'MaxSlew',800,'SlewUnit','T/m/s',...
    'rfRingdownTime', 20e-6, 'rfDeadtime', 100e-6,...
    'adcDeadTime',20e-6, 'B0', 2.89,  ... % this is Siemens' 3T
    'gradRasterTime', res.dt ...
);  

seq=mr.Sequence(sys);      % Create a new sequence object

seq.addBlock(mr.makeArbitraryGrad('x', res.out*sys.gamma, 'system', sys))

seq.plot('stacked', true)

%%
seq.calcPNS('F:/data/2025_0704_ascfiles/MP_GradSys_P034_X60.asc', true, true)