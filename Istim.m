 function tsform=Istim(stim)
%   =======================================================================
%   Istim.m
%   Created by Alan Horsager
%   04-15-2005
%   Designed to create a biphasic stimulus input patttern
%  revised Ione Fine 3/2017
%   =======================================================================

dt = stim.tsample;               % step size, ms
Tx = ceil(stim.pulsedur./dt);    % duration of first phase
Ty = ceil(stim.pulsedur./dt);    % duration of interphase
Tz = ceil(stim.pulsedur./dt);    % duration of second phase
Tp=Tx+Ty+Tz;
np=(stim.dur/1000)*stim.freq; % number of pulses
Tdur=ceil((1000./stim.freq)./dt); % total duration of the actual pulses
Trem=Tdur-Tp; % remainder for each pulse

% creation of a single phase of the stimulus input
sform(1:Tx,1) = -1;
sform(end+1:Tx+Ty, 1) = 0;
sform(end+1:Tx+Ty+Tz, 1) = 1;

tsform=[];
for n=1:np
    tsform=cat(1, tsform, zeros(Trem, 1), sform);
end
tsform(end+1:floor(stim.dur/dt))=0;
tsform=tsform(1:floor(stim.dur/dt));

