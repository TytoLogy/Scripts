%------------------------------------------------------------------------
% OutputInput.m
%------------------------------------------------------------------------
% 
% Script that demonstrates input and output on RosenLab Rig
% 
%------------------------------------------------------------------------

%------------------------------------------------------------------------
%  Sharad Shanbhag
%	sshanbhag@neomed.edu
%------------------------------------------------------------------------
% Created: 7 January, 2013
%
% Revisions:
%------------------------------------------------------------------------
% To Do:
%------------------------------------------------------------------------

%---------------------------------------------------------------
%---------------------------------------------------------------
%% Constants
%---------------------------------------------------------------
%---------------------------------------------------------------
% max possible PA5 attenuation
MAX_ATTEN = 120;
% medusa channel number for spike data
SPIKECHAN = 1;

%---------------------------------------------------------------
%---------------------------------------------------------------
%% Setup paths
%---------------------------------------------------------------
%---------------------------------------------------------------
disp([mfilename ': checking paths'])
% directory when using installed version:
pdir = ['C:\TytoLogy\TytoLogySettings\' getenv('USERNAME')];

if isempty(which('RPload'))
	% could not find the RPload.m function (which is in TytoLogy
	% toolbox) which suggests that the paths are not set or are 
	% incorrect for this setup.  load the paths using the tytopaths program.
	disp([mfilename ': loading paths using ' pdir '\tytopaths.m']);
	run(fullfile(pdir, 'tytopaths'));
	% now recheck
	if isempty(which('RPload'))
		error('%s: tried setting paths via %s, but failed.  sorry.', ...
					mfilename, fullfile(pdir, 'tytopaths.m'));
	end
else
	% seems okay, so continue
	disp([mfilename ': paths ok, launching programn'])
end

%---------------------------------------------------------------
%---------------------------------------------------------------
%% load calibration
%---------------------------------------------------------------
%---------------------------------------------------------------
% get a fake cal structure
caldata = fake_caldata;

%---------------------------------------------------------------
%---------------------------------------------------------------
%% Stimulus settings
%---------------------------------------------------------------
%---------------------------------------------------------------

%------------------------------------------------------
% noise settings
%------------------------------------------------------
% noise low freq
noise.Flo = 200;
% noise high freq
noise.Fhi = 10000;
% noise ITD (interaural time difference, usec)
noise.ITD = 0;
% noise binaural correlation
noise.BC = 100;
% on/off ramp time in msec
noise.Ramp = 5;
% stim intensity (db SPL)
noise.SPL = 50;

%------------------------------------------------------
% tone settings
%------------------------------------------------------
% tone freq
tone.F = 440;
% tone ITD (interaural time difference, usec)
tone.ITD = 0;
% vary onset phase? (0 = no, 1 = yes)
tone.RadVary = 0;
% on/off ramp time
tone.Ramp = 5;
% stim intensity (db SPL)
tone.SPL = 50;

%------------------------------------------------------
% SAM settings
%------------------------------------------------------
% SAM carrier freq range
sam.NoiseF = [200 10000];
% SAM ITD
sam.ITD = 0;
% SAM binaural corr
sam.BC = 100;
% Modulation depth (percent)
sam.sAMPercent = 100;
% SAM modulation freq (Hz)
sam.sAMFreq = 20;
% on/off ramp time
sam.Ramp = 5;
% stim intensity (db SPL)
sam.SPL = 50;

% Left output only
LRenable = [1 0];

%---------------------------------------------------------------
%---------------------------------------------------------------
%% Hardware settings
%---------------------------------------------------------------
%---------------------------------------------------------------

%------------------------------------------------------
% TDT - these are general processing settings
%------------------------------------------------------
% ISI
tdt.StimInterval = 500;
tdt.StimDuration = 250;
tdt.AcqDuration = 300;
% sweep period is usually AcqDuration plus a small factor (10 msec)
tdt.SweepPeriod = tdt.AcqDuration + 10;
tdt.StimDelay = 20;
tdt.HeadstageGain = 1000;			% gain for headstage
tdt.MonitorChannel = 1;				% monitor channel on Rz5 (from medusa)
tdt.MonitorGain = 1000;				% monitor channel gain (for display and monitor channel only)
tdt.decifactor = 1;					% factor to reduce input data sample rate
tdt.HPEnable = 1;						% enable HP filter
tdt.HPFreq = 200;						% HP frequency
tdt.LPEnable = 1;						% enable LP filter
tdt.LPFreq = 10000;					% LP frequency
% # of input (spike) channels
tdt.nChannels = 1;
tdt.InputChannel = zeros(tdt.nChannels, 1);
tdt.OutputChannel = [1 2];
%TTL pulse duration (msec)
tdt.TTLPulseDur = 1;

%------------------------------------------------------
% indev - these are settings for input device (RZ5)
%------------------------------------------------------
% input (spike) device is Medusa on RZ5, 1 channel input
indev.Fs = 25000;
% set this to wherever the circuits are stored
indev.Circuit_Path = 'C:\TytoLogy\Toolbox\TDTToolbox\Circuits\RZ5';
% for recording from 16 Channels
indev.Circuit_Name = 'RZ5_1ChannelAcquire_zBus';
% Dnum = device number - this is for RZ5
indev.Dnum=1;
% place holder for ActiveX control
indev.C = [];
% set status to 0
indev.status = 0;

%------------------------------------------------------
% outdev - these are settings for output device (RZ6)
%------------------------------------------------------
outdev.Fs = 50000;
% set this to wherever the RZ6 circuits are stored
outdev.Circuit_Path = 'C:\TytoLogy\Toolbox\TDTToolbox\Circuits\RZ6\';
outdev.Circuit_Name = 'RZ6_SpeakerOutput_zBus';
% Dnum = device number - this is for RZ6, device 1
outdev.Dnum=1;
outdev.C = [];
outdev.status = 0;

%------------------------------------------------------
% initialize some variables for TDT control objects
%------------------------------------------------------
% zBus is used for triggering
zBUS = [];
% attenuator controls
PA5L = [];
PA5R = [];

%---------------------------------------------------------------
%---------------------------------------------------------------
%% Initialize Hardware
%---------------------------------------------------------------
%---------------------------------------------------------------
disp('...starting TDT hardware...');

%------------------------------------------------------
% Initialize zBus control
%------------------------------------------------------
try
	disp('...starting zBUS...')
	tmpdev = zBUSinit('GB');
	zBUS.C = tmpdev.C;
	zBUS.handle = tmpdev.handle;
	zBUS.status = tmpdev.status;
catch err
	disp([mfilename ': error initializing zBUS'])
	disp(err.message);
	disp(err.identifier);
	disp(err.stack);
end
%------------------------------------------------------
% Initialize RZ5/Medusa
%------------------------------------------------------
try
	disp('...starting Medusa attached to RZ5...')
	tmpdev = RZ5init('GB');
	indev.C = tmpdev.C;
	indev.handle = tmpdev.handle;
	indev.status = tmpdev.status;
catch err
	disp([mfilename ': error initializing zBUS'])
	disp(err.message);
	disp(err.identifier);
	disp(err.stack);
end
%------------------------------------------------------
% Initialize RZ6
%------------------------------------------------------
try
	disp('...starting RZ6 for loudspeaker output...')
	tmpdev = RZ6init('GB', outdev.Dnum);
	outdev.C = tmpdev.C;
	outdev.handle = tmpdev.handle;
	outdev.status = tmpdev.status;
catch err
	disp([mfilename ': error initializing zBUS'])
	disp(err.message);
	disp(err.identifier);
	disp(err.stack);
end
%------------------------------------------------------
% Initialize Attenuators
%------------------------------------------------------
try
	PA5L = PA5init('GB', 1, MAX_ATTEN);
	PA5R = PA5init('GB', 2, MAX_ATTEN);
catch err
	disp([mfilename ': error initializing zBUS'])
	disp(err.message);
	disp(err.identifier);
	disp(err.stack);
end

%------------------------------------------------------
% Loads circuits
%------------------------------------------------------
try
	indev.status = RPload(indev);
catch err
	disp([mfilename ': error loading indev'])
	disp(err.message);
	disp(err.identifier);
	disp(err.stack);
end
try 
	outdev.status = RPload(outdev);
catch err
	disp([mfilename ': error initializing zBUS'])
	disp(err.message);
	disp(err.identifier);
	disp(err.stack);
end

%------------------------------------------------------
% Starts Circuits
%------------------------------------------------------
inStatus = RPrun(indev);
outStatus = RPrun(outdev);

%------------------------------------------------------
% Send zBus A and B triggers to initialize and check enable
% status
%------------------------------------------------------
zBUStrigA_PULSE(zBUS);
zBUStrigB_PULSE(zBUS);
tmp = RPgettagval(indev, 'Enable');
fprintf('indev enable: %d\n', tmp);
tmp = RPgettagval(outdev, 'Enable');
fprintf('outdev enable: %d\n', tmp);

%------------------------------------------------------
% get the input and output sampling rates
%------------------------------------------------------
outdev.Fs = RPsamplefreq(outdev);
indev.Fs = RPsamplefreq(indev);
% get the tags and values for the circuits
tmptags = RPtagnames(outdev);
outdev.TagName = tmptags;				
tmptags = RPtagnames(indev);
indev.TagName = tmptags;


%------------------------------------------------------
% set up circuits
%------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input Settings
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set the length of time to acquire data (in samples)
RPsettag(indev, 'AcqDur', ms2samples(tdt.AcqDuration, indev.Fs));
% Set the total sweep period time (samples)
RPsettag(indev, 'SwPeriod', ms2samples(tdt.SweepPeriod, indev.Fs));
% set the HP filter
if tdt.HPEnable == 1
	RPsettag(indev, 'HPEnable', 1);
	RPsettag(indev, 'HPFreq', tdt.HPFreq);
else
	RPsettag(indev, 'HPEnable', 0);
end
% set the LP filter
if tdt.LPEnable == 1
	RPsettag(indev, 'LPEnable', 1);
	RPsettag(indev, 'LPFreq', tdt.LPFreq);
else
	RPsettag(indev, 'LPEnable', 0);
end
% set the HeadstageGain
status = RPsettag(indev, 'mcGain', tdt.HeadstageGain);
% get the buffer index
index_in = RPgettag(indev, 'mcIndex');
% Set the sweep count to 1
RPsettag(indev, 'SwCount', 1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Output Settings
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set up some of the buffer/stimulus parameters
% Set the total sweep period time (in samples)
RPsettag(outdev, 'SwPeriod', ms2samples(tdt.SweepPeriod, outdev.Fs));
% Set the sweep count to 1
RPsettag(outdev, 'SwCount', 1);
% Set the Stimulus Delay (in samples)
RPsettag(outdev, 'StimDelay', ms2samples(tdt.StimDelay, outdev.Fs));
% Set the Stimulus Duration (in samples)
RPsettag(outdev, 'StimDur', ms2samples(tdt.StimDuration, outdev.Fs));
% Set the monitor D/A channel on RX5 and monitor gain
RPsettag(indev, 'MonChan', tdt.MonitorChannel);
RPsettag(indev, 'MonGain', tdt.MonitorGain);	
% Turn on the monitor channel
status = RPtrig(indev, 1);

%------------------------------------------------------
% scope/audio monitor for spike data
%------------------------------------------------------
RPsettag(indev, 'MonChan', SPIKECHAN);
RPsettag(indev, 'MonitorEnable', 1);
RPsettag(indev, 'MonGain', tdt.MonitorGain);

%------------------------------------------------------
% get the number of points to acquire - needed to know sampling rate to
% do this, which is why it's done here
%------------------------------------------------------
AcquirePoints = tdt.nChannels * ms2samples(tdt.AcqDuration, indev.Fs);

%---------------------------------------------------------------
%---------------------------------------------------------------
%% synthesize stimuli
%---------------------------------------------------------------
%---------------------------------------------------------------
%------------------------------------------------------
% synthesize noise stimulus
%------------------------------------------------------
[noise.S, noise.rms_val] = syn_headphone_noise(tdt.StimDuration, outdev.Fs, ...
												noise.Flo, noise.Fhi, ...
												noise.ITD, noise.BC, caldata);
% ramp stimulus on and off
noise.S = sin2array(noise.S, noise.Ramp, outdev.Fs);
%------------------------------------------------------
% synthesize tone stim
%------------------------------------------------------
[tone.S, tone.rms_val] = syn_headphone_tone(tdt.StimDuration, outdev.Fs, ...
											tone.F, tone.ITD, tone.RadVary, caldata);
% ramp stimulus on and off
tone.S = sin2array(tone.S, tone.Ramp, outdev.Fs);
%------------------------------------------------------
% synthesize S.A.M. stimulus
%------------------------------------------------------
sam.NoiseF = [noise.Flo noise.Fhi];
[sam.S, sam.rms_val, sam.rms_mod, sam.modPhi]  = ...
					syn_headphone_amnoise(tdt.StimDuration, outdev.Fs, sam.NoiseF, ...
													sam.ITD, sam.BC, ...
													sam.sAMPercent, sam.sAMFreq, [], ...
													caldata);
% ramp stimulus on and off
sam.S = sin2array(sam.S, sam.Ramp, outdev.Fs);

%---------------------------------------------------------------
%---------------------------------------------------------------
%% pause
%---------------------------------------------------------------
%---------------------------------------------------------------
fprintf('Ready to go!\npress key to continue...\n')
pause


%---------------------------------------------------------------
%---------------------------------------------------------------
%% Play/record noise
%---------------------------------------------------------------
%---------------------------------------------------------------
%------------------------------------------------------
% compute attenuator settings
%------------------------------------------------------
spl_val = (noise.SPL * [1 1]);
[atten, spl_val] = figure_headphone_atten(spl_val, noise.rms_val, ...
													caldata, LRenable);
%------------------------------------------------------
% set the attenuators
%------------------------------------------------------
PA5setatten(PA5L, atten(1));
PA5setatten(PA5R, atten(2));
%------------------------------------------------------
% play the sound and return the response
%------------------------------------------------------
[resp, rate] = speakerstim_medusarec_1chan(noise.S, AcquirePoints, indev, outdev, zBUS);
% de-multiplex the data if more than 1 channel was collected
% mcDeMux returns an array that is [nChannels, nPoints]
if tdt.nChannels > 1
	resp = mcDeMux(resp, tdt.nChannels);
else
	resp = resp';
end

%------------------------------------------------------
% Plotting/display setup
%------------------------------------------------------
% get the # of points to send out and to collect
acqpts = ms2samples(tdt.AcqDuration, indev.Fs);
outpts = ms2samples(tdt.StimDuration, outdev.Fs);
% stimulus start and end bins used for analyzing input data
stim_start = ms2samples(tdt.StimDelay, indev.Fs);
stim_end = stim_start + ms2samples(tdt.StimDuration, indev.Fs);
% time vector for plots
dt = 1/indev.Fs;
tvec = 1000*dt*(0:(acqpts-1));
figure(10)
subplot(211)
plot(noise.S(1, :));
subplot(212)
plot(tvec, resp);

%---------------------------------------------------------------
%---------------------------------------------------------------
%% pause
%---------------------------------------------------------------
%---------------------------------------------------------------
disp('Noise output complete.  press key to continue...')
pause

%---------------------------------------------------------------
%---------------------------------------------------------------
%% Play/record tone
%---------------------------------------------------------------
%---------------------------------------------------------------
%------------------------------------------------------
% compute attenuator settings
%------------------------------------------------------
spl_val = (tone.SPL * [1 1]);
[atten, spl_val] = figure_headphone_atten(spl_val, tone.rms_val, ...
															caldata, LRenable);
%------------------------------------------------------
% set the attenuators
%------------------------------------------------------
PA5setatten(PA5L, atten(1));
PA5setatten(PA5R, atten(2));
%------------------------------------------------------
% play the sound and return the response
%------------------------------------------------------
[resp, rate] = speakerstim_medusarec_1chan(tone.S, AcquirePoints, indev, outdev, zBUS);
% de-multiplex the data if more than 1 channel was collected
% mcDeMux returns an array that is [nChannels, nPoints]
if tdt.nChannels > 1
	resp = mcDeMux(resp, tdt.nChannels);
else
	resp = resp';
end
% plot
figure(11)
subplot(211)
plot(tone.S(1, :));
subplot(212)
plot(resp);

%---------------------------------------------------------------
%---------------------------------------------------------------
%% pause
%---------------------------------------------------------------
%---------------------------------------------------------------
disp('Tone output complete.  press key to continue...')
pause


%---------------------------------------------------------------
%---------------------------------------------------------------
%% Play/record SAM
%---------------------------------------------------------------
%---------------------------------------------------------------
%------------------------------------------------------
% compute attenuator settings
%------------------------------------------------------
spl_val = (sam.SPL * [1 1]);
[atten, spl_val] = figure_headphone_atten(spl_val, sam.rms_mod, ...
															caldata, LRenable);
%------------------------------------------------------
% set the attenuators
%------------------------------------------------------
PA5setatten(PA5L, atten(1));
PA5setatten(PA5R, atten(2));
%------------------------------------------------------
% play the sound and return the response
%------------------------------------------------------
[resp, rate] = speakerstim_medusarec_1chan(sam.S, AcquirePoints, indev, outdev, zBUS);
% de-multiplex the data if more than 1 channel was collected
% mcDeMux returns an array that is [nChannels, nPoints]
if tdt.nChannels > 1
	resp = mcDeMux(resp, tdt.nChannels);
else
	resp = resp';
end
% plot
figure(12)
subplot(211)
plot(sam.S(1, :));
subplot(212)
plot(resp);

%---------------------------------------------------------------
%---------------------------------------------------------------
%% pause
%---------------------------------------------------------------
%---------------------------------------------------------------
disp('sAM output complete.  press key to continue...')
pause

%---------------------------------------------------------------
%---------------------------------------------------------------
%% close TDT objects, etc
%---------------------------------------------------------------
%---------------------------------------------------------------

% turn off the monitor via software trigger 2
RPtrig(indev, 2);
disp('...closing PA5L')
PA5L.status = PA5close(PA5L);
disp('...closing PA5R')
PA5R.status = PA5close(PA5R);
disp('...closing indev')
indev.status = RPclose(indev);
disp('...closing outdev')
outdev.status = RPclose(outdev);
disp('...closing zBUS')
zBUS.status = zBUSclose(zBUS);


