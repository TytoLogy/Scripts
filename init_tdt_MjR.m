function tdt = init_tdt_MjR    % NB: this is automatically run by a line at the end of the function makegui

gui = get(gcf,'userdata');   % references data previously set so it can be passed into this function
tdt = gui.handles.tdt;

%------------------------------------------------------
% Constants
%------------------------------------------------------
% max possible PA5 attenuation
MAX_ATTEN = 120;
% medusa channel number for spike data
SPIKECHAN = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Device settings
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize Hardware
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('...starting TDT hardware...');

% create dummy figure
hDUMMY = figure;

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

set(hDUMMY, 'Visible', 'off');

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

keyboard
% WEIRD HERE

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

% MjR added
gui.handles.PA5L = PA5L;
gui.handles.PA5R = PA5R;
gui.handles.zBUS = zBUS;
gui.handles.outdev = outdev;
gui.handles.indev = indev;

set(gui.fig,'userdata',gui);

% end function init_tdt_MjR

