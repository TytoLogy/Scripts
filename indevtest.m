%% outdev settings
outdev.Fs = [];
% set this to wherever the circuits are stored
outdev.Circuit_Path = 'C:\TytoLogy\Toolbox\TDTToolbox\Circuits\RZ6\';
outdev.Circuit_Name = 'RZ6_SpeakerOutput_zBus';
% Dnum = device number - this is for RX6, device 1
outdev.Dnum=1;
outdev.C = [];
outdev.status = 0;

%% other settings
interface = 'GB';
device_num = 1;

%% Initialize zBus control
disp('...starting zBUS...')
tmpdev = zBUSinit('GB');
zBUS.C = tmpdev.C;
zBUS.handle = tmpdev.handle;
zBUS.status = tmpdev.status;

%% initialize RZ6
disp('...starting RZ6 for lautspracher output...')
tmpdev = RZ6init('GB', outdev.Dnum);
outdev.C = tmpdev.C;
outdev.handle = tmpdev.handle;
outdev.status = tmpdev.status;

%% Loads circuits
outdev.status = RPload(outdev);

%% Starts Circuits
outStatus = RPrun(outdev);

%% Get circuit information
% get the output sampling rates
outdev.Fs = RPsamplefreq(outdev);
% get the tags and values for the circuits
tmptags = RPtagnames(outdev);
outdev.TagName = tmptags;

RPgettagval(outdev, 'Enable')
outdev.Fs = RPsamplefreq(outdev);
outdev.Fs

% zBUStrigA_PULSE(zBUS)
% outdev.Fs = RPsamplefreq(outdev);
% outdev.Fs
% pause(1)
% zBUStrigB_PULSE(zBUS)
% outdev.Fs = RPsamplefreq(outdev);
% outdev.Fs

keyboard

%% clean up
disp('...closing outdev')
outdev.status = RPclose(outdev);
disp('...closing zBUS')
zBUS.status = zBUSclose(zBUS);