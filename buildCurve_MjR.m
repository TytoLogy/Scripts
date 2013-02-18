function buildCurve_MjR
gui=get(gcf,'userdata');

% initialize tdt here after user has time to choose stimulus and
% acquisition duration settings that are needed for proper initialization
init_tdt_MjR;   % initialize TDT hardware

%------------------------------------------------------
%------------------------------------------------------
% load calibration

% get a fake cal structure
caldata = fake_caldata;


%------------------------------------------------------
%------------------------------------------------------
%------------------------------------------------------
% Hardware settings

%------------------------------------------------------
% TDT - these are general processing settings
%------------------------------------------------------
tdt = gui.handles.tdt;
% ISI - these are specified in GUI
tdt.StimInterval = (str2num(get(gui.StimInt,'string')));       % stimulus interval (ms)
tdt.StimDuration = (str2num(get(gui.StimDur,'string')));       % stimulus duration (ms)
tdt.AcqDuration = (str2num(get(gui.AcqDur,'string')));         % acquisition duration (ms)
tdt.StimDelay = (str2num(get(gui.StimDelay,'string')));        % stimulus delay (ms)
% sweep period is usually AcqDuration plus a small factor (10 msec)
tdt.SweepPeriod = tdt.AcqDuration + 10;

% % sweep period is usually AcqDuration plus a small factor (10 msec)
% tdt.SweepPeriod = tdt.AcqDuration + 10;
% tdt.StimDelay = 20;
% tdt.HeadstageGain = 1000;	% gain for headstage
% tdt.MonitorChannel = 1;	% monitor channel on Rz5 (from medusa)
% tdt.MonitorGain = 1000;	% monitor channel gain (for display and monitor channel only)
% tdt.decifactor = 1;	% factor to reduce input data sample rate
% tdt.HPEnable = 1;	% enable HP filter
% tdt.HPFreq = 200;	% HP frequency
% tdt.LPEnable = 1;	% enable LP filter
% tdt.LPFreq = 10000;	% LP frequency
% % # of input (spike) channels
% tdt.nChannels = 1;
% tdt.InputChannel = zeros(tdt.nChannels, 1);
% tdt.OutputChannel = [1 2];
% %TTL pulse duration (msec)
% tdt.TTLPulseDur = 1;

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
%------------------------------------------------------
% Stimulus settings

%------------------------------------------------------
% tAM settings
%------------------------------------------------------
tAM.NoiseF = (str2num(get(gui.NoiseF,'string')));       % tAM carrier freq range
tAM.RiseFall = (str2num(get(gui.RiseFall,'string')));   % rise/fall time for each period (ms)
tAM.tAMDepth = (str2num(get(gui.tAMDepth,'string')));   % Modulation depth (percent)
tAM.tAMFreq = (str2num(get(gui.tAMFreq,'string')));     % tAM modulation freq (Hz)
tAM.SPL = (str2num(get(gui.SPL,'string')));             % stim intensity (db SPL)

% tAM ITD
tAM.ITD = 0;
% tAM ILD
tAM.ILD = 0;
% tAM binaural corr
tAM.BC = 100;
% on/off ramp time for entire signal (ms)
% NB: whether to do this is an issue that needs resolving
tAM.Ramp = 1;


%------------------------------------------------------
%------------------------------------------------------
% set/get parameters
c.nreps =       (str2num(get(gui.nReps,'string')));         % number of repetitions
c.freezeStim =  (str2num(get(gui.freezeStim,'string')));    % 1 or 0 indicating FROZEN NOISE
% get string for type of stimulus - noise or tone
if length(tAM.NoiseF) ~= 1
    c.stimtype = 'noise';
else
    c.stimtype = 'tone';
end

allcurves = get(gui.curvetype,'string');    % syntax can also be read_ui_str(gui.curvetype)
curveidx = get(gui.curvetype,'value');      % syntax can also be read_ui_val(gui.curvetype)
curvetype = allcurves(curveidx,:);

% REMEMBER: allcurves will have dimensions matching the longest string, so
% need to pad shorter strings with spaces

% ntrials == # of stim values (risetimes, moddepths, etc.)
% set number of trials, determined by range of varied parameter
if strcmp(curvetype,'tAM_RISEFALL')
    if length(tAM.RiseFall) > 1
        c.ntrials = length(tAM.RiseFall);
        if length(tAM.tAMDepth) > 1 | length(tAM.tAMFreq) > 1
            warndlg('Make sure other fields have one parameter each.')
        end
    else
        warndlg('Hey you!  Put a range of values into RiseFall!   e.g: [2 5 10 20 50]')
    end
elseif strcmp(curvetype,'tAM_MODDEPTH')
    if length(tAM.tAMDepth) > 1
        c.ntrials = length(tAM.tAMDepth);
        if length(tAM.RiseFall) > 1 | length(tAM.tAMFreq) > 1
            warndlg('Make sure other fields have one parameter each.')
        end
    else
        warndlg('Hey you!  Put a range of values into ModDepth!   e.g: [100 75 50 40 30 20 10 5 0]')
    end
elseif strcmp(curvetype,'tAM_MODFREQ ')
    if length(tAM.tAMFreq) > 1
        c.ntrials = length(tAM.tAMFreq);
        if length(tAM.RiseFall) > 1 | length(tAM.tAMDepth) > 1
            warndlg('Make sure other fields have one parameter each.')
        end
    else
        warndlg('Hey you!  Put a range of values into ModFreq!   e.g: [2 5 10 20 50 100 200]')
    end
end

% allocate some arrays for storage
c.nstims = c.nreps * c.ntrials;
c.repnum = zeros(c.nstims, 1);
c.trialnum = zeros(c.nstims, 1);
sindex = 0;
for rep = 1:c.nreps
	for trial = 1:c.ntrials
		sindex = sindex + 1;
		c.repnum(sindex) = rep;
		c.trialnum(sindex) = trial;
	end
end


%------------------------------------------------------
%------------------------------------------------------
% build curves

disp([mfilename ' is building stimuli for ' curvetype ' curve...'])
switch curvetype
    
    
%% tAM risefall modulation Curve

    case 'tAM_RISEFALL'
        varpar = [];
        gui=get(gcbf,'userdata');   % references data previously set so it can be passed into this function

        % Stimulus parameter to vary (varName) and the range (stimvar)
		c.vname = upper(curvetype);
		c.vrange = tAM.RiseFall;
		
		% for tAM_RISEFALL curves, these parameters are fixed:
		ITD = tAM.ITD;
        ILD = tAM.ILD;
		BC = tAM.BC;
		SPL = tAM.SPL;  % this was ABI - see if you can substitute SPL below
		tAMFreq = tAM.tAMFreq;
        tAMDepth = tAM.tAMDepth;
        CarrierFREQ = tAM.NoiseF;
        if tAM.NoiseF(end) >= caldata.freq(end) | tAM.NoiseF(1) <= caldata.freq(1)
            error('CAN''T DO THIS! Your carrier frequency is outside the range of frequencies in the calibration file.')
        end
        
		% If noise is frozen, generate zero ITD spectrum or tone, to be replicated later
		if c.freezeStim & length(CarrierFREQ)>1
			% get ITD = 0 Smag and Sphase
            [c.S0, c.rms_val0, c.rms_mod0, c.modPhi0, c.Smag0, c.Sphase0] = ...
                syn_headphone_tAM(5, tdt.StimDuration, outdev.Fs, CarrierFREQ, ...
                                  tAM.ITD, tAM.BC, ...
                                  tAMDepth, tAMFreq, [], ...
                                  caldata);
        end
        
        % Randomize trial presentations
        stimseq = HPCurve_randomSequence(c.nreps, c.ntrials);
        c.trialRandomSequence = stimseq;

		sindex = 0;
        plotsignals = figure(2);
		% now loop through the randomized trials
		for rep = 1:c.nreps
			for trial = 1:c.ntrials
				sindex = sindex + 1;

				% Get the randomized stimulus variable value from c.stimvar 
				% indices stored in c.trialRandomSequence
				tAMRiseFall = c.vrange(c.trialRandomSequence(rep, trial));

				% spl_val sets the L and R channel db levels, and the ILD
				spl_val = computeLRspl(ILD, SPL);  % substituted SPL for ABI

				% Synthesize noise or tone, frozed or unfrozed and get rms values for setting attenuator
                if ~c.freezeStim | length(CarrierFREQ)==1 % stimulus is unfrozen or tone
                    [Sn, rms_val, rms_mod, modPhi] = ...
                        syn_headphone_tAM(tAMRiseFall, tdt.StimDuration, outdev.Fs, CarrierFREQ, ...
                                          tAM.ITD, tAM.BC, ...
                                          tAMDepth, tAMFreq, [], ...
                                          caldata);

                else	% stimulus is frozen
                    [Sn, rms_val, rms_mod, modPhi] = ...
                        syn_headphone_tAM(tAMRiseFall, tdt.StimDuration, outdev.Fs, CarrierFREQ, ...
                                          tAM.ITD, tAM.BC, ...
                                          tAMDepth, tAMFreq, [], ...
                                          caldata, c.Smag0, c.Sphase0);
				end

                % ramp stimulus on and off
                Sn = sin2array(Sn, tAM.Ramp, outdev.Fs);


				% get the attenuator settings for the desired SPL
                atten = figure_headphone_atten(spl_val, rms_mod, caldata);

				% Store the parameters in the stimulus cache struct
				c.stimvar{sindex} = tAMRiseFall;
				c.Sn{sindex} = Sn;
				c.splval{sindex} = spl_val;
				c.rmsval{sindex} = rms_mod;
				c.atten{sindex} = atten;
				c.ITD(sindex) = ITD;
				c.ILD(sindex) = ILD;
				c.BC(sindex) = BC;
                c.tAMRiseFall = tAMRiseFall;
				c.CarrierFREQ{sindex} = CarrierFREQ;
				c.tAMDepth(sindex) = tAMDepth;
				c.tAMFreq(sindex) = tAMFreq;
                
                varpar=strvcat(varpar, ['RiseFall (ms): ', num2str(tAMRiseFall)]);
                
                % plot each different type of signal
                if rep == 1
                    figure(plotsignals);  subplot(c.ntrials,1,trial);
                    plot(c.Sn{sindex}(1,:),'r-');  title(['RiseFall (ms) = ' num2str(tAMRiseFall)])
                    if trial ~= c.ntrials
                        set(gca,'xticklabel','')
                    end
                end
            end	%%% End of TRIAL LOOP
     
		end %%% End of REPS LOOP
        
        % 1) passes struct c to userdata within gui.stimulus
        % 2) updates listbox with list of parameter that was varied for each stimulus
        set(gui.stimulus, 'string', varpar, 'userdata', c, 'value', 1);

        c.FREQ = c.CarrierFREQ;
        c.ABI = c.splval;
        c.curvetype = curvetype;

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% tAM ModDepth Curve  

    case 'tAM_MODDEPTH'
        varpar = [];
        gui=get(gcbf,'userdata');   % references data previously set so it can be passed into this function

        % Stimulus parameter to vary (varName) and the range (stimvar)
		c.vname = upper(curvetype);
		c.vrange = tAM.tAMDepth;
		
		% for tAM_MODDEPTH curves, these parameters are fixed:
		ITD = tAM.ITD;
        ILD = tAM.ILD;
		BC = tAM.BC;
		SPL = tAM.SPL;  % this was ABI - see if you can substitute SPL below
		tAMFreq = tAM.tAMFreq;
        tAMRiseFall = tAM.RiseFall;
        CarrierFREQ = tAM.NoiseF;
        if tAM.NoiseF(end) >= caldata.freq(end) | tAM.NoiseF(1) <= caldata.freq(1)
            error('CAN''T DO THIS! Your carrier frequency is outside the range of frequencies in the calibration file.')
        end
        
		% If noise is frozen, generate zero ITD spectrum or tone, to be replicated later
		if c.freezeStim & length(CarrierFREQ)>1
			% get ITD = 0 Smag and Sphase
            [c.S0, c.rms_val0, c.rms_mod0, c.modPhi0, c.Smag0, c.Sphase0] = ...
                syn_headphone_tAM(tAMRiseFall, tdt.StimDuration, outdev.Fs, CarrierFREQ, ...
                                  tAM.ITD, tAM.BC, ...
                                  100, tAMFreq, [], ...
                                  caldata);
        end
        
        % Randomize trial presentations
        stimseq = HPCurve_randomSequence(c.nreps, c.nTrials);
        c.trialRandomSequence = stimseq;

		sindex = 0;
        plotsignals = figure(2);
		% now loop through the randomized trials
		for rep = 1:c.nreps
			for trial = 1:c.nTrials
				sindex = sindex + 1;

				% Get the randomized stimulus variable value from c.stimvar 
				% indices stored in c.trialRandomSequence
				tAMDepth = c.vrange(c.trialRandomSequence(rep, trial));

				% spl_val sets the L and R channel db levels, and the ILD
				spl_val = computeLRspl(ILD, SPL);  % substituted SPL for ABI

				% Synthesize noise or tone, frozed or unfrozed and get rms values for setting attenuator
                if ~c.freezeStim | length(CarrierFREQ)==1 % stimulus is unfrozen or tone
                    [Sn, rms_val, rms_mod, modPhi] = ...
                        syn_headphone_tAM(tAMRiseFall, tdt.StimDuration, outdev.Fs, CarrierFREQ, ...
                                          tAM.ITD, tAM.BC, ...
                                          tAMDepth, tAMFreq, [], ...
                                          caldata);

                else	% stimulus is frozen
                    [Sn, rms_val, rms_mod, modPhi] = ...
                        syn_headphone_tAM(tAMRiseFall, tdt.StimDuration, outdev.Fs, CarrierFREQ, ...
                                          tAM.ITD, tAM.BC, ...
                                          tAMDepth, tAMFreq, [], ...
                                          caldata, c.Smag0, c.Sphase0);
				end

                % ramp stimulus on and off
                Sn = sin2array(Sn, tAM.Ramp, outdev.Fs);


				% get the attenuator settings for the desired SPL
                atten = figure_headphone_atten(spl_val, rms_mod, caldata);

				% Store the parameters in the stimulus cache struct
				c.stimvar{sindex} = tAMDepth;
				c.Sn{sindex} = Sn;
				c.splval{sindex} = spl_val;
				c.rmsval{sindex} = rms_mod;
				c.atten{sindex} = atten;
				c.ITD(sindex) = ITD;
				c.ILD(sindex) = ILD;
				c.BC(sindex) = BC;
                c.tAMRiseFall = tAMRiseFall;
				c.CarrierFREQ{sindex} = CarrierFREQ;
				c.tAMDepth(sindex) = tAMDepth;
				c.tAMFreq(sindex) = tAMFreq;
                
                varpar=strvcat(varpar, ['ModDepth (%): ', num2str(tAMDepth)]);
                
                % plot each different type of signal
                if rep == 1
                    figure(plotsignals);  subplot(c.nTrials,1,trial);
                    plot(c.Sn{sindex}(1,:),'r-');  title(['ModDepth (%) = ' num2str(tAMDepth)])
                    if trial ~= c.nTrials
                        set(gca,'xticklabel','')
                    end
                end
            end	%%% End of TRIAL LOOP
     
		end %%% End of REPS LOOP
        
        % 1) passes struct c to userdata within gui.stimulus
        % 2) updates listbox with list of parameter that was varied for each stimulus
        set(gui.stimulus, 'string', varpar, 'userdata', c, 'value', 1);

        c.FREQ = c.CarrierFREQ;
        c.ABI = c.splval;
        c.curvetype = curvetype;
        
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% tAM ModFreq Curve  

    case 'tAM_MODFREQ '
        varpar = [];
        gui=get(gcbf,'userdata');   % references data previously set so it can be passed into this function

        % Stimulus parameter to vary (varName) and the range (stimvar)
		c.vname = upper(curvetype);
		c.vrange = tAM.tAMFreq;
		
		% for tAM_RISEFALL curves, these parameters are fixed:
		ITD = tAM.ITD;
        ILD = tAM.ILD;
		BC = tAM.BC;
		SPL = tAM.SPL;  % this was ABI - see if you can substitute SPL below
        tAMRiseFall = tAM.RiseFall;
        tAMDepth = tAM.tAMDepth;
        CarrierFREQ = tAM.NoiseF;
        if tAM.NoiseF(end) >= caldata.freq(end) | tAM.NoiseF(1) <= caldata.freq(1)
            error('CAN''T DO THIS! Your carrier frequency is outside the range of frequencies in the calibration file.')
        end
        
		% If noise is frozen, generate zero ITD spectrum or tone, to be replicated later
		if c.freezeStim & length(CarrierFREQ)>1
			% get ITD = 0 Smag and Sphase
            [c.S0, c.rms_val0, c.rms_mod0, c.modPhi0, c.Smag0, c.Sphase0] = ...
                syn_headphone_tAM(tAMRiseFall, tdt.StimDuration, outdev.Fs, CarrierFREQ, ...
                                  tAM.ITD, tAM.BC, ...
                                  tAMDepth, 5, [], ...
                                  caldata);
        end
        
        % Randomize trial presentations
        stimseq = HPCurve_randomSequence(c.nreps, c.nTrials);
        c.trialRandomSequence = stimseq;

		sindex = 0;
        plotsignals = figure(2);
		% now loop through the randomized trials
		for rep = 1:c.nreps
			for trial = 1:c.nTrials
				sindex = sindex + 1;

				% Get the randomized stimulus variable value from c.stimvar 
				% indices stored in c.trialRandomSequence
				tAMFreq = c.vrange(c.trialRandomSequence(rep, trial));

				% spl_val sets the L and R channel db levels, and the ILD
				spl_val = computeLRspl(ILD, SPL);  % substituted SPL for ABI

				% Synthesize noise or tone, frozed or unfrozed and get rms values for setting attenuator
                if ~c.freezeStim | length(CarrierFREQ)==1 % stimulus is unfrozen or tone
                    [Sn, rms_val, rms_mod, modPhi] = ...
                        syn_headphone_tAM(tAMRiseFall, tdt.StimDuration, outdev.Fs, CarrierFREQ, ...
                                          tAM.ITD, tAM.BC, ...
                                          tAMDepth, tAMFreq, [], ...
                                          caldata);

                else	% stimulus is frozen
                    [Sn, rms_val, rms_mod, modPhi] = ...
                        syn_headphone_tAM(tAMRiseFall, tdt.StimDuration, outdev.Fs, CarrierFREQ, ...
                                          tAM.ITD, tAM.BC, ...
                                          tAMDepth, tAMFreq, [], ...
                                          caldata, c.Smag0, c.Sphase0);
				end

                % ramp stimulus on and off
                Sn = sin2array(Sn, tAM.Ramp, outdev.Fs);


				% get the attenuator settings for the desired SPL
                atten = figure_headphone_atten(spl_val, rms_mod, caldata);

				% Store the parameters in the stimulus cache struct
				c.stimvar{sindex} = tAMFreq;
				c.Sn{sindex} = Sn;
				c.splval{sindex} = spl_val;
				c.rmsval{sindex} = rms_mod;
				c.atten{sindex} = atten;
				c.ITD(sindex) = ITD;
				c.ILD(sindex) = ILD;
				c.BC(sindex) = BC;
                c.tAMRiseFall = tAMRiseFall;
				c.CarrierFREQ{sindex} = CarrierFREQ;
				c.tAMDepth(sindex) = tAMDepth;
				c.tAMFreq(sindex) = tAMFreq;
                
                varpar=strvcat(varpar, ['tAM Modfreq (Hz): ', num2str(tAMFreq)]);
                
                % plot each different type of signal
                if rep == 1
                    figure(plotsignals);  subplot(c.nTrials,1,trial);
                    plot(c.Sn{sindex}(1,:),'r-');  title(['tAM ModFreq (Hz) = ' num2str(tAMFreq)])
                    if trial ~= c.nTrials
                        set(gca,'xticklabel','')
                    end
                end
            end	%%% End of TRIAL LOOP
     
		end %%% End of REPS LOOP
        
        % 1) passes struct c to userdata within gui.stimulus
        % 2) updates listbox with list of parameter that was varied for each stimulus
        set(gui.stimulus, 'string', varpar, 'userdata', c, 'value', 1);

        c.FREQ = c.CarrierFREQ;
        c.ABI = c.splval;
        c.curvetype = curvetype;
        
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

gui.handles.tdt = tdt;

set(gui.fig,'userdata',gui);

% end of function buildCurve
