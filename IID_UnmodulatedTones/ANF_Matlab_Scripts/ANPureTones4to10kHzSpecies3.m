% ANPureTones4to10kHzSpecies3.m
% Andrew Brughera 2021, Apr 14
% Modified from ANTrTonesRepsSp3.m
% Applying the Zilany et al. 2014 model with pure tones,
% ? responses to multiple presentations written to sudfolders?
% Modeling ANF responses to pure-tone bursts from 0 to 80 dB SPL
% To be applied for model-LSO IID responses

rng('shuffle'); % Seed random number generator rng, based on present time
rng_state = rng; % Save the rng Type, Seed, and State
% model fiber parameters
% (CF varied) CF in Hz;
cohc  = 1.0;   % normal ohc function
cihc  = 1.0;   % normal ihc function
%fiberType = 3; % spontaneous rate (in spikes/s) of the fiber BEFORE refractory effects; "1" = Low, 0 spk/s; "2" = Medium, 5 spk/s; "3" = High, 100 spk/s
implnt = 1;    % "0" for approximate or "1" for actual implementation of the power-law functions in the Synapse

% stimulus parameters
Fs_Hz = 100e3;  % sampling rate of 100,000 Hz;
Fs_kHzStr = num2str(round(Fs_Hz/1000));
Ts_s = 1/Fs_Hz; % sample period
% Zilany et al. 2009, 2014, adaptation coefficients designed for
% sample rate of 100 kHz, "must be 100 kHz".

% (Zilany ANF) PSTH parameters -> make these SpikeTime parameters instead
% We might make a psth later, ...
% but now we want SPIKE TIMES for each stimulus presentation, so nrep = 1:
nrep = 1; % # of stimulus repetitions for psth (e.g., 50) from Zilany model
%psthbinwidth_s = 1/(80 * F0_Hz) % 0.5e-3; % binwidth in seconds;
psthbinwidth_s = Ts_s; % binwidth in seconds at full temporal resolution;

% fiberTypeStrings = {'Low Spont'; 'Medium Spont'; 'High Spont'};
% fiberTypeStringsShort = {'Lo'; 'Md'; 'Hi'};
fiberTypeStrings = {'Medium Spont'};
fiberTypeStringsShort = {'Md'};
nFiberTypes = length(fiberTypeStringsShort);

CarrierFreqs_Hz = [4000 6000 10000];
% ModFreqs_Hz = [32 64 128 256 512 800];

% Bernstein & Trahiotis, 2002: Starting phases random for each component
% *** Probably reduce the number of starting phases
StartPhases_rad = 2*pi*rand(1,49); % carriers: 49 random phases
% On/Off ramps will be the only envelope for pure tones
% EnvStartPhases_deg = 0:-15:-720; % envelopes: 24+24+1 start-phase decrements
minimumZilanyCF_Hz = 56;
CFs_Hz = max(CarrierFreqs_Hz,minimumZilanyCF_Hz);
% Minimum CF allowed by Zilany et al. model is 91.2 Hz
% In c code by Zilany et al. 2009, 2014:
% /** Calculate the location on basilar membrane from CF */
% bmplace = 11.9 * log10(0.80 + cf / 456.0);

% Tone-burst stimuli (Sanes & Rubel, 1988)
% Duration (dur)
durSilence_s = 0.15;
durBurst_s = 0.05;
nBursts = 40; % with leading and trailing silence
durTotal_s = nBursts * (durSilence_s + durBurst_s) + durSilence_s;
t_s = 0:Ts_s:durTotal_s-Ts_s;
% rise/fall time in seconds
durRamp_s = 0.004; % ramps 4-ms (ours cosine-squared)
nPtsSilence = durSilence_s * Fs_Hz;
silence = zeros(1,nPtsSilence);
% Plateau at maximum amplitude for each burst
durPlateau_s = durBurst_s - (2 * durRamp_s);
nPtsPlateau = durPlateau_s * Fs_Hz;
nPtsPerRamp = durRamp_s * Fs_Hz;
Plateau = ones(1,nPtsPlateau);
n = 2; % raised sine power
rt_s = t_s(1:(nPtsPerRamp));
onRamp = sin(2*pi*(1/(4*durRamp_s)) * rt_s).^n;
offRamp = cos(2*pi*(1/(4*durRamp_s)) * rt_s).^n;
envBurst = [onRamp Plateau offRamp];
% The envelope contains 41 silences & 40 bursts
% (alternating, starting and ending in silence)
% There are 5 silences and bursts per second
env_1s = [silence envBurst ...
    silence envBurst ...
    silence envBurst ...
    silence envBurst ...
    silence envBurst];
% The final burst ends at 8 seconds,
% and is followed by 150 ms of silence
env = [env_1s env_1s ...
    env_1s env_1s ...
    env_1s env_1s ...
    env_1s env_1s ...
    silence];

dBstep = 2;
stimIs_dBSPLrms = 0:dBstep:80; % rms stimulus intensity in dB SPL

%nANFsPerFiberType = 2; % 2 for testing
nANFsPerFiberType = 40;
% 40 excitatory, 8 inhibitory LSO inputs (Gjoni et al. 2018)
%
% (NO ITDs)
% 
% For IID we apply a single stimulus
nStimReps = 1; % (If ever returning to multiple stimulus presentations, restore StimRep in the filename)
% spikeTimes_s = cell([length(CFs_Hz) length(ModFreqs_Hz) length(EnvStartPhases_deg) nStimReps length(fiberTypeStrings) nANFsPerFiberType]);
spikeTimes_s = cell([length(CFs_Hz) length(stimIs_dBSPLrms) nStimReps length(fiberTypeStrings) nANFsPerFiberType]);

% (ear can be used to set the lead and lag frequencies)
ears = cell(1,2);
ears{1} = 'LeadEar'; % left, right, ispi, contra, (often contra or left)
ears{2} = 'LagEar'; % left, right, ispi, contra, (often ipsi or right)

% +1 for contra, -1 for ipsi
LeadLag = [1 -1];

% species
% 1: Cat
% 2: Human (Q10 from Shera et al. 2002)
% 3: Human (Q10 from Glasberg & Brown 1990)
species = 3;

% noiseType
% 0: fixed fGn: (noise will be same in every simulation)
% 1: variable fGn:
noiseType = 1;


for ear = 1:length(ears)
    
    for idxCF = 1:length(CFs_Hz)
        
        % T_s = T_s_vector(idxCF);
        CF_Hz = CFs_Hz(idxCF); % Characteristic Frequency in Hz;
        
        % stimI_dBSPLrms = 0:dBstep:60;
        %for idxModF = 1:length(ModFreqs_Hz)
        for idxStimI = 1:length(stimIs_dBSPLrms)
            
            stimI_dBSPLrms = stimIs_dBSPLrms(idxStimI);
            % Stimulate at CF for both ears
            F0_Hz = CarrierFreqs_Hz(idxCF); % carrier frequency in Hz
            
            % Scale according to actual rms
            tone_unscaled = sin(2*pi*F0_Hz*t_s + 2*pi*rand); % carrier has random phase
            toneRMS1 = tone_unscaled/(rms(tone_unscaled)); % RMS set prior to onset/offset ramps
            unmodSoundPressure_Pa = 20e-6 * 10^(stimI_dBSPLrms/20) * toneRMS1; % scale for 0 dBSPL = 20e-6 Pascal
            soundPressure_Pa = env .* unmodSoundPressure_Pa; % apply the envelope
            
            if ear==1
                pltStr = 'b-';
            elseif ear==2
                pltStr = 'r-';
            end
            %figure; plot(t_s,pressureInput_Pa,pltStr)
            %title(['Carrier ' num2str(F0_Hz) ' Hz, TTRate' num2str(ModFreqHz) ' Hz, EnvStPhs ' num2str(EnvStartPhase_deg)])
            
            % Inner hair cell response from stimulus, not fiber type
            vihc = model2014_IHC_30Hz(soundPressure_Pa,CF_Hz,nrep,Ts_s,durTotal_s*2,cohc,cihc,species);
            
            for stimRep = 1:nStimReps
                for fiberType = 1:nFiberTypes
                    spikeCountThisFiberType = 0;
                    for ANFiber = 1:nANFsPerFiberType
                        % nrep = 1 for spike times from 1 stimulus
                        [meanrate,varrate,psth] = model2014_Synapse_30Hz(vihc,CF_Hz,nrep,Ts_s,fiberType,noiseType,implnt);
                        timeout_s = (1:length(psth))*Ts_s;
                        spikeTimes_s{idxCF,stimRep,fiberType,ANFiber} = timeout_s(psth==1);
                        spikeCountThisFiberType = spikeCountThisFiberType + length(spikeTimes_s{idxCF,stimRep,fiberType,ANFiber});
                    end  % ANFiber
                    
                    % (for computation efficiency in Matlab
                    % do not repeatedly concatenate spike times)
                    % Now that all ANFs are completed for this CF,
                    % loop again on ANFs.
                    % Sort fiber number & spike times for this CF by time,
                    % and write result to text file.
                    
                    ANFNumber_SpikeTimes = -1*ones(spikeCountThisFiberType,2);
                    spikeIdx = 0;
                    
                    for ANFiber = 1:nANFsPerFiberType
                        newSpikeTimes_s = spikeTimes_s{idxCF,stimRep,fiberType,ANFiber};
                        newSpikeCount = length(newSpikeTimes_s);
                        % EarLab, 1st column: fiber number, 2nd column: spiketimes
                        ANFNumber_SpikeTimes((spikeIdx + 1):(spikeIdx + newSpikeCount),1) = (ANFiber - 1);
                        ANFNumber_SpikeTimes((spikeIdx + 1):(spikeIdx + newSpikeCount),2) = newSpikeTimes_s;
                        % prepare spikeIdx for next AN fiber
                        spikeIdx = spikeIdx + newSpikeCount;
                    end
                    
                    % Sort by times, saving the indexes
                    [sortedSpikeTimes,sortTimesIdxs] = sort(ANFNumber_SpikeTimes(:,2));
                    ANFNumber_SpikeTimes_sortedByTime = [ANFNumber_SpikeTimes(sortTimesIdxs,1) sortedSpikeTimes];
                    
                    % For many stimReps, renameAN.py is unwieldy;
                    % prepare AN filename for Python, starting w strings
                    CFStr = num2str(CF_Hz/1000);
                    CFStr(find(CFStr == '.')) = 'p';
                    CFStrNoPad = CFStr;
                    roundCFStr = num2str(round(CF_Hz/1000));
                    while (length(roundCFStr) < 2)
                        roundCFStr = ['0' roundCFStr];
                        CFStr = ['0' CFStr];
                    end
                    
                    CaStr = num2str(F0_Hz/1000);
                    roundCaStr = num2str(round(F0_Hz/1000));
                    while (length(roundCaStr) < 2)
                        roundCaStr = ['0' roundCaStr];
                        CaStr = ['0' CaStr];
                    end
                    CaStr(find(CaStr == '.')) = 'p';
                    if (F0_Hz == CF_Hz)
                        CFCaStr = ['CaCF' CFStr 'kHz'];
                    else
                        CFCaStr = ['CF' CFStr 'kHzCa' CaStr 'kHz'];
                    end
                    stimRepStr = num2str(stimRep);
                    while (length(stimRepStr) < 2)
                        stimRepStr = ['0' stimRepStr];
                    end
                    %                         modFreqStr = num2str(ModFreqHz);
                    %                         while (length(modFreqStr) < 3)
                    %                             modFreqStr = ['0' modFreqStr];
                    %                         end
                    stimIStr = num2str(stimI_dBSPLrms);
                    while (length(stimIStr) < 2)
                        stimIStr = ['0' stimIStr];
                    end
                    
%                     filename = ['ANSpTs' 'Tone' CFCaStr 'Ear' num2str(ear) 'dBSPL' stimIStr 'Rep' stimRepStr fiberTypeStringsShort{fiberType} 'Sp' num2str(species)];
                    % As nStimReps equals 1, remove 'Rep' stimRepStr from the filename
                    filename = ['ANSpTs' 'Tone' CFCaStr 'Ear' num2str(ear) 'dBSPL' stimIStr fiberTypeStringsShort{fiberType} 'Sp' num2str(species)];
                    minusIndx = find(filename == '-');
                    for m = 1:length(minusIndx)
                        filename(minusIndx(m)) = 'N';
                    end
                    pointIndx = find(filename == '.');
                    for m = 1:length(pointIndx)
                        filename(pointIndx(m)) = 'p';
                    end
                    filename = [filename '.txt'];
                    % *** ADD cd according stimRep
                    % Update as necessary:
                    fpath = ['C:\Users\Andrew Brughera\lso\IID\' CFStrNoPad 'kHz\'];
                    % fpath = ['C:\Users\Andrew Brughera\lso\IID\' CFStrNoPad 'kHz\stimReps\rep' stimRepStr '\'];
                    fid0 = fopen([fpath filename],'w');
                    
                    for m = 1:size(ANFNumber_SpikeTimes_sortedByTime,1)
                        count_fpr = fprintf(fid0,'%2d %11.9f\n', ANFNumber_SpikeTimes_sortedByTime(m,:));
                    end
                    status_fcl = fclose(fid0);
                    
                end  % fiberType
                
            end % stimRep
            
        end % idxStimI
        
    end  % idxCF
    
end % ear
