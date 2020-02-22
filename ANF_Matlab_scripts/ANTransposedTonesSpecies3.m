% ANTransposedTones.m
% Andrew Brughera 2019, Jan 8
% Modified from ANSR100kHzAMwFSPhsSpecies3.m (of 2018 Aug 9)
% which lists ANSR100kHzAMStaticIPD1Species3.m as a source.
% (AM and Static IPD, was a control condition for 
% amplitude modulated binaural beats)
% (also consider 2 control stimuli:
% binaural beats with no AM, AM tones with a fixed, static IPD)
% Andrew Brughera 8 June, 2017
% 2019, Jan 8, Now Applying the Zilany et al. 2014 model 
% with transposed tones:
%
% Fs_Hz to 100 kHz, compatible with Zilany et al. 2014 model
% Checking the effects of CF on response to 32-Hz pure tone at 60 dBSPL rms
% CFs_Hz = [4000 6000 10000];
% low-, medium-, and high-spontaneous-rate AN fibers simulated
% only 2 fibers per CF for quick preliminary results
%
rng('shuffle'); % Seed random number generator rng, based on present time
rng_state = rng; % Save the rng Type, Seed, and State 
% model fiber parameters
%clear all;
% (CF varied) CF in Hz;
cohc  = 1.0;   % normal ohc function
cihc  = 1.0;   % normal ihc function
%fiberType = 3; % spontaneous rate (in spikes/s) of the fiber BEFORE refractory effects; "1" = Low, 0 spk/s; "2" = Medium, 5 spk/s; "3" = High, 100 spk/s
implnt = 1;    % "0" for approximate or "1" for actual implementation of the power-law functions in the Synapse

% stimulus parameters
Fs_Hz = 100e3;  % sampling rate of 100,000 Hz;
Fs_kHzStr = num2str(round(Fs_Hz/1000));
Ts_s = 1/Fs_Hz;
% Zilany et al. 2009, 2014, adaptation coefficients designed for
% sample rate of 100 kHz, "must be 100 kHz".

% PSTH parameters -> make these SpikeTime parameters instead
% We can make a psth later, but
nrep = 1;              % number of stimulus repetitions (e.g., 50);
%psthbinwidth_s = 1/(80 * F0_Hz) % 0.5e-3; % binwidth in seconds;
psthbinwidth_s = Ts_s; % binwidth in seconds at full temporal resolution;

% fiberTypeStrings = {'Low Spont'; 'Medium Spont'; 'High Spont'};
% fiberTypeStringsShort = {'LoSp'; 'MdSp'; 'HiSp'};
fiberTypeStrings = {'Medium Spont'};
fiberTypeStringsShort = {'MdSp'};
nFiberTypes = length(fiberTypeStringsShort);

CarrierFreqs_Hz = [4000 6000 10000];
%BBeatFreqs_Hz = [4 8 16 32 64];
ModFreqs_Hz = [32 64 128 256 512 600 700 800 900];
%ModFreqs_Hz = [64 128 256 512 800];
% (Generate stimuli and model AN spike responses in Matlab.)
% (Generate bushy cell spike responses in EarLab)
% (Apply the beat direction at the model MSO inputs in NEURON)
% BBeatDirections = [1 -1];
% Dietz et al.:
% "The polarity of the IPD was defined with respect to the AMBB beat direction:
% positive IPD means leading on the side with the higher carrier
% frequency."

% Bernstein & Rtahiotis, 2002: Starting phases random for each component
%StartPhases_deg = [0 -45 -90 -135 -180 -225 -270 -315]; % phase lags
%StartPhases_deg = 2*pi*rand(1,24); % carriers: 24 random phases 
StartPhases_rad = 2*pi*rand(1,49); % carriers: 49 random phases 
EnvStartPhases_deg = 0:-15:-720; % envelopes: 24+24+1 start-phase decrements
% allowing 25 instances of either ear trailing by up to 360 deg
% Also present AM stimuli with same modulation freqs and static IPD
minimumZilanyCF_Hz = 56;
CFs_Hz = max(CarrierFreqs_Hz,minimumZilanyCF_Hz);
% CFs_Hz = [56];
% Minimum CF allowed by Zilany et al. model is 91.2 Hz
% In c code by Zilany et al. 2009, 2014:
% /** Calculate the location on basilar membrane from CF */
% bmplace = 11.9 * log10(0.80 + cf / 456.0);

% Stimulus duration 300 ms, T_s = 0.3
Tdur_s = 0.3;
% rise/fall time in seconds
rTdur_s = 0.02; % Bernstein & Trahiotis, 2002: ramps 20-ms cosine-squared
t_s = 0:Ts_s:Tdur_s-Ts_s;
mxpts = length(t_s);
irpts = rTdur_s*Fs_Hz;
n = 2; % raised sine power
rt_s = t_s(1:(irpts));
onRamp = sin(2*pi*(1/(4*rTdur_s)) * rt_s).^n;
offRamp = cos(2*pi*(1/(4*rTdur_s)) * rt_s).^n;

stimI_dBSPLrms = 75; % rms stimulus intensity in dB SPL


% setup for  maximum number of MSO inputs, Bushy Cells (BCs), AN fibers
% maxLSOInputsContraI = 4;
% maxLSOInputsIpsiI = maxMSOInputsContraI;
% maxLSOInputsContraE = 4;
% maxLSOInputsIpsiE = maxMSOInputsContraE;
% nBCs = maxMSOInputsContraI + maxMSOInputsIpsiI + maxMSOInputsContraE + maxMSOInputsIpsiE;
% nANFsPerSBC = 3;
% nANFsPerGBC = 30;
% On each side (left and right)
%nANFsPerFiberType = 132; % The EarLab BushyCell modules select M fibers at random from nANFsPerFiberType
%nANFsPerFiberType = 1; % 1 for testing
nANFsPerFiberType = 40; % 20 contraI, 20 ipsiE direct inputs to model LSO neuron,
% the default values in Wang and Colburn, 2012.
% 40 excitatory, 8 inhibitory inputs estimated by Gjoni et al. 2018
% ITD will be generated from shuffled phases of envelope, so don't need
% % independent spikes for each ITD tested
% ITDsPerCarrierPeriod = 20;
% nITDs = 2*ITDsPerCarrierPeriod + 1;

spikeTimes_s = cell([length(CFs_Hz) length(ModFreqs_Hz) length(EnvStartPhases_deg) length(fiberTypeStrings) nANFsPerFiberType]);
% Synchrony Index (help synccalc)
%syncIdx   = -1 * ones([length(CFs_Hz) length(ModFreqs_Hz) length(EnvStartPhases_deg) length(fiberTypeStrings) nANFsPerFiberType]);
% Rayleigh Criterion for significance: 2*number_of_spikes*(syncIdx^2) > 13.8
%two_n_Rsq  = -1 * ones([length(CFs_Hz) length(ModFreqs_Hz) length(EnvStartPhases_deg) length(fiberTypeStrings) nANFsPerFiberType]);

% (ear can be used to set the lead and lag frequencies)
ears = cell(1,2);
ears{1} = 'LeadEar'; % left, right, ispi, contra, (often contra or left)
ears{2} = 'LagEar'; % left, right, ispi, contra, (often ipsi or right)

% +1 for contra, -1 for ipsi
LeadLag = [1 -1];

% %better to use a structure?
% earStruct = struct('side',{'Contra','Ipsi'},'LeadLag',{1,-1});

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
        
        for idxModF = 1:length(ModFreqs_Hz)
            
            ModFreqHz = ModFreqs_Hz(idxModF);
            % For AM with static IPD, stimulate at CF for both ears
            F0_Hz = CarrierFreqs_Hz(idxCF); % carrier frequency in Hz
            %F0_Hz = CarrierFreqs_Hz(idxCF) + (LeadLag(ear)*(ModFreqHz/2)); % carrier frequency in Hz
            
            for idxEnvStartPhase = 1:length(EnvStartPhases_deg)
                
                EnvStartPhase_deg = EnvStartPhases_deg(idxEnvStartPhase)
                EnvStartPhase_rad = (2*pi/360)*EnvStartPhase_deg;
                % SETTLE
                %EnvStartPhase_rad = 0
                EnvStPhsDegString = num2str(EnvStartPhase_deg);
                leadingZeros = '';
                if (EnvStartPhase_deg == 0)
                    EnvStPhsDegString = '000';
                elseif (EnvStartPhase_deg < 0)
                    EnvStPhsDegString(1) = 'N';
                    if (EnvStartPhase_deg > -10)
                        leadingZeros = '00';
                    elseif (EnvStartPhase_deg > -100)
                        leadingZeros = '0';
                    end
                    EnvStPhsDegString =[EnvStPhsDegString(1) leadingZeros EnvStPhsDegString(2:end)];
                elseif (EnvStartPhase_deg > 0)
                    if (EnvStartPhase_deg < 10)
                        leadingZeros = '00';
                    elseif (EnvStartPhase_deg < 100)
                        leadingZeros = '0';
                    end
                    EnvStPhsDegString =[leadingZeros EnvStPhsDegString];
                end
                
                % Transposed tone
                % Time domain: half-wave rectified envelope
                AM = sin(2*pi*ModFreqHz*t_s + EnvStartPhase_rad);
                % (the negative phases are delays)
                env = AM .* (AM >= 0);
                % Frequency domain: remove envelope components above 2 kHz
                ENV_ALLF = fft(env);
                % fres_Hz = Fs_Hz/(Tdur_s * Fs_Hz)
                fres_Hz = 1/Tdur_s; 
                % Truncate components above 2 kHz
                freq_Hz = 0:fres_Hz:(Fs_Hz - fres_Hz);
                fLim_Hz = 2000;
                freq_mask = ((freq_Hz <= fLim_Hz) | (freq_Hz >= (Fs_Hz - fLim_Hz)));
                ENV_2KHZ = ENV_ALLF .* freq_mask;
                env_2kHz = ifft(ENV_2KHZ);
                % figure; plot(t_s, env_2kHz)
                % axis([0 6/ModFreqHz -0.1 1.1])
                
                % Scale according to actual rms
                ttone_unscaled = env_2kHz.*sin(2*pi*F0_Hz*t_s + 2*pi*rand); % carrier has random phase
                ttoneRMS1 = ttone_unscaled/(rms(ttone_unscaled)); % RMS set prior to onset/offset ramps
                pressureInput_Pa = 20e-6 * 10^(stimI_dBSPLrms/20) * ttoneRMS1; % scale for 0 dBSPL = 20e-6 Pascal
                
                % Apply on/off ramps (sin^2 / cos^2)
                pressureInput_Pa(1:(irpts))=pressureInput_Pa(1:irpts).*onRamp;
                pressureInput_Pa((mxpts-irpts+1):mxpts)=pressureInput_Pa((mxpts-irpts+1):mxpts).*offRamp;
                if ear==1
                    pltStr = 'b-';
                elseif ear==2
                    pltStr = 'r-';
                end
                %figure; plot(t_s,pressureInput_Pa,pltStr)
                %title(['Carrier ' num2str(F0_Hz) ' Hz, TTRate' num2str(ModFreqHz) ' Hz, EnvStPhs ' num2str(EnvStartPhase_deg)])
                
                for fiberType = 1:nFiberTypes
                    
                    spikeCountThisFiberType = 0;
                    
                    for ANFiber = 1:nANFsPerFiberType
                        %for ANFiber = 1:2   % short test
                        
                        vihc = model2014_IHC_30Hz(pressureInput_Pa,CF_Hz,nrep,Ts_s,Tdur_s*2,cohc,cihc,species);
                        [meanrate,varrate,psth] = model2014_Synapse_30Hz(vihc,CF_Hz,nrep,Ts_s,fiberType,noiseType,implnt);
                        
                        timeout_s = (1:length(psth))*Ts_s;
                        %     psth_samples_per_bin = round(psthbinwidth_s*Fs_Hz);  % number of psth bins per psth bin
                        %     psthtime_s = timeout_s(1:psth_samples_per_bin:end); % time vector for psth
                        %     pr = sum(reshape(psth,psth_samples_per_bin,length(psth)/psth_samples_per_bin))/nrep; % pr of spike in each bin
                        %     Psth = pr/psthbinwidth_s; % psth in units of spikes/s
                        
                        spikeTimes_s{idxCF,fiberType,ANFiber} = timeout_s(psth==1);
                        spikeCountThisFiberType = spikeCountThisFiberType + length(spikeTimes_s{idxCF,fiberType,ANFiber});
                        %[syncIdx(idxCF,fiberType,ANFiber),two_n_Rsq(idxCF,fiberType,ANFiber)] = synccalc(spikeTimes_s{idxCF,fiberType,ANFiber},1/ModFreqHz);
                        
                    end  % ANFiber
                    
                    % (for computation efficiency in Matlab
                    % do not repeatedly concatenate spike times)
                    % Now that all ANFs are completed for this CF,
                    % loop again on ANFs.
                    % Sort fiber number & spike times for this CF by time,
                    % and write result to text file.
                    
                    % filename
                    % THE F0 AND AM FREQUENCY:
                    % 'AM' num2str(ModFreq_Hz) 'Hz')
                    % For AM with static IPD, insert the ear number after the carrier frequency 
                    filename = ['ANSpTs' 'CF' num2str(CF_Hz/1000) 'kHz' 'Ca' num2str(F0_Hz/1000) 'kHz' num2str(ear) 'EPhs' EnvStPhsDegString 'TT' num2str(ModFreqHz) 'Hz' num2str(stimI_dBSPLrms) 'dBspl' fiberTypeStringsShort{fiberType} 'Spec' num2str(species)];
                    minusIndx = find(filename == '-');
                    for m = 1:length(minusIndx)
                        filename(minusIndx(m)) = 'N';
                    end
                    pointIndx = find(filename == '.');
                    for m = 1:length(pointIndx)
                        filename(pointIndx(m)) = 'p';
                    end
                    filename = [filename '.txt'];
                    fid0 = fopen(filename,'w');
                    
                    ANFNumber_SpikeTimes = -1*ones(spikeCountThisFiberType,2);
                    spikeIdx = 0;
                    
                    for ANFiber = 1:nANFsPerFiberType
                        newSpikeTimes_s = (spikeTimes_s{idxCF,fiberType,ANFiber});
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
                    
                    for m = 1:size(ANFNumber_SpikeTimes_sortedByTime,1)
                        count_fpr = fprintf(fid0,'%2d %11.9f\n', ANFNumber_SpikeTimes_sortedByTime(m,:));
                    end
                    status_fcl = fclose(fid0);
                    
                end  % fiberType
                
            end % idxEnvStartPhase
            
        end % idxModF
        
    end  % idxCF
    
end % ear

% List CFs, Fiber Types, synchrony indexes, Rayleigh criterion value
% CFs_Hz
% fiberTypeStrings
% syncIdx
% two_n_Rsq
