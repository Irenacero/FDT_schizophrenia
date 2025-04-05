
% Aim: This computes the average frequency of the signal for each ROI (region of interest) for each subject
% Input: 
% - raw BOLD fMRI of dataset (ts), matrix of node x time points
% Output: frequencies (f_diff), matrix of 1 x node

% Irene Acero & Paulina Clara Dagnino, Upf, April 2023 
% Original code sent by Gustavo Deco

clear all;

%%%%%%%%%%%%%%%%%%%%%%%%%% Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%

N=80; % number of brain nodes to analyze 
indexN=1:N; % indexes of nodes to analyze (e.g., if subcorticals want to be excluded: [1:31 50:80])
NSUB_groups=[50,50]; % number of subjects in each group 
groups={'ucla_schizophrenia_dbs80','ucla_subsetcontrols_dbs80'}; % names of groups to analyze 

% parameters of the data
TR=2;  % repetition Time (seconds)


system='linux';

filePath = matlab.desktop.editor.getActiveFilename;
fprintf('%s\n',filePath);
if strcmp(system,'linux')
    myFolders = split(filePath,"/");
else
    myFolders = split(filePath,"\");
end

generalPathScripts = join(myFolders(1:length(myFolders)-3),"\");
pathResults=join([generalPathScripts{1},'\Output\Results\']);
pathFigures=join([generalPathScripts{1},'\Output\Figures\']);
pathRawData=join([generalPathScripts{1},'\Data']);
pathDependencies=join([generalPathScripts{1},'\Scripts\Dependencies\']);

if strcmp(system,'linux')
    generalPathScripts= replace(generalPathScripts,'\','/');
    pathResults= replace(pathResults,'\','/');
    pathFigures= replace(pathFigures,'\','/');
    pathRawData= replace(pathRawData,'\','/');
    pathDependencies= replace(pathDependencies,'\','/');
end

addpath(pathResults)
addpath(pathFigures)
addpath(pathRawData)
addpath(pathDependencies)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% in each group, for each subject, access the BOLD signal and calculate the frequency
for g=1:length(groups)

    disp(join(['Group:', groups(g)]))

    ts_g=load(join([pathRawData, groups{g}, '.mat'])); % load BOLD fMRI signals as a cell each row a subject
    name=fieldnames(ts_g); % name of the fieldnames in the structure
    name=name(1); % name of the first fieldname in the structure
    ts_g=ts_g.(string(name)); % BOLD fMRI
    NSUB=NSUB_groups(g); % look for number of subjects in this group


    for sub=1:NSUB

        disp(join(['Subject: ', num2str(sub)]))

        clear signal_filt Power_Areas;
        ts=ts_g{1,sub}.dbs80ts; % access the time signal of a particular subject
        ts=ts(indexN,:); % obtain the brain nodes of interest
    
        disp('Preprocessing BOLD fMRI signal')

        % filter and clean the data
        for seed=1:N % for each voxel (region of interest)
            ts(seed,:)=detrend(ts(seed,:)-mean(ts(seed,:))); % remove the mean value from the timeseries and detrend to remove any linear trend 
        end
        signal_filt=ts(:,10:end-10); % remove filter effects in extremities
        [Ns, Tmaxred]=size(signal_filt); % obtain the number of brain nodes and maximum time duration of signal (as points)
        TT=Tmaxred;
        Ts = TT*TR; % convert the time of signals to seconds
        freq = (0:TT/2-1)/Ts; % convert to frequency as 1/time(s)
        nfreqs=length(freq); % look for the amount of frequency points
    
       disp('Computing frequencies for each brain node')

        % for each brain node, look for the frequencies
        for seed=1:N
            pw = abs(fft(zscore(signal_filt(seed,:)))); % absolute value of the Fast Fourier Transform of the z-scored signal
            PowSpect = pw(1:floor(TT/2)).^2/(TT/TR); % calculate the power spectrum
            Power_Areas=gaussfilt(freq,PowSpect,0.005); % apply Gaussian filter to smooth the power spectrum and reduce noise 
            [~,index]=max(Power_Areas); % look for index of the the maximum power
            index=squeeze(index);
            f_diff_sub(sub,seed)=freq(index); % save the frequency which has highest power 
        end

    end
    
    f_diff = mean(f_diff_sub,1); % average frequencies for all subjects 
    
    save(join([pathResults,sprintf('empirical_%s.mat',groups{g})],'/'), 'f_diff', 'f_diff_sub');

end 
