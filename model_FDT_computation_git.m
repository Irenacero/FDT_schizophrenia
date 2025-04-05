% Aim: Compute FDT violation and perturbablitiy maps 
% Input: 
% - raw BOLD fMRI of dataset (ts), matrix of node x time points
% - structural connectivity of dataset (SC), matrix node x node 
% - empirical frequencies computed in empirical_freq.m (f_diff),  matrix of 1 x node
% Output: 
% - FDT deviation values (FDTm_subjects), matrix of 1 x subjects 
% - perturbability maps of FDT (perFDT_subjects), matrix of subjects x
% nodes
% - Effective Connectivity (perCeff_subjects), matrix of subjects x nodes
% - Functional Connectivity (perFC_subjects), matrix of subjects x nodes

% Irene Acero & Paulina Clara Dagnino, Upf, April 2023 
% Original code sent by Gustavo Deco

clear all;

%%%%%%%%%%%%%%%%%%%%%%%%%% Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%

rng shuffle;
N=80; % number of brain nodes to analyze 
indexN=1:N; % indexes of nodes to analyze (e.g., if subcorticals want to be excluded: [1:31 50:80])
Tau=2; %; % lag of time (in seconds)
groups={'ucla_schizophrenia_dbs80'}; % names of groups to analyze 

% parameters of the data
TR=2;  % repetition Time (in seconds)

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


sigma=0.01;
epsFC=0.004;
epsFCtau=0.01;
maxC=0.2;

sc = join([pathRawData, 'SC_dbs80HARDIFULL.mat']); % load structural connectivity matrix
load (sc)
C = SC_dbs80HARDI(indexN,indexN);  % select the structrual connectivity matrix and the nodes of interest
C = C/max(max(C))*maxC; 

Isubdiag = find(tril(ones(N),-1));

%% FDT computation
for i = 1:length(groups)
    clear COVtauemp sigratiosim  FDTm_subjects FC_subjects COVtau_subjects Ceff_subjects fittFC_subjects fittCVtau_subjects perCeff_subjects perFC_subjects perFDT_subjects;
    group_name = groups{i}; % name of the group

   
    path = join([pathRawData, group_name, '.mat']); % path to the time series data of the group
    path_empirical = join([pathResults, 'empirical_', group_name, '.mat']); 

    load(path) % load the group time series
    load(path_empirical) % load the group empirical data

    
    disp(join(['Computing FDT of group ', group_name]))
    data_group = subject;
    NSUB = size(data_group, 2);
    
    f_diff = f_diff; % empirical frequencies
    
    % % Group analysis
    for nsub=1:NSUB
        ts=data_group{1, nsub}.dbs80ts;
        clear tss
        for seed=1:N 
            ts(seed,:)=detrend(ts(seed,:)-nanmean(ts(seed,:))); % remove the mean value from the timeseries and detrend to remove any linear trend 
            tss = ts;
        end
        % FC
        ts2=tss(indexN,10:end-10); % remove filter effects in extremities
        Tm=size(ts2,2); % number of time points
        FCemp=corrcoef(ts2'); % empirical functional connectivity
        FC_subjects(nsub,:,:)=FCemp; 
        COVemp=cov(ts2'); % empirical covariance matrix
        % COV(tau)
        tst=ts2';
        for i=1:N
            for j=1:N
                sigratio(i,j)=1/sqrt(COVemp(i,i))/sqrt(COVemp(j,j)); % scaling factor based on the empirical covariance matrix
                [clag lags] = xcov(tst(:,i),tst(:,j),Tau); %  empirical cross-covariance signals at time lag Tau 
                indx=find(lags==Tau); % index of the time lag of interest (Tau)
                COVtauemp(i,j)=clag(indx)/size(tst,1);  % empirical covariance at time lag Tau 
            end
        end
        COVtauemp=COVtauemp.*sigratio; % empirical covariance at time lag Tau scaled
        COVtau_subjects(nsub,:,:)=COVtauemp; % store data for each subject
    end

    FCemp=squeeze(mean(FC_subjects)); % mean of empirical functional connectivity of all subjects
    COVtauemp=squeeze(mean(COVtau_subjects)); % mean of empirical covariance at time lag Tau of all subjects

    Cnew=C; % initial structural connectivity matrix 
    olderror=100000; % initial error
    for iter=1:500
        [FCsim,COVsim,COVsimtotal,A]=hopf_int(Cnew,f_diff,sigma); % simulated functional connectivity matrix (FCsim), covariance matrix (COVsim), total covariance matrix (COVsimtotal), Jacobian matrix (A)
        COVtausim=expm((Tau*TR)*A)*COVsimtotal; % total simulated covariance at time lag Tau
        COVtausim=COVtausim(1:N,1:N); % simulated covariance at time lag Tau (nodes of interest)
        for i=1:N
            for j=1:N
                sigratiosim(i,j)=1/sqrt(COVsim(i,i))/sqrt(COVsim(j,j)); % scaling factor based on the simulated covariance matrix
            end
        end
        COVtausim=COVtausim.*sigratiosim; % simulated covariance at time lag Tau scaled
        errorFC(iter)=mean(mean((FCemp-FCsim).^2)); % difference between empirical and simulated functional connectivity
        errorCOVtau(iter)=mean(mean((COVtauemp-COVtausim).^2)); % difference between empirical and simulated covariance matrix

        if mod(iter,100)<0.1
            errornow=mean(mean((FCemp-FCsim).^2))+mean(mean((COVtauemp-COVtausim).^2));
            if  (olderror-errornow)/errornow<0.001 % if the curent error is smaller than 0.001 from last iteration
                break;
            end
            if  olderror<errornow % if the current error is larger than the one from last iteration
                break;
            end
            olderror=errornow; % update old error by current error
        end

        for i=1:N  
            for j=1:N
                if (C(i,j)>0 || j==N-i+1)
                    Cnew(i,j)=Cnew(i,j)+epsFC*(FCemp(i,j)-FCsim(i,j)) ...
                        +epsFCtau*(COVtauemp(i,j)-COVtausim(i,j)); % update SC matrix 
                    if Cnew(i,j)<0 % put to 0 negative values
                        Cnew(i,j)=0;
                    end
                end
            end
        end
        Cnew = Cnew/max(max(Cnew))*maxC;
    end
    Ceffgroup=Cnew;

    % Individual analysis: fit an indiviualized coupling matrix C for each
    % participant
    for nsub=1:NSUB
        disp(nsub)
        ts=data_group{1, nsub}.dbs80ts;
        clear tss
        for seed=1:N
            ts(seed,:)=detrend(ts(seed,:)-nanmean(ts(seed,:)));
            tss = ts;
        end
        % FC(0)
        ts2=tss(indexN,10:end-10);
        Tm=size(ts2,2);
        FCemp=corrcoef(ts2');
        FC_subjects(nsub,:,:)=FCemp;
        COVemp=cov(ts2');
        % COV(tau)
        tst=ts2';
        for i=1:N
            for j=1:N
                sigratio(i,j)=1/sqrt(COVemp(i,i))/sqrt(COVemp(j,j));
                [clag lags] = xcov(tst(:,i),tst(:,j),Tau);
                indx=find(lags==Tau);
                COVtauemp(i,j)=clag(indx)/size(tst,1);
            end
        end
        COVtauemp=COVtauemp.*sigratio;
        Cnew=C;%Ceffgroup;
        olderror=100000;
        for iter=1:500
            % Obtain simulated data
            [FCsim,COVsim,COVsimtotal,A]=hopf_int(Cnew,f_diff,sigma);
            COVtausim=expm((Tau*TR)*A)*COVsimtotal;
            COVtausim=COVtausim(1:N,1:N);
            for i=1:N
                for j=1:N
                    sigratiosim(i,j)=1/sqrt(COVsim(i,i))/sqrt(COVsim(j,j));
                end
            end
            COVtausim=COVtausim.*sigratiosim;
            errorFC_subjects(nsub,iter)=mean(mean((FCemp-FCsim).^2));
            errorCOVtau_subjects(nsub,iter)=mean(mean((COVtauemp-COVtausim).^2));
    
            if mod(iter,50)<0.1
                errornow=mean(mean((FCemp-FCsim).^2))+mean(mean((COVtauemp-COVtausim).^2));
                if  (olderror-errornow)/errornow<0.001
                    break;
                end
                if  olderror<errornow
                    break;
                end
                olderror=errornow;
            end
    
            for i=1:N 
                for j=1:N
                    if (C(i,j)>0 || j==N-i+1)
                        Cnew(i,j)=Cnew(i,j)+epsFC*(FCemp(i,j)-FCsim(i,j)) ...
                            +epsFCtau*(COVtauemp(i,j)-COVtausim(i,j));
                        if Cnew(i,j)<0
                            Cnew(i,j)=0;
                        end
                    end
                end
            end
            Cnew = Cnew/max(max(Cnew))*maxC;
        end
        Ceff=Cnew;
        Ceff_subjects(nsub,:,:)=Ceff;
        [FCsim,COVsim,COVsimtotal,A]=hopf_int(Ceff,f_diff,sigma);
        fittFC_subjects(nsub)=corr2(FCemp(Isubdiag),FCsim(Isubdiag));
        COVtausim=expm((Tau*TR)*A)*COVsimtotal;
        COVtausim=COVtausim(1:N,1:N);
        for i=1:N
            for j=1:N
                sigratiosim(i,j)=1/sqrt(COVsim(i,i))/sqrt(COVsim(j,j));
            end
        end
        COVtausim=COVtausim.*sigratiosim;
        fittCVtau_subjects(nsub)=corr2(COVtauemp(Isubdiag),COVtausim(Isubdiag));
    end

    % FDT violation: complute the expectation values when a perturbation is
    % applied

    for nsub=1:NSUB
        Ceff=squeeze(Ceff_subjects(nsub,:,:));
        [FCsim,COVsim,COVsimtotal,A]=hopf_int(Ceff,f_diff,sigma); 
        invA=inv(A); % inverse of the Jacobian matrix 
        for i=1:2*N
            for j=1:2*N
                hh=zeros(2*N,1); % vector with all components equal to 0 exccelt the j component, which equals the value of the perturbation
                hh(j)=1; % the perturbation equals 1
                xepsilon=-invA*hh; 
                chi(i,j)=abs((2*COVsimtotal(i,j)/sigma^2)-xepsilon(i)); % final deviation from the FDT
                chi2(i,j)=abs(xepsilon(i));
            end
        end
        chij=mean(chi(1:N,1:N))./mean(chi2(1:N,1:N)); % average over regions
        FDTm_subjects(nsub)=mean(chij); % level of non-equilibrium
        perFDT_subjects(nsub,:)=chij; % perturbability map over all brain regions
        perCeff_subjects(nsub,:)=mean(Ceff);
        perFC_subjects(nsub,:)=mean(squeeze(FC_subjects(nsub,:,:)));

       
    end
    
    group_name = join([group_name, '_FDT_results_Tau2_nofiltfilt']);
    save (fullfile(pathResults, group_name), 'FDTm_subjects', 'perFDT_subjects', 'perCeff_subjects', 'perFC_subjects', 'Ceff_subjects', 'fittFC_subjects', 'fittCVtau_subjects', 'errorFC_subjects', 'errorCOVtau_subjects')
end

