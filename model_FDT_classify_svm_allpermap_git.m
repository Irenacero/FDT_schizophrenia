% Aim: Classification between two groups of data with a linear
% support vector machine (SVM) model using cross-validation
% Input: FDT calculation of two groups of data
% Output: accuracy and confusion matrices of classification

% Irene Acero & Paulina Clara Dagnino, Upf, April 2023 
% Original code sent by Gustavo Deco 

%%%%%%%%%%%%%%%%%%%%%%%%%%%% Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;

N=80; % number of brain nodes to analyze 
name='UCLA'; % name of dataset, string, for storing workspace and figures
groups={'ucla_schizophrenia_dbs80', 'ucla_subsetcontrols_dbs80'};

% paths
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
pathDependencies=join([generalPathScripts{1},'\Scripts\Dependencies\']);

if strcmp(system,'linux')
    generalPathScripts= replace(generalPathScripts,'\','/');
    pathResults= replace(pathResults,'\','/');
    pathFigures= replace(pathFigures,'\','/');
    pathDependencies= replace(pathDependencies,'\','/');
end

addpath(pathResults)
addpath(pathFigures)
addpath(pathDependencies)


thres=0.01; % threshold for stats
perm=10; %10000; % iterations for stats
kfold=1000; % number of iterations for SVM cross-validation
perc=80; % percentage for training svm
test = 'ranksum'; % statistical test

normalization='zscore'; % 'zscore' if normalization with zscore
nodes='all'; % 'all' if analysis is done with all nodes, 'top' if analysis is done with top nodes according to SSND 
flag_zscore=1; % If zscored values, 1. If not zscored values, 0.
flag_stats=1; % If top nodes, 1. If all nodes, 0.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i = 1:length(groups)
    group_name = groups{i}; % name of the group
    
    path = join([pathResults, group_name, '_FDT_results_Tau2_nofiltfilt.mat']); % path to the data of the group
    load(path) % load the data   

    FDTvarname = sprintf('perFDT_subjects%d', i);
    eval([FDTvarname ' = perFDT_subjects;']);

    FDTvarname = sprintf('perCeff_subjects%d', i);
    eval([FDTvarname ' = perCeff_subjects;']);

    FDTvarname = sprintf('perFC_subjects%d', i);
    eval([FDTvarname ' = perFC_subjects;']);

end

n_groups = floor(1:size(groups, 2));
combinations = nchoosek(n_groups, 2);

for i=1:size(combinations, 1)

    i_comb = combinations(i, :);
    group1_FDT =  eval(join(['perFDT_subjects', num2str(i_comb(1))]));
    group2_FDT =  eval(join(['perFDT_subjects', num2str(i_comb(2))]));
    group1_Ceff =  eval(join(['perCeff_subjects', num2str(i_comb(1))]));
    group2_Ceff =  eval(join(['perCeff_subjects', num2str(i_comb(2))]));
    group1_FC =  eval(join(['perFC_subjects', num2str(i_comb(1))]));
    group2_FC =  eval(join(['perFC_subjects', num2str(i_comb(2))]));

        
    NSUB1=size(group1_FDT,1); % number of subjects in group 1
    NSUB2=size(group2_FDT,1); % number of subjects in group 2


    % loop in all (FC, Ceff, FDT)

    group1_all={group1_FC, group1_Ceff, group1_FDT};
    group2_all={group2_FC, group2_Ceff, group2_FDT};

    for j=1:length(group1_all)
        clear DataAll1 DataAll2
        group1=group1_all{j};
        group2=group2_all{j};

        if strcmp(nodes,'all') 
            xxdata=[group1;group2]; 
        else
            [sopm, indpp, pp] = function_statistics_ssnd(group1,group2,N,test,thres, perm);
            [aux top]=min(pp);
            xxdata=[group1(:,indpp(1:top));group2(:,indpp(1:top))]; 
            top_all(j,:)=top;
            indpp_all(j,:)=indpp;
        end
        
        if strcmp(normalization,'zscore') 
            xxdata=zscore(xxdata); % zscore normalizes every value in a dataset such that 
            % the mean of all of the values is 0 and the standard deviation is 1. 
        end
        [coeff, score, latent, tsquared, explained, mu] = pca(xxdata);
        DataAll1=score(1:NSUB1,1:2);
        DataAll2=score(NSUB1+1:NSUB1+NSUB2,1:2);
        [pcmat(i,j,:,:),acc(i,j), acc_all{i,j},top_10_indices{i,j}]=function_svm(DataAll1, DataAll2, kfold, perc);
        disp(pcmat)
    end
    

    %%%% concatenated (FDT and Ceff) %%%%function
    
    if strcmp(nodes,'all') 
        xxdata=[group1_FDT group1_Ceff; group2_FDT group2_Ceff];

    else
        xxdata=[group1_FDT(:,indpp_all(3,(1:top_all(3,:)))) group1_Ceff(:,indpp_all(2,(1:top_all(2,:))));group2_FDT(:,indpp_all(3,(1:top_all(3,:)))) group2_Ceff(:,indpp_all(2,(1:top_all(2,:))))];
    end

    if strcmp(normalization,'zscore') 
        xxdata=zscore(xxdata);
    end

    [coeff, score, latent, tsquared, explained, mu] = pca(xxdata);
    DataAll1=score(1:NSUB1,1:2);
    DataAll2=score(NSUB1+1:NSUB1+NSUB2,1:2);
    [pcmat(i,4,:,:),acc(i,4), acc_all{i,4}, top_10_indices{i,4}]=function_svm(DataAll1, DataAll2, kfold, perc);
    disp(pcmat)
end

% save workspace
save(join([pathResults,sprintf('results_FDT_classification_%s.mat',name)]), 'acc', 'pcmat');

%EOF