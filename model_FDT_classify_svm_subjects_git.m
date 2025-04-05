% Aim: Classification between two groups of data with a linear
% support vector machine (SVM) model using cross-validation
% Input: FDT calculation of two groups of data
% Output: accuracy and confusion matrices of classification

% Irene Acero & Paulina Clara Dagnino, Upf, April 2023 
% Original code sent by Gustavo Deco 

%%%%%%%%%%%%%%%%%%%%%%%%%%%% Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;

N=80; % number of brain nodes to analyze 
groups={'ucla_schizophrenia_dbs80','ucla_subsetcontrols_dbs80'}; % names of groups to analyze 

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


thres=0.05; % threshold for stats
perm=10000; %10000; % iterations for stats
kfold=1000; %10000; % number of iterations for SVM cross-validation
perc=75; % percentage for training svm
test = 'ranksum'; % statistical test

normalization=0; % 'zscore' if normalization with zscore
nodes='all'; % 'all' if analysis is done with all nodes, 'top' if analysis is done with top nodes according to SSND 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i = 1:length(groups)
    group_name = groups{i}; % name of the group
    
    path = join([pathResults, group_name, '_FDT_results_Tau2_nofiltfilt.mat']); % path to the data of the group
    load(path) % load the data   

    FDTvarname = sprintf('perFDT_subjects%d', i);
    eval([FDTvarname ' = perFDT_subjects;']);

end

n_groups = floor(1:size(groups, 2));
combinations = nchoosek(n_groups, 2);

for i=1:size(combinations, 1)
    i_comb = combinations(i, :);
    group1_FDT =  eval(join(['perFDT_subjects', num2str(i_comb(1))]));
    group2_FDT =  eval(join(['perFDT_subjects', num2str(i_comb(2))]));

        
    NSUB1=size(group1_FDT,1); % number of subjects in group 1
    NSUB2=size(group2_FDT,1); % number of subjects in group 2


    % loop in all (FC, Ceff, FDT)

    group1_all={group1_FDT};
    group2_all={group2_FDT};

    for j=1:length(group1_all)
        clear DataAll1 DataAll2 scores predictedLabels positive_scores trueLabels X Y
        group1=group1_all{j};%mean(group1_all{j},2);
        group2=group2_all{j};%mean(group2_all{j},2);

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
        
        DataAll1=xxdata(1:NSUB1,:);
        DataAll2=xxdata(NSUB1+1:NSUB1+NSUB2,:);
        [pcmat(i,j,:,:),acc(i,j), svm_model, acc_all]=function_svm_subject(DataAll1, DataAll2, kfold,perc);
         
    end
   
end

