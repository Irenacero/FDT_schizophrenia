% Aim: Perform the statistical analysis of FDT deviation using a selected
% statistical test and create boxplots
% Input (obtained in FDT_computation.m): 
% - FDT deviation values (FDTm), matrix of 1 x subjects
% Output: Statistical output (p-values) and boxplots

% Irene Acero & Paulina Clara Dagnino, Upf, April 2023 
% Original code sent by Gustavo Deco



clear all;

%%%%%%%%%%%%%%%%%%%%%%%%%% Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%

perm = 10000; % number of permuatations
thres = 0.05; % threshold for the statistical test
groups={'ucla_schizophrenia_dbs80','ucla_subsetcontrols_dbs80'}; % names of groups to analyze 

test = 'ranksum';
%Parametric statistical test (assumes normality): 
%-Paired t-test (‘ttest’): dependent groups, reduces variability in the data by using each pair of observations as a single unit of analysis.
%-Unpaired t-test or two-sample t-test (‘ttest2’): independent groups. The same matlab code is used but in reality the two-sample t-test is a type of unpaired t-test in which the former assumes equal variance and the latter assumes similar but not exactly equal variances.

%Non-parametric statistical test:
%-Wilxocon sign-rank test / Paired Wilcoxon test: no normality or no equal variance, paired/matched samples (e.g., same group before and after),  
%-Wilxocon rank-sum test / Mann-Whitney U test: no normality and no equal variance, independent groups and unequal variances.

%Montecarlo is used when the assumptions of traditional parametric or non-parametric tests are not met.

system='linux';

% paths 
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
pathDependencies=join([generalPathScripts{1},'\Scripts\Dependencies']);

if strcmp(system,'linux')
    generalPathScripts= replace(generalPathScripts,'\','/');
    pathResults= replace(pathResults,'\','/');
    pathFigures= replace(pathFigures,'\','/');
    pathDependencies= replace(pathDependencies,'\','/');
end

addpath(pathResults)
addpath(pathFigures)
addpath(pathDependencies)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load the variable of interest
for i = 1:length(groups)
    group_name = groups{i}; % name of the group
    
    path = join([pathResults, group_name, '_FDT_results_Tau2_nofiltfilt.mat']); % path to the data of the group
    load(path) % load the data   

    FDTvarname = sprintf('sFDT%d', i); % temporal name of the variable
    eval([FDTvarname ' = FDTm_subjects;']);

end

n_groups = floor(1:size(groups, 2)); % number of groups to analyse
combinations = nchoosek(n_groups, 2); % all possible pairwise combinations

% perform statistical test for all pairwise analysis
for i=1:size(combinations, 1)
    i_comb = combinations(i, :);
    a = eval(join(['sFDT', num2str(i_comb(1))]));
    b = eval(join(['sFDT', num2str(i_comb(2))]));

    %remove outliers
    a = rmoutliers(a,"mean");
    b = rmoutliers(b,"mean");

    varname = sprintf('stats_pval_%d', i);
    if strcmp(test,'signrank') 
        stats=permutation_htest2_np([a,b],[ones(1,numel(a)) 2*ones(1,numel(b))],perm,thres,'signrank');
        pvalue=min(stats.pvals); 
    elseif strcmp(test,'ranksum') 
        stats=permutation_htest2_np([a,b],[ones(1,numel(a)) 2*ones(1,numel(b))],perm,thres,'ranksum');
        pvalue=min(stats.pvals);
    elseif strcmp(test,'montecarlo') 
        pvalue=PERM_Montecarlo_permtest_multiple(a,b,perm);
    elseif strcmp(test,'ttest') 
        stats=permutation_htest2_np([a,b],[ones(1,numel(a)) 2*ones(1,numel(b))],perm,thres,'ttest');
        pvalue=min(stats.pvals);
    elseif strcmp(test,'ttest2') 
        stats=permutation_htest2_np([a,b],[ones(1,numel(a)) 2*ones(1,numel(b))],perm,thres,'ttest2');
        pvalue=min(stats.pvals);
    end

%     [s] = PERM_Montecarlo_permtest_multiple(a,b,nPerm);
    eval([varname ' = pvalue;']);
    v{i} = i_comb;
    disp(join(['The p-value obtained in the pairwise comparison of ', string(groups(i_comb(1))), ' and ', string(groups(i_comb(2))), ' is ', string(pvalue)]))
end

% create a matrix with all the values of FDT of all groups
allvars = who;
pattern = 'sFDT'; 
C = allvars(startsWith(allvars, pattern));
C = C';
varData{1}=a;
varData{2}=b;


% create a matrix with all the stats of all groups
allvars = who;
pattern = 'stats_pval_'; 
st = allvars(startsWith(allvars, pattern));
st = st';
for i = 1:length(st)
    v_stats(i) = eval(st{i});
end

% create boxplots of the statistical test
figure
maxNumEl = max(cellfun(@numel,varData));
Cpad = cellfun(@(x){padarray(x(:),[maxNumEl-numel(x),0],NaN,'post')}, varData);
Cmat = cell2mat(Cpad);
boxplot(Cmat,'Labels',groups)
title ("FDT deviation"); 
H=sigstar(v,v_stats); % add stars to significantly different boxplots
f = gcf;


