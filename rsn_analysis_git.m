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

groups={'ucla_schizophrenia_dbs80','ucla_subsetcontrols_dbs80'}; % names of groups to analyze 

dk_nets = load('DK80_nets.mat');

for i = 1:length(groups)
    load(join([pathResults, groups{i}, '_FDT_results_Tau2_nofiltfilt.mat'])) % path to the data of the group
    data = perFDT_subjects;
    g = string(groups(i));
    total_subjects = cell(1, size(data, 1));
    total_subjects_net = cell(size(data, 1), 8);
    for j = 1:size(data, 1)
        data_subject = data(j,:);
        rsn = dk_nets.DK_nets.*data_subject;
        total_subjects{j} = rsn;
        for k = 1:8
            net_data = nonzeros(rsn(k, :));
            total_subjects_net{j, k} = net_data;
        end
    end
    rsn_subjects.(g) = total_subjects;
    rsn_subjects_nets.(g) = total_subjects_net; %matrix of 4 fields (one for each group) which consist of all different rsn of each subject 
    
    savename = join([pathResults, groups{i}, '_rsn_results.mat']);
    save (savename, 'rsn_subjects', 'rsn_subjects_nets')
end


%% STATISTICAL ANALYSIS

nPerm = 10000;
rsn_names = {'Vis', 'SomMot', 'DorsAttn', 'SalVentAttn', 'Limbic', 'Cont', 'Default', 'Subcortical'};
variable = struct();

% Uncomment for the analysis by subject
% for i = 1:length(groups)
%     group = groups{i};
%     load(join([pathResults, groups{i}, '_rsn_results.mat']))
%     data = rsn_subjects_nets.(string(group));
%     for j = 1:length(rsn_names)
%         a = data(:,j);
%         a = cellfun(@(x) mean(x, 'all'), a);
%         net_data = a;
%         net = rsn_names(j);
%         temp_var = strcat(string(net), '_', string(groups(i)));
%         variable.(char(temp_var)) = net_data;  
%         disp(variable)
%     end
% 
% end

% Uncomment for the analysis by node
for i = 1:length(groups)
    group = groups{i};
    load(join([pathResults, groups{i}, '_rsn_results.mat']))
    data = rsn_subjects_nets.(string(group));
    for j = 1:length(rsn_names)
        a = data(:,j);
        a = cell2mat(a'); 
        net_data = mean(a, 2); 
        net = rsn_names(j);
        temp_var = strcat(string(net), '_', string(groups(i)));
        variable.(char(temp_var)) = net_data;  
        disp(variable)
    end

end


%create boxplots
N = 8;
pval = 0.05;
p = nan(N, 1);

figure
for y = 1:length(rsn_names)

    disp(y)
    subplot(2,4,y)

    temp_a = strcat(rsn_names(y), "_ucla_schizophrenia_dbs80");
    temp_b = strcat(rsn_names(y), "_ucla_subsetcontrols_dbs80");
    
    a=variable.(temp_a);
    b=variable.(temp_b);

    a = a';
    b = b';

    %remove outliers
    a = rmoutliers(a,"mean");
    b = rmoutliers(b,"mean");
    disp(median(a))
    disp(iqr(a))
    disp(median(b))
    disp(iqr(b))
    clear pvalue H
    stats=permutation_htest2_np([a,b],[ones(1,numel(a)) 2*ones(1,numel(b))],1000,0.05,'ranksum');
    pvalue=min(stats.pvals);
    ps(y) = pvalue;
    clear STATS p 

    C = {a, b};
    maxNumEl = max(cellfun(@numel,C));
    Cpad = cellfun(@(x){padarray(x(:),[maxNumEl-numel(x),0],NaN,'post')}, C);
    Cmat = cell2mat(Cpad);
    boxplot(Cmat,'Labels',groups)
    title (rsn_names(y)); 
    H=sigstar({[1,2]},[pvalue]);
    disp([pvalue])
    f = gcf;
end

q_values = mafdr(ps, 'BHFDR', true);

%% radar plot

a = load(join([pathResults, 'ucla_schizophrenia_dbs80_rsn_results.mat']));
b = load(join([pathResults, 'ucla_subsetcontrols_dbs80_rsn_results.mat']));
a = a.rsn_subjects_nets.ucla_schizophrenia_dbs80;
b = b.rsn_subjects_nets.ucla_subsetcontrols_dbs80;


for i = 1:8
   
    a_rsn = a(:,i);
    a_rsn = cellfun(@(x) mean(x, 'all'), a_rsn);
    disp(i)
    disp(mean(a_rsn))
    disp(std(a_rsn))
    c = sum(a_rsn)/sum(a_rsn ~= 0);
    schizophrenia_total_radar(i) = c;
    b_rsn = b(:,i);
    b_rsn = cellfun(@(x) mean(x, 'all'), b_rsn);
    disp(mean(b_rsn))
    disp(std(b_rsn))
    d = sum(b_rsn)/sum(b_rsn ~= 0);
    controls_total_radar(i) = d;

    mean1 = mean(a_rsn);
    mean2 = mean(b_rsn);
    std1 = std(a_rsn);
    std2 = std(b_rsn);
    n1 = length(a_rsn);
    n2 = length(b_rsn);
    
    % Compute pooled standard deviation
    sp = sqrt(((n1 - 1) * std1^2 + (n2 - 1) * std2^2) / (n1 + n2 - 2));
    
    % Compute Cohen's d
    d = (mean1 - mean2) / sp;
    smd(i) = d;
end

figure
P = [schizophrenia_total_radar; controls_total_radar];

% Spider plot
s = spider_plot_class(P);
% Legend properties
s.LegendLabels = {'Schizophrenia', 'Controls'};
s.AxesLabels = {'Vis', 'SomMot', 'DorsAttn', 'SalVentAttn', 'Limbic', 'Cont', 'Default', 'Subcortical'};
s.LegendHandle.Location = 'northeastoutside';
max = 70*ones(1,8);
min = 25*ones(1,8);
s.AxesLimits = [min;max];

