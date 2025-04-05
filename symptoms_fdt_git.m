clear all

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

path = join([pathResults, 'ucla_schizophrenia_dbs80_FDT_results_tau2_nofiltfilt.mat']); 
load(path)


symptoms_data = readmatrix('scores_03.csv');
symptoms_data = symptoms_data(:,2:end);
%% CORRELATIONS

clear coefficients pvals
for symptom_index = 1:3
    symptom_data = symptoms_data(:, symptom_index);
    figure(); scatter(symptom_data,std(perFDT_subjects,0,2))
    [correlation_coefficients, p_values] = corrcoef([std(perFDT_subjects,0,2), symptom_data]);
    [a,b] = corr(std(perFDT_subjects,0,2), symptom_data, 'type', 'Spearman')
    correlation = correlation_coefficients(1:end-1, end);
    p_values = p_values(1:end-1, end);
    ps(symptom_index) = b;
    coefficients(symptom_index) = correlation;
    pvals(symptom_index) = p_values;
    figure(); scatter(symptom_data,std(perFDT_subjects,0,2));
    hold on;
    x = symptom_data;
    y = std(perFDT_subjects,0,2);
    coeffs = polyfit(x, y, 1); 
    x_fit = linspace(min(x), max(x), 100);
    y_fit = polyval(coeffs, x_fit);
    plot(x_fit, y_fit, 'r-', 'LineWidth', 2);

end

q_values = mafdr(ps, 'BHFDR', true);


%% Nodes characterizing each factor
symptom_data = symptoms_data(:, 1);
med = median(symptom_data);
group1 = perFDT_subjects(find(symptom_data > med),:);
group2 = perFDT_subjects(find(symptom_data < med),:);
group1_mean = mean(group1);
group2_mean = mean(group2);
distinct_nodes1 = group1_mean-group2_mean;
[sorted_values, sorted_indices] = sort(distinct_nodes1, 'descend');
top_distinct_nodes_1 = sorted_indices(1:10);


symptom_data = symptoms_data(:, 2);
med = median(symptom_data);
group1 = perFDT_subjects(find(symptom_data > med),:);
group2 = perFDT_subjects(find(symptom_data < med),:);
group1_mean = mean(group1);
group2_mean = mean(group2);
distinct_nodes2 = group1_mean-group2_mean;
[sorted_values, sorted_indices] = sort(distinct_nodes2, 'descend');
top_distinct_nodes_2 = sorted_indices(1:10);


symptom_data = symptoms_data(:, 3);
med = median(symptom_data);
group1 = perFDT_subjects(find(symptom_data > med),:);
group2 = perFDT_subjects(find(symptom_data < med),:);
group1_mean = mean(group1);
group2_mean = mean(group2);
distinct_nodes3 = group1_mean-group2_mean;
[sorted_values, sorted_indices] = sort(distinct_nodes3, 'descend');
top_distinct_nodes_3 = sorted_indices(1:10);

