function [val_ssnd, ind_ssnd, pvals] = function_statistics_ssnd(group1,group2,num,test,threshold,permutation)
    % Aim: calculate statistics between groups based on ssnd 

    % Input
    % group1: group 1, matrix of subjects x brain node
    % group2: group 2, matrix of subjects x brain node
    % num: number of brain nodes, int
    % test: statistical test, string ('signrank', 'ranksum','montecarlo', 'ttest','ttest2',
    % does not assume normality: signrank, ranksum, montecarlo
    % asumes normality: ttest, ttest2
    % threshold: threshold for statistical test, float (e.g., 0.01)
    % permutation: amount of permutations for statistical test, int (e.g., 1000)

    % Output
    % val_ssnd: SSNR values in descending order
    % ind_ssnd: index of brain nodes in descending order of the SSNR values 
    % pvals: array of all pvals
  
    % Irene Acero & Paulina Clara Dagnino, Upf, April 2023 
    % Original code sent by Gustavo Deco

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    
    % calculate SSNR (Standardized Squared Norm Difference), standardized difference between 
    % the means of all brain nodes of groups considering their size and variance.
    for n=1:num
        a=group1(:,n)';
        b=group2(:,n)';
        SSND(n)=abs(mean(a)-mean(b))./sqrt(var(a)+var(b)); 
    end
    
    % sort SSNR in descending order values stored in sopm, indexes stored in ind
    [val_ssnd ind_ssnd]=sort(SSND,'descend'); 
    
    % for n = 1 to the total number of brain nodes, in descending order based on the SSNR between
    % groups, obtain the mean of the first n brain nodes, compute statistics between groups and 
    % save the p values in pvals
    for n=1:num
        a=mean(group1(:,ind_ssnd(1:n)),2)';
        b=mean(group2(:,ind_ssnd(1:n)),2)';
        if strcmp(test,'signrank') 
            stats=permutation_htest2_np([a,b],[ones(1,numel(a)) 2*ones(1,numel(b))],permutation,threshold,'signrank');
            pvals(n)=min(stats.pvals); 
        elseif strcmp(test,'ranksum') 
            stats=permutation_htest2_np([a,b],[ones(1,numel(a)) 2*ones(1,numel(b))],permutation,threshold,'ranksum');
            pvals(n)=min(stats.pvals);
        elseif strcmp(test,'montecarlo') 
            pvals(n)=PERM_Montecarlo_permtest(a,b,permutation);
        elseif strcmp(test,'ttest') 
            stats=permutation_htest2_np([a,b],[ones(1,numel(a)) 2*ones(1,numel(b))],permutation,threshold,'ttest');
            pvals(n)=min(stats.pvals);
        elseif strcmp(test,'ttest2') 
            stats=permutation_htest2_np([a,b],[ones(1,numel(a)) 2*ones(1,numel(b))],permutation,threshold,'ttest2');
            pvals(n)=min(stats.pvals);
        end
    end


end

