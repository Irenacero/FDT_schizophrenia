function [pc, ac, svmmodel, acc_all, feature_importance, features_total, all_importance] = function_svm(DataAll1, DataAll2, kfold, percentage)

    % Aim: Classification between two groups of data with a linear
    % support vector machine (SVM) model using cross-validation

    % Input 
    % DataAll1: metrics from group 1, array, size xxxxx
    % DataAll2: metrics from group 2, array, size yyyy 
    % kfold: amount of permutations
    % percentage: percentage of subjects to train svm (based on minimum amount of
    % subject from both groups). Value in percentage (e.g,. 80).

    % Output
    % ac: accuracy of classification, variable, size
    % pcmat: confusion matrix of classification, variable, size 

    % Irene Acero & Paulina Clara Dagnino, Upf, April 2023 
    % Original code sent by Gustavo Deco 
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    NSUB1=size(DataAll1,1);
    NSUB2=size(DataAll2,1);
    NSUB=min(NSUB1,NSUB2); % look for minimum number of subjects between both groups
    NSUB_train=floor(percentage/100*NSUB); % calculate number of subjects for  training svm
    % the test set will be NSUB - NSUB_train

    cl=1:2; % define classes 1 and 2
    pc=zeros(2,2); % create zero matrix of size 2 x 2 for results

    %fi=zeros(80,1); 


    for nfold=1:kfold % loop kfold number of times for cross-validation

        % Obtain the train and test set for group 1
        shuffling=randperm(NSUB1); % shuffle indexes of subjects 
        Data=DataAll1; % acess data from group 1
        Data=Data(shuffling(1:NSUB),:); % access data of shuffled subject until the minimum number of subjects between groups 1 and 2 
        TrainData1=Data(1:NSUB_train,:); % split into train set
        XValidation1=Data(NSUB_train:NSUB,:); % split into test set 
        Responses1=categorical(ones(NSUB_train,1),cl); % create responses for group 1 and convert to categorical information
       

        % Obtain the train and test set for group 2
        shuffling=randperm(NSUB2); 
        Data=DataAll2;
        Data=Data(shuffling(1:NSUB),:);
        TrainData2=Data(1:NSUB_train,:);
        XValidation2=Data(NSUB_train:NSUB,:);
        Responses2=categorical(2*ones(NSUB_train,1),cl);

        % combine data from both groups
        TrainData=vertcat(TrainData1,TrainData2); % combine training data from both groups 
        Responses=vertcat(Responses1,Responses2); % combine response data from both groups 
        
        % train the SVM model
        t = templateSVM('KernelFunction', 'linear',  'KernelScale', 'auto', 'BoxConstraint', 2.5); % template for SVM model with linear kernel function (linear decision boundary for classification)
        svmmodel=fitcecoc(TrainData,Responses,'Learners',t); % train model with fitcecoc function,
        % "fitting an error-correcting output codes (ECOC) classifier", which encodes each class into a binary code used to train multiple binary SVM classifiers. %NEW

        % predict the class of the validation set from each group with the previously trained SVM model

        con=zeros(2,2); % matrix with zeros to store confusion matrices of each fold of the cross-validation
        % the confusion matrix con provides a summary of the performance of the model by counting the number
        % of true positives, false positives, true negatives, and false negatives for each class. 
        % row 1 is group 1, row 2 is group 2
        % column 1 is assignment to group 1, column 2 is assignment to group 2

        test1=predict(svmmodel,XValidation1); 
        for i=1:(NSUB-NSUB_train)
            winclass=test1(i);
            con(1,winclass)=con(1,winclass)+1;
        end
        con(1,:) = con(1,:)/(NSUB-NSUB_train)
        test2=predict(svmmodel,XValidation2);
        for i=1:(NSUB-NSUB_train)
            winclass=test2(i);
            con(2,winclass)=con(2,winclass)+1;
        end
        con(2,:) = con(2,:)/(NSUB-NSUB_train)
        
        pc=pc+con; % add the confusion matrix in matrix pc
        accdist(nfold,:)=diag(con); % add the diagonal elements of the confusion matrix to calculate the accuracy of each fold and save it in accdist 
        acc_all(nfold)=sum(diag(con))/2;

        % feature importance
%         coefficients = svmmodel.BinaryLearners{1}.Beta;
%         features = [1:length(coefficients)];
%         feature_importances_table = table(features', abs(coefficients), 'VariableNames', {'Feature', 'Importance'});
%         feature_importances_table_s = sortrows(feature_importances_table, 'Importance', 'descend');
%         top_n = 10;
%         features = feature_importances_table_s(1:top_n, :).Feature;
%         coeffs = feature_importances_table_s(1:top_n, :).Importance;
%         quotient = floor(features / 80);
%         remainder = mod(features, 80);
%         result = features - (80 * quotient);
%         result(remainder == 0) = quotient(remainder == 0);
%         feature_importance{nfold} = result;
%         features_total{nfold} = features;
%         all_importance{nfold} = feature_importances_table;


    end

    pc=pc/kfold; % obtain the average confusion matrix pc
    ac=sum(diag(pc))/2; % add the diagional elements of the averaged vconfusion matrix to calculate the accuracy and divide it by 2 becausevthere are 2 classes

end