clc; clear all; close all

%Values from Excel file (file import could be added, but is not currently done)
weights = [0.37, 0.16, 0.16, 0.11, 0.07, 0.07, 0.05, 0.02];
var = diag([0.136, 0.052, 0.033, 0.065, 0.033, 0.036, 0.027, 0.011]);
eval = [4.2,4.2,3.1,1.9,2.8,1.4,
    3.1,3.1,2.3,1.0,4.0,3.1;
    2.9,2.9,1.0,1.5,4.1,2.8;
    2.1,1.9,2.1,2.2,4.4,4.7;
    1.6,1.7,2.8,4.3,4.5,5.5;
    1.1,1.3,2.6,3.5,4.1,5.3;
    1.2,1.9,3.0,3.9,4.3,5.5;
    2.0,2.1,2.0,1.8,4.2,3.8];

%baseline evaluation
blscore = weights*eval;
[~,blrank] = sort(blscore);

%Memory allocationç
scoreMem = zeros(2*length(weights)+1,size(eval,2));
rankMem = zeros(2*length(weights)+1,size(eval,2));
ranktemp = zeros(1,size(eval,2));

%initialising with baseline
scoreMem(1,:) =  blscore;
[~,rankMem(1,:)] = sort(blscore);
indx = 2;

for i = 1:length(weights)
    %upsetting scores with +var at weight i
    weightsvar = weights+var(i,:);
    %storing score and rank
    scoreMem(indx,:) =  weightsvar*eval;
    [~,rankMem(indx,:)] = sort(scoreMem(indx,:));
    indx = indx+1;

    %upsetting scores with -var at weight i
    weightsvar = weights-var(i,:);
    %storing score and rank
    scoreMem(indx,:) =  weightsvar*eval;
    [~,sortindx] = sort(scoreMem(indx,:));
    ranktemp(sortindx) = 1:length(sortindx);
    rankMem(indx,:) = ranktemp;
    indx = indx+1;
end

%Average score 
meanScore = mean(scoreMem);
stdScore = std(scoreMem);
[~,sortindx] = sort(meanScore);
ranktemp(sortindx) = 1:length(sortindx);
rankwmeanscore = ranktemp;
%Average rank
meanRank = mean(rankMem);
stdRank = std(rankMem);
[~,sortindx] = sort(meanRank);
ranktemp(sortindx) = 1:length(sortindx);
rankwmeanrank = ranktemp;


%% Printing out the results
fprintf('The baseline results are:\n')
T = array2table([blscore;blrank],'VariableNames',{'O/L','O/L+S/M','NSO+SL','SO+SL','EM','ML'},'RowName',{'Score','Rank'}); 
% Display table
disp(T) 

fprintf('The results from the mean scores are:\n')
T = array2table([meanScore;rankwmeanscore;stdScore],'VariableNames',{'O/L','O/L+S/M','NSO+SL','SO+SL','EM','ML'},'RowName',{'Score','Rank','std'}); 
% Display table
disp(T) 

fprintf('The results from the mean ranks are:\n')
T = array2table([meanRank;rankwmeanrank;stdRank],'VariableNames',{'O/L','O/L+S/M','NSO+SL','SO+SL','EM','ML'},'RowName',{'Score','Rank','std'}); 
% Display table
disp(T) 

%% Plotting out the results
figure
hold on
barid = bar(0:5,rankwmeanscore);
set(gca,'XTickLabel',{'O/L','O/L+S/M','NSO+SL','SO+SL','EM','ML'});
er = errorbar(0:5,rankwmeanscore,stdRank,stdRank);  
er.Color = [0 0 0];                            
er.LineStyle = 'none'; 

