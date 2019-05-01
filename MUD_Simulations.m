% code for multivariate-univariate dependence (MUD) analysis simulations
%
% looks at the relationship between univariate activity in a simulated
% brain region, pattern similarity in that region, and the relationship
% between voxels' contributions to pattern similarity and univariate
% activity
%
% variables of interest are
% 'patternActivationCorrelation','meanOverallActivation','patternSimilarity',
% which have the univariate activity, pattern similarity, and 
% MUD (patternActivationCorrelation) for each of the n iterations
%
% for more details (and results), see Aly M & Turk-Browne NB. (2016). Attention 
% stabilizes representations in the human hippocampus. Cerebral Cortex, 26, 783-796.
%
% these are simulations only; adapt for your own needs
% code very inelegant, because 2014 Mariam was not good at writing code
%
% Mariam Aly, August 2014
%
% ==========================================================================


clear all

fileName = 'Act-MUD';

for z = 1:1000 % iterations

% -a,+p,-m (simulation 6 in paper)
% parameters for negative overall activity, positive pattern similarity, and positive MUD 
% (see Aly & Turk-Browne, 2016, for all parameters)

    patt1(1:500,1)=-0.25+(randn(500,1));
    patt2(1:500,1)=patt1(1:500,1);
    patt3(1:500,1)=patt1(1:500,1);
    patt4(1:500,1)=patt1(1:500,1);
    patt5(1:500,1)=patt1(1:500,1);
    patt6(1:500,1)=patt1(1:500,1);


    patt1(501:1000,1)=(randn(500,1));
    patt2(501:1000,1)=(randn(500,1));
    patt3(501:1000,1)=(randn(500,1));
    patt4(501:1000,1)=(randn(500,1));
    patt5(501:1000,1)=(randn(500,1));
    patt6(501:1000,1)=(randn(500,1));


% put all patterns in a single variable, where rows are voxels and columns
% are different patterns

    patt = [patt1,patt2,patt3,patt4,patt5,patt6];
    pattsim = atanh(corrcoef(patt)); % Fisher transform the correlations

% mean subtract
     for i = 1:size(patt,2) % for all columns / patterns
         patt_mean_sub(:,i) = patt(:,i) - mean([patt(:,i)]); % remove the mean of that column / pattern
     end
 
 % divide by root sum of squares
     for i = 1:size(patt_mean_sub,2) 
         patt_norm(:,i) = patt_mean_sub(:,i)/(sqrt(sum(patt_mean_sub(:,i).^2)));
     end
 
 % calculate correlations to make sure didn't make mistake
 
     for i = 1:size(patt_norm,2) 
         for j = 1:size(patt_norm,2)
             r(i,j) = sum(patt_norm(:,i).*patt_norm(:,j));

         end
     end
 
    check_answer = round(pattsim - r); % should be all 0s if no mistake made
 
% now get voxel contributions to pattern similarity
 
     for i = 1:size(patt_norm,2)
         for j = 1:size(patt_norm,2)
            patt_prod{i,j} = patt_norm(:,i).*patt_norm(:,j);
         end
     end

% get mean voxel contributions to each of the off-diagonal correlations
    voxelContribs = [patt_prod{1,(2:6)}];
    voxelContribs = [patt_prod{1,(2:6)},patt_prod{2,(3:6)},patt_prod{3,(4:6)},patt_prod{4,(5:6)},patt_prod{5,6}];
    meanVoxelContribs = mean(voxelContribs,2);

% get mean activation
    meanAct = mean(patt,2); % per voxel
    meanActivation = mean(meanAct); % across all voxels

% now calculate the correlation between activation and voxels'
% contributions to pattern similarity

    pattAct = [meanAct, meanVoxelContribs];
    pattActCorr = atanh(corrcoef(pattAct)); % Fisher transform the correlations



    patternActivationCorrelation(z) = pattActCorr(1,2);
    meanOverallActivation(z) = meanActivation;
    patternSimilarity(z) = mean([pattsim(1,(2:6)),pattsim(2,(3:6)),pattsim(3,(4:6)),pattsim(4,(5:6)),pattsim(5,6)]);


end

save(fileName, 'patternActivationCorrelation','meanOverallActivation','patternSimilarity');

