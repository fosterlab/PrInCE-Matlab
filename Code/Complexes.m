


%% 0. Initialize
tic
fprintf('\n    0. Initialize')

minrep = user.minrep; % minimum number of replicates an interaction has to be in

% Define folders, i.e. define where everything lives.
datadir = [user.maindir 'Data/']; % where data files live
figdir = [user.maindir 'Figures/Enrichment/']; % where figures live
% Make folders if necessary
if ~exist(datadir, 'dir'); mkdir(datadir); end
if ~exist([datadir '/Enrichment/'], 'dir'); mkdir([datadir '/Enrichment/']); end
if ~exist(figdir, 'dir'); mkdir(figdir); end

InteractionIn = cell(length(user.desiredPrecision),1);
countPrec = 0;
for ii = 1:length(user.desiredPrecision)
  s1 = ['Final_Interactions_list_' num2str(user.desiredPrecision(ii)*100) '_precision.csv'];
  s = [user.maindir 'Data/ROC/CombinedResults/' s1];
  
  if ~exist(s,'file')
    fprintf('\n    Error: Enrichment: Following interaction file not found:')
    fprintf('\n        %s\n',s1)
  else
    countPrec = countPrec+1;
    InteractionIn{countPrec} = s;
  end
end
if countPrec == 0
  fprintf('\n    Error: Enrichment: No interaction files found!')
end
InteractionIn = InteractionIn(1:countPrec);


tt = toc;
fprintf('  ...  %.2f seconds\n',tt)




%% 1. Read input data
tic
fprintf('    1. Read input data')

% i) Predicted pairwise nteractions

% Load interactions detected by ROC_PCPSILAC
interactionPairs = cell(countPrec,1);
for ii = 1:countPrec
  tmp = importdata(InteractionIn{ii});
  
  % remove header from textdata if necessary
  if size(tmp.textdata,1)==size(tmp.data,1)+1
    tmp.textdata = tmp.textdata(2:end,:);
  elseif size(tmp.textdata,1)==size(tmp.data,1)
  else
    disp('Error: Complexes: Mismatch between imported numerical and text data sizes')
  end
  
  % Need 4 columns: Protein-A, Protein-B, replicates, scores
  % Hard code these to textdata(:,[2 3 6 13])
  interactionPairs{ii} = tmp.textdata(:,[2:3 6 13]);
end

% List of N unique protein names, made from predicted interactions
uniqueProteins = [];
for ii = 1:countPrec
  tmp = interactionPairs{ii}(:,1:2);
  uniqueProteins = unique([uniqueProteins; tmp(:)]);
end
Nprot_pred = length(uniqueProteins);

% Turn interactionPairs into indices of uniqueProteins
% (numerical pairs are easier to work with!)
% In interaction index pairs, also include how many replicates it was seen in.
interactionPairs2 = cell(countPrec,1);
for ii = 1:countPrec
  interactionPairs2{ii} = nan(size(interactionPairs{ii},1),3);

  for jj = 1:size(interactionPairs{ii},1)
    I1 = find(ismember(uniqueProteins,interactionPairs{ii}{jj,1}));
    I2 = find(ismember(uniqueProteins,interactionPairs{ii}{jj,2}));
    tmp = interactionPairs{ii}{jj,3};
    nrep = length(unique(tmp(isstrprop(tmp,'digit'))));
        
    interactionPairs2{ii}(jj,:) = [I1 I2 nrep];
  end
end


% ii) Reference (CORUM) pairwise nteractions

corumPairs_tmp = importdata(user.corumpairwisefile);
corumPairs_tmp = unique(corumPairs_tmp);

% Check that corumPairs is an mx1 cell
if ~iscell(corumPairs_tmp) || size(corumPairs_tmp,2)~=1
  error('Error: Complexes: Incorrectly formatted CORUM pairwise file')
end

% split corumPairs into a mx2 cell
% this method is surprisingly faster than cellfun(@strsplit)
corumPairs = cell(length(corumPairs_tmp),2);
for ii = 1:length(corumPairs_tmp)
  fn = corumPairs_tmp{ii};
  for jj = 1:length(fn)
    if fn(jj) == ','
      corumPairs{ii,1} = fn(1:jj-1);
      corumPairs{ii,2} = fn(jj+1:end);
    end
  end
end
clear corumPairs_tmp

% augment uniqueProteins with ones from CORUM pairwise
% Goal: uniqueProteins(1:Nprot) are predicted proteins, the rest are CORUM.
% this is a little awkward. can likely be improved.
for ii = 1:size(corumPairs,1)
  for jj = 1:2
    if ~ismember(corumPairs{ii,jj},uniqueProteins)
      I = length(uniqueProteins);
      uniqueProteins{I+1} = corumPairs{ii,jj};
    end
  end
end

% Turn corumPairs into indices of uniqueProteins
corumPairs2 = nan(size(corumPairs));
for ii = 1:size(corumPairs2,1)
  I1 = find(ismember(uniqueProteins, corumPairs{ii,1}));
  I2 = find(ismember(uniqueProteins, corumPairs{ii,2}));
  corumPairs2(ii,:) = [I1 I2];
end

% corumPairs has both A-B and B-A interactions.
% Remove the redundant entries.
corumPairs2 = sort(corumPairs2,2);
[corumPairs2,I1,I2] = unique(corumPairs2,'rows');
% Use the same indices to remove redundant entries from the string list.
corumPairs = corumPairs(I1,:);


% iii) Reference (CORUM) complexes

% Load CORUM complexes
corumComplex = importdata(user.corumcomplexfile, ',');
if ~isvector(corumComplex)
  error('Error: Complexes: Incorrectly formatted CORUM complex file')
end

% augment uniqueProteins with ones from CORUM complex
% Goal: uniqueProteins(1:Nprot) are predicted proteins, the rest are CORUM.
% this is a little awkward. can likely be improved.
for ii = 1:size(corumComplex,1)
  cmplx = corumComplex{ii};

  Idelim = [0 strfind(cmplx, ',') length(cmplx)+1];
  Nprot = length(Idelim) - 1;
  if Nprot<2
    error('Error: Complexes: Complex of size 1 detected.')
  end
  for jj = 1:Nprot
    prot1 = cmplx(Idelim(jj)+1 : Idelim(jj+1)-1);
  
    if ~ismember(prot1,uniqueProteins)
      I = length(uniqueProteins);
      uniqueProteins{I+1} = prot1;
    end
  end
end

% Turn corumComplex into indices of uniqueProteins
% Reject complexes that have zero overlap with predicted-interaction proteins.
corumComplex2 = cell(size(corumComplex));
for ii = 1:length(corumComplex)
  cmplx = corumComplex{ii};
  
  Idelim = [0 strfind(cmplx, ',') length(cmplx)+1];
  Nprot = length(Idelim) - 1;
  
  I = nan(1, Nprot);
  for jj = 1:Nprot
    prot1 = cmplx(Idelim(jj)+1 : Idelim(jj+1)-1);
    I(jj) = find(ismember(uniqueProteins, prot1));
  end
  
  if sum(isempty(I))>0 || sum(isnan(I))>0
    error('Error: Complexes: CORUM complex protein not found in pairwise proteins.')
  end
  
  % Check that complex includes at least one predicted-interaction protein.
  if sum(I<=Nprot_pred)>0
    corumComplex2{ii} = I;
  end
end
corumComplex2(cellfun('isempty',corumComplex2)) = [];


tt = toc;
fprintf('  ...  %.2f seconds\n',tt)




%% 2. Optimize complex-building parameters, wide grid search
%   i) p, the uncertainty parameter in clusterONE
%   ii) density_threshold, how loose can a complex be? (penalizes large complexes in practice)
tic
fprintf('    2. Build complexes')

pRange = [0 1 10 100 1000 10000 100000];
densRange = linspace(0,5,6);

clear ComplexList
for ii = 1:1%countPrec
  
  % i) Make interaction matrix
  intMatrix = zeros(Nprot_pred,Nprot_pred);
  for jj = 1:size(interactionPairs2{ii},1)
    x = interactionPairs2{ii}(jj,1:2);
    nrep = interactionPairs2{ii}(jj,3);
    intMatrix(x(1),x(2)) = intMatrix(x(1),x(2)) + nrep;  % add nrep
    intMatrix(x(2),x(1)) = intMatrix(x(2),x(1)) + nrep;
  end
  intMatrix(intMatrix<minrep) = 0;
  
  % ii) Explore parameters
  %mr = nan(length(pRange),length(densRange));
  ga = nan(length(pRange),length(densRange));
  nc = nan(length(pRange),length(densRange));
  Members = cell(length(pRange),length(densRange));
  for jj = 1:length(pRange)
    params.p = pRange(jj);
    for kk = 1:length(densRange)
      params.density_threshold = densRange(kk);
      
      disp([num2str(jj) ', ' num2str(kk)])
      
      % iii) make complexes
      Members{jj,kk} = myclusterone(intMatrix,params);
      nc(jj,kk) = length(Members{jj,kk});
      
      % iva) calculate Matching Ratio
      % mr(jj,kk) = matchingratio(Members{jj,kk},corumComplex2);
      
      % ivb) calculate Geometric Accuracy
      ga(jj,kk) = geomacc(Members{jj,kk},corumComplex2);
      
    end
  end
  
  % v) choose parameters
end


tt = toc;
fprintf('  ...  %.2f seconds\n',tt)


%% 3. Optimize complex-building parameters, local hill-climbing

for ii = 1:1%countPrec

  
end


%% 4. Build final complex list

[ComplexList(ii).Members, ComplexList(ii).Connections] = myclusterone(intMatrix2);
Ncomplex(ii) = length(ComplexList(ii).Members);
