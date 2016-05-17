
% To fix:
% * You remake the interaction matrix a number of times. Since there you want to build an identical
% one each time, and since there are a few parameters, it's better to make a simple function for
% this. Then you know you're making it the same each time!
% * Collapse exploration (2) and intensification (3) into a single section.


%% 0. Initialize
tic
fprintf('\n    0. Initialize')

minrep = user.minrep; % minimum number of replicates an interaction has to be in

% Define folders, i.e. define where everything lives.
datadir = [user.maindir 'Data/']; % where data files live
% Make folders if necessary
if ~exist(datadir, 'dir'); mkdir(datadir); end
if ~exist([datadir '/Complexes/'], 'dir'); mkdir([datadir '/Complexes/']); end

% pairwise interactions, found in Final_Interactions_list
InteractionIn = cell(length(user.desiredPrecision),1);
countPrec = 0;
for ii = 1:length(user.desiredPrecision)
  s1 = ['Final_Interactions_list_' num2str(user.desiredPrecision(ii)*100) '_precision.csv'];
  s = [user.maindir 'Data/ROC/CombinedResults/' s1];
  
  if ~exist(s,'file')
    fprintf('\n    Error: Complexes: Following interaction file not found:')
    fprintf('\n        %s\n',s1)
  else
    countPrec = countPrec+1;
    InteractionIn{countPrec} = s;
  end
end
if countPrec == 0
  fprintf('\n    Error: Complexes: No interaction files found!')
end
InteractionIn = InteractionIn(1:countPrec);

% scoreMatrix+Protein, found in score_repx.mat
scoreIn = cell(length(user.Nreplicate),1);
countScore = 0;
for ii = 1:user.Nreplicate
  fn = ['score_rep' num2str(ii) '.mat'];
  s = [user.maindir 'Data/ROC/tmp/' fn];
  
  if ~exist(s,'file')
    fprintf('\n    Error: Complexes: Following scoreMatrix file not found:')
    fprintf('\n        %s\n',s1)
  else
    countScore = countScore+1;
    scoreIn{ii} = s;
  end
end
if countScore < user.Nreplicate
  fprintf('\n    Error: Complexes: Missing a scoreMatrix file!')
end

tt = toc;
fprintf('  ...  %.2f seconds\n',tt)




%% 1. Read input data
% Final_Interactions_list
% scoreMatrix+Protein
% corum
tic
fprintf('    1. Read input data')

% Final_Interactions_list
interactionPairs = cell(countPrec,1);
for ii = 1:countPrec
  %tmp = importdata(InteractionIn{ii});
  tmp = readFinalInteractionList(InteractionIn{ii});
  
  % Need 4 columns: Protein-A, Protein-B, replicates, precision-cutoff
  interactionPairs{ii} = tmp(:,[2:3 6 15]);
end

% scoreMatrix+Protein
Proteins = cell(user.Nreplicate,1);
for ii = 1:user.Nreplicate
  load(scoreIn{ii})
  Nprot_rep = sqrt(length(scoreMatrix));
  % sanity check
  if mod(Nprot_rep,1)~=0 || length(Protein.Isoform)~=Nprot_rep
    error('Complexes: Badly formatted scoreMatrix')
  end
  Proteins{ii} = Protein;
end
clear scoreMatrix TP_Matrix inMatrix possibleInts inverse_self Protein Dist

% corum
corumComplex = importdata(user.corumcomplexfile, ',');
if ~isvector(corumComplex)
  error('Error: Complexes: Incorrectly formatted CORUM complex file')
end


tt = toc;
fprintf('  ...  %.2f seconds\n',tt)



%% 2. Make uniqueProteins
tic
fprintf('    2. Make uniqueProteins')

uniqueProteins = [];

% Start with interactionPairs
for ii = 1:countPrec
  tmp = interactionPairs{ii}(:,1:2);
  uniqueProteins = unique([uniqueProteins; tmp(:)]);
end

% Add interactionScore (Proteins)
for ii = 1:user.Nreplicate
  uniqueProteins = unique([uniqueProteins; Proteins{ii}.Isoform]);
end

% uniqueProteins(1:Nprot_pred) are predicted proteins, the rest are CORUM.
Nprot_pred = length(uniqueProteins);

% Add corum
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

tt = toc;
fprintf('  ...  %.2f seconds\n',tt)



%% 3. Convert protein IDs to indices
% (numerical values are easier to work with!)
tic
fprintf('    3. Convert IDs to indices')

% interactionPairs --> interactionPairs2, indices of uniqueProteins
% also include how many replicates it was seen in.
interactionPairs2 = cell(countPrec,1);
for ii = 1:countPrec
  interactionPairs2{ii} = nan(size(interactionPairs{ii},1),7);
  
  for jj = 1:size(interactionPairs{ii},1)
    I1 = find(ismember(uniqueProteins,interactionPairs{ii}{jj,1}));
    I2 = find(ismember(uniqueProteins,interactionPairs{ii}{jj,2}));
    tmp = interactionPairs{ii}{jj,3};
    nrep = length(unique(tmp(isstrprop(tmp,'digit'))));
    
    tmp2 = zeros(1,3);
    for kk = 1:3
      tmp2(kk) = double(ismember(num2str(kk),interactionPairs{ii}{jj,3}));
    end
    
    precdrop = str2num(interactionPairs{ii}{jj,4});
    
    interactionPairs2{ii}(jj,:) = [I1 I2 nrep precdrop tmp2];
  end
end

% Proteins{ii}.Isoform --> Proteins2, indices of uniqueProteins
% also include how many replicates it was seen in.
Proteins2 = cell(user.Nreplicate,1);
for ii = 1:user.Nreplicate
  Proteins2{ii} = nan(size(Proteins{ii}.Isoform,1),1);
  
  for jj = 1:size(Proteins2{ii},1)
    I = find(ismember(uniqueProteins,Proteins{ii}.Isoform{jj}));
    Proteins2{ii}(jj) = I;
  end
end

% corumComplex --> corumComplex2, indices of uniqueProteins
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



%% 4. Optimize complex-building parameters, wide grid search
% For each replicate, optimize:
%   i) score cutoff
%   ii) p, the uncertainty parameter in clusterONE
%   iii) density threshold, how "loose" is too loose for clusterONE?
tic
fprintf('    4. Build complexes, exploration')

% % Legacy code.
% % This makes intMatrix from interactionPairs2.
% intMatrix = zeros(Nprot_pred,Nprot_pred);
% for jj = 1:size(interactionPairs2{ii},1)
%   x = interactionPairs2{ii}(jj,1:2);
%   nrep = interactionPairs2{ii}(jj,3);
%   intMatrix(x(1),x(2)) = intMatrix(x(1),x(2)) + nrep;  % add nrep
%   intMatrix(x(2),x(1)) = intMatrix(x(2),x(1)) + nrep;
% end
% intMatrix(intMatrix<minrep) = 0;

precRange = [0.5 0.7 0.9];
pRange = [0 100 1000];
densRange = [ 0 1 2];

ga = nan(user.Nreplicate,length(precRange),length(pRange), length(densRange));
nc = nan(user.Nreplicate,length(precRange),length(pRange), length(densRange));
avgsize = nan(user.Nreplicate,length(precRange),length(pRange), length(densRange));
mr = nan(user.Nreplicate,length(precRange),length(pRange), length(densRange));
Members2 = cell(user.Nreplicate,length(precRange),length(pRange), length(densRange));
Members3 = cell(user.Nreplicate,length(precRange),length(pRange), length(densRange));
best_params = nan(user.Nreplicate,2);
for repi = 1:user.Nreplicate
  for ii = 2:length(precRange)
    
    % i) make intMatrix from scoreMatrix, i.e. optimize xcutoff
    intMatrix = zeros(Nprot_pred,Nprot_pred);
    for jj = 1:size(interactionPairs2{1},1)
      if interactionPairs2{1}(jj,end - repi + 1)
        x = interactionPairs2{1}(jj,1:2);
        precdrop = interactionPairs2{1}(jj,4);
        intMatrix(x(1),x(2)) = precdrop;
        intMatrix(x(2),x(1)) = precdrop;
      end
    end
    intMatrix(intMatrix < precRange(ii)) = 0;
    intMatrix(intMatrix >= precRange(ii)) = 1;
    
    for jj = 1:length(pRange)
      
      for mm = 1:length(densRange)
        
        disp([num2str(repi) ', ' num2str(ii) ', ' num2str(jj) ', ' num2str(mm)])
        
        % ii) make complexes
        Members2{repi,ii,jj,mm} = myclusterone(intMatrix, pRange(jj), densRange(mm));
        if isempty(Members2{repi,ii,jj,mm}); continue; end
                
        % iv) assess how good the complexes are
        ga(repi,ii,jj,mm) = geomacc(Members2{repi,ii,jj,mm},corumComplex2);
        mr(repi,ii,jj,mm) = matchingratio(Members2{repi,ii,jj,mm},corumComplex2);
        nc(repi,ii,jj,mm) = length(Members2{repi,ii,jj,mm});
        tmp = zeros(nc(repi,ii,jj,mm),1);
        for kk = 1:nc(repi,ii,jj,mm)
          tmp(kk) = length(Members2{repi,ii,jj,mm}{kk});
        end
        avgsize(repi,ii,jj,mm) = prctile(tmp,75);
        
      end
      
    end
  end
  
  % v) choose parameters
  tmp = sq(ga(repi,:,:,:));
  [~,I] = max(tmp(:));
  [i1, i2, i3] = ind2sub(size(tmp),I);
  best_params(repi,:) = [precRange(i1) pRange(i2) densRange(i3)];
end

tt = toc;
fprintf('  ...  %.2f seconds\n',tt)



%% 5. Build final complex list

tic
fprintf('    5. Build final complex list')

Ncomplex = zeros(countPrec,1);
for ii = 1:countPrec
  % i) Make interaction matrix
  intMatrix = zeros(Nprot_pred,Nprot_pred);
  for jj = 1:size(interactionPairs2{ii},1)
    x = interactionPairs2{ii}(jj,1:2);
    nrep = interactionPairs2{ii}(jj,3);
    intMatrix(x(1),x(2)) = intMatrix(x(1),x(2)) + nrep;  % add nrep
    intMatrix(x(2),x(1)) = intMatrix(x(2),x(1)) + nrep;
  end
  intMatrix(intMatrix<minrep) = 0;
  
  %[ComplexList(ii).Members, ComplexList(ii).Connections] = myclusterone(intMatrix, best_params(ii,1), best_params(ii,2));
  [ComplexList(ii).Members, ComplexList(ii).Connections] = myclusterone(intMatrix, 1, 1);
  Ncomplex(ii) = length(ComplexList(ii).Members);
end

tt = toc;
fprintf('  ...  %.2f seconds\n',tt)



%% 6. Match each predicted complex to a CORUM complex
% for each predicted complex
%   calculate the overlap with each CORUM complex
%   report the CORUM complex with the highest overlap

tic
fprintf('    6. Match each complex to a CORUM complex')

corumMatches = cell(countPrec,1);
for ii = 1:countPrec
  
  % 5 columns: predicted complex number (Members), corum complex number (corumComplex2),
  %            match rank (e.g. best, 2nd best), N overlap, overlap ratio (N_overlap/size_of_corum)
  corumMatches{ii} = zeros(Ncomplex(ii)*5, 2);
  cc = 0;
  
  corumComplex2 = corumComplex3{ii};
  
  overlap = zeros(Ncomplex(ii),length(corumComplex2));
  for jj = 1:Ncomplex(ii)
    for kk = 1:length(corumComplex2)
      overlap(jj,kk) = length(intersect(ComplexList(ii).Members{jj},corumComplex2{kk})) / length(corumComplex2{kk});
    end
    
    Nover = sort(unique(overlap(jj,:)), 'descend');
    Nover = Nover(Nover>0);
    for mm = 1:min(2,length(Nover))
      matches = find(overlap(jj,:) == Nover(mm));
      for nn = 1:length(matches)
        cc = cc+1;
        corumMatches{ii}(cc,1) = jj;
        corumMatches{ii}(cc,2) = matches(nn);
        corumMatches{ii}(cc,3) = mm;
        corumMatches{ii}(cc,4) = length(intersect(ComplexList(ii).Members{jj},corumComplex2{matches(nn)}));
        corumMatches{ii}(cc,5) = overlap(jj,matches(nn));
      end
    end
  end
  corumMatches{ii} = corumMatches{ii}(1:cc,:);
  
end


tt = toc;
fprintf('  ...  %.2f seconds\n',tt)



% %% optimize the sparse interaction matrix for corum
%
% scoreRange = [0.9 0.95 0.99 0.99 0.999 0.9999 0.99999];
% pRange = [0 10 100 1000 10000];
%
% Members2 = cell(user.Nreplicate,length(scoreRange),length(pRange));
% ga = nan(user.Nreplicate,length(scoreRange),length(pRange));
% nc = nan(user.Nreplicate,length(scoreRange),length(pRange));
% avgsize = nan(user.Nreplicate,length(scoreRange),length(pRange));
% mr = nan(user.Nreplicate,length(scoreRange),length(pRange));
% for repi = 1:user.Nreplicate
%   for ii = 1:length(scoreRange)
%
%     % i) Make interaction matrix
%     xcutoff = scoreRange(ii);
%     I = interactionScore{repi}>xcutoff;
%     intMatrix = interactionScore{repi};
%     intMatrix(~I) = 0;
%
%     % ii) Explore parameters
%     for jj = 1:length(pRange)
%       disp([num2str(ii) ', ' num2str(jj)])
%
%       % iii) make complexes
%       Members2{repi,ii,jj} = myclusterone(intMatrix, pRange(jj), 0);
%       if isempty(Members2{repi,ii,jj}); continue; end
%
%       % iv) calculate Geometric Accuracy
%       ga(repi,ii,jj) = geomacc(Members2{repi,ii,jj},corumComplex2);
%
%       % iva) calculate Matching Ratio
%       mr(repi,ii,jj) = matchingratio(Members2{repi,ii,jj},corumComplex2);
%
%       nc(repi,ii,jj) = length(Members2{repi,ii,jj});
%       tmp = zeros(nc(repi,ii,jj),1);
%       for kk = 1:nc(repi,ii,jj)
%         tmp(kk) = length(Members2{repi,ii,jj}{kk});
%       end
%       avgsize(repi,ii,jj) = mean(tmp);
%
%     end
%   end
% end





%% x. Write output

tic
fprintf('    5. Match each complex to a CORUM complex')

writeOutput_complexes

tt = toc;
fprintf('  ...  %.2f seconds\n',tt)


