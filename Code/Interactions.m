%INTERACTIONS Predict protein-protein interactions using co-fractionation.
%   INTERACTIONS classifies every pair of proteins as either interacting or 
%   non-interacting. Within a replicate and condition, INTERACTIONS
%   calculates the co-fractionation dissimilarity of every pair of proteins.
%   Five 'distance measures' quantify how dissimilar co-fractionation
%   profiles are. Distance measures and gold-standard reference
%   interactions are used by the classifier as features and training
%   labels, respectively.
%
%   Master script PRINCE and module GAUSSBUILD must be run before 
%   INTERACTIONS.
%
%   INTERACTIONS produces two output folders: Output/Data/Interactions,
%   which contains all output csv tables, and Output/Figures/Interactions,
%   which contains all figures.
%
%   Workflow:
%   for every replicate:
%       1. Read input
%       2. Calculate distance measures between every pair of
%          co-fractionation profiles
%       3. Collate
%       4. Incorporate reference information
%       5. Calculate interaction probability
%   6. Find probability threshold that gives desired precision
%   7. Find and concatenate interactions at desired precision
%   8. Create a list of all the unique interactions
%   9. Calculate final TP, FP, FN and TN
%   10. Find treatment-specific interactions
%   11. Write output tables
%   12. Make figures
%
%   See also PRINCE, GAUSSBUILD, SCORENB, CALCPPITHRESHOLD.


diary([user.maindir 'logfile.txt'])
disp('Interactions.m')


%% 0. Initialize
tic
fprintf('\n    0. Initialize')

% Load user settings
maindir = user.maindir;
Experimental_channels = user.silacratios;
desiredPrecision = user.desiredPrecision;
number_of_channels = length(user.silacratios);

% Define folders, i.e. define where everything lives.
codedir = [maindir 'Code/']; % where this script lives
funcdir = [maindir 'Code/Functions/']; % where small pieces of code live
datadir = [maindir 'Output/Data/Interactions/']; % where data files live
figdir = [maindir 'Output/Figures/Interactions/']; % where figures live
% Make folders if necessary
if ~exist([maindir '/Output'], 'dir'); mkdir([maindir '/Output']); end
if ~exist([maindir '/Output/Data'], 'dir'); mkdir([maindir '/Output/Data']); end
if ~exist([maindir '/Output/Figures'], 'dir'); mkdir([maindir '/Output/Figures']); end
if ~exist([maindir '/Output/tmp'], 'dir'); mkdir([maindir '/Output/tmp']); end
if ~exist(figdir, 'dir'); mkdir(figdir); end
if ~exist(datadir, 'dir'); mkdir(datadir); end

dd = dir([maindir '/Output/tmp/Adjusted_*_Combined_OutputGaus.csv']);

% Define Raw SILAC ratios data. Dynamically find filenames
if user.skipalignment==1 || isempty(dd)
  % If Alignment was skipped, use raw data + Gauss_Build output
  
  dd = dir([maindir 'Output/tmp/*_Raw_data_maxquant_rep*.csv']);
  isgood = zeros(size(dd));
  for ii = 1:length(dd)
    for jj = 1:length(user.silacratios)
      isgood(ii) = ~isempty(strfind(dd(ii).name,user.silacratios{jj})) | isgood(ii)==1;
    end
  end
  dd = dd(isgood==1);
  ChromatogramIn = cell(size(dd));
  for di = 1:length(dd)
    ChromatogramIn{di} = [maindir 'Output/tmp/' dd(di).name];
  end
  
  dd = dir([maindir 'Output/tmp/*Combined_OutputGaus_rep*csv']);
  isgood = zeros(size(dd));
  for ii = 1:length(dd)
    for jj = 1:length(user.silacratios)
      isgood(ii) = ~isempty(strfind(dd(ii).name,user.silacratios{jj})) | isgood(ii)==1;
    end
  end
  dd = dd(isgood==1);
  GaussIn = cell(size(dd));
  for di = 1:length(dd)
    GaussIn{di} = [maindir 'Output/tmp/' dd(di).name];
  end
else
  % If Alignment was not skipped, use Alignment output
  
  dd = dir([maindir 'Output/tmp/Adjusted*Raw_for_ROC_analysis*rep*csv']);
  isgood = zeros(size(dd));
  for ii = 1:length(dd)
    for jj = 1:length(user.silacratios)
      isgood(ii) = ~isempty(strfind(dd(ii).name,user.silacratios{jj})) | isgood(ii)==1;
    end
  end
  dd = dd(isgood==1);
  ChromatogramIn = cell(size(dd));
  for di = 1:length(dd)
    ChromatogramIn{di} = [maindir 'Output/tmp/' dd(di).name];
  end
  
  GaussIn = cell(size(dd));
  dd = dir([maindir 'Output/tmp/Adjusted_Combined_OutputGaus*rep*csv']);
  isgood = zeros(size(dd));
  for ii = 1:length(dd)
    for jj = 1:length(user.silacratios)
      isgood(ii) = ~isempty(strfind(dd(ii).name,user.silacratios{jj})) | isgood(ii)==1;
    end
  end
  dd = dd(isgood==1);
  for di = 1:length(dd)
    GaussIn{di} = [maindir 'Output/tmp/' dd(di).name];
  end
end

number_of_replicates = length(ChromatogramIn) / number_of_channels;

% Map each rep*channel to a channel
rep2channel = zeros(1,length(GaussIn));
for ii = 1:length(GaussIn)
  ratioInFilename = zeros(length(user.silacratios),1);
  for jj = 1:length(user.silacratios)
    if ~isempty(strfind(GaussIn{ii},user.silacratios{jj}))
      ratioInFilename(jj) = 1;
    end
  end
  
  if sum(ratioInFilename)~=1
    error('Gaussian fitting filenames are badly formatted')
  end
  
  rep2channel(ii) = find(ratioInFilename);
end

% pre-allocate final results
Recall = nan(number_of_replicates,number_of_channels);
Precision = nan((size(Recall)));
TPR = nan((size(Recall)));
FPR = nan((size(Recall)));
No_Int = nan((size(Recall)));
notes = cell(0,1);

tt = toc;
fprintf('  ...  %.2f seconds\n',tt)


%% 
interaction_count = 0;
binary_interaction_list = cell(10^6,15);

for channel_counter = 1:number_of_channels
  s = ['\n    Channel ' num2str(user.silacratios{channel_counter})];
  fprintf(s)
  
  replicatesThisChannel = find(rep2channel == channel_counter);
  for replicate_counter = replicatesThisChannel
    
    s = ['\n        Replicate ' num2str(replicate_counter)];
    fprintf(s)
    
    %% 1. Read input data
    tic
    s = '\n        1. Read input data';
    fprintf(s)
    
    % Corum binary interactions
    fid=fopen(user.corumpairwisefile, 'rt');    %Input corum data base.
    Corum_Import= textscan(fid,'%s\t', 'Delimiter',',');
    fclose(fid);
    No=length(Corum_Import{1})/2;
    Corum_Protein_names = (reshape(Corum_Import{1,1},2,No)');
    Unique_Corum = unique(Corum_Protein_names);
    
    
    % Raw SILAC ratios
    tmp = importdata(ChromatogramIn{replicate_counter});
    if isfield(tmp.data,'Sheet1')
      Chromatograms_raw = tmp.data.Sheet1;
    else
      Chromatograms_raw = tmp.data;
    end
    % if the first column is Replicate, remove it
    tmp = length(unique(Chromatograms_raw(:,1)));
    if tmp < user.Nreplicate*2+1 % i.e. the column has unique elements in the ballpark of user.Nreplicate
      Chromatograms_raw = Chromatograms_raw(:,2:end);
    end
    
    % Gaussian parameters
    tmp = importdata(GaussIn{replicate_counter}, ',');
    if isfield(tmp.data,'Sheet1')
      tmp1 = tmp.data.Sheet1;
      tmp2 = tmp.textdata.Sheet1;
    else
      tmp1 = tmp.data;
      tmp2 = tmp.textdata;
    end
    % Ensure tmp2 is a single column of protein names
    if size(tmp2,2)~=1
      disp('Comparison: Error: violated assumption, textdata is not single column of protein names')
    end
    tmp2 = tmp2(:,1);
    % Ensure protein names are formatted correctly
    for ii = 2:size(tmp2,1)
      tmp2{ii} = strrep(tmp2{ii},',','');
    end
    % Correct bug found for Nick's tissue data.
    % Isoform number (should be textdata) becomes the first column of data.
    if size(tmp1,2)==7
      I = find(~isnan(tmp1(:,7)))';
      if ~isempty(I)
        for jj = I
          tmp2{jj+1} = [tmp2{jj+1} num2str(tmp1(jj,1))];
          tmp1(jj,1:6) = tmp1(jj,2:end);
          tmp1(jj,end) = nan;
        end
      end
    end
    tmp1 = tmp1(:,1:6);
    Gaus_import = [tmp2 num2cell( cat(1,zeros(1,size(tmp1,2)),tmp1) )];
    
    % do a bit of housekeeping
    clear Maxquant_raw Raw_data_reshaped Corum_Import tmp
    
    tt = toc;
    fprintf('  ...  %.2f seconds\n',tt)
    
    
    %% 2. Pre-process, make Protein.
    % - In chromatograms, replaces nans with 0.05
    % - Remove gaussians and chromatograms with C<5
    % - Make normalized gaussians, i.e. set height to 1
    tic
    fprintf('        2. Pre-process, make Protein.')
    
    % Pre-process data to remove Gaussian above fraction five
    H_raw = cell2mat(Gaus_import(2:end,2));
    C_raw = cell2mat(Gaus_import(2:end,3));
    W_raw = cell2mat(Gaus_import(2:end,4));
    
    % The data is nan-padded. Find where the real data starts and stops.
    tmp = find(sum(isnan(Chromatograms_raw))==size(Chromatograms_raw,1));
    if isempty(tmp)
      tmp = -1;
    end
    frac1 = max([1 tmp(find(tmp<user.Nfraction/2,1,'last'))]); % start of real data
    frac2 = min([size(Chromatograms_raw,2) tmp(find(tmp>user.Nfraction/2,1,'first'))]); % end of real data
    
    % Replace NaN values with 0.05
    Chromatograms = Chromatograms_raw;
    I = find(isnan(Chromatograms));
    Chromatograms(I)= 0.05 * rand(size(I));
    
    % remove gaussian with centers below frac1
    Ibad = zeros(size(C_raw));
    H = H_raw(~Ibad);
    C = C_raw(~Ibad);
    W = W_raw(~Ibad);
    Chromatograms = Chromatograms(~Ibad,:);
    Chromatograms_raw = Chromatograms_raw(~Ibad,:);
    
    % How many Gaussians?
    I = find(~Ibad);
    unique_names = Gaus_import(I+1,1);
    X = repmat(unique_names,1,length(unique_names));
    Ngauss = sum(strcmp(X,X'),1)';
    clear X unique_names
    
    % Make Protein structure
    % This summarizes the proteins in our sample
    
    clear Protein
    
    % Combine Majority Protein ID with Proteins Identified in replicate
    Protein.Isoform=Gaus_import(find(~Ibad)+1,1);
    Protein.Height = H;
    Protein.Width = W;
    Protein.Center = C;
    
    % Read Major protein groups file
    % If it doesn't exist, set it equal to the protein IDs
    if ~isempty(user.majorproteingroupsfile)
      Protein_IDs = readproteingroupsfile(user);
    else
      Protein_IDs = Protein.Isoform;
    end
    nn = min(size(Protein_IDs,2),23);
    Protein_IDs = Protein_IDs(:,1:nn);
    
    Dimensions_Gaus_import = length(Protein.Isoform);
    Dimension_of_Protein_IDs = size(Protein_IDs);
    
    % Create Protein.NoIsoform by removing the isoform tag from Protein.Isoform
    for Gaus_import_counter1 = 1:Dimensions_Gaus_import(1)
      Find_Gaus_isoforms=strfind(Protein.Isoform(Gaus_import_counter1), '-');
      if ~isempty(Find_Gaus_isoforms{1}) == 1
        length_of_Find_Gaus_Protein=length(Protein.Isoform{Gaus_import_counter1});
        number_of_characters_to_remove2=Find_Gaus_isoforms{1}-1;
        Protein.NoIsoform(Gaus_import_counter1,1)={Protein.Isoform{Gaus_import_counter1}(1:number_of_characters_to_remove2)};
      elseif isempty(Find_Gaus_isoforms{1}) == 1
        Protein.NoIsoform(Gaus_import_counter1,1)=Protein.Isoform(Gaus_import_counter1);
      end
    end
    
    % Create Protein_IDS_no_isoform by removing the isoform tag from Protein_IDs
    Protein_IDS_no_isoform = cell(size(Protein_IDs));
    for Protein_ID_counter1 = 2:Dimension_of_Protein_IDs(1)
      for Protein_ID_counter2 = 1:Dimension_of_Protein_IDs(2)-sum(cellfun('isempty',Protein_IDs(Protein_ID_counter1,:)))
        Find_isoforms=strfind(Protein_IDs(Protein_ID_counter1,Protein_ID_counter2), '-');
        if ~isempty(Find_isoforms{1}) == 1
          length_of_Protein_ID=length(Protein_IDs{Protein_ID_counter1,Protein_ID_counter2});
          number_of_characters_to_remove=Find_isoforms{1}-1;
          Protein_IDS_no_isoform(Protein_ID_counter1,Protein_ID_counter2)=cellstr(Protein_IDs{Protein_ID_counter1,Protein_ID_counter2}(1:number_of_characters_to_remove));
        elseif isempty(Find_isoforms{1}) == 1
          Protein_IDS_no_isoform(Protein_ID_counter1,Protein_ID_counter2)=Protein_IDs(Protein_ID_counter1,Protein_ID_counter2);
        end
      end
    end
    
    % Remove non-unique entries from each row of Protein_IDS_no_isoform
    Protein_IDS_no_isoform_no_dup_processed = cell(size(Protein_IDs));
    for ii = 2:Dimension_of_Protein_IDs(1)
      tmp = Protein_IDS_no_isoform(ii,:);
      I = ~cellfun('isempty',tmp);
      uniqueCells = unique(tmp(I(:)));
      Protein_IDS_no_isoform_no_dup_processed(ii,1:length(uniqueCells)) = uniqueCells;
    end
    
    %Look up Major Protein group Information
    % For each protein that was fitted by a Gaussian
    %   Find the rows containing that protein name in Protein_IDs
    %   Store all protein names from that row of Protein_IDs
    Protein.MajorIDs = cell(Dimensions_Gaus_import,Dimension_of_Protein_IDs(2));
    Protein.MajorID_NoIsoforms = cell(Dimensions_Gaus_import,Dimension_of_Protein_IDs(2));
    for ii = 1:Dimensions_Gaus_import(1)
      
      %disp(Lookup_protein_name)
      protName = Protein.Isoform(ii);
      
      [ia, ib] = find(strcmp(protName,Protein_IDs));
      
      %copy names to MajorID from Protein_IDs raw data
      tmp = Protein_IDs(ia,:);
      I = ~cellfun('isempty',tmp);
      tmp = unique(tmp(I(:)));
      Protein.MajorIDs(ii,1:length(tmp)) = tmp;
      
      %copy names to MajorID_NoIsoforms from Protein_IDS_no_isoform_no_dup_processed
      tmp = Protein_IDS_no_isoform_no_dup_processed(ia,:);
      I = ~cellfun('isempty',tmp);
      tmp = unique(tmp(I(:)));
      Protein.MajorID_NoIsoforms(ii,1:length(tmp)) = tmp;
    end
        
    tt = toc;
    fprintf('  ...  %.2f seconds\n',tt)
    
    
    
    %% 3. Make chromatogram-chromatogram distance matrices, Dist
    tic
    fprintf('        3. Reduce to protein level, make Dist.')
    
    % Co-Apex 1: minimum distance between Gaussian parameter doublets (C, W)
    gaussParamDist = squareform(pdist([C W],'euclidean')); % distance between every (C,W) pair
    chromNames = Protein.Isoform; % save the protein ID of each Gaussian
    
    % Reduce variables from gaussian-level to protein-level
    [~,Ireduce] = unique(Protein.Isoform);
    Chromatograms = Chromatograms(Ireduce,:);
    Chromatograms_raw = Chromatograms_raw(Ireduce,:);
    C = C(Ireduce,:);
    W = W(Ireduce,:);
    Ngauss = Ngauss(Ireduce);
    fn = fieldnames(Protein);
    for ii = 1:length(fn)
      Protein.(fn{ii}) = Protein.(fn{ii})(Ireduce,:);
    end
    Dimensions_Gaus_import = length(Ireduce);
    
    % Area under the chromatogram
    auc = sum(Chromatograms,2);
    
    % Calculate distance matrices
    clear Dist
    Dist.Euc = squareform(pdist(Chromatograms,'euclidean'));              % 1. Euclidean distance
    Dist.R = -1*corr(Chromatograms');                                     % 2. Cleaned chromatogram R^2
    [R,p] = corrcoef(Chromatograms_raw','rows','pairwise');
    Dist.Rraw = 1 - R;                                                    % 3. Raw chromatogram R^2
    Dist.Rpraw = p;                                                       % 4. Raw chromatogram correlation p-value
    Dist.CoApex = zeros(length(Protein.Isoform),length(Protein.Isoform));
    I = cell(size(Protein.Isoform));
    for ii = 1:length(Protein.Isoform)
      I{ii} = ismember(chromNames,Protein.Isoform{ii});
    end
    for ii = 1:length(Protein.Isoform)
      for jj = 1:length(Protein.Isoform)
        if ii == jj; continue; end
        tmp = gaussParamDist(I{ii},I{jj});
        Dist.CoApex(ii,jj) = min(tmp(:));                                 % 5. Co-Apex score 1
      end
    end
    [~,mx] = max(Chromatograms,[],2);
    Dist.CoApex2 = squareform(pdist(mx,'euclidean'));                     % 6. Co-Apex score 2
    fn = fieldnames(Dist);
    NDistFields = length(fn);
    
    clear gaussParamDist
    
    tt = toc;
    fprintf('  ...  %.2f seconds\n',tt)
    
    
    
    %% 4. Make TP, FP matrices
    tic
    fprintf('        4. Make TP, FP, FN matrices')
    
    % Make Pos_Corum_proteins, which is all individual proteins that are in corum and our dataset
    try
      Unique_MajorProteinIDs = unique(vertcat(Protein.MajorID_NoIsoforms{:}),'rows');
    catch
      tmp = Protein.MajorID_NoIsoforms(:);
      I = cellfun('isempty',tmp);
      tmp(I) = [];
      Unique_MajorProteinIDs = unique(vertcat(tmp));    
    end
    Pos_Corum_proteins = Unique_Corum(ismember(Unique_Corum, Unique_MajorProteinIDs));
    
    % Make Corum_Dataset, which is all interactions that are in corum and our dataset
    cc = 0;
    Corum_Dataset = cell(1000,1);
    for jj = 1:length(Corum_Protein_names); %Write out values as i
      Prot1 = Corum_Protein_names(jj,1);
      Prot2 = Corum_Protein_names(jj,2);
      
      % find these protein names in our sample
      tmp1 = sum(strcmp(Prot1,Pos_Corum_proteins));
      tmp2 = sum(strcmp(Prot2,Pos_Corum_proteins));
      
      if tmp1>0 && tmp2>0
        cc = cc+1;
        Corum_Dataset{cc} = [Prot1{1},'-',Prot2{1}];
      end
    end
    
    % Expand Corum_Dataset to deal with protein groups
    % This all has to do with Major_protein_groups.xlsx.
    % In effect, some proteins have multiple names.
    % e.g. Interaction A-B is not in corum, but A-C is, and B and C are in the same protein group, so
    %   treat A-B as a true positive interaction.
    
    TP_Matrix= zeros(Dimensions_Gaus_import,Dimensions_Gaus_import);
    %Corum_Dataset_expanded = cell(size(Corum_Dataset));
    ii = 0;
    for cc=1:size(Corum_Dataset,1)
      I_hyphen = strfind(Corum_Dataset{cc},'-');
      Prot1 = Corum_Dataset{cc}(1:I_hyphen-1);
      Prot2 = Corum_Dataset{cc}(I_hyphen+1:end);
      
      % find where Prot1 and Prot2 are in Protein.MajorID_NoIsoforms
      ind1 = find(strcmp(Prot1,Protein.MajorID_NoIsoforms));
      [i1,j1] = ind2sub(size(Protein.MajorID_NoIsoforms),ind1);
      ind2 = find(strcmp(Prot2,Protein.MajorID_NoIsoforms));
      [i2,j2] = ind2sub(size(Protein.MajorID_NoIsoforms),ind2);
      
      TP_Matrix(i1,i2) = 1;
      TP_Matrix(i2,i1) = 1;
    end
    
    
    % Make other matrices
    % - Possible interaction matrix
    % - Self interaction matrix
    
    % check which individual proteins are in corum
    % Important! Also check whether any member of that protein's group is in corum
    Ngroup = size(Protein.MajorID_NoIsoforms,2);
    inCorum = zeros(Ngroup,Dimensions_Gaus_import);
    for ni = 1:Ngroup
      tmp = Protein.MajorID_NoIsoforms(:,ni);
      % ismember can't handle empty cells, so as a workaround fill with a string that won't be in Corum
      isemptyCellArray = cellfun('isempty',tmp);
      tmp(isemptyCellArray) = {'--fake'};
      inCorum(ni,:) = ismember(tmp,Unique_Corum);
    end
    inCorum = sum(inCorum,1);
    
    % Calculate possible interaction matrix, asks are both proteins in corum
    Int_matrix = inCorum'*inCorum;
    
    % Create self interaction matrix
    Self_Int_matrix = zeros(Dimensions_Gaus_import,Dimensions_Gaus_import);
    for i=1:Dimensions_Gaus_import
      for ii=1:Dimensions_Gaus_import
        if isequal(Protein.NoIsoform{i},Protein.NoIsoform{ii})
          Self_Int_matrix(i,ii)=1;
        end
      end
    end
    inverse_self = ~Self_Int_matrix; %creation of inverse of self interactions
    clear Self_Int_matrix
    
    % Possible interactions, i.e. both in corum and not a self interaction
    possibleInts = inverse_self & Int_matrix;
    
    % FP matrix
    Inverse_TP_matrix = ~TP_Matrix;
    %FP_Matrix = (Int_matrix & Inverse_TP_matrix);
    
    % Save and clear replicate-specific variables
    sf = [maindir '/Output/tmp/' 'data_rep' num2str(replicate_counter) '_chan' num2str(channel_counter) '.mat'];
    save(sf,'TP_Matrix','possibleInts','Protein','inverse_self','Chromatograms_raw','Chromatograms','Dist')
    clear TP_Matrix possibleInts Protein inverse_self Chromatograms Dist
    
    tt = toc;
    fprintf('  ...  %.2f seconds\n',tt)
    
  end
  
  
  %% 5. Calculate interaction score
  tic
  fprintf('\n    5. Score each interaction via Naive Bayes classifier')
  
  % combine Dist, TP_Matrix, possibleInts from this channel's replicates
  proteinAll = {};
  indexList = nan(0,2);
  classList = nan(0,1);
  possList = nan(0,1);
  DistList = nan(10^7,NDistFields * length(replicatesThisChannel));
  ff = rand(1,15);
  kk = 0;
  for ii = 1:length(replicatesThisChannel)
    replicate_counter = replicatesThisChannel(ii);
    sf = [maindir '/Output/tmp/' 'data_rep' num2str(replicate_counter) '_chan' num2str(channel_counter) '.mat'];
    load(sf,'Dist','TP_Matrix','possibleInts','Protein','inverse_self')
    
    fn = fieldnames(Dist);
    NDistFields = length(fn);
    Icolumn = (ii-1)*NDistFields+1 : ii*NDistFields;
    
    for jj = 1:NDistFields
      N = size(Dist.(fn{jj}),1);
      Dist.(fn{jj})(~inverse_self==1) = nan;
    end
    
    % Read indicis of upper-triangular interaction matrices
    a = find(triu(~isnan(Dist.(fn{1}))));
    [prot1i, prot2i] = ind2sub(size(possibleInts),a);
    
    % Grow IsoformAll, list of all proteins in this channel
    proteinAll0 = proteinAll;
    I = ismember(Protein.Isoform,proteinAll);
    proteinAll(length(proteinAll)+1 : length(proteinAll)+sum(~I)) = Protein.Isoform(~I);
    [~,Iprot] = ismember(Protein.Isoform,proteinAll);
    
    % Where are this replicate's proteins in proteinAll?
    clear tmpprotlist
    tmpprotlist(:,1) = Iprot(prot1i);
    tmpprotlist(:,2) = Iprot(prot2i);
    tmpprotlist = sort(tmpprotlist,2);
    
    % Find interactions that overlap with previous
    [Ia,Ib] = ismember(tmpprotlist,indexList,'rows');
    % Add overlapping interactions to existing rows
    Iinsert = zeros(length(tmpprotlist),1);
    Iinsert(Ia) = Ib(Ia);
    % Append new interactions to the end of the list
    Inew = Iinsert==0;
    Iinsert(Inew) = length(indexList)+1 : length(indexList)+length(tmpprotlist)-sum(~Inew);
    
    indexList0 = indexList;
    classList0 = classList;
    possList0 = possList;
    indexList(Iinsert,:) = tmpprotlist;
    classList(Iinsert(Inew)) = TP_Matrix(a(Inew));           % Class label
    possList(Iinsert(Inew)) = possibleInts(a(Inew));           % interactions in the reference label
    for jj = 1:NDistFields
      DistList(Iinsert, Icolumn(jj)) = Dist.(fn{jj})(a);
    end
    kk = max([kk max(Iinsert)]);
    
    % sanity check: does indexList, possList, classList, proteinAll change?
    if ii>1
      tmp = indexList(1:size(indexList0,1),:);
      x1 = sum(indexList0(:) - tmp(:));
      x2 = sum(classList0(:) - classList(1:size(indexList0,1)));
      x3 = sum(possList0(:) - possList(1:size(indexList0,1)));
      x4 = sum(~ismember(proteinAll(1:length(proteinAll0)), proteinAll0));
      if x1~=0 || x2~=0 || x3~=0 || x4~=0
        error('Error combining distance measures across replicates')
      end
    end
    
    % housekeeping
    clear indexList0 classList0 possList0 Dist tmpprotlist tmp Ia Ib
  end
  DistList = DistList(1:kk,:);
  
  % Classifier method
  [scoreMatrix, feats_new] = scorenb(DistList,possList,classList);
  score = nanmedian(scoreMatrix,2);
  
  tt = toc;
  fprintf('  ...  %.2f seconds\n',tt)
  
  
  
  %% 6. Calculate score threshold for desired precision
  tic
  fprintf('\n    6. Calculate score threshold')
  
  [xcutoff, tmp] = calcPPIthreshold(score(possList==1),classList(possList==1),desiredPrecision);
  
  % save score for Complexes.m
  sf = [maindir '/Output/tmp/' 'score_chan' num2str(channel_counter) '.mat'];
  save(sf,'score','proteinAll','possList','classList','indexList','feats_new','xcutoff')
  
  precRange = tmp.precRange;
  recRange = tmp.recRange;
  
  % Check that desiredPrecision is within [min(precRange) max(precRange)]
  bad_desiredPrecision1 = desiredPrecision-.001>max(precRange);
  bad_desiredPrecision2 = desiredPrecision-.001<min(precRange);
  bad_desiredPrecision = bad_desiredPrecision1 | bad_desiredPrecision2;
  % If desiredPrecision is too high, lower it
  if ~isempty(find(bad_desiredPrecision1, 1))
    max_achievable_precision = max(precRange) * 0.99;
    fprintf('\n\n -------------------------- WARNING --------------------------')
    fprintf('\n Desired precision %6.4f is not achievable with current settings.', desiredPrecision(bad_desiredPrecision1))
    fprintf('\n Lowering it to %6.4f and recalculating score threshold... \n\n', max_achievable_precision)
    desiredPrecision(bad_desiredPrecision1) = floor(max_achievable_precision*100)/100;
    user.desiredPrecision = desiredPrecision;
    notei = length(notes);
    notes{notei+1} = ['Desired precision ' num2str(desiredPrecision(bad_desiredPrecision1)) ' was not achievable. It was lowered to ' num2str(max_achievable_precision) '.'];
  end
  % If desiredPrecision is too low, raise it
  if ~isempty(find(bad_desiredPrecision2, 1))
    min_achievable_precision = min(precRange);
    fprintf('\n\n -------------------------- WARNING --------------------------')
    fprintf('\n Desired precision %6.4f is not achievable with current settings.', desiredPrecision(bad_desiredPrecision1))
    fprintf('\n Raising it to %6.4f and recalculating score threshold... \n\n', min_achievable_precision)
    desiredPrecision(bad_desiredPrecision1) = round(min_achievable_precision*100)/100;
    user.desiredPrecision = desiredPrecision;
    notei = length(notes);
    notes{notei+1} = ['Desired precision ' num2str(desiredPrecision(bad_desiredPrecision1)) ' was not achievable. It was raised to ' num2str(min_achievable_precision) '.'];
  end
  
  tt = toc;
  fprintf('  ...  %.2f seconds\n',tt)
    
  
  %% 7. Find and concatenate interactions
  
  %for di = 1:length(desiredPrecision)
  di = 1;
  tic
  s = ['    7. Find and concatenate interactions at precision = ' num2str(desiredPrecision(di)) '\n'];
  fprintf(s)
  
  FN = sum(score(:)<xcutoff & possList==1 & classList==1);
  
  % Find binary interactions
  Final_interactions = find(score > xcutoff);
  for ii = 1:length(Final_interactions);
    interaction_count = interaction_count+1;
    ia = indexList(Final_interactions(ii),1);
    ib = indexList(Final_interactions(ii),2);
    protA = proteinAll{ia};
    protB = proteinAll{ib};
    clear tmp2
    tmp2{1} = protA;
    tmp2{2} = protB;
    tmp2 = sort(tmp2);
    binary_interaction_list{interaction_count,1} = strcat(tmp2{1},'-',tmp2{1}); % Interactions
    binary_interaction_list{interaction_count,2} = -1; %Dist.Euc(vii,viii); % Euc
    binary_interaction_list{interaction_count,3} = -1; %Dist.CoApex2(vii,viii); % Co-apex
    binary_interaction_list{interaction_count,4} = -1; %Dist.R2(vii,viii); % R^2
    binary_interaction_list{interaction_count,5} = -1; %0; % empty for now
    binary_interaction_list{interaction_count,6} = possList(Final_interactions(ii)); % Both_proteins_Corum
    binary_interaction_list{interaction_count,7} = classList(Final_interactions(ii)); % Interaction_in_Corum
    binary_interaction_list{interaction_count,8} = protA; % Protein_interactions1
    binary_interaction_list{interaction_count,9} = protB; % Protein_interactions2
    binary_interaction_list{interaction_count,10} = -1; %Protein.Center(vii); % Protein_interactions_center1
    binary_interaction_list{interaction_count,11} = -1; %Protein.Center(viii); % Protein_interactions_center2
    binary_interaction_list{interaction_count,12} = channel_counter;
    binary_interaction_list{interaction_count,13} = score(Final_interactions(ii)); % Interaction score
    binary_interaction_list{interaction_count,14} = -1; %scoreMatrix(vii,viii)>xcutoff(di); % Global interaction?
    binary_interaction_list{interaction_count,15} = -1; %scoreMatrix(vii,viii)>xcutoff_rep(replicate_counter,di); % Replicate-specific interaction?
  end
  binary_interaction_list = binary_interaction_list(1:interaction_count,:);
  
  tt = toc;
  fprintf('... %d interactions  ...  %.2f seconds\n',interaction_count,tt)
  
end

clear Dist TP_Matrix scoreMatrix possibleInts inverse_self


%% 8. Build the final interaction list
tic
fprintf('    8. Create a list of all the unique interactions')

interaction_pairs = cell(size(binary_interaction_list,1),1);
for ri=1:size(binary_interaction_list,1)
  interaction_pairs{ri} = [binary_interaction_list{ri,8} '_' binary_interaction_list{ri,9}];
end
interaction_pairs(cellfun('isempty',interaction_pairs)) = [];
[unique_interaction,x1,x2] = unique(interaction_pairs);
length_unique_inter = length(unique_interaction);

%Create array to store values
interaction_final.unique_interactions=cell(length_unique_inter,0);
interaction_final.replicate_numbers=nan(length_unique_inter,0);
interaction_final.proteinA=cell(length_unique_inter,0);
interaction_final.proteinB=cell(length_unique_inter,0);
interaction_final.CenterA=cell(length_unique_inter,0);
interaction_final.CenterB=cell(length_unique_inter,0);
interaction_final.DeltaCenter=cell(length_unique_inter,0);
interaction_final.R2=cell(length_unique_inter,0);
interaction_final.DeltaEucDist=cell(length_unique_inter,0);
interaction_final.proteinInCorum = nan(length_unique_inter,0);
interaction_final.interactionInCorum = nan(length_unique_inter,0);
interaction_final.score = nan(length_unique_inter,number_of_channels);
interaction_final.global=nan(length_unique_inter,1);
interaction_final.replicatespecific=nan(length_unique_inter,1);
cc = 0;
for ii = 1:length_unique_inter
  if mod(ii,1000)==0;disp(num2str(ii));end
  
  %Find location of unique protein interaction
  %Ibin = find(strcmp(unique_interaction(ii), interaction_pairs));
  Ibin = find(x2 == ii);
  
  %find which replicate the interactions were seen in
  diffC = abs(cellfun(@minus,binary_interaction_list(Ibin,10),...
    binary_interaction_list(Ibin,11)));
  goodC = find(diffC<2);
  
  if ~nnz(goodC)==0
    
    cc=cc+1;
    
    interaction_final.unique_interactions(cc,1) = unique_interaction(ii);
    interaction_final.proteinInCorum(cc,1) = binary_interaction_list{Ibin(goodC(1)),6};
    interaction_final.interactionInCorum(cc,1) = binary_interaction_list{Ibin(goodC(1)),7};
    
    for row = 1:length(goodC)
      %check size of replicate_numbers
      interaction_final.replicate_numbers(cc,row) = binary_interaction_list{Ibin(goodC(row)),12};
      
      % copy the values for the first unique interaction to the list
      interaction_final.proteinA(cc,row) = binary_interaction_list(Ibin(goodC(row)),8);
      interaction_final.proteinB(cc,row) = binary_interaction_list(Ibin(goodC(row)),9);
      interaction_final.CenterA(cc,row) = binary_interaction_list(Ibin(goodC(row)),10);
      interaction_final.CenterB(cc,row) = binary_interaction_list(Ibin(goodC(row)),11);
      interaction_final.DeltaCenter(cc,row) = binary_interaction_list(Ibin(goodC(row)),3);
      interaction_final.R2(cc,row) = binary_interaction_list(Ibin(goodC(row)),4);
      interaction_final.DeltaEucDist(cc,row) = binary_interaction_list(Ibin(goodC(row)),2);
      interaction_final.score(cc,row) = binary_interaction_list{Ibin(goodC(row)),13};
    end
    % At least one Gaussian meets global score threshold?
    interaction_final.global(cc) = sum(cell2mat(binary_interaction_list(Ibin(goodC),14)))>0;
    interaction_final.replicatespecific(cc) = sum(cell2mat(binary_interaction_list(Ibin(goodC),15)))>0;
    
    %If center not within two fractions
  elseif nnz(goodC)==0
    cc=cc+1;
    
    interaction_final.unique_interactions(cc,1)= unique_interaction(ii);
    %check size of replicate_numbers
    interaction_final.replicate_numbers(cc,1) = binary_interaction_list{Ibin(1),12};
    
    %copy the values for the first unique interaction to the list
    interaction_final.proteinA(cc,1) = binary_interaction_list(Ibin(1),8);
    interaction_final.proteinB(cc,1) = binary_interaction_list(Ibin(1),9);
    interaction_final.CenterA(cc,1) = binary_interaction_list(Ibin(1),10);
    interaction_final.CenterB(cc,1) = binary_interaction_list(Ibin(1),11);
    interaction_final.DeltaCenter(cc,1) = binary_interaction_list(Ibin(1),3);
    interaction_final.R2(cc,1) = binary_interaction_list(Ibin(1),4);
    interaction_final.DeltaEucDist(cc,1) = binary_interaction_list(Ibin(1),2);
    interaction_final.proteinInCorum(cc,1) = binary_interaction_list{Ibin(1),6};
    interaction_final.interactionInCorum(cc,1) = binary_interaction_list{Ibin(1),7};
    interaction_final.score(cc,1) = binary_interaction_list{Ibin(1),13};
    interaction_final.global(cc) = binary_interaction_list{Ibin(1),14};
    interaction_final.replicatespecific(cc) = binary_interaction_list{Ibin(1),15};
  end
end

% Calculate the average score for each pairwise interaction
%interaction_final.scoreavg = nanmean(interaction_final.score,2);
interaction_final.scoreavg = max(interaction_final.score,[],2);

%clear binary_interaction_list

tt = toc;
fprintf('  ...  %.2f seconds\n',tt)



%% Calculate final TP, FP, Recall, Precision
tic
fprintf('        9. Calculate final TP, FP, FN and TN')

% True interactions Observation across replicate
% Count how many times a true positive interactions was observed
Total_unique_interactions=length(interaction_final.unique_interactions);
number_observation_pre_interaction=zeros(Total_unique_interactions,1);
number_unique_interaction=zeros(Total_unique_interactions,1);
%Record the number of times a interaction was detected in irrespective of replicate number
interactions_prep_replicate=zeros(Total_unique_interactions,replicate_counter);
for ii = 1:Total_unique_interactions
  
  number_observation_pre_interaction(ii) = sum(interaction_final.replicate_numbers(ii,:)>0);
  
  %Determine if multiple gaussians within the same replicate were indentified
  testA = interaction_final.replicate_numbers(ii,:);
  testA(testA==0) = [];
  unique_testA = unique(testA);
  number_of_unique = length(unique_testA);
  number_unique_interaction(ii) = number_of_unique;
end

%Copy number of unique_interaction across replicate to interaction array
interaction_final.times_observed_across_reps = number_unique_interaction;

True_proteinInCorum = interaction_final.proteinInCorum > 0;
True_interactionInCorum = interaction_final.interactionInCorum > 0;
all_positives = sum(interaction_final.proteinInCorum & interaction_final.interactionInCorum) + sum(FN);

final_FP = sum(True_proteinInCorum & ~True_interactionInCorum);
final_TP = sum(True_proteinInCorum & True_interactionInCorum);

% final results
final_Recall = final_TP / all_positives;
final_Precision = final_TP/(final_TP+final_FP);
clearvars Neg_binary_interaction_list Ineg Ipos

% Determine channel-specific interactions
interaction_final.channel = cell(Total_unique_interactions,1);
interaction_final.channel_formatted = cell(Total_unique_interactions,1);
for ii = 1:Total_unique_interactions
  tmp = interaction_final.replicate_numbers(ii,:);
  tmp(tmp==0) = [];
  %interaction_final.channel{ii} = unique(rep2channel(tmp))';
  interaction_final.channel{ii} = tmp;
  tmp = cell(1,length(interaction_final.channel{ii}));
  for jj = 1:length(interaction_final.channel{ii})
    tmp{jj} = user.silacratios{interaction_final.channel{ii}(jj)};
  end
  if length(tmp)>1
    interaction_final.channel_formatted{ii} = strjoin(tmp, ' ; ');
  else
    interaction_final.channel_formatted{ii} = tmp{1};
  end
end

%Create list of scores, note ensure strjoin function is avalible
% interaction_final.scores_formated = cell(size(interaction_final.score,1),1);
% for format_loop=1:Total_unique_interactions
%   tmp = interaction_final.score(format_loop,:);
%   tmp(isnan(tmp)) = [];
%   stringout = cell(length(tmp),1);
%   for jj= 1:length(tmp)
%     stringout{jj} = num2str(tmp(jj));
%   end
%   interaction_final.scores_formated{format_loop} = strjoin(stringout,' ; ');
% end

% Create scoreRank, scorePrctile
%[~, ~, rankDescend] = unique(interaction_final.score);
%interaction_final.scoreRank = (1-rankDescend) + max(rankDescend);
%interaction_final.scorePrctile = interaction_final.scoreRank/length(interaction_final.scoreRank);

% Create precision-dropout
interaction_final.precisionDropout = nan(size(interaction_final.score));
interaction_final.precisionDropoutavg = nan(Total_unique_interactions,1);
interaction_final.precisionDropout_formatted = cell(Total_unique_interactions,1);
Recall_plot = nan(Total_unique_interactions,1);
FPR_plot = nan(Total_unique_interactions,1);
for ii = 1:Total_unique_interactions
  % condition-specific precision
  tmp1 = interaction_final.score(ii,:) - 10e-20;
  tmp2 = interaction_final.replicate_numbers(ii,:);
  stringout = cell(nnz(tmp2),1);
  for jj = 1:nnz(tmp2)
    cutoff = tmp1(jj);
    chan = tmp2(jj);
    Ithischan = interaction_final.replicate_numbers == chan;
    Igoodscore = interaction_final.score >= cutoff;
    Iinteract = sum(Ithischan,2)>0 & sum(Igoodscore,2)>0;
    TP = sum(interaction_final.proteinInCorum(Iinteract)==1 & interaction_final.interactionInCorum(Iinteract)==1);
    FP = sum(interaction_final.proteinInCorum(Iinteract)==1 & interaction_final.interactionInCorum(Iinteract)==0);
    if TP==0 && FP==0
      prec = 1;
    else
      prec = TP / (TP + FP);
    end
    interaction_final.precisionDropout(ii,jj) = prec;
    stringout{jj} = num2str(prec);
  end
  interaction_final.precisionDropout_formatted{ii} = strjoin(stringout,' ; ');
  
  % global precision
  cutoff = interaction_final.scoreavg(ii) - 10e-20;
  I = interaction_final.scoreavg >= cutoff;
  TP = sum(interaction_final.proteinInCorum(I)==1 & interaction_final.interactionInCorum(I)==1);
  FP = sum(interaction_final.proteinInCorum(I)==1 & interaction_final.interactionInCorum(I)==0);
  if TP==0 && FP==0
    prec = 1;
  else
    prec = TP / (TP + FP);
  end
  interaction_final.precisionDropoutavg(ii) = prec;
  Recall_plot(ii) = TP / all_positives;
end

% Count TP, FP, and novel interactions
I = (1:Total_unique_interactions)';
Precision_array = zeros(number_of_channels,3);
% 1 - FP, 2 - TP, 3 - novel
for jj = 1:number_of_channels
  % Find position of interaction observed within the required number of replicates
  position_to_sum = find(interaction_final.times_observed_across_reps>=jj & I);
  
  intInCor = interaction_final.interactionInCorum(position_to_sum);
  protInCor = interaction_final.proteinInCorum(position_to_sum);
  
  %Count the number values for non redundant list
  Precision_array(jj,1) = sum(protInCor & ~intInCor);
  Precision_array(jj,2) = sum(protInCor & intInCor);
  Precision_array(jj,3) = length(position_to_sum) - Precision_array(jj,1) - Precision_array(jj,2);
end


tt = toc;
fprintf('  ...  %.2f seconds\n',tt)


%%

tic
fprintf('    11. Write output files')
writeOutput_interactions
tt = toc;
fprintf('  ...  %.2f seconds\n',tt)


%%

tic
fprintf('    12. Make figures')
makeFigures_interactions
tt = toc;
fprintf('  ...  %.2f seconds\n',tt)


diary('off')