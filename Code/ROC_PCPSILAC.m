%%%%%%%%%%%%%%% Instructions for running Gauss_Build.m:
%
%
%
%%%%%%%%%%%%%%% Logic:      
%   0. Initialize
%   for replicate
%       1. Read input
%       2. Pre-process Gaussians, make Chromatograms and Dist.
%       3. Make Protein structure
%       4. Make TP, FP, possibleInts
%       5. Score each interaction (Naive Bayes)
%   6. Calculate score threshold with all replicates
%   for desiredPrecision
%       7. Find and concatenate interactions at desired precision
%       8. Create a list of all the unique interactions
%       9. Calculate final TP, FP, FN and TN
%       10. Find Treatment interactions
%       11. Write output files
%       12. Make figures
%     
%
%%%%%%%%%%%%%%% To do:          
% - Why is this script pre-processing chromatograms? Shouldn't that have been done earlier?
% - What's going on in 8b? We keep interactions with Delta_Center>2... so why check for it?
% - FN_binary_interaction_list can get up to 3 GB. Can we make it smaller?
%
%
%%%%%%%%%%%%%%% Fix the logic:
% - 'Look up Major Protein group Information' can be seriously improved with a few findstr
% - FN and TN are approximated. Instead of looking at all the combinations of Gaussians, it just
%   treats them all the same. This is for speed reasons: there were too many combinations to go
%   through individually. A fix is to treat things exactly how we treat TP, that is go through and
%   assess each Gaussian-Gaussian pair.
% - FN is calculated using only Int_matrix==1, but Int_matrix has values up to 4. Is it better to
%   just use Int_matrix>0 ?
%
%
%%%%%%%%%%%%%%% Bugs in my (this) code:
% - My Corum_Dataset is 11516x1. Nick's Corum_Dataset is 1x16393. Why are they different sizes?


disp('ROC_PCPSILAC.m')



%% 0. Initialize
tic
fprintf('\n    0. Initialize')



% Load user settings
maindir = user.maindir;
Experimental_channels = user.silacratios;
user.treatmentcondition;
desiredPrecision = user.desiredPrecision;
number_of_channels = length(user.silacratios);
InputFile{1} = user.majorproteingroupsfile;
InputFile{2} = user.corumfile;
InputFile{4} = user.mastergaussian;



% Define folders, i.e. define where everything lives.
codedir = [maindir 'Code/']; % where this script lives
funcdir = [maindir 'Code/Functions/']; % where small pieces of code live
datadir = [maindir 'Data/']; % where data files live
datadir1 = [datadir 'ROC/']; % where data files live
datadir2 = [datadir 'ROC/tmp/']; % where data files live
datadir3 = [datadir 'ROC/CombinedResults/']; % where data files live
figdir1 = [maindir 'Figures/ROC/']; % where figures live
tmpdir = '/Users/Mercy/Academics/Foster/NickCodeData/4B_ROC_homologue_DB/Combined Cyto new DB/';
%tmpdir = '/Users/Mercy/Academics/Foster/Jenny_PCPSILAC/PCPSILAC_Analysis/Data/Alignment/';
% Make folders if necessary
if ~exist(codedir, 'dir'); mkdir(codedir); end
if ~exist(funcdir, 'dir'); mkdir(funcdir); end
if ~exist(datadir, 'dir'); mkdir(datadir); end
if ~exist(datadir1, 'dir'); mkdir(datadir1); end
if ~exist(datadir2, 'dir'); mkdir(datadir2); end
if ~exist(datadir3, 'dir'); mkdir(datadir3); end
if ~exist(figdir1, 'dir'); mkdir(figdir1); end

%define Raw SILAC ratios data, This is the output from the alignment script
MvsL_filename_Raw_rep1=[tmpdir 'MvsL_alignment/Realignment/Adjusted_MvsL_Raw_for_ROC_analysis_rep1.csv'];
MvsL_filename_Raw_rep2=[tmpdir 'MvsL_alignment/Realignment/Adjusted_MvsL_Raw_for_ROC_analysis_rep2.csv'];
MvsL_filename_Raw_rep3=[tmpdir 'MvsL_alignment/Realignment/Adjusted_MvsL_Raw_for_ROC_analysis_rep3.csv'];
HvsL_filename_Raw_rep1=[tmpdir 'HvsL_alignment/Realignment/Adjusted_HvsL_Raw_for_ROC_analysis_rep1.csv'];
HvsL_filename_Raw_rep2=[tmpdir 'HvsL_alignment/Realignment/Adjusted_HvsL_Raw_for_ROC_analysis_rep2.csv'];
HvsL_filename_Raw_rep3=[tmpdir 'HvsL_alignment/Realignment/Adjusted_HvsL_Raw_for_ROC_analysis_rep3.csv'];
MvsL_filename_gaus_rep1=[tmpdir 'MvsL_alignment/Realignment/Adjusted_Combined_OutputGaus_rep1.csv'];
MvsL_filename_gaus_rep2=[tmpdir 'MvsL_alignment/Realignment/Adjusted_Combined_OutputGaus_rep2.csv'];
MvsL_filename_gaus_rep3=[tmpdir 'MvsL_alignment/Realignment/Adjusted_Combined_OutputGaus_rep3.csv'];
HvsL_filename_gaus_rep1=[tmpdir 'HvsL_alignment/Realignment/Adjusted_Combined_OutputGaus_rep1.csv'];
HvsL_filename_gaus_rep2=[tmpdir 'HvsL_alignment/Realignment/Adjusted_Combined_OutputGaus_rep2.csv'];
HvsL_filename_gaus_rep3=[tmpdir 'HvsL_alignment/Realignment/Adjusted_Combined_OutputGaus_rep3.csv'];
List_of_Raw_filename={MvsL_filename_Raw_rep1,MvsL_filename_Raw_rep2,MvsL_filename_Raw_rep3,...
  HvsL_filename_Raw_rep1,HvsL_filename_Raw_rep2,HvsL_filename_Raw_rep3};
List_of_Gaus_filename={MvsL_filename_gaus_rep1,MvsL_filename_gaus_rep2,MvsL_filename_gaus_rep3,...
  HvsL_filename_gaus_rep1,HvsL_filename_gaus_rep2,HvsL_filename_gaus_rep3};

number_of_replicates = length(List_of_Raw_filename) / number_of_channels;
I = find(strcmp(user.silacratios,user.treatmentcondition));
treatment_replicates = (1:number_of_replicates)+(I-1)*number_of_replicates;
I = 1:number_of_replicates*number_of_channels;
untreatment_replicates = I(~ismember(I,treatment_replicates));


% %define Raw SILAC ratios data, This is the output from the alignment script
% % Do this dynamically. Find filenames
% dd = dir([tmpdir 'Adjusted*Raw_for_ROC_analysis*rep*csv']);
% List_of_Raw_filename = cell(size(dd));
% for di = 1:length(dd)
%   List_of_Raw_filename{di} = dd(di).name;
% end
% dd = dir([tmpdir 'Adjusted_Combined_OutputGaus*rep*csv']);
% List_of_Gaus_filename = cell(size(dd));
% for di = 1:length(dd)
%   List_of_Gaus_filename{di} = dd(di).name;
% end

% pre-allocate final results
Recall = nan(number_of_replicates,number_of_channels);
Precision = nan((size(Recall)));
TPR = nan((size(Recall)));
FPR = nan((size(Recall)));
No_Int = nan((size(Recall)));

tt = toc;
fprintf('  ...  %.2f seconds\n',tt)


%%%%% Replicate counter starts
for replicate_counter = 1:number_of_replicates*number_of_channels
  
  s = ['\n        Replicate ' num2str(replicate_counter)];
  fprintf(s)
  
  %% 1. Read input data
  tic
  s = '\n        1. Read input data';
  fprintf(s)
  
  [~, Protein_IDs] = xlsread(InputFile{1});
  
  % Corum binary interactions
  fid=fopen(InputFile{2}, 'rt');    %Input corum data base.
  Corum_Import= textscan(fid,'%s\t', 'Delimiter',',');
  fclose(fid);
  No=length(Corum_Import{1})/2;
  Corum_Protein_names = (reshape(Corum_Import{1,1},2,No)');
  Unique_Corum = unique(Corum_Protein_names);
  
  %Master_Gaussian_list=importdata(InputFile{4},',');
  
  %Raw data for analysis of euclidean distance
  %import data files summary data
  Maxquant = fopen (List_of_Raw_filename{replicate_counter});
  Maxquant_raw = textscan(Maxquant, '%s', 'Delimiter',',');
  fclose(Maxquant);
  
  %Reshape data for use in analysis
  Raw_Dimension=size(Maxquant_raw{:});
  Raw_data_reshaped=reshape(Maxquant_raw{:},68,(Raw_Dimension(1)/68))';
  Raw_data_reshaped_size=size(Raw_data_reshaped);
  
  %divide up reshaped data in correct form
  Chromatograms_raw = Raw_data_reshaped(2:end,3:end);
  Chromatograms_raw = cellfun(@str2num,Chromatograms_raw);
  
  %Gaus data data for analysis of center analysis
  Gaus_import_name = fopen (List_of_Gaus_filename{replicate_counter});
  Gaus_import = textscan(Gaus_import_name, '%s', 'Delimiter',',');
  fclose(Gaus_import_name);
  
  %Reshape data for use in analysis
  Gaus_import_Dimension=size(Gaus_import{:});
  Gaus_import_reshaped=reshape(Gaus_import{:},7,(Gaus_import_Dimension(1)/7))';
  
  % do a bit of housekeeping
  clear Maxquant_raw Raw_data_reshaped Corum_Import
  
  
  tt = toc;
  fprintf('  ...  %.2f seconds\n',tt)
  
  
  
  %% 2. Pre-process gaussians, chromatograms. Get Dist
  % - In chromatograms, replaces nans with 0.05
  % - Remove gaussians and chromatograms with C<5
  % - Make normalized gaussians, i.e. set height to 1
  tic
  fprintf('        2. Clean gaussians and chromatograms, calculate Dist.')
  
  % Pre-process data to remove Gaussian above fraction five
  % Create array to use for comparsions
  H_raw = cellfun(@str2num,Gaus_import_reshaped(2:end,2));
  C_raw = cellfun(@str2num,Gaus_import_reshaped(2:end,3));
  W_raw = cellfun(@str2num,Gaus_import_reshaped(2:end,4));
  
  % Replace NaN values with 0.05
  Chromatograms = Chromatograms_raw;
  Chromatograms(isnan(Chromatograms))= 0.05;
  
  % remove gaussian with centers below five
  Ibad = C_raw<5;
  H = H_raw(~Ibad);
  C = C_raw(~Ibad);
  W = W_raw(~Ibad);
  Chromatograms = Chromatograms(~Ibad,:);
  Chromatograms_raw = Chromatograms_raw(~Ibad,:);
  
  % Make normalized gaussians
  %   x = 0.25:0.25:100;
  %   gaussian_fit=zeros(length(H),length(x));
  %   for ii = 1:length(H)
  %     gaussian_fit(ii,:) = 1*exp(-(x- C(ii)).^2 /2/(W(ii).^2));
  %   end
  
  % How many Gaussians?
  I = find(~Ibad);
  unique_names = Gaus_import_reshaped(I+1,1);
  X = repmat(unique_names,1,length(unique_names));
  Ngauss = sum(strcmp(X,X'),1)';
  clear X unique_names
  
  % Co-Apex score = norm(C) / sqrt(length(C))
  I = find(abs(diff(Ngauss))>0 | Ngauss(2:end)==1);
  I = [1; I+1];
  I0 = 0;
  CoApex = zeros(size(Ngauss));
  for ii = 1:length(I)
    I2 = I0+1 : I0 + Ngauss(I(ii));
    CoApex(I2) = norm(C(I2)) / sqrt(length(I2));
    I0 = max(I2);
  end
  
  % Area under the chromatogram
  auc = sum(Chromatograms,2);
  
  % Reduce variables from gaussian-level to protein-level
  %   Ireduce = find(abs(diff(Ngauss))>0 | Ngauss(2:end)==1);
  %   Ireduce = [1; Ireduce+1];
  %   Chromatograms = Chromatograms(Ireduce,:);
  %   C = C(Ireduce,:);
  %   Ngauss = Ngauss(Ireduce);
  %   CoApex = CoApex(Ireduce);
  %   auc = auc(Ireduce,:);
  %   Gaus_import_reshaped(2:end,:) = Gaus_import_reshaped(Ireduce+1,:);
  
  % Calculate distance matrices
  Dist.Euc = squareform(pdist(Chromatograms,'euclidean'));
  Dist.Center = squareform(pdist(C,'euclidean'));
  Dist.R2 = 1 - corr(Chromatograms').^2; % one minus R squared
  Dist.Ngauss = squareform(pdist(Ngauss));
  Dist.CoApex = squareform(pdist(CoApex));
  Dist.AUC = squareform(pdist(auc));
  
  % GRAVEYARD OF RELEGATED DIST MATRICES
  %   Dist.RawOverlap = nan(size(Chromatograms_raw,1),size(Chromatograms_raw,1));
  %   for ii = 1:size(Chromatograms_raw,1)
  %     for jj = 1:size(Chromatograms_raw,1)
  %       Dist.RawOverlap(ii,jj) = sum(~isnan(Chromatograms_raw(ii,:)) & ~isnan(Chromatograms_raw(jj,:)));
  %     end
  %   end
  %   Dist.R2raw = corr(Chromatograms_raw','rows','pairwise').^2; % one minus R squared
  %   Dist.R2raw2 = Dist.R2raw.*Noverlap;
  %   Dist.Height = squareform(pdist(H,'euclidean'));
  %   Dist.Width = squareform(pdist(W,'euclidean'));
  %   Dist.Gaussian_fits = squareform(pdist(gaussian_fit,'euclidean'));
  %   Dist.dtw = ones(size(Dist.R2))*100;
  %    for ii = 1:size(Chromatograms,1)
  %      if mod(ii,100)==0;disp(num2str(ii));end
  %      for jj = 1:size(Chromatograms,1)
  %        if Dist.R2(ii,jj)<0.4
  %          Dist.dtw(ii,jj) = dtw(Chromatograms(ii,:)',Chromatograms(jj,:)');
  %          %Dist.dtw(ii,jj) = dtw_old(Chromatograms(ii,:),Chromatograms(jj,:));
  %        end
  %      end
  %    end
  
  tt = toc;
  fprintf('  ...  %.2f seconds\n',tt)
  
  
  
  %% 3. Make Protein structure
  % This summarizes the proteins in our sample
  tic
  fprintf('        3. Make Protein structure.')
  
  % Combine Majority Protein ID with Proteins Identified in replicate
  Protein.Isoform=Gaus_import_reshaped(2:end,1);
  Protein.Height=cellfun(@str2num,Gaus_import_reshaped(2:end,2));
  Protein.Width=cellfun(@str2num,Gaus_import_reshaped(2:end,3));
  Protein.Center=cellfun(@str2num,Gaus_import_reshaped(2:end,4));
  
  % remove protein name with center less then 5
  Protein.Isoform(Ibad)=[];
  Protein.Height(Ibad)=[];
  Protein.Width(Ibad)=[];
  Protein.Center(Ibad)=[];
  
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
  for Protein_ID_counter3 = 2:Dimension_of_Protein_IDs(1)
    Array_to_check1=Protein_IDS_no_isoform(Protein_ID_counter3,:);
    Array_to_check_empty=cellfun(@isempty,Array_to_check1);
    Array_to_check1(Array_to_check_empty)=[]; % Replace empty cells with NaN
    uniqueCells=unique(Array_to_check1);
    for Protein_ID_counter4 = 1:length(uniqueCells)
      Protein_IDS_no_isoform_no_dup_processed{Protein_ID_counter3,Protein_ID_counter4}=uniqueCells(Protein_ID_counter4);
    end
  end
  
  %Look up Major Protein group Information
  % For each protein that was fitted by a Gaussian
  %   Find the row containing that protein name in Protein_IDs
  %   Store all protein names from that row of Protein_IDs
  Protein.MajorIDs=cell(Dimensions_Gaus_import,Dimension_of_Protein_IDs(2));
  Protein.MajorID_NoIsoforms=cell(Dimensions_Gaus_import,Dimension_of_Protein_IDs(2));
  f = 0;
  for Gaus_import_counter1 = 1:Dimensions_Gaus_import(1)
    
    %disp(Lookup_protein_name)
    Lookup_protein_name = Protein.Isoform(Gaus_import_counter1);
    
    tmp = find(strcmp(Lookup_protein_name,Protein_IDs));
    [Check_counter1, Check_counter2] = ind2sub(size(Protein_IDs),tmp);
    
    %copy names to MajorID from Protein_IDs raw data
    Array_to_check1 = Protein_IDs(Check_counter1,:);
    Array_to_check_empty1 = cellfun(@isempty,Array_to_check1);
    for Check_counter3 = 1:nnz(~Array_to_check_empty1)
      Protein.MajorIDs{Gaus_import_counter1,Check_counter3} = {Protein_IDs{Check_counter1,Check_counter3}};
    end
    
    %copy names to MajorID_NoIsoforms from Protein_IDS_no_isoform_no_dup_processed
    Array_to_check2=Protein_IDS_no_isoform_no_dup_processed(Check_counter1,:);
    Array_to_check_empty2=cellfun(@isempty,Array_to_check2);
    for Check_counter4 = 1:nnz(~Array_to_check_empty2)
      a = Protein_IDS_no_isoform_no_dup_processed{Check_counter1,Check_counter4};
      Protein.MajorID_NoIsoforms{Gaus_import_counter1,Check_counter4}=a{1};
    end
  end
  
  dimension_Protein_MajorID_NoIsoforms=size(Protein.MajorID_NoIsoforms);
  
  tt = toc;
  fprintf('  ...  %.2f seconds\n',tt)
  
  
  
  %% 4. Make TP, FP matrices
  tic
  fprintf('        4. Make TP, FP, FN matrices')
  
  % Make Pos_Corum_proteins, which is all individual proteins that are in corum and our dataset
  Unique_MajorProteinIDs = unique(vertcat(Protein.MajorID_NoIsoforms{:}),'rows');
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
  Corum_Dataset_expanded = cell(size(Corum_Dataset));
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
  end
  %I = unique(cell2mat(Corum_Dataset_expanded),'rows');
  %Corum_Dataset_expanded = unique(Corum_Dataset_expanded{ii});
  
  
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
  inCorum = sum(inCorum);
  
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
  FP_Matrix = (Int_matrix & Inverse_TP_matrix);
  
  tt = toc;
  fprintf('  ...  %.2f seconds\n',tt)
  
  
  
  %% 5. Calculate score, the likelihood that a protein pair is interacting
  tic
  fprintf('        5. Score each interaciton via Naive Bayes classifier')
  % score approx_eq the likelihood that a protein pair is interacting
  
  % Predict interactions through Euc
  %scoreMatrix = 1 - Dist.Euc;
  %score0 = scoreMatrix(possibleInts(:));
  
  % Predict interactions through R^2
  %scoreMatrix = 1 - Dist.R2;
  
  % Predict interactions with a classifier
  %scoreMatrix = scoresvm(Dist,possibleInts,TP_Matrix);
  %scoreMatrix = nanmean(scoreMatrix,2);
  %scoreMatrix = reshape(scoreMatrix,size(Dist.R2,1),size(Dist.R2,1));
  scoreMatrix = scorenb(Dist,possibleInts,TP_Matrix);
  scoreMatrix = median(scoreMatrix,2);
  %scoreMatrix = nanmean(scoreMatrix,2);
  %scoreMatrix = reshape(scoreMatrix,size(Dist.R2,1),size(Dist.R2,1));
  
  sf = [datadir2 'score_rep' num2str(replicate_counter) '.mat'];
  save(sf,'scoreMatrix','TP_Matrix','possibleInts','Protein','inverse_self','Chromatograms')
  clear scoreMatrix TP_Matrix possibleInts Protein inverse_self Chromatograms
  
  tt = toc;
  fprintf('  ...  %.2f seconds\n',tt)  
end
%%%%% Replicate counter ends




% do a bit of housekeeping
clear Dist Int_matrix





%% 6. Find the score thresholds

tic
fprintf('    6. Find the score cutoff across all replicates')

% Concatenate scores across replicates.
% The following code converts Protein1 and Protein2 to unique (hopefully!) numeric values.
% This greatly reduces the memory used, as otherwise allScores is ~2.5GB.
% NB: It's possible that two protein names will be mapped to the same numeric value.
allScores = zeros(10^7,5);
ff = rand(1,15);
kk = 0;
for replicate_counter = 1:(number_of_replicates*number_of_channels)
  sf = [datadir2 'score_rep' num2str(replicate_counter) '.mat'];
  load(sf)
  a = find(triu(possibleInts));
  [prot1i, prot2i] = ind2sub(size(possibleInts),a);
  I = kk+1 : kk+length(a);
  allScores(I,1) = TP_Matrix(a);           % Class label
  allScores(I,4) = replicate_counter;         % Replicate
  allScores(I,5) = median(scoreMatrix(a),2);  % Median score (equivalent to ensemble voting)
  
  % For space reasons, let's try converting Protein1, Protein2 into numeric
  kk0 = kk;
  for ii = 1:length(prot1i)
    kk = kk+1;
    tmp = double(Protein.Isoform{prot1i(ii)});
    allScores(kk,2) = sum(tmp .* ff(1:length(tmp)));
    tmp = double(Protein.Isoform{prot2i(ii)});
    allScores(kk,3) = sum(tmp .* ff(1:length(tmp)));
  end
  
  % Check the name-to-numeric conversion
  % Did we lose any unique names in the conversion?
  N1a = length(unique(Protein.Isoform(prot1i)));
  N2a = length(unique(Protein.Isoform(prot2i)));
  N1b = length(unique(allScores(kk0+1:kk,2)));
  N2b = length(unique(allScores(kk0+1:kk,3)));
  
  if N1a~=N1b || N2a~=N2b
    disp('Two protein names stored as the same number!')
  end
  
  clear scoreMatrix inverse_self Protein possibleInts TP_Matrix
end
allScores = allScores(1:kk,:);

% Remove redundant interactions
% Options: i) max score, ii) mean score, iii) median score
[~,I1,I2] = unique(allScores(:,2:3),'rows');
allScores2 = nan(length(I1),100);
countV = zeros(length(I1),1);
for ii = 1:length(I2)
  countV(I2(ii)) = countV(I2(ii))+1; % how many times have we seen this interaction?
  allScores2(I2(ii),1) = allScores(ii,1);
  allScores2(I2(ii),countV(I2(ii))+1) = allScores(ii,5);
end
class = allScores2(:,1);
score = max(allScores2(:,2:end),[],2);      % i) max score
%score = nanmean(allScores2(:,2:end),2);    % ii) mean score
%score = nanmedian(allScores2(:,2:end),2);  % iii) median score
clear allScores2 allScores TP_Matrix I1 I2


% 6b. Calculate the score threshold that gives the desired precision
% Algorithm: Calculate precision as a function of score very coarsely. Identify at what score
% value precision crosses the desired level. Zoom in on that score value. Iterate.
nn = 25;
Tol = 0.001; % get within 0.1% precision
maxIter = 20; % zoom in 20 times at most

xcutoff = nan(size(desiredPrecision));
calcprec = zeros(size(xcutoff));
calcrec = zeros(size(xcutoff));
scoreRange = nan(nn*maxIter*length(desiredPrecision),1);  % Save 
precRange = nan(size(scoreRange));                        % these
recRange = nan(size(scoreRange));                         % for
tprRange = nan(size(scoreRange));                         % plotting
fprRange = nan(size(scoreRange));                         % later.
for di = 1:length(desiredPrecision)
  ds = linspace(min(score(:)),max(score(:)),nn); % start off with coarse search
  calcTol = 10^10;
  iter = 0;
  deltaPrec = 1;
  prec0 = zeros(nn,1);
  % Stop zooming in when one of three things happens:
  % i) you get close enough to the desired precision (within calcTol)
  % ii) you've been through maxIter iterations
  % iii) zooming in stops being useful (precision changes by less than deltaPrec b/w iterations)
  while calcTol>Tol && iter<maxIter && deltaPrec>1e-3
    iter=iter+1;
    rec = nan(nn,1);
    prec = nan(nn,1);
    fpr = nan(nn,1);
    tpr = nan(nn,1);
    for dd =1:length(ds)
      % ensemble VOTING
      TP = sum(score>ds(dd) & class==1);
      FP = sum(score>ds(dd) & class==0);
      FN = sum(score<ds(dd) & class==1);
      TN = sum(score<ds(dd) & class==0);
      prec(dd) = TP/(TP+FP);
      rec(dd) = TP/(TP+FN);
      fpr(dd) = FP/(FP+TN);
      tpr(dd) = TP/(TP+FN);
    end
    deltaPrec = nanmean(abs(prec - prec0));
    
    % Save vectors for plotting
    I = (iter-1)*nn+1 : iter*nn;
    scoreRange(I) = ds;
    precRange(I) = prec;
    recRange(I) = rec;
    tprRange(I) = tpr;
    fprRange(I) = fpr;
    
    % Calculate how close to desiredPrecision(di) you got
    [calcTol,I] = min(abs(prec - desiredPrecision(di)));
    
    % Zoom in on region of interest
    i1 = find(prec>desiredPrecision(di));
    if isempty(i1);
      mx = max(score(:));
    else
      mx = ds(i1(1));
    end
    i2 = find(prec<desiredPrecision(di));
    if isempty(i2);
      mn = min(score(:));
    else
      mn = ds(i2(end));
    end
    ds = linspace(mn,mx,nn);
    prec0 = prec;
  end
  xcutoff(di) = ds(I);
  calcprec(di) = prec(I);
  calcrec(di) = rec(I);
end
[scoreRange,I] = sort(scoreRange);
precRange = precRange(I);
recRange = recRange(I);
tprRange = tprRange(I);
fprRange = fprRange(I);

tt = toc;
fprintf('  ...  %.2f seconds\n',tt)




%% 7. Find and concatenate interactions at desired precision'

for pri = 1:length(desiredPrecision)
  tic
  s = ['        7. Find and concatenate interactions at precision = ' num2str(desiredPrecision(pri)) '\n'];
  fprintf(s)
  
  % Set these here, since we want to concatenate across replicates
  fn_count = 0; % false negative interaction counter, used in 7b
  interaction_count = 0; % interaction counter, used in 7c
  
  binary_interaction_list = cell(2*10^4,13);
  FN_binary_interaction_list = zeros(10^7,13);
  
  for replicate_counter = 1:number_of_replicates*number_of_channels
    
    % 7b. Calculate recall and precision, negative interactions
    %Determine the false "binary" interaction list
    tic
    fprintf(['            Replicate ' num2str(replicate_counter)])
    
    % Re-load scoreMatrix, TP_Matrix, possibleInts, and Protein
    sf = [datadir2 'score_rep' num2str(replicate_counter) '.mat'];
    load(sf)
    
    % determine final interaction matrices
    pos = scoreMatrix>xcutoff(pri); % positive hits
    pos = reshape(pos,size(inverse_self,1),size(inverse_self,1));
    Final_interactions = pos & inverse_self;

    % Calculate TP, FP, TN, FN
    TP = sum(Final_interactions(:) & TP_Matrix(:) & possibleInts(:));
    FP = sum(Final_interactions(:) & ~TP_Matrix(:) & possibleInts(:));
    TN = sum(~Final_interactions(:) & ~TP_Matrix(:) & possibleInts(:));
    FN = sum(~Final_interactions(:) & TP_Matrix(:) & possibleInts(:));
    
    %Generate negative matrix for final interaction analysis
    Negative2_matrix = ~Final_interactions & possibleInts;
    
    % final results
    Recall(replicate_counter,pri) = TP/(TP+FN);
    Precision(replicate_counter,pri) = TP/(TP+FP);
    TPR(replicate_counter,pri) = TP/(TP+FN);
    FPR(replicate_counter,pri) = FP/(FP+TN);
    No_Int(replicate_counter,pri) = length(find(triu(Final_interactions==1)));
    
    % Predicted negative interactions
    scoreMatrix = median(scoreMatrix,2);
    scoreMatrix = reshape(scoreMatrix,size(Final_interactions,1),size(Final_interactions,1));
    Protein_elution = strcat(num2str(Protein.Center),'*',Protein.Isoform);
    U1 = triu(Negative2_matrix); % Take half of matrix therefor unique interactions
    [ia,ib] = find(U1 == 1); % protein1, protein2 indices in false negatives
    for ri=1:length(ia)
      ri1 = ia(ri);
      ri2 = ib(ri);
      Int1=Protein_elution(ri1);
      Int2=Protein_elution(ri2);
      if Protein.Center(ri1)==Protein.Center(ri2) && ri1==ri2 % make sure it's not a self interaction
      else
        fn_count = fn_count+1;
        tmp = double(Protein.Isoform{ri1});
        tmpname1 = sum(tmp .* ff(1:length(tmp)));
        tmp = double(Protein.Isoform{ri2});
        tmpname2 = sum(tmp .* ff(1:length(tmp)));
        FN_binary_interaction_list(fn_count,1) = tmpname1; % Protein_interactions1
        FN_binary_interaction_list(fn_count,2) = tmpname2; % Protein_interactions2
        FN_binary_interaction_list(fn_count,3) = TP_Matrix(ri1,ri2); % Interaction_in_Corum
      end
      % sanity check
      if ~strcmp(Int1, Int2) && Protein.Center(ri1)==Protein.Center(ri2) && ri1==ri2
        disp('uh oh!')
      end
    end
    
    
    % 7c. Find "Binary" interactions
    Protein_elution= strcat(num2str(Protein.Center),'*',Protein.Isoform);
    [n,m]=size(Final_interactions);
    U = triu(Final_interactions);    %Take half of matrix therefor unique interactions
    for vii=1:n
      for viii= 1:m
        if U(vii,viii)==1
          Int1=Protein_elution{vii};
          Int2=Protein_elution{viii};
          if strcmp(Int1, Int2)
          else
            interaction_count=1+interaction_count;
            binary_interaction_list{interaction_count,1} = strcat(Int1,'-',Int2); % Interactions
            binary_interaction_list{interaction_count,2} = 0; % Delta_EucDist
            %binary_interaction_list{interaction_count,2} = norm(Chromatograms(vii,:)-Chromatograms(viii,:)); % Delta_EucDist
            binary_interaction_list{interaction_count,3} = abs(Protein.Center(vii)-Protein.Center(viii)); % Delta_Center
            binary_interaction_list{interaction_count,4} = abs(Protein.Height(vii)-Protein.Height(viii)); % Delta_Height
            binary_interaction_list{interaction_count,5} = abs(Protein.Width(vii)-Protein.Width(viii)); % Delta_Width
            binary_interaction_list{interaction_count,6} = possibleInts(vii,viii); % Both_proteins_Corum
            binary_interaction_list{interaction_count,7} = TP_Matrix(vii,viii); % Interaction_in_Corum
            binary_interaction_list{interaction_count,8} = Protein.Isoform{vii}; % Protein_interactions1
            binary_interaction_list{interaction_count,9} = Protein.Isoform{viii}; % Protein_interactions2
            binary_interaction_list{interaction_count,10} = Protein.Center(vii); % Protein_interactions_center1
            binary_interaction_list{interaction_count,11} = Protein.Center(viii); % Protein_interactions_center2
            binary_interaction_list{interaction_count,12} = replicate_counter; % Protein_interactions_center2
            binary_interaction_list{interaction_count,13} = scoreMatrix(vii,viii); % Interaction score
          end
        else
        end
      end
    end
        
    tt = toc;
    fprintf('  ...  %.2f seconds\n',tt)
    
  end
  
  binary_interaction_list = binary_interaction_list(1:interaction_count,:);
  FN_binary_interaction_list = FN_binary_interaction_list(1:fn_count,:);

  
  
  tic
  fprintf('        8. Create a list of all the unique interactions')
  
  interaction_pairs = cell(size(binary_interaction_list,1),1);
  for ri=1:size(binary_interaction_list,1)
    interaction_pairs{ri} = [binary_interaction_list{ri,8} '_' binary_interaction_list{ri,9}];
  end
  unique_interaction = unique(interaction_pairs);
  length_unique_inter = length(unique_interaction);
  
  %Create array to store values
  interaction_final.unique_interactions=cell(0,0);
  interaction_final.replicate_numbers=cell(0,0);
  interaction_final.proteinA=cell(0,0);
  interaction_final.proteinB=cell(0,0);
  interaction_final.CenterA=cell(0,0);
  interaction_final.CenterB=cell(0,0);
  interaction_final.DeltaHeight=cell(0,0);
  interaction_final.DeltaCenter=cell(0,0);
  interaction_final.DeltaWidth=cell(0,0);
  interaction_final.DeltaEucDist=cell(0,0);
  interaction_final.proteinInCorum=cell(0,0);
  interaction_final.interactionInCorum=cell(0,0);
  interaction_final.score=nan(length_unique_inter,60);
  Unique_interaction_counter=0;
  for interaction_counter2A = 1:length_unique_inter
    %Find location of unique protein interaction
    location_interaction_pairs = find(strcmp(unique_interaction(interaction_counter2A), interaction_pairs));
    
    %find which replicate the interactions were seen in
    diffC = abs(cellfun(@minus,binary_interaction_list(location_interaction_pairs,10),...
      binary_interaction_list(location_interaction_pairs,11)));
    position_within_two = find(diffC<100);
    
    %Check if multiple centers are withing two fractions of each other
    if ~nnz(position_within_two)==0
      %Add one to counters and remove location already checked
      %row_counter=row_counter+1;
      Unique_interaction_counter=Unique_interaction_counter+1;
       
      interaction_final.unique_interactions(Unique_interaction_counter,1) = unique_interaction(interaction_counter2A);
      interaction_final.proteinInCorum(Unique_interaction_counter,1) = binary_interaction_list(location_interaction_pairs(position_within_two(1)),6);
      interaction_final.interactionInCorum(Unique_interaction_counter,1) = binary_interaction_list(location_interaction_pairs(position_within_two(1)),7);

      for row = 1:length(position_within_two)
        %check size of replicate_numbers
        interaction_final.replicate_numbers{Unique_interaction_counter,row} = ...
          binary_interaction_list{location_interaction_pairs(position_within_two(row)),12};
        
        % copy the values for the first unique interaction to the list
        interaction_final.proteinA(Unique_interaction_counter,row) = binary_interaction_list(location_interaction_pairs(position_within_two(row)),8);
        interaction_final.proteinB(Unique_interaction_counter,row) = binary_interaction_list(location_interaction_pairs(position_within_two(row)),9);
        interaction_final.CenterA(Unique_interaction_counter,row) = binary_interaction_list(location_interaction_pairs(position_within_two(row)),10);
        interaction_final.CenterB(Unique_interaction_counter,row) = binary_interaction_list(location_interaction_pairs(position_within_two(row)),11);
        interaction_final.DeltaHeight(Unique_interaction_counter,row) = binary_interaction_list(location_interaction_pairs(position_within_two(row)),4);
        interaction_final.DeltaCenter(Unique_interaction_counter,row) = binary_interaction_list(location_interaction_pairs(position_within_two(row)),3);
        interaction_final.DeltaWidth(Unique_interaction_counter,row) = binary_interaction_list(location_interaction_pairs(position_within_two(row)),5);
        interaction_final.DeltaEucDist(Unique_interaction_counter,row) = binary_interaction_list(location_interaction_pairs(position_within_two(row)),2);
        interaction_final.score(Unique_interaction_counter,row) = binary_interaction_list{location_interaction_pairs(position_within_two(row)),13};
        
      end
      %location_interaction_pairs(position_within_two)=[];
      
      %If center not within two fractions
    elseif nnz(position_within_two)==0
      %Add one to counters and remove location already checked
      %location_interaction_pairs(1)=[];
      Unique_interaction_counter=Unique_interaction_counter+1;
      
      interaction_final.unique_interactions(Unique_interaction_counter,1)= unique_interaction(interaction_counter2A);
      %check size of replicate_numbers
      interaction_final.replicate_numbers{Unique_interaction_counter,1}=...
        binary_interaction_list{location_interaction_pairs(1),12};
      
      %copy the values for the first unique interaction to the list
      interaction_final.proteinA(Unique_interaction_counter,1) = binary_interaction_list(location_interaction_pairs(1),8);
      interaction_final.proteinB(Unique_interaction_counter,1) = binary_interaction_list(location_interaction_pairs(1),9);
      interaction_final.CenterA(Unique_interaction_counter,1) = binary_interaction_list(location_interaction_pairs(1),10);
      interaction_final.CenterB(Unique_interaction_counter,1) = binary_interaction_list(location_interaction_pairs(1),11);
      interaction_final.DeltaHeight(Unique_interaction_counter,1) = binary_interaction_list(location_interaction_pairs(1),4);
      interaction_final.DeltaCenter(Unique_interaction_counter,1) = binary_interaction_list(location_interaction_pairs(1),3);
      interaction_final.DeltaWidth(Unique_interaction_counter,1) = binary_interaction_list(location_interaction_pairs(1),5);
      interaction_final.DeltaEucDist(Unique_interaction_counter,1) = binary_interaction_list(location_interaction_pairs(1),2);
      interaction_final.proteinInCorum(Unique_interaction_counter,1) = binary_interaction_list(location_interaction_pairs(1),6);
      interaction_final.interactionInCorum(Unique_interaction_counter,1) = binary_interaction_list(location_interaction_pairs(1),7);
      interaction_final.score(Unique_interaction_counter,row) = binary_interaction_list{location_interaction_pairs(1),13};
    end
  end
  
  tt = toc;
  fprintf('  ...  %.2f seconds\n',tt)
  
  
  
  % Calculate final TP, FP, FN and TN
  tic
  fprintf('        9. Calculate final TP, FP, FN and TN')
  % True interactions Observation across replicate
  % Count how many times a true positive interactions was observed
  Total_unique_interactions=length(interaction_final.unique_interactions);
  number_observation_pre_interaction=zeros(Total_unique_interactions,1);
  number_unique_interaction=zeros(Total_unique_interactions,1);
  %Record the number of times a interaction was detected in irrespective of replicate number
  interactions_prep_replicate=zeros(Total_unique_interactions,replicate_counter);
  for counter=1:Total_unique_interactions
    
    %convert cell array to double
    a = interaction_final.replicate_numbers(counter,:);
    ix = ~cellfun('isempty',a);
    testA = [a{ix}];
    %Save values
    number_observation_pre_interaction(counter) = nnz(testA);
    
    %Determine if multiple gaussians within the same replicate were indentified
    unique_testA = unique(testA);
    number_of_unique = length(unique_testA);
    number_unique_interaction(counter) = number_of_unique;
    
    %Which replicate is the interaction in
    for replicate_counter = 1:(number_of_replicates*number_of_channels)
      interactions_prep_replicate(counter,replicate_counter) = nnz(find(testA(:)==replicate_counter));
    end
  end
  
  %Copy number of unique_interaction across replicate to interaction array
  interaction_final.times_observed_across_reps = number_unique_interaction;
  
  %total number of interactions of unique replicates
  total_observation=sum(number_observation_pre_interaction);
  %total number of interactions of unique replicates
  total_unique_observation=sum(number_unique_interaction);
  
  
  % create unique identifier for false interactions
  % ##### THIS CODE IS AN APPROXIMATION!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  % ##### If protein A and B are not interacting in corum, this sets ALL THE GAUSSIAN COMBINATIONS
  %       between A and B to not interacting. See Dec 3 notes for more info.
  %       The proper fix is to individually treat all Gaussian combinations exactly like we do for
  %       the positive interactions.
  interaction_FN_pairs = FN_binary_interaction_list(:,1) + FN_binary_interaction_list(:,2);
  length_unique_FN_inter = length(unique(interaction_FN_pairs));
  
  %Determine the number of unique false negatives
  interaction_as_number = cat(1,FN_binary_interaction_list(:,3));
  interaction_FN_incorum_pairs = interaction_FN_pairs(interaction_as_number==1);
  length_unique_FN_incorum = length(unique(interaction_FN_incorum_pairs));
    
  %clear varibles to speed up processing
  clearvars interaction_FN_incorum_pairs;
  clearvars interaction_FN_pairs;
  clearvars FN_binary_interaction_list;
  clear interaction_as_number
  
  % Caliculate TP, FP, FN and TN
  %Caliculate values from true interaction data set
  True_proteinInCorum=sum([cell2mat(interaction_final.proteinInCorum)]);
  True_interactionInCorum=sum([cell2mat(interaction_final.interactionInCorum)]);
  
  %Caliculate values from true interaction data set
  False_proteinInCorum=length_unique_FN_inter;
  False_interactionInCorum=length_unique_FN_incorum;
  
  %caliculate TP, FP, FN and TN
  final_TN=False_proteinInCorum-False_interactionInCorum;
  final_FP=True_proteinInCorum-True_interactionInCorum;
  final_FN=False_interactionInCorum;
  final_TP=True_interactionInCorum;
  
  % final results
  final_Recall = final_TP/(final_TP+final_TN);
  final_Precision = final_TP/(final_TP+final_FP);
  final_TPR = final_TP/(final_TP+final_FN);
  final_FPR =final_FP/(final_FP+final_TN);
  
  
  %Format data to be written out
  %Crete list of centerA, note ensure strjoin function is avalible
  interaction_final.centerA_formated=cell(Total_unique_interactions,1);
  for format_loop=1:Total_unique_interactions
    centerA_2bformated1=interaction_final.CenterA(format_loop,:);
    centerA_2bformated1(cellfun('isempty',centerA_2bformated1)) = [];
    length_centerA=length(centerA_2bformated1);
    
    %Test if the array is longer then one entry
    if length_centerA<2
      interaction_final.centerA_formated(format_loop)=centerA_2bformated1(1);
    elseif length_centerA>=2
      for jj= 1:length(centerA_2bformated1)
        centerA_2bformated1{jj} = num2str(centerA_2bformated1{jj});
      end
      interaction_final.centerA_formated(format_loop)= {strjoin(centerA_2bformated1,' ; ')};
    end
  end
  
  %Crete list of centerB, note ensure strjoin function is avalible
  interaction_final.centerB_formated=cell(Total_unique_interactions,1);
  for format_loop=1:Total_unique_interactions
    centerB_2bformated1=interaction_final.CenterB(format_loop,:);
    centerB_2bformated1(cellfun('isempty',centerB_2bformated1)) = [];
    length_centerB=length(centerB_2bformated1);
    
    %Test if the array is longer then one entry
    if length_centerB<2
      interaction_final.centerB_formated(format_loop)=centerB_2bformated1(1);
    elseif length_centerB>=2
      for jj= 1:length(centerB_2bformated1)
        centerB_2bformated1{jj} = num2str(centerB_2bformated1{jj});
      end
      interaction_final.centerB_formated(format_loop)= {strjoin(centerB_2bformated1,' ; ')};
    end
  end
  
  %Crete list of DeltaEuc, note ensure strjoin function is avalible
  interaction_final.DeltaEuc_formated=cell(Total_unique_interactions,1);
  for format_loop=1:Total_unique_interactions
    DeltaEuc_2bformated1=interaction_final.DeltaEucDist(format_loop,:);
    DeltaEuc_2bformated1(cellfun('isempty',DeltaEuc_2bformated1)) = [];
    length_DeltaEuc=length(DeltaEuc_2bformated1);
    %Test if the array is longer then one entry
    if length_DeltaEuc<2
      interaction_final.DeltaEuc_formated(format_loop)=DeltaEuc_2bformated1(1);
    elseif length_DeltaEuc>=2
      for jj= 1:length(DeltaEuc_2bformated1)
        DeltaEuc_2bformated1{jj} = num2str(DeltaEuc_2bformated1{jj});
      end
      interaction_final.DeltaEuc_formated(format_loop)= {strjoin(DeltaEuc_2bformated1,' ; ')};
    end
  end
  
  %Crete list of DeltaCenter, note ensure strjoin function is avalible
  interaction_final.DeltaCenter_formated=cell(Total_unique_interactions,1);
  for format_loop=1:Total_unique_interactions
    DeltaCenter_2bformated1=interaction_final.DeltaCenter(format_loop,:);
    DeltaCenter_2bformated1(cellfun('isempty',DeltaCenter_2bformated1)) = [];
    length_DeltaCenter=length(DeltaCenter_2bformated1);
    
    %Test if the array is longer then one entry
    if length_DeltaCenter<2
      interaction_final.DeltaCenter_formated(format_loop)=DeltaCenter_2bformated1(1);
    elseif length_DeltaCenter>=2
      for jj= 1:length(DeltaCenter_2bformated1)
        DeltaCenter_2bformated1{jj} = num2str(DeltaCenter_2bformated1{jj});
      end
      interaction_final.DeltaCenter_formated(format_loop)= {strjoin(DeltaCenter_2bformated1,' ; ')};
    end
  end
  
  %Crete list of DeltaHeight, note ensure strjoin function is avalible
  interaction_final.DeltaHeight_formated=cell(Total_unique_interactions,1);
  
  for format_loop=1:Total_unique_interactions
    DeltaHeight_2bformated1=interaction_final.DeltaHeight(format_loop,:);
    DeltaHeight_2bformated1(cellfun('isempty',DeltaHeight_2bformated1)) = [];
    length_DeltaHeight=length(DeltaHeight_2bformated1);
    
    %Test if the array is longer then one entry
    if length_DeltaHeight<2
      interaction_final.DeltaHeight_formated(format_loop)=DeltaHeight_2bformated1(1);
    elseif length_DeltaHeight>=2
      for jj= 1:length(DeltaHeight_2bformated1)
        DeltaHeight_2bformated1{jj} = num2str(DeltaHeight_2bformated1{jj});
      end
      interaction_final.DeltaHeight_formated(format_loop)= {strjoin(DeltaHeight_2bformated1,' ; ')};
    end
  end
  
  %Crete list of Deltawidth, note ensure strjoin function is avalible
  interaction_final.Deltawidth_formated=cell(Total_unique_interactions,1);
  for format_loop=1:Total_unique_interactions
    Deltawidth_2bformated1=interaction_final.DeltaWidth(format_loop,:);
    Deltawidth_2bformated1(cellfun('isempty',Deltawidth_2bformated1)) = [];
    length_Deltawidth=length(Deltawidth_2bformated1);
    %Test if the array is longer then one entry
    if length_Deltawidth<2
      interaction_final.Deltawidth_formated(format_loop)=Deltawidth_2bformated1(1);
    elseif length_Deltawidth>=2
      for jj= 1:length(Deltawidth_2bformated1)
        Deltawidth_2bformated1{jj} = num2str(Deltawidth_2bformated1{jj});
      end
      interaction_final.Deltawidth_formated(format_loop)= {strjoin(Deltawidth_2bformated1,' ; ')};
    end
  end
  
  %Crete list of replicates, note ensure strjoin function is avalible
  interaction_final.Replicates_formated=cell(Total_unique_interactions,1);
  for format_loop=1:Total_unique_interactions
    Replicates_2bformated1=interaction_final.replicate_numbers(format_loop,:);
    Replicates_2bformated1(cellfun('isempty',Replicates_2bformated1)) = [];
    length_Replicates=length(Replicates_2bformated1);
    %Test if the array is longer then one entry
    if length_Replicates<2
      interaction_final.Replicates_formated{format_loop}=mat2str(Replicates_2bformated1{1});
    elseif length_Replicates>=2
      for jj= 1:length(Replicates_2bformated1)
        Replicates_2bformated1{jj} = num2str(Replicates_2bformated1{jj});
      end
      interaction_final.Replicates_formated{format_loop}= strjoin(Replicates_2bformated1,' ; ');
    end
  end
  
  %Crete list of scores, note ensure strjoin function is avalible
  interaction_final.scores_formated=cell(Total_unique_interactions,1);
  for format_loop=1:Total_unique_interactions
    tmp = interaction_final.score(format_loop,:);
    tmp(isnan(tmp)) = [];
    length_Replicates=length(tmp);
    %Test if the array is longer then one entry
    if length_Replicates<2
      interaction_final.scores_formated{format_loop}=mat2str(tmp(1));
    elseif length_Replicates>=2
      Replicates_2bformated1 = cell(length_Replicates,1);
      for jj= 1:length(tmp)
        Replicates_2bformated1{jj} = num2str(tmp(jj));
      end
      interaction_final.scores_formated{format_loop} = strjoin(Replicates_2bformated1,' ; ');
    end
  end
  
  for ii = 1:Total_unique_interactions
    tmp = interaction_final.score(ii,:);
    tmp = tmp(~isnan(tmp));
    interaction_final.scores_formated{ii} = strjoin(cellstr(num2str(tmp(:)))',' ; ');
  end
  
  
  %create figure of precision replicates over replicates
  %create array to fill with number of interaction with both proteins in
  %corum and the number of interactions for these proteins in corum
  Precision_array=zeros((number_of_replicates*number_of_channels),3);
  Precision_redundant=zeros((number_of_replicates*number_of_channels),3);
  
  %Determine the precision for interactions observed across replicates
  % Precision_redundant(:,1) = protein in corum
  % Precision_redundant(:,2) = interaction in corum
  for precision_counter = 1:(number_of_replicates*number_of_channels)
    %Find position of interaction observed within the required number of
    %replicates
    position_to_sum=find(interaction_final.times_observed_across_reps==precision_counter);
    
    %Count the number values for non redundant list
    Precision_array(precision_counter,1) = sum(cell2mat(interaction_final.proteinInCorum(position_to_sum)));
    Precision_array(precision_counter,2) = sum(cell2mat(interaction_final.interactionInCorum(position_to_sum)));
    Precision_array(precision_counter,3) = length(position_to_sum);
    %Count the number values for redundant list
    Precision_redundant(precision_counter,1) = (sum(cell2mat(interaction_final.proteinInCorum(position_to_sum)))*precision_counter);
    Precision_redundant(precision_counter,2) = (sum(cell2mat(interaction_final.interactionInCorum(position_to_sum)))*precision_counter);
    Precision_redundant(precision_counter,3) = (length(position_to_sum)*precision_counter);
  end
  
  %Generation precision as a precentage for non redundant interactions
  Precent_precision_R=zeros((number_of_replicates*number_of_channels),1);
  for precision_counter=1:(number_of_replicates*number_of_channels)
    Precent_precision_R(precision_counter)=Precision_array(precision_counter,2)/Precision_array(precision_counter,1);
  end
  
  %Generate precision as a precentage for redundant interactions
  Precent_precision_NR=sum(Precision_redundant(:,2))/sum(Precision_redundant(:,1));
  
  tt = toc;
  fprintf('  ...  %.2f seconds\n',tt)
  
  
  
  tic
  fprintf('        10. Interactions observed only under treatment')
  Observed_treatment=zeros(Total_unique_interactions,1);
  Observed_untreatment=zeros(Total_unique_interactions,1);
  treatment_specific_interaction=zeros(Total_unique_interactions,1);
  untreatment_specific_interaction=zeros(Total_unique_interactions,1);
  
  %Determine If there are condition specific interactions detected
  for treatment_counter=1:Total_unique_interactions
    
    Observed_treatment(treatment_counter,1)=nnz(interactions_prep_replicate(treatment_counter,treatment_replicates));
    Observed_untreatment(treatment_counter,1)=nnz(interactions_prep_replicate(treatment_counter,untreatment_replicates));
    
    %Determine if there are treatment specific interactions detected
    if Observed_treatment(treatment_counter,1)>=2 & Observed_untreatment(treatment_counter,1)==0
      treatment_specific_interaction(treatment_counter,1)=1;
    end
    
    %Determine if there are untreatment specific interactions detected
    if Observed_untreatment(treatment_counter,1)>=2 & Observed_treatment(treatment_counter,1)==0
      untreatment_specific_interaction(treatment_counter,1)=1;
    end
    
  end
  
  %Copy treatment specific interaction to list to write out
  treatment_specific.unique_interactions=cell(0,0);
  treatment_specific.proteinA=cell(0,0);
  treatment_specific.proteinB=cell(0,0);
  treatment_specific.centerA_formated=cell(0,0);
  treatment_specific.centerB_formated=cell(0,0);
  treatment_specific.Replicates_formated=cell(0,0);
  treatment_specific.DeltaHeight_formated=cell(0,0);
  treatment_specific.DeltaCenter_formated=cell(0,0);
  treatment_specific.Deltawidth_formated=cell(0,0);
  treatment_specific.DeltaEuc_formated=cell(0,0);
  treatment_specific.proteinInCorum=cell(0,0);
  treatment_specific.interactionInCorum=cell(0,0);
  
  %Create counter
  treatment_row_counter=1;
  for treatment_counter=1:Total_unique_interactions
    if  treatment_specific_interaction(treatment_counter,1)==1;
      treatment_specific.unique_interactions(treatment_row_counter,1)=interaction_final.unique_interactions(treatment_counter,1);
      treatment_specific.proteinA(treatment_row_counter,1)=interaction_final.proteinA(treatment_counter,1);
      treatment_specific.proteinB(treatment_row_counter,1)=interaction_final.proteinB(treatment_counter,1);
      treatment_specific.centerA_formated(treatment_row_counter,1)=interaction_final.centerA_formated(treatment_counter,1);
      treatment_specific.centerB_formated(treatment_row_counter,1)=interaction_final.centerB_formated(treatment_counter,1);
      treatment_specific.Replicates_formated(treatment_row_counter,1)=interaction_final.Replicates_formated(treatment_counter,1);
      treatment_specific.DeltaHeight_formated(treatment_row_counter,1)=interaction_final.DeltaHeight_formated(treatment_counter,1);
      treatment_specific.DeltaCenter_formated(treatment_row_counter,1)=interaction_final.DeltaCenter_formated(treatment_counter,1);
      treatment_specific.Deltawidth_formated(treatment_row_counter,1)=interaction_final.Deltawidth_formated(treatment_counter,1);
      treatment_specific.DeltaEuc_formated(treatment_row_counter,1)=interaction_final.DeltaEuc_formated(treatment_counter,1);
      treatment_specific.proteinInCorum(treatment_row_counter,1)=interaction_final.proteinInCorum(treatment_counter,1);
      treatment_specific.interactionInCorum(treatment_row_counter,1)=interaction_final.interactionInCorum(treatment_counter,1);
      treatment_row_counter=1+treatment_row_counter;
    end
  end
  
  %Determine the number of interactions
  treament_specific_interactionsnum=length(treatment_specific.unique_interactions);
  
  
  %Copy treatment specific interaction to list to write out
  untreatment_specific.unique_interactions=cell(0,0);
  untreatment_specific.proteinA=cell(0,0);
  untreatment_specific.proteinB=cell(0,0);
  untreatment_specific.centerA_formated=cell(0,0);
  untreatment_specific.centerB_formated=cell(0,0);
  untreatment_specific.Replicates_formated=cell(0,0);
  untreatment_specific.DeltaHeight_formated=cell(0,0);
  untreatment_specific.DeltaCenter_formated=cell(0,0);
  untreatment_specific.Deltawidth_formated=cell(0,0);
  untreatment_specific.DeltaEuc_formated=cell(0,0);
  untreatment_specific.proteinInCorum=cell(0,0);
  untreatment_specific.interactionInCorum=cell(0,0);
  
  %Create counter
  untreatment_row_counter=1;
  for treatment_counter=1:Total_unique_interactions
    if  untreatment_specific_interaction(treatment_counter,1)==1;
      untreatment_specific.unique_interactions(untreatment_row_counter,1)=interaction_final.unique_interactions(treatment_counter,1);
      untreatment_specific.proteinA(untreatment_row_counter,1)=interaction_final.proteinA(treatment_counter,1);
      untreatment_specific.proteinB(untreatment_row_counter,1)=interaction_final.proteinB(treatment_counter,1);
      untreatment_specific.centerA_formated(untreatment_row_counter,1)=interaction_final.centerA_formated(treatment_counter,1);
      untreatment_specific.centerB_formated(untreatment_row_counter,1)=interaction_final.centerB_formated(treatment_counter,1);
      untreatment_specific.Replicates_formated(untreatment_row_counter,1)=interaction_final.Replicates_formated(treatment_counter,1);
      untreatment_specific.DeltaHeight_formated(untreatment_row_counter,1)=interaction_final.DeltaHeight_formated(treatment_counter,1);
      untreatment_specific.DeltaCenter_formated(untreatment_row_counter,1)=interaction_final.DeltaCenter_formated(treatment_counter,1);
      untreatment_specific.Deltawidth_formated(untreatment_row_counter,1)=interaction_final.Deltawidth_formated(treatment_counter,1);
      untreatment_specific.DeltaEuc_formated(untreatment_row_counter,1)=interaction_final.DeltaEuc_formated(treatment_counter,1);
      untreatment_specific.proteinInCorum(untreatment_row_counter,1)=interaction_final.proteinInCorum(treatment_counter,1);
      untreatment_specific.interactionInCorum(untreatment_row_counter,1)=interaction_final.interactionInCorum(treatment_counter,1);
      untreatment_row_counter=1+untreatment_row_counter;
    end
  end
  
  %Determine the number of interactions
  untreament_specific_interactionsnum=length(untreatment_specific.unique_interactions);
  
  tt = toc;
  fprintf('  ...  %.2f seconds\n',tt)
  
  
  
  tic
  fprintf('        11. Write output files')
  writeOutput_ROC
  tt = toc;
  fprintf('  ...  %.2f seconds\n',tt)
  
  tic
  fprintf('        12. Make figures')
  makeFigures_ROC
  tt = toc;
  fprintf('  ...  %.2f seconds\n',tt)
  
end




