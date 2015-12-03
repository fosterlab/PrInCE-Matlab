%%%%%%%%%%%%%%% Instructions for running Gauss_Build.m:
%
%
%
%%%%%%%%%%%%%%% Logic:
% 0. Initalize
% 1. Read input data
% 2. Pre-process gaussians, chromatograms. Make Protein structure.
% 3. Make Protein struct.
% 4. Make list of interactions common b/w sample and CORUM, expand for protein groups.
% 5. Make TP matrix.
% 6.
%
%
%%%%%%%%%%%%%%% Custom functions called by Gauss_build.m:
% -
%
%
%%%%%%%%%%%%%%% To do:
% - Get all the input/output files sorted out.
% - Why is this script pre-processing chromatograms? Shouldn't that have been done earlier?
%
%
%%%%%%%%%%%%%%% Fix the logic:
% - Normalized Gaussians were wrong. Missing the 2 in the denominator.
% - 'Look up Major Protein group Information' can be seriously improved with a few findstr
%
%
%%%%%%%%%%%%%%% Bugs in my (this) code:
% - My Corum_Dataset is 11516x1. Nick's Corum_Dataset is 1x16393. Why are they different sizes?



%% 0. Initialize

% user defined settings
number_of_replicates=3;
number_of_channels=2;
replicate_counter = 4; % Do a single replicate. This is Nick's loop index.

% Define folders, i.e. define where everything lives.
maindir = '/Users/Mercy/Academics/Foster/NickCodeData/GregPCP-SILAC/'; % where everything lives
codedir = [maindir 'Code/']; % where this script lives
funcdir = [maindir 'Code/Functions/']; % where small pieces of code live
datadir = [maindir 'DataFiles/']; % where data files live
figdir = [maindir 'Figures/']; % where figures live
tmpdir = '/Users/Mercy/Academics/Foster/NickCodeData/4B_ROC_homologue_DB/Combined Cyto new DB/';

% List all input files. These contain data that will be read by this script.
InputFile{1} = [datadir 'Input/Major_protein_groups.xlsx'];
InputFile{2} = [datadir 'Input/Corum_correctly_formated_Uniprot_IDs.csv'];
%InputFile{3} = [datadir 'Input/Corum_2012_human.xlsx'];
InputFile{4} = [datadir 'Input/Master_guassian_list.csv'];

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



%% 1. Read input data

[~, Protein_IDs] = xlsread(InputFile{1});

% Corum binary interactions
fid=fopen(InputFile{2}, 'rt');    %Input corum data base.
Corum_Import= textscan(fid,'%s\t', 'Delimiter',',');
fclose(fid);
No=length(Corum_Import{1})/2;
Corum_Protein_names = (reshape(Corum_Import{1,1},2,No)');
Unique_Corum = unique(Corum_Protein_names);

Master_Gaussian_list=importdata(InputFile{4},',');

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



%% 2. Pre-process gaussians, chromatograms. Get Dist
% - In chromatograms, replaces nans with 0.05
% - Remove gaussians and chromatograms with C<5
% - Make normalized gaussians, i.e. set height to 1

% Pre-process data to remove Gaussian above fraction five
% Create array to use for comparsions
H_raw = cellfun(@str2num,Gaus_import_reshaped(2:end,2));
C_raw = cellfun(@str2num,Gaus_import_reshaped(2:end,3));
W_raw = cellfun(@str2num,Gaus_import_reshaped(2:end,4));

% Replace NaN values with 0.05
Chromatograms_raw(isnan(Chromatograms_raw))= 0.05;

% remove gaussian with centers below five
Ibad = C_raw<5;
H = H_raw(~Ibad);
C = C_raw(~Ibad);
W = W_raw(~Ibad);
Chromatograms = Chromatograms_raw(~Ibad,:);

% Make normalized gaussians
x = 0.25:0.25:100;
gaussian_fit=zeros(length(H),length(x));
for ii = 1:length(H)
  gaussian_fit(ii,:) = 1*exp(-(x- C(ii)).^2 /2/(W(ii).^2));
end

%Calculate distances
Dist.Euc = squareform(pdist(Chromatograms,'euclidean'));
Dist.Center = squareform(pdist(C,'euclidean'));
Dist.Height = squareform(pdist(H,'euclidean'));
Dist.Width = squareform(pdist(W,'euclidean'));
Dist.Gaussian_fits= squareform(pdist(gaussian_fit,'euclidean'));

% Make Dist.R^2
Dist.R2 = zeros(Nprot,Nprot);
for ii=1:length(H)
  for jj=1:length(H)
    R = corrcoef(Chromatograms(ii,:),Chromatograms(jj,:));
    Dist.R2(ii,jj) = 1-R(1,2)^2;
  end
end



%% 3. Make Protein structure
% This summarizes the proteins in our sample

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
for Gaus_import_counter1 = 1:Dimensions_Gaus_import(1)
  Lookup_protein_name=Protein.Isoform(Gaus_import_counter1);
  disp(Lookup_protein_name)
  
  tmp = find(strcmp(Lookup_protein_name,Protein_IDs));
  [Check_counter1, Check_counter2] = ind2sub(size(Protein_IDs),tmp);
  
  %copy names to MajorID from Protein_IDs raw data
  Array_to_check1 = Protein_IDs(Check_counter1,:);
  Array_to_check_empty1 = cellfun(@isempty,Array_to_check1);
  for Check_counter3 = 1:nnz(~Array_to_check_empty1)
    Protein.MajorIDs{Gaus_import_counter1,Check_counter3}={Protein_IDs{Check_counter1,Check_counter3}};
  end
  
  %copy names to MajorID_NoIsoforms from Protein_IDS_no_isoform_no_dup_processed
  Array_to_check2=Protein_IDS_no_isoform_no_dup_processed(Check_counter1,:);
  Array_to_check_empty2=cellfun(@isempty,Array_to_check2);
  for Check_counter4 = 1:nnz(~Array_to_check_empty2)
    a = Protein_IDS_no_isoform_no_dup_processed{Check_counter1,Check_counter4};
    Protein.MajorID_NoIsoforms{Gaus_import_counter1,Check_counter4}=a{1};
  end
end




%% 4. Make true positive matrix

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
  
  Prot1_group = Protein.MajorID_NoIsoforms(i1(1),:); % all proteins in the same group as Prot1
  Prot1_group = Prot1_group(~cellfun('isempty',Prot1_group)); % remove empty cells
  
  Prot2_group = Protein.MajorID_NoIsoforms(i2(1),:); % all proteins in the same group as Prot2
  Prot2_group = Prot2_group(~cellfun('isempty',Prot2_group)); % remove empty cells
  
  TP_Matrix(i1,i2) = TP_Matrix(i1(1),i2(1))+1;
  
  for jj = 1:length(Prot1_group)
    for kk = 1:length(Prot2_group)
      ii = ii+1;
      a = Prot1_group{jj};
      b = Prot2_group{kk};
      Corum_Dataset_expanded{ii} = [a,'-',b];
      
      %disp(Corum_Dataset_expanded{ii})
    end
  end
end
%I = unique(cell2mat(Corum_Dataset_expanded),'rows');
%Corum_Dataset_expanded = unique(Corum_Dataset_expanded{ii});



%% 5. Make True Positive matrix

% check which individual proteins are in corum
% Important! Also check whether any member of that protein's group is in corum
inCorum = zeros(1,Dimensions_Gaus_import);
isemptyCellArray = ~cellfun('isempty',Protein.MajorID_NoIsoforms);
for ii = 1:Dimensions_Gaus_import
  Ngroup = sum(isemptyCellArray(ii,:));
  for ni = 1:Ngroup
    if ismember(Protein.MajorID_NoIsoforms{ii,ni},Unique_Corum)
      inCorum(ii) = inCorum(ii)+1;
    end
  end
end
inCorum = find(inCorum>0);
%
% TP_Matrix2 = zeros(Dimensions_Gaus_import,Dimensions_Gaus_import);
% %Multiple_known_interaction_in_proteingroup_Matrix= zeros(Dimensions_Gaus_import,Dimensions_Gaus_import);
%
% for ii = inCorum
%   ii
%   for jj = inCorum
%     Interaction1 = strcat(Protein.MajorID_NoIsoforms{ii,1},'-',Protein.MajorID_NoIsoforms{jj,1});
%     Interaction2 = strcat(Protein.MajorID_NoIsoforms{jj,1},'-',Protein.MajorID_NoIsoforms{ii,1});
%
%     % is this interaction in the expanded Corum_dataset list?
%     b1 = sum(strcmp(Interaction1,Corum_Dataset_expanded));
%     b2 = sum(strcmp(Interaction1,Corum_Dataset_expanded));
%
%     if b1>0 || b2>0
%       TP_Matrix2(ii,jj) = 1;
%     end
%   end
% end



%% 6. Make other matrices
% - Possible interaction matrix
% - Self interaction matrix

% Calculate possible interaction matrix, asks are both proteins in corum
Int_matrix3 = zeros(Dimensions_Gaus_import,Dimensions_Gaus_import);
% Make matrix of protein group with multiple possible positive interactions
Multiple_matches_within_protein_group_Matrix = zeros(Dimensions_Gaus_import,Dimensions_Gaus_import);
for i=inCorum
  i %Display protein number the loop is up to to enable users to monitor progress
  for i_counter=1:Number_MajorIDs(i)
    if ismember(Protein.MajorID_NoIsoforms{i,i_counter}, Pos_Corum_proteins)
      for ii=inCorum
        for ii_counter=1:Number_MajorIDs(ii)
          if ismember(Protein.MajorID_NoIsoforms{ii,ii_counter}, Pos_Corum_proteins)
            Int_matrix3(i,ii) = Int_matrix3(i,ii)+1;
          end
        end
      end
    end
  end
end

% Create self interaction matrix
Self_Int_matrix= zeros(Dimensions_Gaus_import,Dimensions_Gaus_import);
for i=1:Dimensions_Gaus_import
  i
  for ii=1:Dimensions_Gaus_import
    Protein.NoIsoform;
    NS=isequal(Protein.NoIsoform{i},Protein.NoIsoform{ii});
    if NS==1
      Self_Int_matrix(i,ii)=1;
    end
  end
end
inverse_self=~Self_Int_matrix; %creation of inverse of self interactions

% FP matrix
Inverse_TP_matrix=~TP_Matrix;
FP_Matrix= (Int_matrix & Inverse_TP_matrix);


%% 7. Make ROC curves for different Dist fieldnames

fn = fieldnames(Dist);


%cd('ROC analysis related Figures');
for i=1:length(fn)
  i
  
  count=0;
  nr=50;
  Recall_out.(fn{i})=zeros(nr,1);
  Precision_out.(fn{i})=zeros(nr,1);
  TPR_out.(fn{i})=zeros(nr,1);
  FPR_out.(fn{i})=zeros(nr,1);
  
  xDist = linspace(0,max(Dist.(fn{i})(:)),nr);
  
  for Distance = xDist %Loop to calculate performance of variables  %for Distance= 0.1:0.1:50 %Loop to calculate performance of variables
    MatrixGaussian = Dist.(fn{i}) < Distance;
    count=1+count;
    
    %As Euclidean distance represents a measurement across the total
    %graident filter out measurements of guassian< 2 fractions away
    if isequal(fn(i), fn(1))
      MatrixGaussian= (MatrixGaussian & (Dist.(fn{2}) < 2));
    else
      MatrixGaussian = Dist.(fn{i}) < Distance;
    end
    
    %Total number of interaction within distance
    Total_Int = (MatrixGaussian & inverse_self);
    Total_Int2_Triu =triu(Total_Int); %Only consider interactions in the top right of the matrix
    Total_Int_number= length(find(Total_Int2_Triu==1));
    
    %Protein interactions of protein in Corum within distance
    Possible_Int1 = (MatrixGaussian & Int_matrix);
    Possible_Int2 = (Possible_Int1 & inverse_self);
    Possible_Int2_Triu = triu(Possible_Int2); %Only consider interactions in the top right of the matrix
    Possible_Int_number= length(find(Possible_Int2_Triu==1));
    
    % Calculate True possitives
    TP_Matrix_Int1 = (Possible_Int1 & TP_Matrix);
    TP_Matrix_Int2 = (TP_Matrix_Int1 & inverse_self);
    TP_Matrix_Int2_Triu = triu(TP_Matrix_Int2); %Only consider interactions in the top right of the matrix
    TP= length(find(TP_Matrix_Int2_Triu==1));
    
    % Calculate false negatives
    Inverse_TP_matrix1=~TP_Matrix_Int2;
    MaxTP_matrix1 = (TP_Matrix & Inverse_TP_matrix1);
    MaxTP_matrix2 = (MaxTP_matrix1 & inverse_self);
    MaxTP_matrix2_Triu = triu(MaxTP_matrix2); %Only consider interactions in the top right of the matrix
    FN= length(find(MaxTP_matrix2_Triu==1));
    
    % Calculate false positives
    FP_Matrix_Int1= (Possible_Int2 & FP_Matrix);
    FP_Matrix_Int2= (FP_Matrix_Int1 & inverse_self);
    FP_Matrix_Int2_Triu = triu(FP_Matrix_Int2); %Only consider interactions in the top right of the matrix
    FP= length(find(FP_Matrix_Int2_Triu==1));
    
    % Find True negatives
    Inverse_false_postive_matrix=~FP_Matrix_Int1;
    Int_matrix_minus_TP = (Int_matrix & Inverse_TP_matrix);
    Int5_matrix= (Int_matrix_minus_TP & Inverse_false_postive_matrix);
    Int6_matrix = (Int5_matrix & inverse_self);
    Int6_matrix_Triu = triu(Int6_matrix); %Only consider interactions in the top right of the matrix
    TN= (length(find(Int6_matrix_Triu==1)));
    
    Recall = TP/(TP+FN);
    Precision = TP/(TP+FP);
    TPR = TP/(TP+FN);
    FPR =FP/(FP+TN);
    
    Recall_out.(fn{i})(count,1)= Recall;
    Precision_out.(fn{i})(count,1)=Precision;
    TPR_out.(fn{i})(count,1)=TPR;
    FPR_out.(fn{i})(count,1)=FPR;
  end
end


% plot all this
if 1
  cols = {'b' 'g' 'c' [.5 .5 .5] 'r' 'm' [1 .5 0]};
  
  figure,hold on
  for i=1:length(fn)
    plot(FPR_out.(fn{i}),TPR_out.(fn{i}),'color',cols{i})
  end
  plot([0 1],[0 1],'--k')
  axis([0 1 0 1])
  legend('Euc','Center','Height','Width','Gauss fits','R2','location','northwest')
  xlabel('FPR')
  ylabel('TPR')
  
  figure,hold on
  for i=1:length(fn)
    plot(Recall_out.(fn{i}),Precision_out.(fn{i}),'color',cols{i})
  end
  %axis([0 1 0 1])
  legend('Euc','Center','Height','Width','Gauss fits','R2','location','northeast')
  xlabel('Recall')
  ylabel('Precision')
end



%% 8.

xDist = linspace(0,max(Dist.Euc(:)),nr);

% Make interactions_based_only_Ecldian_distance at 80% precision
prec_thresh = 0.80;
Precision = Precision_out.Euc; % Precision_out.Euc

%Test if the Euclidian distance ever gives a precision above 0.8
[~, idx] = min(abs(Precision - prec_thresh));
if max(Precision) < prec_thresh;idx=1;end

%Determine interactions within precision of >0.80 based on Eu distance
goodDist = xDist(idx); %determines the Eucldian distance which gives a precision closest to 0.8
Interactions_detected_due_to_High_Eucldian_distance = Dist.Euc < goodDist;
Interactions_Based_only_Eucldian_distance = (Interactions_detected_due_to_High_Eucldian_distance & inverse_self & (Dist.Center < 2));

%create matrix for Eu dis interactions with precision >0.8
Final_Interactions_Based_only_Eucldian_distance = Interactions_Based_only_Eucldian_distance;


% Make optimisation_matrix_precision
if 0
  Height_Distance = 20.0; % As height has very little effect on interaction determination, this has been set very large
  Width_Distangce = 20.0; % As width has very little effect on interaction determination, this has been set very large
  MatrixHeight = Dist.Height < Height_Distance;
  MatrixWidth = Dist.Width < Width_Distangce;
  %Determine interactions for a FPR of 0.5%
  [~, idx] = min(abs(FPR_out.Euc - 0.05)); %determines the postion which gives a FPR of 0.5%
  Eucldian_distance_with_0_05FPR = ceil(xDist(idx)); %determines the Eucldian distance which gives a FPR of 0.5%
  %Define number of data points
  Optimisation_Dimension = Eucldian_distance_with_0_05FPR/0.025;
  optimisation_matrix_number_interaction = zeros(Optimisation_Dimension,Optimisation_Dimension);
  optimisation_matrix_precision = zeros(Optimisation_Dimension,Optimisation_Dimension);
  Euclidian_Distance_optimisation = 0.025:0.025:Eucldian_distance_with_0_05FPR;
  Gaussian_fit_Distance_optimisation = 0.025:0.025:Eucldian_distance_with_0_05FPR;
  for iii=1:Optimisation_Dimension
    for iv=1:Optimisation_Dimension
      NS1 = Euclidian_Distance_optimisation(iii); % e-6
      NS2 = Gaussian_fit_Distance_optimisation(iv); % e-6
      MatrixGaussian = Dist.Gaussian_fits < NS1; % e-2
      Matrix_Distance = Dist.Euc < NS2; % e-2
      Interactions_optimisation_counter = (MatrixGaussian & MatrixHeight & MatrixWidth & Matrix_Distance & inverse_self); % e-2
      total_interaction_number = length(find(triu(Interactions_optimisation_counter==1))); % e-2
      optimisation_matrix_number_interaction(iii,iv) = total_interaction_number;
      %Calculate true positives
      TP_Matrix_Int_optimisation = (Interactions_optimisation_counter & TP_Matrix);
      TP_optimisation = length(find(TP_Matrix_Int_optimisation==1));
      %Calculate false positives
      Int2_matrix_optimisation = (Interactions_optimisation_counter & Int_matrix);
      Int3_matrix_optimisation = (Int2_matrix_optimisation & inverse_self);
      FP_optimisation = length(find(Int3_matrix_optimisation==1))-TP_optimisation;
      %Precision
      Precision = TP_optimisation/(TP_optimisation+FP_optimisation);
      optimisation_matrix_precision(iii,iv) = Precision;
    end
  end
else
  % Takes a while! Just load it for now
  b = '/Users/Mercy/Academics/Foster/NickCodeData/4B_ROC_homologue_DB/Combined Cyto new DB/Replicate4/ROC analysis data/';

end
