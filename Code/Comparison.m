%%%%%%%%%%%%%%% Logic:
% This script tests whether there are any changes in protein interactions.
% To do this it compares chromatograms (right?) via:
%   i. t-test
%   ii. Mann-Whitney test, i.e. a rank-sum test
%   iii. visualization, i.e. making figures.
% 0. Initialize
% 1. Read input
% 2. Make HvsL_Gaussians, MvsL_Gaussians structures
% WITHIN REPLICATES (?)
%   3. Determine which Gaussians are unique to each channel
%   4. Detect 2-fold changes
%   5. Determine the trends of guassian curves
% ACROSS REPLICATES
%   6. Shared gaussian across replicate.
%   7. Use fold-change > 2 to determine if complexes change
%   8. Use mwwtest and ttest to determine change
%   9. Determine the coverage (number of datapoint observed) for reach protein of each replicate
% 10. Write out tables
% 11. Write out txt for Go enrichment
% 12. Make figures
%
%
%%%%%%%%%%%%%%% Small changes:
% 1. Make all file references and directories absolute. Remove 'cd'. Make this directory-independent!
% 2. Removed the section for Alignment_binary==0. This just adds 5 nans to either side... but this
%    is already done by cleanChromatograms.m in GaussBuild. Don't need it here.
% 3. Removed the missing value imputation. We're already using cleaned chromatograms!
% 4. Replace Adjusted_HvsL_Raw_data_maxquant_modified.xlsx with Adjusted_HvsL_Raw_data_maxquant.csv
% 5. Turn variables like Gaus_MvsL into matrices. Remove the top header.
% 6. Make section 10 (ttest + mwwtest) handle replicates better. Automate.
% 7. Fix the fprintf displays.
% 8. Move all write-to-file output to the end.
%
%
%%%%%%%%%%%%%%% Big changes:
% 1. Make the code able to handle any number of replicates. Automation!
%
%
%%%%%%%%%%%%%%% Questions:
% 1. What is 'Diltuion_factor_master_mix', i.e. "the percentage of standard mixed into samples".
%   - reference dilution factor, since the reference can be too strong
%   - reference is all 50 fractions combined into one sample
% 2. Where does Adjusted_HvsL_Raw_data_maxquant_modified.xlsx come from?
%   - This is a duplicate of the csv file (without the _modified tag). The csv file was just saved
%     as excel table, Nick says because he was playing around with how to read in data. I could
%     probably just ignore the excel and use the original csv file..?
% 3. What is 2408? 7443?
% 4. What's shown in the top Gaussian plot? How are those Gaussian picked?
% 5. What Gaussians "represent" all the replicates? How are these picked?
% 6. What happens to a Gaussian that's only in one replicate?


disp('Comparison.m')
tic


%% 0. Initialize
fprintf('\n    0. Initialize')

% User defined variables
replicate_num=3; %Number of replicates to compare
position_fraction1=6; %Denote the position of the first fraction in fraction_to_plot of aligned samples
fraction_to_plot=55; %For figure define fraction to plot and to considered, all Gaussian under this value will be ignored
%Compare the Areas of the Gaussian curves to use for downstream analysis
Diltuion_factor_master_mix=0.70; %Define the precentage of standard mixed into samples
User_Window=2; %User defined window to consider Guassian the same
SILAC_ratio_combinations={'MvsL','HvsL','HvsM'};


% Plotting variables
myC= [30/255 144/255 255/255
  255/255 215/255 0/255
  178/255 34/255 34/255
  193/255 205/255 193/255];
colour_to_use=[0.254 0.411 0.882 %Colour: Royal Blue
  135/255 206/255 250/255 %Colour: Sky Blue
  0.28 0.82 0.8 %Colour: Medium Turquoise
  205/255 92/255 92/255 %Colour: Indian Red
  178/255 34/255 34/255 %Colour: Firebrick
  255/255 69/255 0 %Colour: Orange Red
  244/255 238/255 224/255 %Colour: Honeydew 2
  106/255 90/255 205/255 %Colour: Slate Blue
  34/255 139/255 34/255 %Colour: Forest Green
  222/255 184/255 135/255 %Colour: Burlywood
  186/255 85/255 211/255 %Colour: Medium Orchid
  219/255 112/255 147/255]; %Colour: Pale Violet Red


% Define folders, i.e. define where everything lives.
maindir = '/Users/Mercy/Academics/Foster/NickCodeData/GregPCP-SILAC/'; % where everything lives
codedir = [maindir 'Code/']; % where this script lives
funcdir = [maindir 'Code/Functions/']; % where small pieces of code live
datadir = [maindir 'Data/']; % where data files live
datadir1 = [datadir 'Comparison/']; % where data files live
datadir2 = [datadir 'Comparison/forROC/']; % where data files live
figdir1 = [maindir 'Figures/Comparison/']; % where figures live
figdir2 = [maindir 'Figures/Comparison/ProteinGaussianMaps/']; % where figures live
% Make folders if necessary
if ~exist(codedir, 'dir'); mkdir(codedir); end
if ~exist(funcdir, 'dir'); mkdir(funcdir); end
if ~exist(datadir, 'dir'); mkdir(datadir); end
if ~exist(datadir1, 'dir'); mkdir(datadir1); end
if ~exist(datadir2, 'dir'); mkdir(datadir2); end
if ~exist(figdir1, 'dir'); mkdir(figdir1); end
if ~exist(figdir2, 'dir'); mkdir(figdir2); end

% List all input files. These contain data that will be read by this script.
InputFile{1}=[datadir 'Adjusted_MvsL_Raw_data_maxquant_modified.xlsx']; % from Alignment
InputFile{2}=[datadir 'Adjusted_HvsL_Raw_data_maxquant_modified.xlsx']; % from Alignment
InputFile{3}=[datadir 'Adjusted_HvsM_Raw_data_maxquant_modified.xlsx']; % from Alignment
InputFile{4}=[datadir 'Alignment/Adjusted_MvsL_Combined_OutputGaus.csv']; % from Alignment
InputFile{5}=[datadir 'Alignment/Adjusted_HvsL_Combined_OutputGaus.csv']; % from Alignment

tt = toc;
fprintf('  ...  %.2f seconds\n',tt)


%% 1. Read input
fprintf('    1. Read input')


% Import data from input files
[num_val_MvsL,txt_MvsL] = xlsread(InputFile{1}); %Import file aligned output, NB 1st column is replicate
[num_val_HvsL,txt_HvsL] = xlsread(InputFile{2}); %Import file aligned output, ''
[num_val_HvsM,txt_HvsM] = xlsread(InputFile{3}); %Import file aligned output, ''
Proteins=length(num_val_MvsL);
%Remove data point with low values
num_val_MvsL(num_val_MvsL<0.20) = nan;
num_val_HvsL(num_val_HvsL<0.20) = nan;

%load results of MvsL protein interaction script
MvsL = fopen(InputFile{4}); %This corresponds to the output from the Gaus script
Gaus_import_MvsL = textscan(MvsL, '%s', 'Delimiter',',');
fclose(MvsL);
%Reshape data for use in analysis
Dimension_MvsL=size(Gaus_import_MvsL{:});
Gaus_MvsL=reshape(Gaus_import_MvsL{:},10,(Dimension_MvsL(1)/10))';

%load results of MvsL protein interaction script
HvsL = fopen(InputFile{5}); %This corresponds to the output from the Gaus script
Gaus_import_HvsL = textscan(HvsL, '%s', 'Delimiter',',');
fclose(HvsL);
%Reshape data for use in analysis
Dimension_HvsL=size(Gaus_import_HvsL{:});
Gaus_HvsL=reshape(Gaus_import_HvsL{:},10,(Dimension_HvsL(1)/10))';

%Number of fractions
fraction_number=size(num_val_MvsL);
fraction_number(2)=fraction_number(2)-1;

%Calculate number of guassians
Gaussian_number_MvsL=length(Gaus_MvsL)-1; %minus one to remove header, corresponds to the gaussian fitted
Gaussian_number_HvsL=length(Gaus_HvsL)-1;


% Add unique identifiers to the Gaussian lists
Unique_indentifer_MvsL=cell(1,1);
Unique_indentifer_MvsL(1,1)={'Unique_identifier'}; % replicate_GaussianNumber_channel_ProteinName
Unique_indentifer_MvsL(1,2)={'Replicate_Protein_identifier'}; % replicate_ProteinName
for unique_identifer_counter1= 1:Gaussian_number_MvsL
  %convert to number from string
  Replicate_number=str2num(Gaus_MvsL{unique_identifer_counter1+1,3});
  %convert to number from string
  Index_Gaussian_number=str2num(Gaus_MvsL{unique_identifer_counter1+1,1});
  Unique_indentifer_MvsL{unique_identifer_counter1+1,1}=strcat(mat2str(Replicate_number),'_',mat2str(Index_Gaussian_number),'_MvsL_',Gaus_MvsL{unique_identifer_counter1+1,4});
  Unique_indentifer_MvsL{unique_identifer_counter1+1,2}=strcat(mat2str(Replicate_number),'_',Gaus_MvsL{unique_identifer_counter1+1,4});
end
Gaus_MvsL=[Unique_indentifer_MvsL, Gaus_MvsL]; %Add unique protein_replicate identifier

Unique_indentifer_HvsL=cell(1,1);
Unique_indentifer_HvsL(1,1)={'Unique_identifier'};
Unique_indentifer_HvsL(1,2)={'Replicate_Protein_identifier'};
for unique_identifer_counter1= 1:Gaussian_number_HvsL
  %convert to number from string
  Replicate_number=str2num(Gaus_HvsL{unique_identifer_counter1+1,3});
  %convert to number from string
  Index_Gaussian_number=str2num(Gaus_HvsL{unique_identifer_counter1+1,1});
  Unique_indentifer_HvsL{unique_identifer_counter1+1,1}=strcat(mat2str(Replicate_number),'_',mat2str(Index_Gaussian_number),'_HvsL_',Gaus_HvsL{unique_identifer_counter1+1,4});
  Unique_indentifer_HvsL{unique_identifer_counter1+1,2}=strcat(mat2str(Replicate_number),'_',Gaus_HvsL{unique_identifer_counter1+1,4});
end
Gaus_HvsL=[Unique_indentifer_HvsL, Gaus_HvsL]; %Add unique protein_replicate identifier


% Add unique identifiers to the chromatograms
%MvsL
Unique_indentifer_MvsL_maxqaunt(1,1)={'Unique_identifier'}; % replicate_ProteinName
for unique_identifer_counter2= 1:Proteins
  %convert to number from string
  Replicate_number=num_val_MvsL(unique_identifer_counter2,1);
  %convert to number from string
  Unique_indentifer_MvsL_maxqaunt{unique_identifer_counter2+1,1}=strcat(mat2str(Replicate_number),'_',txt_MvsL{unique_identifer_counter2+1,1});
end
%HvsL
Unique_indentifer_HvsL_maxqaunt(1,1)={'Unique_identifier'};
for unique_identifer_counter2= 1:Proteins
  %convert to number from string
  Replicate_number=num_val_HvsL(unique_identifer_counter2,1);
  %convert to number from string
  Unique_indentifer_HvsL_maxqaunt{unique_identifer_counter2+1,1}=strcat(mat2str(Replicate_number),'_',txt_HvsL{unique_identifer_counter2+1,1});
end
%HvsM
Unique_indentifer_HvsM_maxqaunt(1,1)={'Unique_identifier'};
for unique_identifer_counter2= 1:Proteins
  %convert to number from string
  Replicate_number=num_val_HvsM(unique_identifer_counter2,1);
  %convert to number from string
  Unique_indentifer_HvsM_maxqaunt{unique_identifer_counter2+1,1}=strcat(mat2str(Replicate_number),'_',txt_HvsM{unique_identifer_counter2+1,1});
end
%Add unique protein_replicate identifier
txt_MvsL=[Unique_indentifer_MvsL_maxqaunt, txt_MvsL];
txt_HvsL=[Unique_indentifer_HvsL_maxqaunt, txt_HvsL];
txt_HvsM=[Unique_indentifer_HvsM_maxqaunt, txt_HvsM];


%expected amount of proteins, note 1 is the height
Standard_area=fraction_number(2)*1*(1/Diltuion_factor_master_mix);


%Copy data to be used for plotting
num_val_MvsL_for_figures = num_val_MvsL;
num_val_HvsL_for_figures = num_val_HvsL;

% Clean chromatograms
for ri = 1:Proteins % loop over proteins
  num_val_MvsL(ri,:) = cleanChromatogram(num_val_MvsL(ri,:),[1 3]);
  num_val_HvsL(ri,:) = cleanChromatogram(num_val_HvsL(ri,:),[1 3]);
  num_val_HvsM(ri,:) = cleanChromatogram(num_val_HvsM(ri,:),[1 3]);
end

tt = toc;
fprintf('  ...  %.2f seconds\n',tt)



%% 2. Make HvsL_Gaussians, MvsL_Gaussians structures
fprintf('    2. Make HvsL_Gaussians, MvsL_Gaussians structures')


%Varible related to gaussians
unique_MvsL_guassians=unique(Gaus_MvsL(2:Gaussian_number_MvsL+1,2));
unique_HvsL_guassians=unique(Gaus_HvsL(2:Gaussian_number_HvsL+1,2));
All_guassians_identified_redunent_list=[unique_MvsL_guassians;unique_HvsL_guassians];
All_guassians_identified_nonredunent_list=unique(All_guassians_identified_redunent_list);
Number_All_guassians_identified_nonredunent_list=length(All_guassians_identified_nonredunent_list);

%create counter to empty number into
Number_guassian_pre_protein_MvsL=zeros(Number_All_guassians_identified_nonredunent_list,1);
Number_guassian_pre_protein_HvsL=zeros(Number_All_guassians_identified_nonredunent_list,1);
%Set to five as maximum number of fitted gaussian curves is five
Location_of_guassians_MvsL=zeros(Number_All_guassians_identified_nonredunent_list,5);
Location_of_guassians_HvsL=zeros(Number_All_guassians_identified_nonredunent_list,5);

% Create array of values grouped by Protein and replicate
for guassian_counter1 = 1:Number_All_guassians_identified_nonredunent_list
  %define protein to look for
  Protein_name=All_guassians_identified_nonredunent_list(guassian_counter1);
  
  %Determine if Protein was observed in MvsL or HvsL channel
  L1 = ismember(Gaus_MvsL(:,2),Protein_name);
  L2 = ismember(Gaus_HvsL(:,2),Protein_name);
  
  %Sum the number of Gaussians
  Number_guassian_pre_protein_MvsL(guassian_counter1,1)= sum(L1);
  Number_guassian_pre_protein_HvsL(guassian_counter1,1)= sum(L2);
  
  %find the position in MvsL
  [internal_location_guassians_MvsL]=ind2sub(size(Gaus_MvsL), strmatch(Protein_name, Gaus_MvsL(:,2), 'exact'));
  number_counter1=length(internal_location_guassians_MvsL);
  counter_balance1=0;
  for internal_counter1 = 1:number_counter1
    if str2num(Gaus_MvsL{internal_location_guassians_MvsL(internal_counter1),8})>fraction_to_plot
      counter_balance1=counter_balance1+1;
    else
      Location_of_guassians_MvsL(guassian_counter1,internal_counter1-counter_balance1)=internal_location_guassians_MvsL(internal_counter1)-1;
      MvsL_Gaussians.Unique_identifier(guassian_counter1,internal_counter1-counter_balance1)=Gaus_MvsL(internal_location_guassians_MvsL(internal_counter1),1);
      MvsL_Gaussians.Guassian_index_number(guassian_counter1,internal_counter1-counter_balance1)=str2num(Gaus_MvsL{internal_location_guassians_MvsL(internal_counter1),3});
      MvsL_Gaussians.Protein_number(guassian_counter1,internal_counter1-counter_balance1)=str2num(Gaus_MvsL{internal_location_guassians_MvsL(internal_counter1),4});
      MvsL_Gaussians.Replicate(guassian_counter1,internal_counter1-counter_balance1)=str2num(Gaus_MvsL{internal_location_guassians_MvsL(internal_counter1),5});
      MvsL_Gaussians.Protein{guassian_counter1,internal_counter1-counter_balance1}=Gaus_MvsL{internal_location_guassians_MvsL(internal_counter1),6};
      MvsL_Gaussians.Height(guassian_counter1,internal_counter1-counter_balance1)=str2num(Gaus_MvsL{internal_location_guassians_MvsL(internal_counter1),7});
      MvsL_Gaussians.Center(guassian_counter1,internal_counter1-counter_balance1)=str2num(Gaus_MvsL{internal_location_guassians_MvsL(internal_counter1),8});
      MvsL_Gaussians.Width(guassian_counter1,internal_counter1-counter_balance1)=str2num(Gaus_MvsL{internal_location_guassians_MvsL(internal_counter1),9});
      MvsL_Gaussians.SSE(guassian_counter1,internal_counter1-counter_balance1)=str2num(Gaus_MvsL{internal_location_guassians_MvsL(internal_counter1),10});
      MvsL_Gaussians.adjrsquare(guassian_counter1,internal_counter1-counter_balance1)=str2num(Gaus_MvsL{internal_location_guassians_MvsL(internal_counter1),11});
      MvsL_Gaussians.Complex_size(guassian_counter1,internal_counter1-counter_balance1)=str2num(Gaus_MvsL{internal_location_guassians_MvsL(internal_counter1),12});
      
      %Calculate Area of Gaussian
      guassian_area_fun =  @(x)( MvsL_Gaussians.Height(guassian_counter1,internal_counter1-counter_balance1)...
        *exp(-((x- MvsL_Gaussians.Center(guassian_counter1,internal_counter1-counter_balance1))...
        /MvsL_Gaussians.Width(guassian_counter1,internal_counter1-counter_balance1)).^2));
      guassian_area= integral(guassian_area_fun,1,55);
      
      %Save Gaussian area for MvsL
      MvsL_Gaussians.Gaussian_area(guassian_counter1,internal_counter1-counter_balance1)= guassian_area;
    end
  end
  
  
  %find the position in HvsL
  [internal_location_guassians_HvsL]=ind2sub(size(Gaus_HvsL), strmatch(Protein_name, Gaus_HvsL(:,2), 'exact'));
  number_counter2=length(internal_location_guassians_HvsL);
  counter_balance2=0;
  for internal_counter2 = 1:number_counter2
    if str2num(Gaus_HvsL{internal_location_guassians_HvsL(internal_counter2),8})>fraction_to_plot
      counter_balance2=counter_balance2+1;
    else
      Location_of_guassians_HvsL(guassian_counter1,internal_counter2-counter_balance2)=internal_location_guassians_HvsL(internal_counter2)-1;
      HvsL_Gaussians.Unique_identifier(guassian_counter1,internal_counter2-counter_balance2)=Gaus_HvsL(internal_location_guassians_HvsL(internal_counter2),1);
      HvsL_Gaussians.Guassian_index_number(guassian_counter1,internal_counter2-counter_balance2)=str2num(Gaus_HvsL{internal_location_guassians_HvsL(internal_counter2),3});
      HvsL_Gaussians.Protein_number(guassian_counter1,internal_counter2-counter_balance2)=str2num(Gaus_HvsL{internal_location_guassians_HvsL(internal_counter2),4});
      HvsL_Gaussians.Replicate(guassian_counter1,internal_counter2-counter_balance2)=str2num(Gaus_HvsL{internal_location_guassians_HvsL(internal_counter2),5});
      HvsL_Gaussians.Protein{guassian_counter1,internal_counter2-counter_balance2}=Gaus_HvsL{internal_location_guassians_HvsL(internal_counter2),6};
      HvsL_Gaussians.Height(guassian_counter1,internal_counter2-counter_balance2)=str2num(Gaus_HvsL{internal_location_guassians_HvsL(internal_counter2),7});
      HvsL_Gaussians.Center(guassian_counter1,internal_counter2-counter_balance2)=str2num(Gaus_HvsL{internal_location_guassians_HvsL(internal_counter2),8});
      HvsL_Gaussians.Width(guassian_counter1,internal_counter2-counter_balance2)=str2num(Gaus_HvsL{internal_location_guassians_HvsL(internal_counter2),9});
      HvsL_Gaussians.SSE(guassian_counter1,internal_counter2-counter_balance2)=str2num(Gaus_HvsL{internal_location_guassians_HvsL(internal_counter2),10});
      HvsL_Gaussians.adjrsquare(guassian_counter1,internal_counter2-counter_balance2)=str2num(Gaus_HvsL{internal_location_guassians_HvsL(internal_counter2),11});
      HvsL_Gaussians.Complex_size(guassian_counter1,internal_counter2-counter_balance2)=str2num(Gaus_HvsL{internal_location_guassians_HvsL(internal_counter2),12});
      
      
      %Calculate Area of Gaussian
      guassian_area_fun =  @(x)( HvsL_Gaussians.Height(guassian_counter1,internal_counter2-counter_balance2)...
        *exp(-((x- HvsL_Gaussians.Center(guassian_counter1,internal_counter2-counter_balance2))...
        /HvsL_Gaussians.Width(guassian_counter1,internal_counter2-counter_balance2)).^2));
      guassian_area= integral(guassian_area_fun,1,55);
      
      %Save Gaussian area for HvsL
      HvsL_Gaussians.Gaussian_area(guassian_counter1,internal_counter2-counter_balance2)= guassian_area;
    end
  end
  
  %Save name of Protein in each row
  if ~isempty(internal_location_guassians_MvsL)
    MvsL_Gaussians.Protein_name(guassian_counter1,1)=Gaus_MvsL(internal_location_guassians_MvsL(1),6);
    MvsL_Gaussians.Replicate_Protein_identifier(guassian_counter1,1)=Gaus_MvsL(internal_location_guassians_MvsL(1),2);
  else
    MvsL_Gaussians.Protein_name(guassian_counter1,1)=Gaus_HvsL(internal_location_guassians_HvsL(1),6);
    MvsL_Gaussians.Replicate_Protein_identifier(guassian_counter1,1)={'NaN'};
  end
  
  if ~isempty(internal_location_guassians_HvsL)
    HvsL_Gaussians.Protein_name(guassian_counter1,1)=Gaus_HvsL(internal_location_guassians_HvsL(1),6);
    HvsL_Gaussians.Replicate_Protein_identifier(guassian_counter1,1)=Gaus_HvsL(internal_location_guassians_HvsL(1),2);
  else
    HvsL_Gaussians.Protein_name(guassian_counter1,1)=Gaus_MvsL(internal_location_guassians_MvsL(1),6);
    HvsL_Gaussians.Replicate_Protein_identifier(guassian_counter1,1)={'NaN'};
  end
  
end


%Determine which row to remove based on if values absent from
rows_to_remove=zeros(Number_All_guassians_identified_nonredunent_list,1);
for guassian_counter1= 1:Number_All_guassians_identified_nonredunent_list
  if sum(nnz(HvsL_Gaussians.Center(guassian_counter1,:))+nnz(MvsL_Gaussians.Center(guassian_counter1,:)))==0
    rows_to_remove(guassian_counter1,1)=1;
  end
end

%remove MvsL
MvsL_Gaussians.Unique_identifier(rows_to_remove(:)==1,:)=[];
MvsL_Gaussians.Guassian_index_number(rows_to_remove==1,:)=[];
MvsL_Gaussians.Protein_number(rows_to_remove==1,:)=[];
MvsL_Gaussians.Replicate(rows_to_remove==1,:)=[];
MvsL_Gaussians.Protein(rows_to_remove==1,:)=[];
MvsL_Gaussians.Height(rows_to_remove==1,:)=[];
MvsL_Gaussians.Center(rows_to_remove==1,:)=[];
MvsL_Gaussians.Width(rows_to_remove==1,:)=[];
MvsL_Gaussians.SSE(rows_to_remove==1,:)=[];
MvsL_Gaussians.adjrsquare(rows_to_remove==1,:)=[];
MvsL_Gaussians.Complex_size(rows_to_remove==1,:)=[];
MvsL_Gaussians.Gaussian_area(rows_to_remove==1,:)=[];
MvsL_Gaussians.Protein_name(rows_to_remove==1,:)=[];
MvsL_Gaussians.Replicate_Protein_identifier(rows_to_remove==1,:)=[];

%remove HvsL
Location_of_guassians_HvsL(rows_to_remove==1,:)=[];
HvsL_Gaussians.Unique_identifier(rows_to_remove==1,:)=[];
HvsL_Gaussians.Guassian_index_number(rows_to_remove==1,:)=[];
HvsL_Gaussians.Protein_number(rows_to_remove==1,:)=[];
HvsL_Gaussians.Replicate(rows_to_remove==1,:)=[];
HvsL_Gaussians.Protein(rows_to_remove==1,:)=[];
HvsL_Gaussians.Height(rows_to_remove==1,:)=[];
HvsL_Gaussians.Center(rows_to_remove==1,:)=[];
HvsL_Gaussians.Width(rows_to_remove==1,:)=[];
HvsL_Gaussians.SSE(rows_to_remove==1,:)=[];
HvsL_Gaussians.adjrsquare(rows_to_remove==1,:)=[];
HvsL_Gaussians.Complex_size(rows_to_remove==1,:)=[];
HvsL_Gaussians.Gaussian_area(rows_to_remove==1,:)=[];
HvsL_Gaussians.Protein_name(rows_to_remove==1,:)=[];
HvsL_Gaussians.Replicate_Protein_identifier(rows_to_remove==1,:)=[];

%Define the size of the MvsL and HvsL arrays
Dimension_of_array=size(MvsL_Gaussians.Center);

tt = toc;
fprintf('  ...  %.2f seconds\n',tt)



%% 3. Determine which Gaussians are unique to each channel
fprintf('    3. Determine which Gaussians are unique to each channel (WITHIN REPLICATE)')


%Determine which gaussians are unique and which are shared between MvsL and
%HvsL of the same replicate
for guassian_counter2 = 1:Dimension_of_array(1)
  
  %determine Center for specific protein, note this is created so as the loop
  %cycle through each guassian the identified gaussians can be removed
  Replicate_MvsL_protein=MvsL_Gaussians.Replicate(guassian_counter2, :);
  Replicate_HvsL_protein=HvsL_Gaussians.Replicate(guassian_counter2, :);
  Centres_MvsL_protein=MvsL_Gaussians.Center(guassian_counter2, :);
  Centres_HvsL_protein=HvsL_Gaussians.Center(guassian_counter2, :);
  Height_MvsL_protein=MvsL_Gaussians.Height(guassian_counter2, :);
  Height_HvsL_protein=HvsL_Gaussians.Height(guassian_counter2, :);
  Width_MvsL_protein=MvsL_Gaussians.Width(guassian_counter2, :);
  Width_HvsL_protein=HvsL_Gaussians.Width(guassian_counter2, :);
  SSE_MvsL_protein=MvsL_Gaussians.SSE(guassian_counter2, :);
  SSE_HvsL_protein=HvsL_Gaussians.SSE(guassian_counter2, :);
  determined_AdjrSequare_of_guassians_MvsL=MvsL_Gaussians.adjrsquare(guassian_counter2,:);
  determined_AdjrSequare_of_guassians_HvsL=HvsL_Gaussians.adjrsquare(guassian_counter2,:);
  Complex_size_MvsL_protein=MvsL_Gaussians.Complex_size(guassian_counter2, :);
  Complex_size_HvsL_protein=HvsL_Gaussians.Complex_size(guassian_counter2, :);
  Protein_name_MvsL=MvsL_Gaussians.Protein_name(guassian_counter2, :);
  Protein_name_HvsL=HvsL_Gaussians.Protein_name(guassian_counter2, :);
  Unique_ID_MvsL=MvsL_Gaussians.Replicate_Protein_identifier(guassian_counter2, :);
  Unique_ID_HvsL=HvsL_Gaussians.Replicate_Protein_identifier(guassian_counter2, :);
  Area_MvsL=MvsL_Gaussians.Gaussian_area(guassian_counter2, :);
  Area_HvsL=HvsL_Gaussians.Gaussian_area(guassian_counter2, :);
  position_counter=1;
  
  %Determine which channel has the best gaussian fit
  AdjrSequare_of_guassians_channels=[determined_AdjrSequare_of_guassians_MvsL(1), determined_AdjrSequare_of_guassians_HvsL(1)];
  Max_AdjrSequare_of_guassians(guassian_counter2,1)=max(AdjrSequare_of_guassians_channels);
  position_of_best_AdjrSequare=find(ismember(AdjrSequare_of_guassians_channels,Max_AdjrSequare_of_guassians(guassian_counter2,1)));
  
  %Determine if the gaussians are the same until all guassians have been
  %compared
  
  ready = false;
  while ~ready
    %count the number of non zero elements
    Number_of_non_zero_elements_MvsL=nnz(Centres_MvsL_protein);
    Number_of_non_zero_elements_HvsL=nnz(Centres_HvsL_protein);
    %create the matrix to be filled
    Matrix_of_interactions=zeros(Number_of_non_zero_elements_MvsL,Number_of_non_zero_elements_HvsL);
    
    
    %create matrix to compare gaussians
    for internal_counter3= 1:Number_of_non_zero_elements_MvsL
      for internal_counter4= 1:Number_of_non_zero_elements_HvsL
        Centres_MvsL_for_matrix=Centres_MvsL_protein(internal_counter3);
        Centres_HvsL_for_matrix=Centres_HvsL_protein(internal_counter4);
        Matrix_of_interactions(internal_counter3,internal_counter4)= abs(Centres_MvsL_for_matrix-Centres_HvsL_for_matrix);
      end
    end
    
    %test matrix to determine which gaussains to keep and which ones to remove
    %test to makre sure the matrix has numbers
    if isempty(Matrix_of_interactions(:)) == 1
      row_to_copy_to_summed_array = nnz(Centres_MvsL_protein);
      col_to_copy_to_summed_array = nnz(Centres_HvsL_protein);
      if row_to_copy_to_summed_array>0
        for copy_counter1=1:row_to_copy_to_summed_array
          Summed_Replicate(guassian_counter2,position_counter)=Replicate_MvsL_protein(copy_counter1);
          Summed_Center_of_guassians(guassian_counter2,position_counter)= Centres_MvsL_protein(copy_counter1);
          Summed_Height_of_guassians(guassian_counter2,position_counter)=Height_MvsL_protein(copy_counter1);
          Summed_Width_of_guassians(guassian_counter2,position_counter)=Width_MvsL_protein(copy_counter1);
          Summed_SSE_of_guassians(guassian_counter2,position_counter)=SSE_MvsL_protein(copy_counter1);
          Summed_adjrsquare_of_guassians(guassian_counter2,position_counter)=determined_AdjrSequare_of_guassians_MvsL(copy_counter1);
          Summed_Complex_size(guassian_counter2,position_counter)=Complex_size_MvsL_protein(copy_counter1);
          Summed_protein_name(guassian_counter2,1)=Protein_name_MvsL(1);
          Summed_Unique_ID(guassian_counter2,1)=Unique_ID_MvsL(1);
          Combined_gaussian_areas(guassian_counter2,position_counter)=Area_MvsL(copy_counter1);
          Area_of_MvsL_gaussian(guassian_counter2,position_counter)=Area_MvsL(copy_counter1);
          Summed_Center_of_guassians_found_in_MvsL(guassian_counter2,position_counter)=1;
          position_counter=position_counter+1;
        end
        ind_to_remove = [1:row_to_copy_to_summed_array] ; % indices to be removed
        Replicate_MvsL_protein(ind_to_remove)=[];
        Centres_MvsL_protein(ind_to_remove)=[];
        Height_MvsL_protein(ind_to_remove)=[];
        Width_MvsL_protein(ind_to_remove)=[];
        SSE_MvsL_protein(ind_to_remove)=[];
        Complex_size_MvsL_protein(ind_to_remove)=[];
        Area_MvsL(ind_to_remove)=[];
        determined_AdjrSequare_of_guassians_MvsL(ind_to_remove)=[];
      end
      if col_to_copy_to_summed_array>0
        for copy_counter2=1:col_to_copy_to_summed_array
          Summed_Replicate(guassian_counter2,position_counter)=Replicate_HvsL_protein(copy_counter2);
          Summed_Center_of_guassians(guassian_counter2,position_counter)= Centres_HvsL_protein(copy_counter2);
          Summed_Height_of_guassians(guassian_counter2,position_counter)=Height_HvsL_protein(copy_counter2);
          Summed_Width_of_guassians(guassian_counter2,position_counter)=Width_HvsL_protein(copy_counter2);
          Summed_SSE_of_guassians(guassian_counter2,position_counter)=SSE_HvsL_protein(copy_counter2);
          Summed_adjrsquare_of_guassians(guassian_counter2,position_counter)=determined_AdjrSequare_of_guassians_HvsL(copy_counter2);
          Summed_Complex_size(guassian_counter2,position_counter)=Complex_size_HvsL_protein(copy_counter2);
          Summed_protein_name(guassian_counter2,1)=Protein_name_HvsL(1);
          Summed_Unique_ID(guassian_counter2,1)=Unique_ID_HvsL(1);
          Combined_gaussian_areas(guassian_counter2,position_counter)=Area_HvsL(copy_counter2);
          Area_of_HvsL_gaussian(guassian_counter2,position_counter)=Area_HvsL(copy_counter2);
          Summed_Center_of_guassians_found_in_HvsL(guassian_counter2,position_counter)=1;
          position_counter=position_counter+1;
        end
        ind_to_remove = [1:col_to_copy_to_summed_array] ; % indices to be removed
        Replicate_HvsL_protein(ind_to_remove)=[];
        Centres_HvsL_protein(ind_to_remove)=[];
        Height_HvsL_protein(ind_to_remove)=[];
        Width_HvsL_protein(ind_to_remove)=[];
        SSE_HvsL_protein(ind_to_remove)=[];
        Complex_size_HvsL_protein(ind_to_remove)=[];
        Area_HvsL(ind_to_remove)=[];
        determined_AdjrSequare_of_guassians_HvsL(ind_to_remove)=[];
      end
      %test to find the best gaussians within user defined number of fractions of each other
    elseif min(Matrix_of_interactions(:)) <=User_Window
      % position of the guassians that are the same between MvsL and HvsL channel
      [row,col] = ind2sub(size(Matrix_of_interactions),find(ismember((Matrix_of_interactions),min(Matrix_of_interactions(:)))));
      %use the Gaussian corresponding to the best AdjrSequare
      
      if position_of_best_AdjrSequare(1) == 1 % the MvsL contains the best adjrsequare
        Summed_Replicate(guassian_counter2,position_counter)=Replicate_MvsL_protein(row(1));
        Summed_Center_of_guassians(guassian_counter2,position_counter)= Centres_MvsL_protein(row(1));
        Summed_Height_of_guassians(guassian_counter2,position_counter)=Height_MvsL_protein(row(1));
        Summed_Width_of_guassians(guassian_counter2,position_counter)=Width_MvsL_protein(row(1));
        Summed_SSE_of_guassians(guassian_counter2,position_counter)=SSE_MvsL_protein(row(1));
        Summed_adjrsquare_of_guassians(guassian_counter2,position_counter)=determined_AdjrSequare_of_guassians_MvsL(row(1));
        Summed_Complex_size(guassian_counter2,position_counter)=Complex_size_MvsL_protein(row(1));
        Combined_gaussian_areas(guassian_counter2,position_counter)=Area_MvsL(row(1));
        Area_of_MvsL_gaussian(guassian_counter2,position_counter)=Area_MvsL(row(1));
        Area_of_HvsL_gaussian(guassian_counter2,position_counter)=Area_HvsL(col(1));
        Summed_protein_name(guassian_counter2,1)=Protein_name_MvsL(1);
        Summed_Unique_ID(guassian_counter2,1)=Unique_ID_MvsL(1);
        Summed_Center_of_guassians_found_in_MvsL(guassian_counter2,position_counter)=1;
        Summed_Center_of_guassians_found_in_HvsL(guassian_counter2,position_counter)=1;
        position_counter=position_counter+1; % add 1 to position
        
      elseif position_of_best_AdjrSequare(1) == 2 % the HvsL contains the best adjrsequare
        Summed_Replicate(guassian_counter2,position_counter)=Replicate_HvsL_protein(col(1));
        Summed_Center_of_guassians(guassian_counter2,position_counter)= Centres_HvsL_protein(col(1));
        Summed_Height_of_guassians(guassian_counter2,position_counter)=Height_HvsL_protein(col(1));
        Summed_Width_of_guassians(guassian_counter2,position_counter)=Width_HvsL_protein(col(1));
        Summed_SSE_of_guassians(guassian_counter2,position_counter)=SSE_HvsL_protein(col(1));
        Summed_adjrsquare_of_guassians(guassian_counter2,position_counter)=determined_AdjrSequare_of_guassians_HvsL(col(1));
        Summed_Complex_size(guassian_counter2,position_counter)=Complex_size_HvsL_protein(col(1));
        Combined_gaussian_areas(guassian_counter2,position_counter)=Area_HvsL(col(1));
        Area_of_MvsL_gaussian(guassian_counter2,position_counter)=Area_MvsL(row(1));
        Area_of_HvsL_gaussian(guassian_counter2,position_counter)=Area_HvsL(col(1));
        Summed_protein_name(guassian_counter2,1)=Protein_name_HvsL(1);
        Summed_Unique_ID(guassian_counter2,1)=Unique_ID_HvsL(1);
        Summed_Center_of_guassians_found_in_MvsL(guassian_counter2,position_counter)=1;
        Summed_Center_of_guassians_found_in_HvsL(guassian_counter2,position_counter)=1;
        position_counter=position_counter+1; % add 1 to position
      end
      
      %remove the most similar gaussains that are within 2 units
      %remove values in MvsL channels
      Replicate_MvsL_protein(:,[row(:)])=[];
      Centres_MvsL_protein(:,[row(:)])=[];
      Height_MvsL_protein(:,[row(:)])=[];
      Width_MvsL_protein(:,[row(:)])=[];
      SSE_MvsL_protein(:,[row(:)])=[];
      Complex_size_MvsL_protein(:,[row(:)])=[];
      Area_MvsL(:,[row(:)])=[];
      determined_AdjrSequare_of_guassians_MvsL(:,[row(:)])=[];
      
      %remove values in HvsL channels
      Replicate_HvsL_protein(:,[col(:)])=[];
      Centres_HvsL_protein(:,[col(:)])=[];
      Height_HvsL_protein(:,[col(:)])=[];
      Width_HvsL_protein(:,[col(:)])=[];
      SSE_HvsL_protein(:,[col(:)])=[];
      Complex_size_HvsL_protein(:,[col(:)])=[];
      Area_HvsL(:,[col(:)])=[];
      determined_AdjrSequare_of_guassians_HvsL(:,[col(:)])=[];
      
      %Consider gaussians not within user defined window of fractions of each other as unique
    elseif min(Matrix_of_interactions(:)) >User_Window
      % position of the guassians that are the same between MvsL and HvsL channel
      [row,col] = ind2sub(size(Matrix_of_interactions),find(ismember((Matrix_of_interactions),min(Matrix_of_interactions(:)))));
      %MvsL
      Summed_Replicate(guassian_counter2,position_counter)=Replicate_MvsL_protein(row(1));
      Summed_Center_of_guassians(guassian_counter2,position_counter)= Centres_MvsL_protein(row(1));
      Summed_Height_of_guassians(guassian_counter2,position_counter)=Height_MvsL_protein(row(1));
      Summed_Width_of_guassians(guassian_counter2,position_counter)=Width_MvsL_protein(row(1));
      Summed_SSE_of_guassians(guassian_counter2,position_counter)=SSE_MvsL_protein(row(1));
      Summed_adjrsquare_of_guassians(guassian_counter2,position_counter)=determined_AdjrSequare_of_guassians_MvsL(row(1));
      Summed_Complex_size(guassian_counter2,position_counter)=Complex_size_MvsL_protein(row(1));
      Area_of_MvsL_gaussian(guassian_counter2,position_counter)=Area_MvsL(row(1));
      Combined_gaussian_areas(guassian_counter2,position_counter)=Area_MvsL(row(1));
      Summed_protein_name(guassian_counter2,1)=Protein_name_MvsL(1);
      Summed_Unique_ID(guassian_counter2,1)=Unique_ID_MvsL(1);
      Summed_Center_of_guassians_found_in_MvsL(guassian_counter2,position_counter)=1;
      %HvsL
      Summed_Replicate(guassian_counter2,position_counter+1)=Replicate_HvsL_protein(col(1));
      Summed_Center_of_guassians(guassian_counter2,position_counter+1)= Centres_HvsL_protein(col(1));
      Summed_Height_of_guassians(guassian_counter2,position_counter+1)=Height_HvsL_protein(col(1));
      Summed_Width_of_guassians(guassian_counter2,position_counter+1)=Width_HvsL_protein(col(1));
      Summed_SSE_of_guassians(guassian_counter2,position_counter+1)=SSE_HvsL_protein(col(1));
      Summed_adjrsquare_of_guassians(guassian_counter2,position_counter+1)=determined_AdjrSequare_of_guassians_HvsL(col(1));
      Summed_Complex_size(guassian_counter2,position_counter+1)=Complex_size_HvsL_protein(col(1));
      Area_of_HvsL_gaussian(guassian_counter2,position_counter+1)=Area_HvsL(col(1));
      Combined_gaussian_areas(guassian_counter2,position_counter+1)=Area_HvsL(col(1));
      Summed_protein_name(guassian_counter2,1)=Protein_name_HvsL(1);
      Summed_Unique_ID(guassian_counter2,1)=Unique_ID_HvsL(1);
      Summed_Center_of_guassians_found_in_HvsL(guassian_counter2,position_counter+1)=1;
      position_counter=position_counter+2; % add 2 to position
      
      %remove the most similar gaussains that are within 2 units
      %remove values in MvsL channels
      Replicate_MvsL_protein(:,[row(:)])=[];
      Centres_MvsL_protein(:,[row(:)])=[];
      Height_MvsL_protein(:,[row(:)])=[];
      Width_MvsL_protein(:,[row(:)])=[];
      SSE_MvsL_protein(:,[row(:)])=[];
      Area_MvsL(:,[row(:)])=[];
      Complex_size_MvsL_protein(:,[row(:)])=[];
      determined_AdjrSequare_of_guassians_MvsL(:,[row(:)])=[];
      
      %remove values in HvsL channels
      Replicate_HvsL_protein(:,[col(:)])=[];
      Centres_HvsL_protein(:,[col(:)])=[];
      Height_HvsL_protein(:,[col(:)])=[];
      Width_HvsL_protein(:,[col(:)])=[];
      SSE_HvsL_protein(:,[col(:)])=[];
      Area_HvsL(:,[col(:)])=[];
      Complex_size_HvsL_protein(:,[col(:)])=[];
      determined_AdjrSequare_of_guassians_HvsL(:,[col(:)])=[];
      
    end
    
    %check the dimension of the matrix to all numbers have been checked
    row_check = nnz(Centres_MvsL_protein);
    col_check = nnz(Centres_HvsL_protein);
    ready = (row_check == 0 && col_check == 0);
  end
  
end

%Total number of Gaussians
Total_number_of_unique_gaussians=nnz(Summed_Center_of_guassians);

%Unique Gaussian in each replicate
Unique_gaussians_in_MvsL=nnz(Summed_Center_of_guassians_found_in_MvsL);
Unique_gaussians_in_HvsL=nnz(Summed_Center_of_guassians_found_in_HvsL);

%Determine th dimension of Summed_Center_of_guassians_found_in_HvsL
[number_of_proteins, number_of_column]=size(Summed_Center_of_guassians);

%determine if Summed_Center_of_guassians and MvsL or HvsL are equal
[number_of_proteins_M, number_of_column_M]=size(Summed_Center_of_guassians_found_in_MvsL);
[number_of_proteins_H, number_of_column_H]=size(Summed_Center_of_guassians_found_in_HvsL);

%test MvsL
if ~(number_of_column == number_of_column_M)
  column_to_add=number_of_column-number_of_column_M;
  require_columns=zeros(number_of_proteins_M,column_to_add);
  %Add values for center measurements
  Summed_Center_of_guassians_found_in_MvsL=[Summed_Center_of_guassians_found_in_MvsL, require_columns];
  %Add vlaues for Area_of_HvsL_gaussian measurements
  Area_of_MvsL_gaussian= [Area_of_MvsL_gaussian, require_columns];
end

%test HvsL
if ~(number_of_column == number_of_column_H)
  column_to_add=number_of_column-number_of_column_H;
  require_columns=zeros(number_of_proteins_H,column_to_add);
  %Add values for center measurements
  Summed_Center_of_guassians_found_in_HvsL=[Summed_Center_of_guassians_found_in_HvsL, require_columns];
  %Add values for Area_of_HvsL_gaussian measurements
  Area_of_HvsL_gaussian= [Area_of_HvsL_gaussian, require_columns];
end

%Shared Gaussians in each replicate
shared_gaussian_counter=0;
shared_gaussian_position=1;
for count_shared_guassians1=1:number_of_proteins
  for count_shared_guassians2=1:number_of_column
    if Summed_Center_of_guassians_found_in_HvsL(count_shared_guassians1,count_shared_guassians2) ==1 && Summed_Center_of_guassians_found_in_MvsL(count_shared_guassians1,count_shared_guassians2) ==1
      shared_gaussian_counter=1+shared_gaussian_counter;
      Comparing_gaussian_in_channels{shared_gaussian_position,1}='Both';
      Comparing_gaussian_in_channels_for_figure(shared_gaussian_position,1)=0.5;
      shared_gaussian_position=1+shared_gaussian_position;
    elseif Summed_Center_of_guassians_found_in_HvsL(count_shared_guassians1,count_shared_guassians2) ==1
      Comparing_gaussian_in_channels{shared_gaussian_position,1}='HvsL';
      Comparing_gaussian_in_channels_for_figure(shared_gaussian_position,1)=1;
      shared_gaussian_position=1+shared_gaussian_position;
    elseif Summed_Center_of_guassians_found_in_MvsL(count_shared_guassians1,count_shared_guassians2) ==1
      Comparing_gaussian_in_channels{shared_gaussian_position,1}='MvsL';
      Comparing_gaussian_in_channels_for_figure(shared_gaussian_position,1)=0.1;
      shared_gaussian_position=1+shared_gaussian_position;
    end
  end
end

%Write out a Combined_OutputGaus table with unique gaussians and if they
%were identified in both channels or only MvsL or HvsL
%Formatting channel information
Formatting_counter1=1;
for writeout_counter1= 1:number_of_proteins
  for writeout_counter2= 1:nnz(Summed_Height_of_guassians(writeout_counter1,:))
    formatted_Comparing_gaussian_in_channels(writeout_counter1,writeout_counter2)=Comparing_gaussian_in_channels(Formatting_counter1);
    Formatting_counter1=1+Formatting_counter1;
  end
end


%Write gaussian found between replicates to structure array
Combined_Gaussians.Unique_identifier=Summed_Unique_ID;
Combined_Gaussians.Protein_name=Summed_protein_name;
Combined_Gaussians.Replicate=Summed_Replicate(:,1);
Combined_Gaussians.Centrer=Summed_Center_of_guassians;
Combined_Gaussians.Height=Summed_Height_of_guassians;
Combined_Gaussians.Width= Summed_Width_of_guassians;
Combined_Gaussians.adjrsquare=Summed_adjrsquare_of_guassians;
Combined_Gaussians.Channels=formatted_Comparing_gaussian_in_channels;
Combined_Gaussians.SSE=Summed_SSE_of_guassians;
Combined_Gaussians.Complex_size=Summed_Complex_size;
Combined_Gaussians.Area_of_combined_gaussian=Combined_gaussian_areas;
Combined_Gaussians.Area_of_HvsL_gaussian=Area_of_HvsL_gaussian;
Combined_Gaussians.Area_of_MvsL_gaussian=Area_of_MvsL_gaussian;
Combined_Gaussians.Summed_combined_area=sum(Combined_gaussian_areas')';
Combined_Gaussians.Summed_HvsL_area=sum(Area_of_HvsL_gaussian')';
Combined_Gaussians.Summed_MvsL_area=sum(Area_of_MvsL_gaussian')';


%Determine the Coverage of complexes within sample based on area
Combined_Gaussians.Precentage_MvsL=Combined_Gaussians.Area_of_MvsL_gaussian/Standard_area;
Combined_Gaussians.Precentage_HvsL=Combined_Gaussians.Area_of_HvsL_gaussian/Standard_area;
Combined_Gaussians.Precentage_combined=Combined_Gaussians.Area_of_combined_gaussian/Standard_area;
Combined_Gaussians.Summed_Precentage_MvsL=Combined_Gaussians.Summed_MvsL_area/Standard_area;
Combined_Gaussians.Summed_Precentage_HvsL=Combined_Gaussians.Summed_HvsL_area/Standard_area;
Combined_Gaussians.Summed_Precentage_combined=Combined_Gaussians.Summed_combined_area/Standard_area;

%Calculate Coverage of the interactome
Coverage_of_complexes_MvsL=sum(mean(Combined_Gaussians.Precentage_MvsL));
Coverage_of_complexes_HvsL=sum(mean(Combined_Gaussians.Precentage_HvsL));

% #output
%Write table of is gaussian were fitted in the HvsL and/or MvsL
fid_combined_gaus_output= fopen([datadir1 'Unique_gaussians_observed_between_replicates.csv'],'wt'); % create the output file with the header infromation
fprintf (fid_combined_gaus_output,'%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n',... %header for OutputGaus output
  'Protein name','Complex Size', 'Height', 'Center', 'Width', 'SSE', 'adjrsquare','Observed in channel', 'Area observed for Unique Gaussians',...
  'Area observed for MvsL Gaussians','Area observed for HvsL Gaussians','Area coverage of unique Gaussians','Area coverage of MvsL Gaussians',...
  'Area coverage of HvsL Gaussians'); %Write Header
for writeout_counter1= 1:number_of_proteins
  for writeout_counter2= 1:nnz(Combined_Gaussians.Centrer(writeout_counter1,:))
    fprintf(fid_combined_gaus_output,'%s,%6.3f,%6.3f,%6.3f,%6.3f,%6.3f,%6.3f,%s,%6.3f,%6.3f,%6.3f,%6.3f,%6.3f,%6.3f,\n',...
      Combined_Gaussians.Unique_identifier{writeout_counter1},...
      Combined_Gaussians.Complex_size(writeout_counter1,writeout_counter2),...
      Combined_Gaussians.Height(writeout_counter1,writeout_counter2),...
      Combined_Gaussians.Centrer(writeout_counter1,writeout_counter2),...
      Combined_Gaussians.Width(writeout_counter1,writeout_counter2),...
      Combined_Gaussians.SSE(writeout_counter1,writeout_counter2),...
      Combined_Gaussians.adjrsquare(writeout_counter1,writeout_counter2),...
      Combined_Gaussians.Channels{writeout_counter1,writeout_counter2},...
      Combined_Gaussians.Area_of_combined_gaussian(writeout_counter1,writeout_counter2),...
      Combined_Gaussians.Area_of_MvsL_gaussian(writeout_counter1,writeout_counter2),...
      Combined_Gaussians.Area_of_HvsL_gaussian(writeout_counter1,writeout_counter2),...
      Combined_Gaussians.Precentage_combined(writeout_counter1,writeout_counter2),...
      Combined_Gaussians.Precentage_MvsL(writeout_counter1,writeout_counter2),...
      Combined_Gaussians.Precentage_HvsL(writeout_counter1,writeout_counter2));
  end
end
fclose(fid_combined_gaus_output);

%write out summary of guassian between replicates
fid_Summary_gaussian= fopen([datadir1 'Summary_gaussian_detected_between_replicates.csv'],'wt'); % create the summary file of the interaction output
fprintf (fid_Summary_gaussian,'%s,%s,%s,%s,%s,%s,\n',... %header for OutputGaus output
  'Total gaussians', 'Gaussians within MvsL channel', 'Gaussians within HvsL channel', 'Shared Gaussians', 'Coverage of complexes observed in MvsL',...
  'Coverage of complexes observed in HvsL'); %Write Header
fprintf (fid_Summary_gaussian,'%6.4f,%6.4f,%6.4f,%6.4f,%6.4f,%6.4f\n',...
  Total_number_of_unique_gaussians, Unique_gaussians_in_MvsL,...
  Unique_gaussians_in_HvsL, shared_gaussian_counter, Coverage_of_complexes_MvsL, Coverage_of_complexes_HvsL);
fprintf (fid_Summary_gaussian,',\n');
fclose(fid_Summary_gaussian);

tt = toc;
fprintf('  ...  %.2f seconds\n',tt)



%% 4. Detect 2-fold changes
% This is done with Gaussian Center
fprintf('    4. Detect 2-fold changes (WITHIN REPLICATE)')


%Unique gene name
All_protein_names=[MvsL_Gaussians.Protein_name; HvsL_Gaussians.Protein_name];
Unique_protein_names=unique(All_protein_names);
%remove NaN from list
test_nan=ismember(Unique_protein_names, {'NaN'});
Unique_protein_names(test_nan==1)=[];
%number of unique names in list
number_of_unique_protein_with_gaussians=length(Unique_protein_names);

%Create arrays to write apex values to
MvsL_Ratio_Apex=zeros(1,5);
HvsL_Ratio_Apex=zeros(1,5);
HvsM_Ratio_Apex=zeros(1,5);
Apex_values=zeros(Total_number_of_unique_gaussians,4);
Increase_at_apex_in_stimulation=zeros(number_of_proteins,1);
Decrease_at_apex_in_stimulation=zeros(number_of_proteins,1);
log_2_fold_comparsion_of_gaussian=zeros(number_of_proteins,1);
Do_gaussian_in_channels_change_for_figure=zeros(number_of_proteins,1);
Do_gaussian_in_channels_change=cell(number_of_proteins,1);
Do_gaussian_in_channels_position=1;
%Counters to detect changes
Increase_counter=0;
Decrease_counter=0;
No_change_counter=0;
Unquantifiable_counter=0;

test=0;
for Gaussian_counter1= 1:Dimension_of_array(1)
  for internal_counter5 = 1:nnz(Combined_Gaussians.Centrer(Gaussian_counter1,:))
    test=1+test;
  end
end

%Determine which Gaussians change by comparing fitted data with raw data
for Gaussian_counter1= 1:Dimension_of_array(1)
  for internal_counter5 = 1:nnz(Combined_Gaussians.Centrer(Gaussian_counter1,:))
    [location_Protein_in_raw_textdata]=ind2sub(size(txt_MvsL(:,1)), strmatch(Combined_Gaussians.Unique_identifier{Gaussian_counter1},txt_MvsL(:,1),'exact'));
    Center_to_test=Combined_Gaussians.Centrer(Gaussian_counter1,internal_counter5); % determine the center of the Gaussian
    rounded_center=round(Center_to_test); % round center of the Gaussian to nearest integer
    
    %Determine if Center is greater the fraction fraction_to_plot
    if rounded_center > fraction_to_plot
      
      %Mark as Unquantifiable
      dyanmic_title=strcat('Unquantifiable_Center_greater_than_',mat2str(fraction_to_plot));
      Do_gaussian_in_channels_change{Gaussian_counter1,internal_counter5}= dyanmic_title;
      Do_gaussian_in_channels_change_for_figure(Do_gaussian_in_channels_position)=0;
      log_2_fold_comparsion_of_gaussian(Do_gaussian_in_channels_position)= 0;
      Unquantifiable_counter=Unquantifiable_counter+1;
      
    elseif rounded_center <= fraction_to_plot
      %retrieve ratios of MvsL
      MvsL_Ratio_Apex(1,1) = num_val_MvsL(location_Protein_in_raw_textdata-1,rounded_center-2 +6 +1); %Add six as Realigned raw data ranges from -5 to 60
      MvsL_Ratio_Apex(1,2) = num_val_MvsL(location_Protein_in_raw_textdata-1,rounded_center-1 +6 +1); %Add six as Realigned raw data ranges from -5 to 60
      MvsL_Ratio_Apex(1,3) = num_val_MvsL(location_Protein_in_raw_textdata-1,rounded_center +6 +1); %Add six as Realigned raw data ranges from -5 to 60
      MvsL_Ratio_Apex(1,4) = num_val_MvsL(location_Protein_in_raw_textdata-1,rounded_center+1 +6 +1); %Add six as Realigned raw data ranges from -5 to 60
      MvsL_Ratio_Apex(1,5) = num_val_MvsL(location_Protein_in_raw_textdata-1,rounded_center+2 +6 +1); %Add six as Realigned raw data ranges from -5 to 60
      
      %retrieve ratios of HvsL
      HvsL_Ratio_Apex(1,1) = num_val_HvsL(location_Protein_in_raw_textdata-1,rounded_center-2+6 +1); %Add six as Realigned raw data ranges from -5 to 60
      HvsL_Ratio_Apex(1,2) = num_val_HvsL(location_Protein_in_raw_textdata-1,rounded_center-1 +6 +1); %Add six as Realigned raw data ranges from -5 to 60
      HvsL_Ratio_Apex(1,3) = num_val_HvsL(location_Protein_in_raw_textdata-1,rounded_center +6 +1); %Add six as Realigned raw data ranges from -5 to 60
      HvsL_Ratio_Apex(1,4) = num_val_HvsL(location_Protein_in_raw_textdata-1,rounded_center+1 +6 +1); %Add six as Realigned raw data ranges from -5 to 60
      HvsL_Ratio_Apex(1,5) = num_val_HvsL(location_Protein_in_raw_textdata-1,rounded_center+2 +6 +1); %Add six as Realigned raw data ranges from -5 to 60
      
      %retrieve ratios of HvsM
      HvsM_Ratio_Apex(1,1) = num_val_HvsM(location_Protein_in_raw_textdata-1,rounded_center-2+6 +1); %Add six as Realigned raw data ranges from -5 to 60
      HvsM_Ratio_Apex(1,2) = num_val_HvsM(location_Protein_in_raw_textdata-1,rounded_center-1+6 +1); %Add six as Realigned raw data ranges from -5 to 60
      HvsM_Ratio_Apex(1,3) = num_val_HvsM(location_Protein_in_raw_textdata-1,rounded_center +6 +1); %Add six as Realigned raw data ranges from -5 to 60
      HvsM_Ratio_Apex(1,4) = num_val_HvsM(location_Protein_in_raw_textdata-1,rounded_center+1 +6 +1); %Add six as Realigned raw data ranges from -5 to 60
      HvsM_Ratio_Apex(1,5) = num_val_HvsM(location_Protein_in_raw_textdata-1,rounded_center+2 +6 +1); %Add six as Realigned raw data ranges from -5 to 60
      
      if (nnz(HvsL_Ratio_Apex) == 0 | nnz(MvsL_Ratio_Apex) == 0 | nnz(HvsM_Ratio_Apex) == 0 |  sum(isnan(HvsL_Ratio_Apex(1,2:4)) >=2 ) | sum(isnan(MvsL_Ratio_Apex(1,2:4)) >= 2))
        %If value for the apex is missing or not values around the apex can be mark as Unquantifiable
        Do_gaussian_in_channels_change{Gaussian_counter1,internal_counter5}= 'Unquantifiable';
        Do_gaussian_in_channels_change_for_figure(Do_gaussian_in_channels_position)=0;
        log_2_fold_comparsion_of_gaussian(Do_gaussian_in_channels_position)= 0;
        Unquantifiable_counter=Unquantifiable_counter+1;
        
      elseif ~(nnz(HvsL_Ratio_Apex) == 0 | nnz(MvsL_Ratio_Apex) == 0 | nnz(HvsM_Ratio_Apex) == 0 |  sum(isnan(HvsL_Ratio_Apex(1,2:4)) >= 2) | sum(isnan(MvsL_Ratio_Apex(1,2:4)) >= 2))
        
        %save values
        MvsL_Ratio_Apex_prefiltered=MvsL_Ratio_Apex;
        HvsL_Ratio_Apex_prefiltered=HvsL_Ratio_Apex;
        
        %Remove values 0.5
        MvsL_Ratio_Apex(MvsL_Ratio_Apex(1,:)<0.5)=[];
        HvsL_Ratio_Apex(HvsL_Ratio_Apex(1,:)<0.5)=[];
        
        %Remove nan values
        HvsM_Ratio_Apex(isnan(HvsM_Ratio_Apex))=[];
        HvsL_Ratio_Apex(isnan(HvsL_Ratio_Apex))=[];
        MvsL_Ratio_Apex(isnan(MvsL_Ratio_Apex))=[];
        
        %Write values to table to output
        Apex_values(Do_gaussian_in_channels_position,1)=mean(MvsL_Ratio_Apex_prefiltered./HvsL_Ratio_Apex_prefiltered);
        Apex_values(Do_gaussian_in_channels_position,2)=mean(MvsL_Ratio_Apex); %save the mean of MvsL ratios
        Apex_values(Do_gaussian_in_channels_position,3)=mean(HvsL_Ratio_Apex); %save the mean of HvsL ratios
        Apex_values(Do_gaussian_in_channels_position,4)=mean(HvsM_Ratio_Apex); %save the mean of HvsM ratios
        
        %Test if apex changes
        
        if (Apex_values(Do_gaussian_in_channels_position,4) <= 2 && Apex_values(Do_gaussian_in_channels_position,4) >= 0.5)
          Increase_at_apex_in_stimulation(Gaussian_counter1,internal_counter5)= 0;
          Decrease_at_apex_in_stimulation(Gaussian_counter1,internal_counter5)= 0;
          Do_gaussian_in_channels_change{Gaussian_counter1,internal_counter5}= 'No change';
          Do_gaussian_in_channels_change_for_figure(Do_gaussian_in_channels_position)=5;
          log_2_fold_comparsion_of_gaussian(Do_gaussian_in_channels_position)= log2(Apex_values(Do_gaussian_in_channels_position,4));
          No_change_counter=No_change_counter+1;
          
        elseif (Apex_values(Do_gaussian_in_channels_position,4) > 2)
          Increase_at_apex_in_stimulation(Gaussian_counter1,internal_counter5)= 1;
          Decrease_at_apex_in_stimulation(Gaussian_counter1,internal_counter5)= 0;
          Do_gaussian_in_channels_change{Gaussian_counter1,internal_counter5}='Increase';
          Do_gaussian_in_channels_change_for_figure(Do_gaussian_in_channels_position)=10;
          log_2_fold_comparsion_of_gaussian(Do_gaussian_in_channels_position)= log2(Apex_values(Do_gaussian_in_channels_position,4));
          Increase_counter=Increase_counter+1;
          
        elseif (Apex_values(Do_gaussian_in_channels_position,4) < 0.5)
          Decrease_at_apex_in_stimulation(Gaussian_counter1,internal_counter5)= 1;
          Increase_at_apex_in_stimulation(Gaussian_counter1,internal_counter5)= 0;
          Do_gaussian_in_channels_change{Gaussian_counter1,internal_counter5}='Decrease';
          Do_gaussian_in_channels_change_for_figure(Do_gaussian_in_channels_position)=1;
          log_2_fold_comparsion_of_gaussian(Do_gaussian_in_channels_position)= log2(Apex_values(Do_gaussian_in_channels_position,4));
          Decrease_counter=Decrease_counter+1;
          
        end
      end
    end
    % add one to position and move to next center to check
    Do_gaussian_in_channels_position=1+Do_gaussian_in_channels_position;
  end
end

%Normalise log2 comparsion
mean_log_2_fold_comparsion_of_gaussian=mean(log_2_fold_comparsion_of_gaussian);
normaliased_log_2_fold_comparsion_of_gaussian=log_2_fold_comparsion_of_gaussian-mean_log_2_fold_comparsion_of_gaussian;

%Copy Log2 values into same format as Height, Center and width
Formatting_counter=1;
for writeout_counter1= 1:Dimension_of_array(1)
  for writeout_counter2= 1:nnz(Combined_Gaussians.Centrer(writeout_counter1,:))
    formated_normaliased_log_2_fold_comparsion_of_gaussian(writeout_counter1,writeout_counter2)=normaliased_log_2_fold_comparsion_of_gaussian(Formatting_counter);
    formated_log_2_fold_comparsion_of_gaussian(writeout_counter1,writeout_counter2)=log_2_fold_comparsion_of_gaussian(Formatting_counter);
    Formatting_counter=1+Formatting_counter;
  end
end

%Write comparsion output to combined structure array
Combined_Gaussians.log2_of_gaussians=formated_log_2_fold_comparsion_of_gaussian;
Combined_Gaussians.log2_normalised_gaussians=formated_normaliased_log_2_fold_comparsion_of_gaussian;
Combined_Gaussians.Increase_at_apex=Increase_at_apex_in_stimulation;
Combined_Gaussians.Decrease_at_apex=Decrease_at_apex_in_stimulation;
Combined_Gaussians.Observed_change=Do_gaussian_in_channels_change;


%%Count how many guassian change in each fraction
Hist_array=zeros(fraction_number(2),3);

for fraction_counter1= 1:number_of_proteins
  for fraction_counter2= 1:nnz(Combined_Gaussians.Centrer(fraction_counter1,:))
    %Define Center of Gaussian
    Center_of_gaussain= floor(Combined_Gaussians.Centrer(fraction_counter1,fraction_counter2));
    %Ensure center is greater then zero
    if  Center_of_gaussain>0
      
      if Combined_Gaussians.log2_normalised_gaussians(fraction_counter1,fraction_counter2) >=1 && Center_of_gaussain <=fraction_to_plot
        
        Hist_array(Center_of_gaussain,1)=Hist_array(Center_of_gaussain,1)+1;
        Hist_array(Center_of_gaussain,2)=Hist_array(Center_of_gaussain,2)+1;
      elseif Combined_Gaussians.log2_normalised_gaussians(fraction_counter1,fraction_counter2) <=-1 && Center_of_gaussain <=fraction_to_plot
        
        Hist_array(Center_of_gaussain,1)=Hist_array(Center_of_gaussain,1)+1;
        Hist_array(Center_of_gaussain,3)=Hist_array(Center_of_gaussain,3)+1;
      elseif Combined_Gaussians.log2_normalised_gaussians(fraction_counter1,fraction_counter2) >=-1....
          && Combined_Gaussians.log2_normalised_gaussians(fraction_counter1,fraction_counter2) <=1....
          && Center_of_gaussain <=fraction_to_plot
        
        Hist_array(Center_of_gaussain,1)=Hist_array(Center_of_gaussain,1)+1;
      end
    end
  end
end

% #output
%Write out comparsion of gaussian in biological replicate
fid_combined_gaus_with_changes_output= fopen([datadir1 'Unique_gaussians_with_changes.csv'],'wt'); % create the output file with the header infromation
fprintf (fid_combined_gaus_with_changes_output,'%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n',... %header for OutputGaus output
  'Protein name', 'Height', 'Center', 'Width', 'SSE', 'adjrsquare','Observed in channel','Changes Observed',...
  'Area observed for Unique Gaussians','Area coverage of unique Gaussians', 'Fold Change (based on raw data)','Normalise Fold Change (based on raw data)'); %Write Header
for writeout_counter1= 1:number_of_proteins
  for writeout_counter2= 1:nnz(Combined_Gaussians.Centrer(writeout_counter1,:))
    fprintf(fid_combined_gaus_with_changes_output,'%s,%6.3f,%6.3f,%6.3f,%6.3f,%6.3f,%s,%s,%6.4f,%6.4f,%6.4f,%6.4f,\n',...
      Combined_Gaussians.Unique_identifier{writeout_counter1},...
      Combined_Gaussians.Height(writeout_counter1,writeout_counter2),...
      Combined_Gaussians.Centrer(writeout_counter1,writeout_counter2),...
      Combined_Gaussians.Width(writeout_counter1,writeout_counter2),...
      Combined_Gaussians.SSE(writeout_counter1,writeout_counter2),...
      Combined_Gaussians.adjrsquare(writeout_counter1,writeout_counter2),...
      Combined_Gaussians.Channels{writeout_counter1,writeout_counter2},...
      Combined_Gaussians.Observed_change{writeout_counter1,writeout_counter2},...
      Combined_Gaussians.Area_of_combined_gaussian(writeout_counter1, writeout_counter2),...
      Combined_Gaussians.Precentage_combined(writeout_counter1,writeout_counter2),...
      Combined_Gaussians.log2_of_gaussians(writeout_counter1,writeout_counter2),...
      Combined_Gaussians.log2_normalised_gaussians(writeout_counter1,writeout_counter2));
  end
end
fclose(fid_combined_gaus_with_changes_output);

%write out summary
fid_Summary_gaussian= fopen([datadir1 'Summary_gaussian_detected.csv'],'wt'); % create the summary file of the interaction output
fprintf (fid_Summary_gaussian,'%s,%s,%s,%s,%s,%s,%s,%s,\n',... %header for OutputGaus output
  'Total gaussians', 'Gaussians within MvsL channel', 'Gaussians within HvsL channel', 'Shared Gaussians',...
  'Number of Gaussians that increase', 'Number of Gaussians that decrease',....
  'Number of Gaussians that do not change','Number of Unquantified' ); %Write Header
fprintf (fid_Summary_gaussian,'%6.4f,%6.4f,%6.4f,%6.4f,%6.4f,%6.4f,%6.4f,%6.4f\n',...
  Total_number_of_unique_gaussians, Unique_gaussians_in_MvsL,...
  Unique_gaussians_in_HvsL, shared_gaussian_counter, Increase_counter, Decrease_counter, No_change_counter, Unquantifiable_counter);
fprintf (fid_Summary_gaussian,',\n');
fclose(fid_Summary_gaussian);

tt = toc;
fprintf('  ...  %.2f seconds\n',tt)



%% 5. Determine the trends of guassian curves
% Aim: look at the determined gaussians within single replicate and compare
% if the all the gaussians change, some change or none of them change
fprintf('    5. Determine the trends of guassian curves (WITHIN REPLICATE)')


%Define Variables
Combined_Gaussians.Trend=cell(Dimension_of_array(1),2);
dont_change=zeros(Dimension_of_array(1),1);
shows_increase=zeros(Dimension_of_array(1),1);
shows_decrease=zeros(Dimension_of_array(1),1);
Combined_Gaussians.Total_number_of_guassian_for_protein=zeros(Dimension_of_array(1),1);
No_change_consistent_across_gaus=0;
Increase_consistent_across_gaus=0;
Decrease_consistent_across_gaus=0;
Increase_inconsistent_across_gaus=0;
Decrease_inconsistent_across_gaus=0;
Increase_and_decrease_across_gaus=0;

%Loop over proteins to assess guassian information
for writeout_counter1= 1:Dimension_of_array(1)
  Protein_to_test_trend= Do_gaussian_in_channels_change(writeout_counter1,:);
  %Count number of proteins
  Combined_Gaussians.Total_number_of_guassian_for_protein(writeout_counter1,1)=sum(~cellfun('isempty',Protein_to_test_trend));
  
  %test if change is uniform 'No change', count
  no_change={'No change'};
  dont_change(writeout_counter1,1)=sum(strcmp(no_change,Protein_to_test_trend));
  %test if change is uniform 'Increase', count
  Increase_change={'Increase'};
  shows_increase(writeout_counter1,1)=sum(strcmp(Increase_change,Protein_to_test_trend));
  %test if change is uniform 'Decrease', count
  Decrease_change={'Decrease'};
  shows_decrease(writeout_counter1,1)=sum(strcmp(Decrease_change,Protein_to_test_trend));
  
  if (dont_change(writeout_counter1,1) >0 && shows_increase(writeout_counter1,1) == 0 && shows_decrease(writeout_counter1,1) == 0)
    Combined_Gaussians.Trend(writeout_counter1,1)={'No change'};
    Combined_Gaussians.Trend(writeout_counter1,2)={'Consitent across gaussians'};
    No_change_consistent_across_gaus=1+ No_change_consistent_across_gaus;
  elseif (dont_change(writeout_counter1,1) ==0 && shows_increase(writeout_counter1,1) > 0 && shows_decrease(writeout_counter1,1) == 0)
    Combined_Gaussians.Trend(writeout_counter1,1)={'Increase'};
    Combined_Gaussians.Trend(writeout_counter1,2)={'Consitent across gaussians'};
    Increase_consistent_across_gaus=1+ Increase_consistent_across_gaus;
  elseif (dont_change(writeout_counter1,1) ==0 && shows_increase(writeout_counter1,1) == 0 && shows_decrease(writeout_counter1,1) > 0)
    Combined_Gaussians.Trend(writeout_counter1,1)={'Decrease'};
    Combined_Gaussians.Trend(writeout_counter1,2)={'Consitent across gaussians'};
    Decrease_consistent_across_gaus=1+ Decrease_consistent_across_gaus;
  elseif (shows_increase(writeout_counter1,1) >0 && dont_change(writeout_counter1,1) > 0)
    Combined_Gaussians.Trend(writeout_counter1,1)={'Increase'};
    Combined_Gaussians.Trend(writeout_counter1,2)={'Inconsitent across gaussians'};
    Increase_inconsistent_across_gaus=1+ Increase_inconsistent_across_gaus;
  elseif (dont_change(writeout_counter1,1) >0 && shows_decrease(writeout_counter1,1) > 0)
    Combined_Gaussians.Trend(writeout_counter1,1)={'Decrease'};
    Combined_Gaussians.Trend(writeout_counter1,2)={'Inconsitent across gaussians'};
    Decrease_inconsistent_across_gaus=1+ Decrease_inconsistent_across_gaus;
  elseif (shows_increase(writeout_counter1,1) >0 && shows_decrease(writeout_counter1,1) > 0)
    Combined_Gaussians.Trend(writeout_counter1,1)={'+/-'};
    Combined_Gaussians.Trend(writeout_counter1,2)={'Inconsitent across gaussians'};
    Increase_and_decrease_across_gaus=1+Increase_and_decrease_across_gaus;
  end
  
end

%Arrange trend data to group by proteins across replicates
Protein_information.Trend_compare_between_replicates=cell(number_of_unique_protein_with_gaussians,6);

for count_shared_guassians4= 1:number_of_unique_protein_with_gaussians
  %find location of proteins within master protein table
  [internal_location_of_protein_of_interest]=ind2sub(Dimension_of_array(1), strmatch(Unique_protein_names(count_shared_guassians4), Combined_Gaussians.Protein_name, 'exact'));
  
  Number_of_times_seen=length(internal_location_of_protein_of_interest);
  for Number_of_times_seen_counter1= 1:Number_of_times_seen
    if Combined_Gaussians.Replicate(internal_location_of_protein_of_interest(Number_of_times_seen_counter1)) == 1
      Protein_information.Trend_compare_between_replicates{count_shared_guassians4,1}= Combined_Gaussians.Trend(internal_location_of_protein_of_interest(Number_of_times_seen_counter1),1);
      Protein_information.Trend_compare_between_replicates{count_shared_guassians4,4}= Combined_Gaussians.Trend(internal_location_of_protein_of_interest(Number_of_times_seen_counter1),2);
    elseif Combined_Gaussians.Replicate(internal_location_of_protein_of_interest(Number_of_times_seen_counter1)) == 2
      Protein_information.Trend_compare_between_replicates{count_shared_guassians4,2}= Combined_Gaussians.Trend(internal_location_of_protein_of_interest(Number_of_times_seen_counter1),1);
      Protein_information.Trend_compare_between_replicates{count_shared_guassians4,5}= Combined_Gaussians.Trend(internal_location_of_protein_of_interest(Number_of_times_seen_counter1),2);
    elseif Combined_Gaussians.Replicate(internal_location_of_protein_of_interest(Number_of_times_seen_counter1)) == 3
      Protein_information.Trend_compare_between_replicates{count_shared_guassians4,3}= Combined_Gaussians.Trend(internal_location_of_protein_of_interest(Number_of_times_seen_counter1),1);
      Protein_information.Trend_compare_between_replicates{count_shared_guassians4,6}= Combined_Gaussians.Trend(internal_location_of_protein_of_interest(Number_of_times_seen_counter1),2);
    end
  end
end


% #output
%Write out trend observed in gaussian grouped by proteins and replicate
fid_summed_area= fopen([datadir1 'Gaussian_trend_analysis_protein_replicate.csv'],'wt'); % create the output file with the header infromation
fprintf (fid_summed_area,'%s,%s,%s,%s,%s,%s,%s,%s\n',... %header for OutputGaus output
  'Protein name', 'Summed Area MvsL', 'Summed Area HvsL', '% of expect area MvsL', '% of expect area HvsL', 'Number of gaussians observed','Observed change',' Where changes consistent across all guassians'); %Write Header
for writeout_counter1= 1:Dimension_of_array(1)
  fprintf(fid_summed_area,'%s,%6.3f,%6.3f,%6.3f,%6.3f,%6.3f,%s,%s,\n',...
    Combined_Gaussians.Unique_identifier{writeout_counter1},...
    Combined_Gaussians.Summed_MvsL_area(writeout_counter1,1),...
    Combined_Gaussians.Summed_HvsL_area(writeout_counter1,1),...
    Combined_Gaussians.Summed_Precentage_MvsL(writeout_counter1,1),...
    Combined_Gaussians.Summed_Precentage_HvsL(writeout_counter1,1),...
    nnz(Combined_Gaussians.Centrer(writeout_counter1,:)),...
    Combined_Gaussians.Trend{writeout_counter1,1},...
    Combined_Gaussians.Trend{writeout_counter1,2});
end
fclose(fid_summed_area);

fid_summed_area1= fopen([datadir1 'Summary_gaussian_trend_analysis_protein_replicate.csv'],'wt'); % create the output file with the header infromation
fprintf (fid_summed_area1,'%s,%s,%s,%s,%s,%s,%s\n',... %header for OutputGaus output
  'Number of Protein observed across replicates', 'No change- consistent across gaussians', 'Increase- consistent across gaussians',...
  'Decrease- consistent across gaussians','Increase- inconsistent across gaussians', 'Decrease- inconsistent across gaussians',...
  'Increase and Decrease observed in Gaussians'); %Write Header
fprintf (fid_summed_area1,'%6.4f,%6.4f,%6.4f,%6.4f,%6.4f,%6.4f,%6.4f\n',...
  Dimension_of_array(1),...
  No_change_consistent_across_gaus, Increase_consistent_across_gaus,...
  Decrease_consistent_across_gaus, Increase_inconsistent_across_gaus,...
  Decrease_inconsistent_across_gaus, Increase_and_decrease_across_gaus);
fclose(fid_summed_area1);

tt = toc;
fprintf('  ...  %.2f seconds\n',tt)



%% 6. Shared gaussian across replicate.
% Comparsion 3 Part A:  Determine which Gaussians are shared across replicates
% Aim: determined gaussians which are shared
fprintf('    6. Shared gaussian across replicate (BETWEEN REPLICATES)')


%Create array of master gaussian list
Master_Gaussian_list.Protein_name=Unique_protein_names;
Master_Gaussian_list.Unique_identifier=cell(number_of_unique_protein_with_gaussians,1);
Master_Gaussian_list.Protein_number=zeros(number_of_unique_protein_with_gaussians,1);
Master_Gaussian_list.Replicate=zeros(number_of_unique_protein_with_gaussians,1);
Master_Gaussian_list.Channel=cell(number_of_unique_protein_with_gaussians,1);
Master_Gaussian_list.Guassian_index_number=zeros(number_of_unique_protein_with_gaussians,1);
Master_Gaussian_list.Center=zeros(number_of_unique_protein_with_gaussians,1);
Master_Gaussian_list.Height=zeros(number_of_unique_protein_with_gaussians,1);
Master_Gaussian_list.Width=zeros(number_of_unique_protein_with_gaussians,1);
Master_Gaussian_list.SSE=zeros(number_of_unique_protein_with_gaussians,1);
Master_Gaussian_list.adjrsquare=zeros(number_of_unique_protein_with_gaussians,1);
Master_Gaussian_list.Complex_size=zeros(number_of_unique_protein_with_gaussians,1);
Master_Gaussian_list.Replicate_Protein_identifier=cell(number_of_unique_protein_with_gaussians,1);
Master_Gaussian_list.Gaussian_area=zeros(number_of_unique_protein_with_gaussians,1);

%Compare MvsL first then HvsL
for Shared_accross_replicate_counter = 1:number_of_unique_protein_with_gaussians
  
  %find location of proteins in MvsL and HvsL arrays
  [internal_location_of_protein_of_interest_MvsL]=ind2sub(Dimension_of_array(1), strmatch(Unique_protein_names(Shared_accross_replicate_counter), MvsL_Gaussians.Protein_name, 'exact'));
  [internal_location_of_protein_of_interest_Master]=ind2sub(Dimension_of_array(1), strmatch(Unique_protein_names(Shared_accross_replicate_counter), Master_Gaussian_list.Protein_name, 'exact'));
  
  %Compare gaussians from each replicate channel
  Number_of_times_seen_MvsL=length(internal_location_of_protein_of_interest_MvsL);
  
  %Compare MvsL Gaussians
  for MvsL_counter=1:Number_of_times_seen_MvsL
    
    %determine Center for specific protein, note this is created so as the loop
    %cycle through each guassian the identified gaussians can be removed
    %Copy replicate data
    Replicate_MvsL_protein=MvsL_Gaussians.Replicate(internal_location_of_protein_of_interest_MvsL(MvsL_counter), :);
    Replicate_Master_list_protein=Master_Gaussian_list.Replicate(internal_location_of_protein_of_interest_Master, :);
    %Copy Unique_identifier
    Unique_identified_MvsL=MvsL_Gaussians.Unique_identifier(internal_location_of_protein_of_interest_MvsL(MvsL_counter), :);
    Unique_identified_Master_list_protein=Master_Gaussian_list.Unique_identifier(internal_location_of_protein_of_interest_Master, :);
    %Copy Protein number
    Protein_number_MvsL_protein=MvsL_Gaussians.Protein_number(internal_location_of_protein_of_interest_MvsL(MvsL_counter), :);
    Protein_number_Master_list_protein=Master_Gaussian_list.Protein_number(internal_location_of_protein_of_interest_Master, :);
    %Copy Guassian_index_number
    Guassian_index_number_MvsL_protein=MvsL_Gaussians.Guassian_index_number(internal_location_of_protein_of_interest_MvsL(MvsL_counter), :);
    Guassian_index_number_Master_list_protein=Master_Gaussian_list.Guassian_index_number(internal_location_of_protein_of_interest_Master, :);
    %Copy Replicate_Protein_identifier
    Replicate_Protein_identifier_MvsL_protein=MvsL_Gaussians.Replicate_Protein_identifier(internal_location_of_protein_of_interest_MvsL(MvsL_counter), :);
    Replicate_Protein_identifier_Master_list_protein=Master_Gaussian_list.Replicate_Protein_identifier(internal_location_of_protein_of_interest_Master, :);
    %Copy Center data
    Centres_MvsL_protein=MvsL_Gaussians.Center(internal_location_of_protein_of_interest_MvsL(MvsL_counter), :);
    Centres_Master_list_protein=Master_Gaussian_list.Center(internal_location_of_protein_of_interest_Master, :);
    %Copy Height data
    Height_MvsL_protein=MvsL_Gaussians.Height(internal_location_of_protein_of_interest_MvsL(MvsL_counter), :);
    Height_Master_list_protein=Master_Gaussian_list.Height(internal_location_of_protein_of_interest_Master, :);
    %Copy Width data
    Width_MvsL_protein=MvsL_Gaussians.Width(internal_location_of_protein_of_interest_MvsL(MvsL_counter), :);
    Width_Master_list_protein=Master_Gaussian_list.Width(internal_location_of_protein_of_interest_Master, :);
    %Copy SSE data
    SSE_MvsL_protein=MvsL_Gaussians.SSE(internal_location_of_protein_of_interest_MvsL(MvsL_counter), :);
    SSE_Master_list_protein=Master_Gaussian_list.SSE(internal_location_of_protein_of_interest_Master, :);
    %Copy adjrsquare data
    determined_AdjrSequare_of_guassians_MvsL=MvsL_Gaussians.adjrsquare(internal_location_of_protein_of_interest_MvsL(MvsL_counter), :);
    determined_AdjrSequare_of_guassians_Master_list=Master_Gaussian_list.adjrsquare(internal_location_of_protein_of_interest_Master, :);
    %Copy Complex size
    Complex_size_MvsL=MvsL_Gaussians.Complex_size(internal_location_of_protein_of_interest_MvsL(MvsL_counter), :);
    Complex_size_Master_list=Master_Gaussian_list.Complex_size(internal_location_of_protein_of_interest_Master, :);
    %Copy Gaussian area
    Gaussian_area_MvsL=MvsL_Gaussians.Gaussian_area(internal_location_of_protein_of_interest_MvsL(MvsL_counter), :);
    Gaussian_area_Master_list=Master_Gaussian_list.Gaussian_area(internal_location_of_protein_of_interest_Master, :);
    
    %Determine what position in Master_Gaussian_list.Center to write to
    if Master_Gaussian_list.Center(internal_location_of_protein_of_interest_Master,1)== 0
      position_counter=0;
    elseif ~Master_Gaussian_list.Center(internal_location_of_protein_of_interest_Master,1) == 0
      position_counter=(nnz(Master_Gaussian_list.Center(internal_location_of_protein_of_interest_Master, :)));
    end
    
    ready = false;
    while ~ready
      %count the number of non zero elements
      Number_of_non_zero_elements_MvsL=nnz(Centres_MvsL_protein);
      Number_of_non_zero_elements_Master_list=nnz(Centres_Master_list_protein);
      %create the matrix to be filled
      Matrix_of_interactions=zeros(Number_of_non_zero_elements_MvsL,Number_of_non_zero_elements_Master_list);
      
      
      %create matrix to compare gaussians
      for internal_counter6= 1:Number_of_non_zero_elements_MvsL
        for internal_counter7= 1:Number_of_non_zero_elements_Master_list
          Centres_MvsL_for_matrix=Centres_MvsL_protein(internal_counter6);
          Centres_Master_list_protein_for_matrix=Centres_Master_list_protein(internal_counter7);
          Matrix_of_interactions(internal_counter6,internal_counter7)= abs(Centres_MvsL_for_matrix-Centres_Master_list_protein_for_matrix);
        end
      end
      
      
      %test matrix to determine which gaussains to keep and which ones to remove
      %test to makre sure the matrix has numbers
      if isempty(Matrix_of_interactions(:)) == 1
        row_to_copy_to_summed_array = nnz(Centres_MvsL_protein);
        col_to_copy_to_summed_array = nnz(Centres_Master_list_protein);
        if row_to_copy_to_summed_array>0
          
          %add one to position counter
          position_counter=position_counter+1;
          
          %populate Master_guassian_list with values
          Master_Gaussian_list.Replicate(Shared_accross_replicate_counter,position_counter)=Replicate_MvsL_protein(1);
          Master_Gaussian_list.Center(Shared_accross_replicate_counter,position_counter)= Centres_MvsL_protein(1);
          Master_Gaussian_list.Height(Shared_accross_replicate_counter,position_counter)= Height_MvsL_protein(1);
          Master_Gaussian_list.Width(Shared_accross_replicate_counter,position_counter)= Width_MvsL_protein(1);
          Master_Gaussian_list.SSE(Shared_accross_replicate_counter,position_counter)= SSE_MvsL_protein(1);
          Master_Gaussian_list.adjrsquare(Shared_accross_replicate_counter,position_counter)= determined_AdjrSequare_of_guassians_MvsL(1);
          Master_Gaussian_list.Unique_identifier(Shared_accross_replicate_counter,position_counter)= Unique_identified_MvsL(1);
          Master_Gaussian_list.Replicate_Protein_identifier(Shared_accross_replicate_counter,position_counter)= Replicate_Protein_identifier_MvsL_protein(1);
          Master_Gaussian_list.Guassian_index_number(Shared_accross_replicate_counter,position_counter)=Guassian_index_number_MvsL_protein(1);
          Master_Gaussian_list.Protein_number(Shared_accross_replicate_counter,position_counter)= Protein_number_MvsL_protein(1);
          Master_Gaussian_list.Complex_size(Shared_accross_replicate_counter,position_counter)=Complex_size_MvsL(1);
          Master_Gaussian_list.Gaussian_area(Shared_accross_replicate_counter,position_counter)=Gaussian_area_MvsL(1);
          Master_Gaussian_list.Channel(Shared_accross_replicate_counter,position_counter)={'MvsL'};
          
          %remove values for arry being tested
          Replicate_MvsL_protein(1)=[];
          Centres_MvsL_protein(1)=[];
          Height_MvsL_protein(1)=[];
          Width_MvsL_protein(1)=[];
          SSE_MvsL_protein(1)=[];
          determined_AdjrSequare_of_guassians_MvsL(1)=[];
          Unique_identified_MvsL(1)=[];
          Guassian_index_number_MvsL_protein(1)=[];
          Protein_number_MvsL_protein(1)=[];
          Complex_size_MvsL(1)=[];
          Gaussian_area_MvsL(1)=[];
        end
        
        %test to find the best gaussians within user defined fractions of each other
      elseif min(Matrix_of_interactions(:)) <=User_Window
        % position of the guassians that are the same between MvsL and HvsL channel
        [row,col] = ind2sub(size(Matrix_of_interactions),find(ismember((Matrix_of_interactions),min(Matrix_of_interactions(:)))));
        %use the Gaussian corresponding to the best AdjrSequare
        
        
        %Determine which channel has the best gaussian fit
        AdjrSequare_of_guassians_channels=[determined_AdjrSequare_of_guassians_MvsL(row(1)), determined_AdjrSequare_of_guassians_Master_list(col(1))];
        Max_AdjrSequare_of_guassians=max(AdjrSequare_of_guassians_channels);
        position_of_best_AdjrSequare=find(ismember(AdjrSequare_of_guassians_channels,Max_AdjrSequare_of_guassians(1)));
        
        
        if position_of_best_AdjrSequare(1) == 1 % the MvsL contains the best adjrsequare, overwrite current value with MvsL values
          
          %as the col varible denotes the position in the master list of
          %the value of interest write the new MvsL to this position
          Master_Gaussian_list.Replicate(Shared_accross_replicate_counter,col)=Replicate_MvsL_protein(row(1));
          Master_Gaussian_list.Center(Shared_accross_replicate_counter,col)= Centres_MvsL_protein(row(1));
          Master_Gaussian_list.Height(Shared_accross_replicate_counter,col)= Height_MvsL_protein(row(1));
          Master_Gaussian_list.Width(Shared_accross_replicate_counter,col)= Width_MvsL_protein(row(1));
          Master_Gaussian_list.SSE(Shared_accross_replicate_counter,col)= SSE_MvsL_protein(row(1));
          Master_Gaussian_list.adjrsquare(Shared_accross_replicate_counter,col)= determined_AdjrSequare_of_guassians_MvsL(row(1));
          Master_Gaussian_list.Unique_identifier(Shared_accross_replicate_counter,col)= Unique_identified_MvsL(row(1));
          Master_Gaussian_list.Replicate_Protein_identifier(Shared_accross_replicate_counter,col)= Replicate_Protein_identifier_MvsL_protein(1);
          Master_Gaussian_list.Guassian_index_number(Shared_accross_replicate_counter,col)=Guassian_index_number_MvsL_protein(row(1));
          Master_Gaussian_list.Protein_number(Shared_accross_replicate_counter,col)= Protein_number_MvsL_protein(row(1));
          Master_Gaussian_list.Complex_size(Shared_accross_replicate_counter,col)=Complex_size_MvsL(row(1));
          Master_Gaussian_list.Gaussian_area(Shared_accross_replicate_counter,col)=Gaussian_area_MvsL(row(1));
          Master_Gaussian_list.Channel(Shared_accross_replicate_counter,col)={'MvsL'};
          
          %if the Master list contains the same guassian, within 2 fractions and
          %contains the best adjrsequare do nothing
          
        end
        
        %remove the most similar gaussains that are within 2 units
        %remove values in MvsL channels
        Replicate_MvsL_protein(:,[row(1)])=[];
        Centres_MvsL_protein(:,[row(1)])=[];
        Height_MvsL_protein(:,[row(1)])=[];
        Width_MvsL_protein(:,[row(1)])=[];
        SSE_MvsL_protein(:,[row(1)])=[];
        determined_AdjrSequare_of_guassians_MvsL(:,[row(1)])=[];
        Unique_identified_MvsL(:,[row(1)])=[];
        Guassian_index_number_MvsL_protein(:,[row(1)])=[];
        Protein_number_MvsL_protein(:,[row(1)])=[];
        Complex_size_MvsL(:,[row(1)])=[];
        Gaussian_area_MvsL(:,[row(1)])=[];
        
        %Consider gaussians not within user defined number of fractions of each other as unique
      elseif min(Matrix_of_interactions(:)) >User_Window
        
        %position of the guassians that are the same between MvsL and HvsL channel
        [row,col] = ind2sub(size(Matrix_of_interactions),find(ismember((Matrix_of_interactions),min(Matrix_of_interactions(:)))));
        
        %add one to position
        position_counter=position_counter+1;
        Master_Gaussian_list.Replicate(Shared_accross_replicate_counter,position_counter)=Replicate_MvsL_protein(row(1));
        Master_Gaussian_list.Center(Shared_accross_replicate_counter,position_counter)= Centres_MvsL_protein(row(1));
        Master_Gaussian_list.Height(Shared_accross_replicate_counter,position_counter)= Height_MvsL_protein(row(1));
        Master_Gaussian_list.Width(Shared_accross_replicate_counter,position_counter)= Width_MvsL_protein(row(1));
        Master_Gaussian_list.SSE(Shared_accross_replicate_counter,position_counter)= SSE_MvsL_protein(row(1));
        Master_Gaussian_list.adjrsquare(Shared_accross_replicate_counter,position_counter)= determined_AdjrSequare_of_guassians_MvsL(row(1));
        Master_Gaussian_list.Unique_identifier(Shared_accross_replicate_counter,position_counter)= Unique_identified_MvsL(row(1));
        Master_Gaussian_list.Replicate_Protein_identifier(Shared_accross_replicate_counter,position_counter)= Replicate_Protein_identifier_MvsL_protein(1);
        Master_Gaussian_list.Guassian_index_number(Shared_accross_replicate_counter,position_counter)=Guassian_index_number_MvsL_protein(row(1));
        Master_Gaussian_list.Protein_number(Shared_accross_replicate_counter,position_counter)= Protein_number_MvsL_protein(row(1));
        Master_Gaussian_list.Complex_size(Shared_accross_replicate_counter,position_counter)=Complex_size_MvsL(row(1));
        Master_Gaussian_list.Gaussian_area(Shared_accross_replicate_counter,position_counter)=Gaussian_area_MvsL(row(1));
        Master_Gaussian_list.Channel(Shared_accross_replicate_counter,position_counter)={'MvsL'};
        
        
        %remove values in MvsL channels
        Replicate_MvsL_protein(:,[row(1)])=[];
        Centres_MvsL_protein(:,[row(1)])=[];
        Height_MvsL_protein(:,[row(1)])=[];
        Width_MvsL_protein(:,[row(1)])=[];
        SSE_MvsL_protein(:,[row(1)])=[];
        determined_AdjrSequare_of_guassians_MvsL(:,[row(1)])=[];
        Unique_identified_MvsL(:,[row(1)])=[];
        Guassian_index_number_MvsL_protein(:,[row(1)])=[];
        Protein_number_MvsL_protein(:,[row(1)])=[];
        Complex_size_MvsL(:,[row(1)])=[];
        Gaussian_area_MvsL(:,[row(1)])=[];
        
      end
      
      %check the dimension of the matrix to all numbers have been checked
      row_check = nnz(Centres_MvsL_protein);
      ready = (row_check == 0);
      
    end
  end
end

for Shared_accross_replicate_counter = 1:number_of_unique_protein_with_gaussians
  
  %find location of proteins in HvsL arrays
  [internal_location_of_protein_of_interest_HvsL]=ind2sub(Dimension_of_array(1), strmatch(Unique_protein_names(Shared_accross_replicate_counter), HvsL_Gaussians.Protein_name, 'exact'));
  [internal_location_of_protein_of_interest_Master]=ind2sub(Dimension_of_array(1), strmatch(Unique_protein_names(Shared_accross_replicate_counter), Master_Gaussian_list.Protein_name, 'exact'));
  
  %Compare gaussians from each replicate channel
  Number_of_times_seen_HvsL=length(internal_location_of_protein_of_interest_HvsL);
  
  %Compare HvsL Gaussians
  for HvsL_counter=1:Number_of_times_seen_HvsL
    
    %determine Center for specific protein, note this is created so as the loop
    %cycle through each guassian the identified gaussians can be removed
    %Copy replicate data
    Replicate_HvsL_protein=HvsL_Gaussians.Replicate(internal_location_of_protein_of_interest_HvsL(HvsL_counter), :);
    Replicate_Master_list_protein=Master_Gaussian_list.Replicate(internal_location_of_protein_of_interest_Master, :);
    %Copy Unique_identifier
    Unique_identified_HvsL=HvsL_Gaussians.Unique_identifier(internal_location_of_protein_of_interest_HvsL(HvsL_counter), :);
    Unique_identified_Master_list_protein=Master_Gaussian_list.Unique_identifier(internal_location_of_protein_of_interest_Master, :);
    %Copy Protein number
    Protein_number_HvsL_protein=HvsL_Gaussians.Protein_number(internal_location_of_protein_of_interest_HvsL(HvsL_counter), :);
    Protein_number_Master_list_protein=Master_Gaussian_list.Protein_number(internal_location_of_protein_of_interest_Master, :);
    %Copy Guassian_index_number
    Guassian_index_number_HvsL_protein=HvsL_Gaussians.Guassian_index_number(internal_location_of_protein_of_interest_HvsL(HvsL_counter), :);
    Guassian_index_number_Master_list_protein=Master_Gaussian_list.Guassian_index_number(internal_location_of_protein_of_interest_Master, :);
    %Copy Replicate_Protein_identifier
    Replicate_Protein_identifier_HvsL_protein=HvsL_Gaussians.Replicate_Protein_identifier(internal_location_of_protein_of_interest_HvsL(HvsL_counter), :);
    Replicate_Protein_identifier_Master_list_protein=Master_Gaussian_list.Replicate_Protein_identifier(internal_location_of_protein_of_interest_Master, :);
    %Copy Center data
    Centres_HvsL_protein=HvsL_Gaussians.Center(internal_location_of_protein_of_interest_HvsL(HvsL_counter), :);
    Centres_Master_list_protein=Master_Gaussian_list.Center(internal_location_of_protein_of_interest_Master, :);
    %Copy Height data
    Height_HvsL_protein=HvsL_Gaussians.Height(internal_location_of_protein_of_interest_HvsL(HvsL_counter), :);
    Height_Master_list_protein=Master_Gaussian_list.Height(internal_location_of_protein_of_interest_Master, :);
    %Copy Width data
    Width_HvsL_protein=HvsL_Gaussians.Width(internal_location_of_protein_of_interest_HvsL(HvsL_counter), :);
    Width_Master_list_protein=Master_Gaussian_list.Width(internal_location_of_protein_of_interest_Master, :);
    %Copy SSE data
    SSE_HvsL_protein=HvsL_Gaussians.SSE(internal_location_of_protein_of_interest_HvsL(HvsL_counter), :);
    SSE_Master_list_protein=Master_Gaussian_list.SSE(internal_location_of_protein_of_interest_Master, :);
    %Copy adjrsquare data
    determined_AdjrSequare_of_guassians_HvsL=HvsL_Gaussians.adjrsquare(internal_location_of_protein_of_interest_HvsL(HvsL_counter), :);
    determined_AdjrSequare_of_guassians_Master_list=Master_Gaussian_list.adjrsquare(internal_location_of_protein_of_interest_Master, :);
    %Copy Complex size
    Complex_size_HvsL=HvsL_Gaussians.Complex_size(internal_location_of_protein_of_interest_HvsL(HvsL_counter), :);
    Complex_size_Master_list=Master_Gaussian_list.Complex_size(internal_location_of_protein_of_interest_Master, :);
    %Copy Gaussian area
    Gaussian_area_HvsL=HvsL_Gaussians.Gaussian_area(internal_location_of_protein_of_interest_HvsL(HvsL_counter), :);
    Gaussian_area_Master_list=Master_Gaussian_list.Gaussian_area(internal_location_of_protein_of_interest_Master, :);
    
    
    %Determine what position in Master_Gaussian_list.Center to write to
    if Master_Gaussian_list.Center(internal_location_of_protein_of_interest_Master,1)== 0
      position_counter=0;
    elseif ~Master_Gaussian_list.Center(internal_location_of_protein_of_interest_Master,1) == 0
      position_counter=(nnz(Master_Gaussian_list.Center(internal_location_of_protein_of_interest_Master, :)));
    end
    
    ready = false;
    while ~ready
      %count the number of non zero elements
      Number_of_non_zero_elements_HvsL=nnz(Centres_HvsL_protein);
      Number_of_non_zero_elements_Master_list=nnz(Centres_Master_list_protein);
      %create the matrix to be filled
      Matrix_of_interactions=zeros(Number_of_non_zero_elements_HvsL,Number_of_non_zero_elements_Master_list);
      
      
      %create matrix to compare gaussians
      for internal_counter6= 1:Number_of_non_zero_elements_HvsL
        for internal_counter7= 1:Number_of_non_zero_elements_Master_list
          Centres_HvsL_for_matrix=Centres_HvsL_protein(internal_counter6);
          Centres_Master_list_protein_for_matrix=Centres_Master_list_protein(internal_counter7);
          Matrix_of_interactions(internal_counter6,internal_counter7)= abs(Centres_HvsL_for_matrix-Centres_Master_list_protein_for_matrix);
        end
      end
      
      
      %test matrix to determine which gaussains to keep and which ones to remove
      %test to makre sure the matrix has numbers
      if isempty(Matrix_of_interactions(:)) == 1
        row_to_copy_to_summed_array = nnz(Centres_HvsL_protein);
        col_to_copy_to_summed_array = nnz(Centres_Master_list_protein);
        if row_to_copy_to_summed_array>0
          
          %add one to position counter
          position_counter=position_counter+1;
          
          %populate Master_guassian_list with values
          Master_Gaussian_list.Replicate(Shared_accross_replicate_counter,position_counter)=Replicate_HvsL_protein(1);
          Master_Gaussian_list.Center(Shared_accross_replicate_counter,position_counter)= Centres_HvsL_protein(1);
          Master_Gaussian_list.Height(Shared_accross_replicate_counter,position_counter)= Height_HvsL_protein(1);
          Master_Gaussian_list.Width(Shared_accross_replicate_counter,position_counter)= Width_HvsL_protein(1);
          Master_Gaussian_list.SSE(Shared_accross_replicate_counter,position_counter)= SSE_HvsL_protein(1);
          Master_Gaussian_list.adjrsquare(Shared_accross_replicate_counter,position_counter)= determined_AdjrSequare_of_guassians_HvsL(1);
          Master_Gaussian_list.Unique_identifier(Shared_accross_replicate_counter,position_counter)= Unique_identified_HvsL(1);
          Master_Gaussian_list.Replicate_Protein_identifier(Shared_accross_replicate_counter,position_counter)= Replicate_Protein_identifier_HvsL_protein(1);
          Master_Gaussian_list.Guassian_index_number(Shared_accross_replicate_counter,position_counter)=Guassian_index_number_HvsL_protein(1);
          Master_Gaussian_list.Protein_number(Shared_accross_replicate_counter,position_counter)= Protein_number_HvsL_protein(1);
          Master_Gaussian_list.Complex_size(Shared_accross_replicate_counter,position_counter)=Complex_size_HvsL(1);
          Master_Gaussian_list.Gaussian_area(Shared_accross_replicate_counter,position_counter)=Gaussian_area_HvsL(1);
          Master_Gaussian_list.Channel(Shared_accross_replicate_counter,position_counter)={'HvsL'};
          
          
          %remove values for arry being tested
          Replicate_HvsL_protein(1)=[];
          Centres_HvsL_protein(1)=[];
          Height_HvsL_protein(1)=[];
          Width_HvsL_protein(1)=[];
          SSE_HvsL_protein(1)=[];
          determined_AdjrSequare_of_guassians_HvsL(1)=[];
          Unique_identified_HvsL(1)=[];
          Guassian_index_number_HvsL_protein(1)=[];
          Protein_number_HvsL_protein(1)=[];
          Complex_size_HvsL(1)=[];
          Gaussian_area_HvsL(1)=[];
          
        end
        
        %test to find the best gaussians within user defined number of fractions of each other
      elseif min(Matrix_of_interactions(:)) <=User_Window
        % position of the guassians that are the same between MvsL and HvsL channel
        [row,col] = ind2sub(size(Matrix_of_interactions),find(ismember((Matrix_of_interactions),min(Matrix_of_interactions(:)))));
        %use the Gaussian corresponding to the best AdjrSequare
        
        
        %Determine which channel has the best gaussian fit
        AdjrSequare_of_guassians_channels=[determined_AdjrSequare_of_guassians_HvsL(row(1)), determined_AdjrSequare_of_guassians_Master_list(col(1))];
        Max_AdjrSequare_of_guassians=max(AdjrSequare_of_guassians_channels);
        position_of_best_AdjrSequare=find(ismember(AdjrSequare_of_guassians_channels,Max_AdjrSequare_of_guassians(1)));
        
        
        if position_of_best_AdjrSequare(1) == 1 % the MvsL contains the best adjrsequare, overwrite current value with MvsL values
          
          %as the col varible denotes the position in the master list of
          %the value of interest write the new MvsL to this position
          Master_Gaussian_list.Replicate(Shared_accross_replicate_counter,col)=Replicate_HvsL_protein(row(1));
          Master_Gaussian_list.Center(Shared_accross_replicate_counter,col)= Centres_HvsL_protein(row(1));
          Master_Gaussian_list.Height(Shared_accross_replicate_counter,col)= Height_HvsL_protein(row(1));
          Master_Gaussian_list.Width(Shared_accross_replicate_counter,col)= Width_HvsL_protein(row(1));
          Master_Gaussian_list.SSE(Shared_accross_replicate_counter,col)= SSE_HvsL_protein(row(1));
          Master_Gaussian_list.adjrsquare(Shared_accross_replicate_counter,col)= determined_AdjrSequare_of_guassians_HvsL(row(1));
          Master_Gaussian_list.Unique_identifier(Shared_accross_replicate_counter,col)= Unique_identified_HvsL(row(1));
          Master_Gaussian_list.Replicate_Protein_identifier(Shared_accross_replicate_counter,col)= Replicate_Protein_identifier_HvsL_protein(1);
          Master_Gaussian_list.Guassian_index_number(Shared_accross_replicate_counter,col)=Guassian_index_number_HvsL_protein(row(1));
          Master_Gaussian_list.Protein_number(Shared_accross_replicate_counter,col)= Protein_number_HvsL_protein(row(1));
          Master_Gaussian_list.Complex_size(Shared_accross_replicate_counter,col)=Complex_size_HvsL(row(1));
          Master_Gaussian_list.Gaussian_area(Shared_accross_replicate_counter,col)=Gaussian_area_HvsL(row(1));
          Master_Gaussian_list.Channel(Shared_accross_replicate_counter,col)={'HvsL'};
          
          
          % if the Master list contains the same guassian, within 2 fractions and
          % contains the best adjrsequare do nothing
          
        end
        
        %remove the most similar gaussains that are within 2 units
        %remove values in HvsL channels
        Replicate_HvsL_protein(:,[row(1)])=[];
        Centres_HvsL_protein(:,[row(1)])=[];
        Height_HvsL_protein(:,[row(1)])=[];
        Width_HvsL_protein(:,[row(1)])=[];
        SSE_HvsL_protein(:,[row(1)])=[];
        determined_AdjrSequare_of_guassians_HvsL(:,[row(1)])=[];
        Unique_identified_HvsL(:,[row(1)])=[];
        Guassian_index_number_HvsL_protein(:,[row(1)])=[];
        Protein_number_HvsL_protein(:,[row(1)])=[];
        Complex_size_HvsL(:,[row(1)])=[];
        Gaussian_area_HvsL(:,[row(1)])=[];
        
        %Consider gaussians not within user defined number of fractions of each other as unique
      elseif min(Matrix_of_interactions(:)) >User_Window
        %position of the guassians that are the same between MvsL and HvsL channel
        [row,col] = ind2sub(size(Matrix_of_interactions),find(ismember((Matrix_of_interactions),min(Matrix_of_interactions(:)))));
        
        %add one to position
        position_counter=position_counter+1;
        
        %populate Master_guassian_list with values
        Master_Gaussian_list.Replicate(Shared_accross_replicate_counter,position_counter)=Replicate_HvsL_protein(row(1));
        Master_Gaussian_list.Center(Shared_accross_replicate_counter,position_counter)= Centres_HvsL_protein(row(1));
        Master_Gaussian_list.Height(Shared_accross_replicate_counter,position_counter)= Height_HvsL_protein(row(1));
        Master_Gaussian_list.Width(Shared_accross_replicate_counter,position_counter)= Width_HvsL_protein(row(1));
        Master_Gaussian_list.SSE(Shared_accross_replicate_counter,position_counter)= SSE_HvsL_protein(row(1));
        Master_Gaussian_list.adjrsquare(Shared_accross_replicate_counter,position_counter)= determined_AdjrSequare_of_guassians_HvsL(row(1));
        Master_Gaussian_list.Unique_identifier(Shared_accross_replicate_counter,position_counter)= Unique_identified_HvsL(row(1));
        Master_Gaussian_list.Replicate_Protein_identifier(Shared_accross_replicate_counter,position_counter)= Replicate_Protein_identifier_HvsL_protein(1);
        Master_Gaussian_list.Guassian_index_number(Shared_accross_replicate_counter,position_counter)=Guassian_index_number_HvsL_protein(row(1));
        Master_Gaussian_list.Protein_number(Shared_accross_replicate_counter,position_counter)= Protein_number_HvsL_protein(row(1));
        Master_Gaussian_list.Complex_size(Shared_accross_replicate_counter,position_counter)=Complex_size_HvsL(row(1));
        Master_Gaussian_list.Gaussian_area(Shared_accross_replicate_counter,position_counter)=Gaussian_area_HvsL(row(1));
        Master_Gaussian_list.Channel(Shared_accross_replicate_counter,position_counter)={'HvsL'};
        
        
        %remove values in HvsL channels
        Replicate_HvsL_protein(:,[row(1)])=[];
        Centres_HvsL_protein(:,[row(1)])=[];
        Height_HvsL_protein(:,[row(1)])=[];
        Width_HvsL_protein(:,[row(1)])=[];
        SSE_HvsL_protein(:,[row(1)])=[];
        determined_AdjrSequare_of_guassians_HvsL(:,[row(1)])=[];
        Unique_identified_HvsL(:,[row(1)])=[];
        Guassian_index_number_HvsL_protein(:,[row(1)])=[];
        Protein_number_HvsL_protein(:,[row(1)])=[];
        Complex_size_HvsL(:,[row(1)])=[];
        Gaussian_area_HvsL(:,[row(1)])=[];
        
      end
      
      %check the dimension of the matrix to all numbers have been checked
      row_check = nnz(Centres_HvsL_protein);
      ready = (row_check == 0);
      
    end
  end
  
end


%Check master_gaussian list to ensure all gaussian measurements are not
%within two fractions of each other

for Shared_accross_replicate_counter1 = 1:number_of_unique_protein_with_gaussians
  
  %Copy Center data
  Centres_Master_list_protein=Master_Gaussian_list.Center(Shared_accross_replicate_counter1, :);
  
  ready = false;
  
  while ~ready
    %count the number of non zero elements
    Number_of_non_zero_elements_Master_list=nnz(Centres_Master_list_protein);
    %create the matrix to be filled
    Matrix_of_interactions=zeros(Number_of_non_zero_elements_Master_list,Number_of_non_zero_elements_Master_list);
    
    
    %create matrix to compare gaussians
    for internal_counter6= 1:Number_of_non_zero_elements_Master_list
      for internal_counter7= 1:Number_of_non_zero_elements_Master_list
        Centres_Master_list_protein_for_matrix1=Centres_Master_list_protein(internal_counter6);
        Centres_Master_list_protein_for_matrix2=Centres_Master_list_protein(internal_counter7);
        Matrix_of_interactions(internal_counter6,internal_counter7)= abs(Centres_Master_list_protein_for_matrix1-Centres_Master_list_protein_for_matrix2);
      end
    end
    
    %Take upper right section and remove diagonal (zero) elements
    Matrix_of_interactions=triu(Matrix_of_interactions);
    Matrix_of_interactions(Matrix_of_interactions==0)=NaN;
    
    %test to find the best gaussians within user defined number of fractions of each other
    if nanmin(Matrix_of_interactions(:)) <=User_Window
      % position of the guassians that are the same between MvsL and HvsL channel
      [row,col] = ind2sub(size(Matrix_of_interactions),find(ismember((Matrix_of_interactions),min(Matrix_of_interactions(:)))));
      %use the Gaussian corresponding to the best AdjrSequare
      
      %Determine which channel has the best gaussian fit
      AdjrSequare_of_guassians_channels=[determined_AdjrSequare_of_guassians_Master_list(row(1)), determined_AdjrSequare_of_guassians_Master_list(col(1))];
      Max_AdjrSequare_of_guassians=max(AdjrSequare_of_guassians_channels);
      position_of_best_AdjrSequare=find(ismember(AdjrSequare_of_guassians_channels,Max_AdjrSequare_of_guassians(1)));
      
      if position_of_best_AdjrSequare(1) == 1 % the row contains the best adjrsequare, remove col
        
        Master_Gaussian_list.Replicate(Shared_accross_replicate_counter1,[col(1)])= 0;
        Master_Gaussian_list.Center(Shared_accross_replicate_counter1,[col(1)])= 0;
        Master_Gaussian_list.Height(Shared_accross_replicate_counter1,[col(1)])= 0;
        Master_Gaussian_list.Width(Shared_accross_replicate_counter1,[col(1)])= 0;
        Master_Gaussian_list.SSE(Shared_accross_replicate_counter1,[col(1)])= 0;
        Master_Gaussian_list.adjrsquare(Shared_accross_replicate_counter1,[col(1)])= 0;
        Master_Gaussian_list.Unique_identifier{Shared_accross_replicate_counter1,[col(1)]}= [];
        Master_Gaussian_list.Replicate_Protein_identifier{Shared_accross_replicate_counter1,[col(1)]}= [];
        Master_Gaussian_list.Guassian_index_number(Shared_accross_replicate_counter1,[col(1)])=0;
        Master_Gaussian_list.Protein_number(Shared_accross_replicate_counter1,[col(1)])= 0;
        Master_Gaussian_list.Complex_size(Shared_accross_replicate_counter1,[col(1)])=0;
        Master_Gaussian_list.Gaussian_area(Shared_accross_replicate_counter1,[col(1)])=0;
        Master_Gaussian_list.Channel{Shared_accross_replicate_counter1,[col(1)]}=[];
        
        Centres_Master_list_protein(:,[col(1)])=[];
        
      elseif position_of_best_AdjrSequare(2) == 1 % the col contains the best adjrsequare, remove row
        
        Master_Gaussian_list.Replicate(Shared_accross_replicate_counter1,[row(1)])= 0;
        Master_Gaussian_list.Center(Shared_accross_replicate_counter1,[row(1)])= 0;
        Master_Gaussian_list.Height(Shared_accross_replicate_counter1,[row(1)])= 0;
        Master_Gaussian_list.Width(Shared_accross_replicate_counter1,[row(1)])= 0;
        Master_Gaussian_list.SSE(Shared_accross_replicate_counter1,[row(1)])= 0;
        Master_Gaussian_list.adjrsquare(Shared_accross_replicate_counter1,[row(1)])= 0;
        Master_Gaussian_list.Unique_identifier{Shared_accross_replicate_counter1,[row(1)]}= [];
        Master_Gaussian_list.Replicate_Protein_identifier{Shared_accross_replicate_counter1,[row(1)]}= [];
        Master_Gaussian_list.Guassian_index_number(Shared_accross_replicate_counter1,[row(1)])=0;
        Master_Gaussian_list.Protein_number(Shared_accross_replicate_counter1,[row(1)])= 0;
        Master_Gaussian_list.Complex_size(Shared_accross_replicate_counter1,[row(1)])=0;
        Master_Gaussian_list.Gaussian_area(Shared_accross_replicate_counter1,[row(1)])=0;
        Master_Gaussian_list.Channel{Shared_accross_replicate_counter1,[row(1)]}=[];
        
        Centres_Master_list_protein(:,[row(1)])=[];
        
      end
      
      %test if the minimum value is greater then user defined window
    elseif min(Matrix_of_interactions(:)) >User_Window
      
      ready = true;
      
      %test if all elements are equal to NaN
    elseif all(isnan(Matrix_of_interactions(:)) == 1)
      
      ready = true;
      
    end
  end
  
  
end

%Clean up Master list for downstream processing
for Shared_accross_replicate_counter1 = 1:number_of_unique_protein_with_gaussians
  %Count how many nnz elements are observed in the protein
  for cleanup_counter=1:nnz(Master_Gaussian_list.Replicate(Shared_accross_replicate_counter1,:))
    %Work out the position of the nnz elements
    [index_cleanUp]= find(Master_Gaussian_list.Replicate(Shared_accross_replicate_counter1,:));
    
    
    Finalised_Master_Gaussian_list.Replicate(Shared_accross_replicate_counter1,cleanup_counter)= Master_Gaussian_list.Replicate(Shared_accross_replicate_counter1,index_cleanUp(cleanup_counter));
    Finalised_Master_Gaussian_list.Center(Shared_accross_replicate_counter1,cleanup_counter)= Master_Gaussian_list.Center(Shared_accross_replicate_counter1,index_cleanUp(cleanup_counter));
    Finalised_Master_Gaussian_list.Height(Shared_accross_replicate_counter1,cleanup_counter)= Master_Gaussian_list.Height(Shared_accross_replicate_counter1,index_cleanUp(cleanup_counter));
    Finalised_Master_Gaussian_list.Width(Shared_accross_replicate_counter1,cleanup_counter)= Master_Gaussian_list.Width(Shared_accross_replicate_counter1,index_cleanUp(cleanup_counter));
    Finalised_Master_Gaussian_list.SSE(Shared_accross_replicate_counter1,cleanup_counter)= Master_Gaussian_list.SSE(Shared_accross_replicate_counter1,index_cleanUp(cleanup_counter));
    Finalised_Master_Gaussian_list.adjrsquare(Shared_accross_replicate_counter1,cleanup_counter)= Master_Gaussian_list.adjrsquare(Shared_accross_replicate_counter1,index_cleanUp(cleanup_counter));
    Finalised_Master_Gaussian_list.Unique_identifier{Shared_accross_replicate_counter1,cleanup_counter}= Master_Gaussian_list.Unique_identifier{Shared_accross_replicate_counter1,index_cleanUp(cleanup_counter)};
    Finalised_Master_Gaussian_list.Replicate_Protein_identifier{Shared_accross_replicate_counter1,cleanup_counter}= Master_Gaussian_list.Replicate_Protein_identifier{Shared_accross_replicate_counter1,index_cleanUp(cleanup_counter)};
    Finalised_Master_Gaussian_list.Guassian_index_number(Shared_accross_replicate_counter1,cleanup_counter)= Master_Gaussian_list.Guassian_index_number(Shared_accross_replicate_counter1,index_cleanUp(cleanup_counter));
    Finalised_Master_Gaussian_list.Protein_number(Shared_accross_replicate_counter1,cleanup_counter)= Master_Gaussian_list.Protein_number(Shared_accross_replicate_counter1,index_cleanUp(cleanup_counter));
    Finalised_Master_Gaussian_list.Complex_size(Shared_accross_replicate_counter1,cleanup_counter)= Master_Gaussian_list.Complex_size(Shared_accross_replicate_counter1,index_cleanUp(cleanup_counter));
    Finalised_Master_Gaussian_list.Gaussian_area(Shared_accross_replicate_counter1,cleanup_counter)= Master_Gaussian_list.Gaussian_area(Shared_accross_replicate_counter1,index_cleanUp(cleanup_counter));
    Finalised_Master_Gaussian_list.Channel{Shared_accross_replicate_counter1,cleanup_counter}= Master_Gaussian_list.Channel{Shared_accross_replicate_counter1,index_cleanUp(cleanup_counter)};
  end
  
end

%copy names to finalised master list
Finalised_Master_Gaussian_list.Protein_name=Master_Gaussian_list.Protein_name;

%Count the number of gassians in the master list
Total_number_of_unique_gaussians_master_list=nnz(Finalised_Master_Gaussian_list.Center);
Dimension_of_master_gaussian_list=size(Finalised_Master_Gaussian_list.Center);

%Create Chromatogram file containing the raw isotoplogue data for each guassian that has been kept
Chromatograms=zeros(1,fraction_number(2));
Chromatogram_counter=1;

for writeout_counter1= 1:number_of_unique_protein_with_gaussians
  for writeout_counter2= 1:nnz(Finalised_Master_Gaussian_list.Center(writeout_counter1,:))
    
    if isequal(Finalised_Master_Gaussian_list.Channel{writeout_counter1,writeout_counter2}, 'MvsL')
      [internal_location_of_protein_of_interest_in_MvsL_raw_data]=strmatch(Finalised_Master_Gaussian_list.Replicate_Protein_identifier{writeout_counter1,writeout_counter2}, txt_MvsL(:,1), 'exact');
      Chromatograms(Chromatogram_counter,:)=num_val_MvsL((internal_location_of_protein_of_interest_in_MvsL_raw_data-1),(2:(fraction_number(2)+1)));
      Chromatogram_counter=1 + Chromatogram_counter;
      
    elseif isequal(Finalised_Master_Gaussian_list.Channel{writeout_counter1,writeout_counter2}, 'HvsL')
      [internal_location_of_protein_of_interest_in_HvsL_raw_data]=strmatch(Finalised_Master_Gaussian_list.Replicate_Protein_identifier{writeout_counter1,writeout_counter2}, txt_HvsL(:,1), 'exact');
      Chromatograms(Chromatogram_counter,:)=num_val_HvsL((internal_location_of_protein_of_interest_in_HvsL_raw_data-1),(2:(fraction_number(2)+1)));
      Chromatogram_counter=1 + Chromatogram_counter;
      
    end
    
  end
end

%Format Name Names, Unique identifier and Guassian index number to write out

Write_out_counter=1;
for writeout_counter1= 1:number_of_unique_protein_with_gaussians
  for writeout_counter2= 1:nnz(Finalised_Master_Gaussian_list.Center(writeout_counter1,:))
    Write_out_unique_identifier{Write_out_counter}=Finalised_Master_Gaussian_list.Unique_identifier{writeout_counter1,writeout_counter2};
    Write_out_Guassian_index_number(Write_out_counter)=Finalised_Master_Gaussian_list.Guassian_index_number(writeout_counter1,writeout_counter2);
    Write_out_Protein_name{Write_out_counter}=Finalised_Master_Gaussian_list.Protein_name{writeout_counter1};
    Write_out_counter=1+Write_out_counter;
  end
end

% #output
% write out table of master guassians list
fid_master_list= fopen([datadir2 'Master_guassian_list.csv'],'wt'); % create the output file with the header infromation
fprintf (fid_master_list,'%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n',... %header for OutputGaus output
  'Protein name', 'Unique_identifier (of best gaussian)','Replicate (of best gaussian)',...
  'Channel (of best gaussian)', 'Guassian_index_number (of best gaussian)',...
  'Center (of best gaussian)', 'Height (of best gaussian)', 'Width (of best gaussian)',...
  'SSE (of best gaussian)', 'adjrsquare (of best gaussian)','Gaussian Area (of best gaussian)','Complex Size (of best gaussian)'); %Write Header
for writeout_counter1= 1:number_of_unique_protein_with_gaussians
  for writeout_counter2= 1:nnz(Finalised_Master_Gaussian_list.Center(writeout_counter1,:))
    fprintf(fid_master_list,'%s,%s,%6.3f,%s,%6.3f,%6.3f,%6.3f,%6.3f,%6.3f,%6.3f,%6.3f,%6.3f,\n',...
      Finalised_Master_Gaussian_list.Protein_name{writeout_counter1},...
      Finalised_Master_Gaussian_list.Unique_identifier{writeout_counter1,writeout_counter2},...
      Finalised_Master_Gaussian_list.Replicate(writeout_counter1,writeout_counter2),...
      Finalised_Master_Gaussian_list.Channel{writeout_counter1,writeout_counter2},...
      Finalised_Master_Gaussian_list.Guassian_index_number(writeout_counter1,writeout_counter2),...
      Finalised_Master_Gaussian_list.Center(writeout_counter1,writeout_counter2),...
      Finalised_Master_Gaussian_list.Height(writeout_counter1,writeout_counter2),...
      Finalised_Master_Gaussian_list.Width(writeout_counter1,writeout_counter2),...
      Finalised_Master_Gaussian_list.SSE(writeout_counter1,writeout_counter2),...
      Finalised_Master_Gaussian_list.adjrsquare(writeout_counter1,writeout_counter2),...
      Finalised_Master_Gaussian_list.Gaussian_area(writeout_counter1,writeout_counter2),...
      Finalised_Master_Gaussian_list.Complex_size(writeout_counter1,writeout_counter2));
  end
end
fclose(fid_master_list);

fid_master_chromatogram_list= fopen([datadir2 'Master_Chromatogram_list.csv'],'wt'); % create the output file with the header infromation
for writeout_counter1= 1:Total_number_of_unique_gaussians_master_list
  
  %write out Protein name and unique identifiers
  fprintf(fid_master_chromatogram_list,'%s,%6.4f,%s,',...
    Write_out_unique_identifier{writeout_counter1},...
    Write_out_Guassian_index_number(writeout_counter1),...
    Write_out_Protein_name{writeout_counter1});
  
  %write out raw data to use for euclean distance
  fprintf(fid_master_chromatogram_list,'%6.4g,', Chromatograms(writeout_counter1,:));
  fprintf(fid_master_chromatogram_list,'\n');
  
end
fclose(fid_master_chromatogram_list);

tt = toc;
fprintf('  ...  %.2f seconds\n',tt)



%% 7. Use fold-change > 2 to determine if complexes change
fprintf('    7. Use fold-change > 2 to determine if complexes change (BETWEEN REPLICATEs)')


% Comparsion 3 Part B: 
% Aim: identify changes by using the master list and comparing if a gaussian is
% observed within two fractions, if so use the center of that Gaussian

%Create arrays to write apex values to
MvsL_Ratio_Apex=zeros(1,5);
HvsL_Ratio_Apex=zeros(1,5);
HvsM_Ratio_Apex=zeros(1,5);
Apex_values_replicate_1=zeros(Total_number_of_unique_gaussians_master_list,4);
Apex_values_replicate_2=zeros(Total_number_of_unique_gaussians_master_list,4);
Apex_values_replicate_3=zeros(Total_number_of_unique_gaussians_master_list,4);
Apex_values_replicate_4=zeros(Total_number_of_unique_gaussians_master_list,4);
Increase_at_apex_in_stimulation_replicate_1=zeros(Dimension_of_master_gaussian_list(1),Dimension_of_master_gaussian_list(2));
Increase_at_apex_in_stimulation_replicate_2=zeros(Dimension_of_master_gaussian_list(1),Dimension_of_master_gaussian_list(2));
Increase_at_apex_in_stimulation_replicate_3=zeros(Dimension_of_master_gaussian_list(1),Dimension_of_master_gaussian_list(2));
Increase_at_apex_in_stimulation_replicate_4=zeros(Dimension_of_master_gaussian_list(1),Dimension_of_master_gaussian_list(2));
Decrease_at_apex_in_stimulation_replicate_1=zeros(Dimension_of_master_gaussian_list(1),Dimension_of_master_gaussian_list(2));
Decrease_at_apex_in_stimulation_replicate_2=zeros(Dimension_of_master_gaussian_list(1),Dimension_of_master_gaussian_list(2));
Decrease_at_apex_in_stimulation_replicate_3=zeros(Dimension_of_master_gaussian_list(1),Dimension_of_master_gaussian_list(2));
Decrease_at_apex_in_stimulation_replicate_4=zeros(Dimension_of_master_gaussian_list(1),Dimension_of_master_gaussian_list(2));
log_2_fold_comparsion_of_gaussian_replicate_1=zeros(Dimension_of_master_gaussian_list(1),Dimension_of_master_gaussian_list(2));
log_2_fold_comparsion_of_gaussian_replicate_2=zeros(Dimension_of_master_gaussian_list(1),Dimension_of_master_gaussian_list(2));
log_2_fold_comparsion_of_gaussian_replicate_3=zeros(Dimension_of_master_gaussian_list(1),Dimension_of_master_gaussian_list(2));
log_2_fold_comparsion_of_gaussian_replicate_4=zeros(Dimension_of_master_gaussian_list(1),Dimension_of_master_gaussian_list(2));
Do_gaussian_in_channels_change_for_figure_replicate_1=zeros(Dimension_of_master_gaussian_list(1),Dimension_of_master_gaussian_list(2));
Do_gaussian_in_channels_change_for_figure_replicate_2=zeros(Dimension_of_master_gaussian_list(1),Dimension_of_master_gaussian_list(2));
Do_gaussian_in_channels_change_for_figure_replicate_3=zeros(Dimension_of_master_gaussian_list(1),Dimension_of_master_gaussian_list(2));
Do_gaussian_in_channels_change_for_figure_replicate_4=zeros(Dimension_of_master_gaussian_list(1),Dimension_of_master_gaussian_list(2));
Do_gaussian_in_channels_change_replicate_1=cell(Dimension_of_master_gaussian_list(1),Dimension_of_master_gaussian_list(2));
Do_gaussian_in_channels_change_replicate_2=cell(Dimension_of_master_gaussian_list(1),Dimension_of_master_gaussian_list(2));
Do_gaussian_in_channels_change_replicate_3=cell(Dimension_of_master_gaussian_list(1),Dimension_of_master_gaussian_list(2));
Do_gaussian_in_channels_change_replicate_4=cell(Dimension_of_master_gaussian_list(1),Dimension_of_master_gaussian_list(2));
Protein_location_replicate1=zeros(Dimension_of_master_gaussian_list(1),Dimension_of_master_gaussian_list(2));
Protein_location_replicate2=zeros(Dimension_of_master_gaussian_list(1),Dimension_of_master_gaussian_list(2));
Protein_location_replicate3=zeros(Dimension_of_master_gaussian_list(1),Dimension_of_master_gaussian_list(2));
Protein_location_replicate4=zeros(Dimension_of_master_gaussian_list(1),Dimension_of_master_gaussian_list(2));
Do_gaussian_in_channels_position=1;

%Counters to detect changes
Increase_counter_replicate_1=0;
Decrease_counter_replicate_1=0;
No_change_counter_replicate_1=0;
Increase_counter_replicate_2=0;
Decrease_counter_replicate_2=0;
No_change_counter_replicate_2=0;
Increase_counter_replicate_3=0;
Decrease_counter_replicate_3=0;
No_change_counter_replicate_3=0;
Increase_counter_replicate_4=0;
Decrease_counter_replicate_4=0;
No_change_counter_replicate_4=0;
Unquantifiable_observed_only_in_MvsL_channel_replicate_1=0;
Unquantifiable_observed_only_in_MvsL_channel_replicate_2=0;
Unquantifiable_observed_only_in_MvsL_channel_replicate_3=0;
Unquantifiable_observed_only_in_MvsL_channel_replicate_4=0;
Unquantifiable_observed_only_in_HvsL_channel_replicate_1=0;
Unquantifiable_observed_only_in_HvsL_channel_replicate_2=0;
Unquantifiable_observed_only_in_HvsL_channel_replicate_3=0;
Unquantifiable_observed_only_in_HvsL_channel_replicate_4=0;
Unquantifiable_Not_observed_in_replicate_replicate_1=0;
Unquantifiable_Not_observed_in_replicate_replicate_2=0;
Unquantifiable_Not_observed_in_replicate_replicate_3=0;
Unquantifiable_Not_observed_in_replicate_replicate_4=0;
MvsL_Ratio_replicate_1=cell(Dimension_of_master_gaussian_list(1),Dimension_of_master_gaussian_list(2));
MvsL_Ratio_replicate_2=cell(Dimension_of_master_gaussian_list(1),Dimension_of_master_gaussian_list(2));
MvsL_Ratio_replicate_3=cell(Dimension_of_master_gaussian_list(1),Dimension_of_master_gaussian_list(2));
MvsL_Ratio_replicate_4=cell(Dimension_of_master_gaussian_list(1),Dimension_of_master_gaussian_list(2));
HvsL_Ratio_replicate_1=cell(Dimension_of_master_gaussian_list(1),Dimension_of_master_gaussian_list(2));
HvsL_Ratio_replicate_2=cell(Dimension_of_master_gaussian_list(1),Dimension_of_master_gaussian_list(2));
HvsL_Ratio_replicate_3=cell(Dimension_of_master_gaussian_list(1),Dimension_of_master_gaussian_list(2));
HvsL_Ratio_replicate_4=cell(Dimension_of_master_gaussian_list(1),Dimension_of_master_gaussian_list(2));
HvsM_Ratio_replicate_1=cell(Dimension_of_master_gaussian_list(1),Dimension_of_master_gaussian_list(2));
HvsM_Ratio_replicate_2=cell(Dimension_of_master_gaussian_list(1),Dimension_of_master_gaussian_list(2));
HvsM_Ratio_replicate_3=cell(Dimension_of_master_gaussian_list(1),Dimension_of_master_gaussian_list(2));
HvsM_Ratio_replicate_4=cell(Dimension_of_master_gaussian_list(1),Dimension_of_master_gaussian_list(2));

%Determine which Gaussians change by comparing fitted data with raw data
for Gaussian_counter1= 1:Dimension_of_master_gaussian_list(1)
  for internal_counter5 =1:nnz(Finalised_Master_Gaussian_list.Center(Gaussian_counter1,:))
    for replicate_counter = 1:replicate_num
      
      %locate protein to compare
      [location_Protein_in_raw_textdata]=ind2sub(size(txt_HvsL(:,2)), strmatch(Finalised_Master_Gaussian_list.Protein_name{Gaussian_counter1},txt_HvsL(:,2),'exact'));
      Center_to_test=Finalised_Master_Gaussian_list.Center(Gaussian_counter1,internal_counter5); % determine the center of the Gaussian
      rounded_center=floor(Center_to_test); % round center of the Gaussian to nearest integer
      
      %Check if a Gaussian peak was detected within two fractions of this master Gaussian value?
      %Find protein in Combined_Gaussian analysis
      [location_Protein_in_combined_Gaus]=ind2sub(length(Combined_Gaussians.Protein_name),...
        strmatch(Finalised_Master_Gaussian_list.Protein_name{Gaussian_counter1},Combined_Gaussians.Protein_name(:),'exact'));
      %Reset rounded_center if Gaussian found in replicate
      for find_rep_counter=1:length(location_Protein_in_combined_Gaus) %As Combined Gaussian is used to create Final will also have atleast one value
        if replicate_counter==Combined_Gaussians.Replicate(location_Protein_in_combined_Gaus(find_rep_counter))
          %check if Center is within 2 fractions of master list center
          temp_value=abs(Combined_Gaussians.Centrer(location_Protein_in_combined_Gaus(find_rep_counter),:)- Center_to_test);
          if  temp_value<=(User_Window+0.49)
            [idx idx] =min(temp_value);
            rounded_center= floor(Combined_Gaussians.Centrer(location_Protein_in_combined_Gaus(find_rep_counter),idx));
          end
        end
      end
      
      %Determine if Center is greater the used define fraction of interest
      if rounded_center > fraction_to_plot && replicate_counter == 1
        
        %Mark as Unquantifiable replicate 1
        Protein_location_replicate1(Gaussian_counter1,internal_counter5)=location_Protein_in_raw_textdata(replicate_counter);
        Do_gaussian_in_channels_change_replicate_1{Gaussian_counter1,internal_counter5}= 'Unquantifiable_Not_observed_in_replicate';
        Do_gaussian_in_channels_change_for_figure_replicate_1(Gaussian_counter1,internal_counter5)=0;
        log_2_fold_comparsion_of_gaussian_replicate_1(Gaussian_counter1,internal_counter5)= NaN;
        Unquantifiable_Not_observed_in_replicate_replicate_1=Unquantifiable_Not_observed_in_replicate_replicate_1+1;
        Apex_values_replicate_1(Do_gaussian_in_channels_position,:)=NaN(1,4);
        
      elseif rounded_center > fraction_to_plot && replicate_counter == 2
        
        %Mark as Unquantifiable replicate 2
        Protein_location_replicate2(Gaussian_counter1,internal_counter5)=location_Protein_in_raw_textdata(replicate_counter);
        Do_gaussian_in_channels_change_replicate_2{Gaussian_counter1,internal_counter5}= 'Unquantifiable_Not_observed_in_replicate';
        Do_gaussian_in_channels_change_for_figure_replicate_2(Gaussian_counter1,internal_counter5)=0;
        log_2_fold_comparsion_of_gaussian_replicate_2(Gaussian_counter1,internal_counter5)= NaN;
        Unquantifiable_Not_observed_in_replicate_replicate_2=Unquantifiable_Not_observed_in_replicate_replicate_2+1;
        Apex_values_replicate_2(Do_gaussian_in_channels_position,:)=NaN(1,4);
        
      elseif rounded_center > fraction_to_plot && replicate_counter == 3
        
        %Mark as Unquantifiable replicate 3
        Protein_location_replicate3(Gaussian_counter1,internal_counter5)=location_Protein_in_raw_textdata(replicate_counter);
        Do_gaussian_in_channels_change_replicate_3{Gaussian_counter1,internal_counter5}= 'Unquantifiable_Not_observed_in_replicate';
        Do_gaussian_in_channels_change_for_figure_replicate_3(Gaussian_counter1,internal_counter5)=0;
        log_2_fold_comparsion_of_gaussian_replicate_3(Gaussian_counter1,internal_counter5)= NaN;
        Unquantifiable_observed_only_in_MvsL_channel_replicate_3=Unquantifiable_observed_only_in_MvsL_channel_replicate_3+1;
        Apex_values_replicate_3(Do_gaussian_in_channels_position,:)=NaN(1,4);
        
      elseif rounded_center > fraction_to_plot && replicate_counter == 4
        
        %Mark as Unquantifiable replicate 4
        Protein_location_replicate4(Gaussian_counter1,internal_counter5)=location_Protein_in_raw_textdata(replicate_counter);
        Do_gaussian_in_channels_change_replicate_4{Gaussian_counter1,internal_counter5}= 'Unquantifiable_Not_observed_in_replicate';
        Do_gaussian_in_channels_change_for_figure_replicate_4(Gaussian_counter1,internal_counter5)=0;
        log_2_fold_comparsion_of_gaussian_replicate_4(Gaussian_counter1,internal_counter5)= NaN;
        Unquantifiable_observed_only_in_MvsL_channel_replicate_4=Unquantifiable_observed_only_in_MvsL_channel_replicate_4+1;
        Apex_values_replicate_4(Do_gaussian_in_channels_position,:)=NaN(1,4);
        
      elseif rounded_center <= fraction_to_plot
        
        %retrieve ratios of MvsL
        MvsL_Ratio_Apex(1,1)= num_val_MvsL((location_Protein_in_raw_textdata(replicate_counter)-1),rounded_center-2 +6 +1); %Add six as Realigned raw data ranges from -5 to 60
        MvsL_Ratio_Apex(1,2)= num_val_MvsL((location_Protein_in_raw_textdata(replicate_counter)-1),rounded_center-1 +6 +1); %Add six as Realigned raw data ranges from -5 to 60
        MvsL_Ratio_Apex(1,3)= num_val_MvsL((location_Protein_in_raw_textdata(replicate_counter)-1),rounded_center +6 +1); %Add six as Realigned raw data ranges from -5 to 60
        MvsL_Ratio_Apex(1,4)= num_val_MvsL((location_Protein_in_raw_textdata(replicate_counter)-1),rounded_center+1 +6 +1); %Add six as Realigned raw data ranges from -5 to 60
        MvsL_Ratio_Apex(1,5)= num_val_MvsL((location_Protein_in_raw_textdata(replicate_counter)-1),rounded_center+2 +6 +1); %Add six as Realigned raw data ranges from -5 to 60
        
        %retrieve ratios of HvsL
        HvsL_Ratio_Apex(1,1)= num_val_HvsL((location_Protein_in_raw_textdata(replicate_counter)-1),rounded_center-2 +6 +1); %Add six as Realigned raw data ranges from -5 to 60
        HvsL_Ratio_Apex(1,2)= num_val_HvsL((location_Protein_in_raw_textdata(replicate_counter)-1),rounded_center-1 +6 +1); %Add six as Realigned raw data ranges from -5 to 60
        HvsL_Ratio_Apex(1,3)= num_val_HvsL((location_Protein_in_raw_textdata(replicate_counter)-1),rounded_center +6 +1); %Add six as Realigned raw data ranges from -5 to 60
        HvsL_Ratio_Apex(1,4)= num_val_HvsL((location_Protein_in_raw_textdata(replicate_counter)-1),rounded_center+1 +6 +1); %Add six as Realigned raw data ranges from -5 to 60
        HvsL_Ratio_Apex(1,5)= num_val_HvsL((location_Protein_in_raw_textdata(replicate_counter)-1),rounded_center+2 +6 +1); %Add six as Realigned raw data ranges from -5 to 60
        
        %retrieve ratios of HvsM
        HvsM_Ratio_Apex(1,1)= num_val_HvsM((location_Protein_in_raw_textdata(replicate_counter)-1),rounded_center-2 +6 +1); %Add six as Realigned raw data ranges from -5 to 60
        HvsM_Ratio_Apex(1,2)= num_val_HvsM((location_Protein_in_raw_textdata(replicate_counter)-1),rounded_center-1 +6 +1); %Add six as Realigned raw data ranges from -5 to 60
        HvsM_Ratio_Apex(1,3)= num_val_HvsM((location_Protein_in_raw_textdata(replicate_counter)-1),rounded_center +6 +1); %Add six as Realigned raw data ranges from -5 to 60
        HvsM_Ratio_Apex(1,4)= num_val_HvsM((location_Protein_in_raw_textdata(replicate_counter)-1),rounded_center+1 +6 +1); %Add six as Realigned raw data ranges from -5 to 60
        HvsM_Ratio_Apex(1,5)= num_val_HvsM((location_Protein_in_raw_textdata(replicate_counter)-1),rounded_center+2+6 +1); %Add six as Realigned raw data ranges from -5 to 60
        
        %If value for the apex is missing or not values around the apex can be mark as Unquantifiable
        if (sum(isnan(MvsL_Ratio_Apex(1,:))) == 5 | (sum(isnan(MvsL_Ratio_Apex(1,2:4))) >= 2) | (sum(isnan(HvsL_Ratio_Apex(1,:))) == 5) | (sum(isnan(HvsL_Ratio_Apex(1,2:4))) >= 2))
          
          if replicate_counter == 1
            Protein_location_replicate1(Gaussian_counter1,internal_counter5)=location_Protein_in_raw_textdata(replicate_counter);
            %Test if observed only in Heavy channels
            if ~(((sum(isnan(HvsL_Ratio_Apex(1,:))) == 5) | (sum(isnan(HvsL_Ratio_Apex(1,2:4))) >= 2)) & ((sum(isnan(MvsL_Ratio_Apex(1,:))) == 5) | (sum(isnan(MvsL_Ratio_Apex(1,2:4))) >= 2)))
              %Write out what channel was observed
              Do_gaussian_in_channels_change_replicate_1{Gaussian_counter1,internal_counter5}= 'Unquantifiable_observed_only_in HvsL channel';
              Do_gaussian_in_channels_change_for_figure_replicate_1(Gaussian_counter1,internal_counter5)=0;
              log_2_fold_comparsion_of_gaussian_replicate_1(Gaussian_counter1,internal_counter5)= NaN;
              Unquantifiable_observed_only_in_HvsL_channel_replicate_1=Unquantifiable_observed_only_in_HvsL_channel_replicate_1+1;
              Apex_values_replicate_1(Do_gaussian_in_channels_position,:)=NaN(1,4);
              
              %Test if observed only in Light channels
            elseif ~(((sum(isnan(MvsL_Ratio_Apex(1,:))) == 5) | (sum(isnan(MvsL_Ratio_Apex(1,2:4))) >= 2)) & ((sum(isnan(HvsL_Ratio_Apex(1,:))) == 5) | (sum(isnan(HvsL_Ratio_Apex(1,2:4))) >= 2)))
              %Write out what channel was observed
              Do_gaussian_in_channels_change_replicate_1{Gaussian_counter1,internal_counter5}= 'Unquantifiable_observed_only_in_MvsL channel';
              Do_gaussian_in_channels_change_for_figure_replicate_1(Gaussian_counter1,internal_counter5)=0;
              log_2_fold_comparsion_of_gaussian_replicate_1(Gaussian_counter1,internal_counter5)= NaN;
              Unquantifiable_observed_only_in_MvsL_channel_replicate_1=Unquantifiable_observed_only_in_MvsL_channel_replicate_1+1;
              Apex_values_replicate_1(Do_gaussian_in_channels_position,:)=NaN(1,4);
              
              %Test if not observed in both channels
            elseif (((sum(isnan(MvsL_Ratio_Apex(1,:))) == 5 | sum(isnan(MvsL_Ratio_Apex(1,2:4))) >= 2)) & ((sum(isnan(HvsL_Ratio_Apex(1,:))) == 5 | sum(isnan(HvsL_Ratio_Apex(1,2:4))) >= 2)))
              %Write out what channel was observed
              Do_gaussian_in_channels_change_replicate_1{Gaussian_counter1,internal_counter5}= 'Unquantifiable_Not_observed_in_replicate';
              Do_gaussian_in_channels_change_for_figure_replicate_1(Gaussian_counter1,internal_counter5)=0;
              log_2_fold_comparsion_of_gaussian_replicate_1(Gaussian_counter1,internal_counter5)= NaN;
              Unquantifiable_Not_observed_in_replicate_replicate_1=Unquantifiable_Not_observed_in_replicate_replicate_1+1;
              Apex_values_replicate_1(Do_gaussian_in_channels_position,:)=NaN(1,4);
              
            end
            
          elseif replicate_counter == 2
            Protein_location_replicate2(Gaussian_counter1,internal_counter5)=location_Protein_in_raw_textdata(replicate_counter);
            %Test if observed only in light channels
            if ~(((sum(isnan(MvsL_Ratio_Apex(1,:))) == 5) | (sum(isnan(MvsL_Ratio_Apex(1,2:4))) >= 2)) & ((sum(isnan(HvsL_Ratio_Apex(1,:))) == 5) | (sum(isnan(HvsL_Ratio_Apex(1,2:4))) >= 2)))
              %Write out what channel was observed
              Do_gaussian_in_channels_change_replicate_2{Gaussian_counter1,internal_counter5}= 'Unquantifiable_observed_only_in MvsL channel';
              Do_gaussian_in_channels_change_for_figure_replicate_2(Gaussian_counter1,internal_counter5)=0;
              log_2_fold_comparsion_of_gaussian_replicate_2(Gaussian_counter1,internal_counter5)= NaN;
              Unquantifiable_observed_only_in_MvsL_channel_replicate_2=Unquantifiable_observed_only_in_MvsL_channel_replicate_2+1;
              Apex_values_replicate_2(Do_gaussian_in_channels_position,:)=NaN(1,4);
              
              %Test if observed only in heavy channels
            elseif ~(((sum(isnan(HvsL_Ratio_Apex(1,:))) == 5) | (sum(isnan(HvsL_Ratio_Apex(1,2:4))) >= 2)) & ((sum(isnan(MvsL_Ratio_Apex(1,:))) == 5) | (sum(isnan(MvsL_Ratio_Apex(1,2:4))) >= 2)))
              %Write out what channel was observed
              Do_gaussian_in_channels_change_replicate_2{Gaussian_counter1,internal_counter5}= 'Unquantifiable_observed_only_in HvsL channel';
              Do_gaussian_in_channels_change_for_figure_replicate_2(Gaussian_counter1,internal_counter5)=0;
              log_2_fold_comparsion_of_gaussian_replicate_2(Gaussian_counter1,internal_counter5)= NaN;
              Unquantifiable_observed_only_in_HvsL_channel_replicate_2=Unquantifiable_observed_only_in_HvsL_channel_replicate_2+1;
              Apex_values_replicate_2(Do_gaussian_in_channels_position,:)=NaN(1,4);
              
              %Test if not observed in both channels
            elseif (((sum(isnan(MvsL_Ratio_Apex(1,:))) == 5 | sum(isnan(MvsL_Ratio_Apex(1,2:4))) >= 2)) & ((sum(isnan(HvsL_Ratio_Apex(1,:))) == 5 | sum(isnan(HvsL_Ratio_Apex(1,2:4))) >= 2)))
              %Write out what channel was observed
              Do_gaussian_in_channels_change_replicate_2{Gaussian_counter1,internal_counter5}= 'Unquantifiable_Not_observed_in_replicate';
              Do_gaussian_in_channels_change_for_figure_replicate_2(Gaussian_counter1,internal_counter5)=0;
              log_2_fold_comparsion_of_gaussian_replicate_2(Gaussian_counter1,internal_counter5)= NaN;
              Unquantifiable_Not_observed_in_replicate_replicate_2=Unquantifiable_Not_observed_in_replicate_replicate_2+1;
              Apex_values_replicate_2(Do_gaussian_in_channels_position,:)=NaN(1,4);
              
            end
            
          elseif replicate_counter == 3
            Protein_location_replicate3(Gaussian_counter1,internal_counter5)=location_Protein_in_raw_textdata(replicate_counter);
            %Test if observed only in light channels
            if ~(((sum(isnan(MvsL_Ratio_Apex(1,:))) == 5) | (sum(isnan(MvsL_Ratio_Apex(1,2:4))) >= 2)) & ((sum(isnan(HvsL_Ratio_Apex(1,:))) == 5) | (sum(isnan(HvsL_Ratio_Apex(1,2:4))) >= 2)))
              %Write out what channel was observed
              Do_gaussian_in_channels_change_replicate_3{Gaussian_counter1,internal_counter5}= 'Unquantifiable_observed_only_in MvsL channel';
              Do_gaussian_in_channels_change_for_figure_replicate_3(Gaussian_counter1,internal_counter5)=0;
              log_2_fold_comparsion_of_gaussian_replicate_3(Gaussian_counter1,internal_counter5)= NaN;
              Unquantifiable_observed_only_in_MvsL_channel_replicate_3=Unquantifiable_observed_only_in_MvsL_channel_replicate_3 +1;
              Apex_values_replicate_3(Do_gaussian_in_channels_position,:)=NaN(1,4);
              
              %Test if observed only in heavy channels
            elseif ~(((sum(isnan(HvsL_Ratio_Apex(1,:))) == 5) | (sum(isnan(HvsL_Ratio_Apex(1,2:4))) >= 2)) & ((sum(isnan(MvsL_Ratio_Apex(1,:))) == 5) | (sum(isnan(MvsL_Ratio_Apex(1,2:4))) >= 2)))
              %Write out what channel was observed
              Do_gaussian_in_channels_change_replicate_3{Gaussian_counter1,internal_counter5}= 'Unquantifiable_observed_only_in HvsL channel';
              Do_gaussian_in_channels_change_for_figure_replicate_3(Gaussian_counter1,internal_counter5)=0;
              log_2_fold_comparsion_of_gaussian_replicate_3(Gaussian_counter1,internal_counter5)= NaN;
              Unquantifiable_observed_only_in_HvsL_channel_replicate_3=Unquantifiable_observed_only_in_HvsL_channel_replicate_3+1;
              Apex_values_replicate_3(Do_gaussian_in_channels_position,:)=NaN(1,4);
              
              %Test if not observed in both channels
            elseif (((sum(isnan(MvsL_Ratio_Apex(1,:))) == 5 | sum(isnan(MvsL_Ratio_Apex(1,2:4))) >= 2)) & ((sum(isnan(HvsL_Ratio_Apex(1,:))) == 5 | sum(isnan(HvsL_Ratio_Apex(1,2:4))) >= 2)))
              %Write out what channel was observed
              Do_gaussian_in_channels_change_replicate_3{Gaussian_counter1,internal_counter5}= 'Unquantifiable_Not_observed_in_replicate';
              Do_gaussian_in_channels_change_for_figure_replicate_3(Gaussian_counter1,internal_counter5)=0;
              log_2_fold_comparsion_of_gaussian_replicate_3(Gaussian_counter1,internal_counter5)= NaN;
              Unquantifiable_Not_observed_in_replicate_replicate_3=Unquantifiable_Not_observed_in_replicate_replicate_3+1;
              Apex_values_replicate_3(Do_gaussian_in_channels_position,:)=NaN(1,4);
            end
            
          elseif replicate_counter == 4
            Protein_location_replicate4(Gaussian_counter1,internal_counter5)=location_Protein_in_raw_textdata(replicate_counter);
            %Test if observed only in light channels
            if ~(((sum(isnan(MvsL_Ratio_Apex(1,:))) == 5) | (sum(isnan(MvsL_Ratio_Apex(1,2:4))) >= 2)) & ((sum(isnan(HvsL_Ratio_Apex(1,:))) == 5) | (sum(isnan(HvsL_Ratio_Apex(1,2:4))) >= 2)))
              %Write out what channel was observed
              Do_gaussian_in_channels_change_replicate_4{Gaussian_counter1,internal_counter5}= 'Unquantifiable_observed_only_in MvsL channel';
              Do_gaussian_in_channels_change_for_figure_replicate_4(Gaussian_counter1,internal_counter5)=0;
              log_2_fold_comparsion_of_gaussian_replicate_4(Gaussian_counter1,internal_counter5)= NaN;
              Unquantifiable_observed_only_in_MvsL_channel_replicate_4=Unquantifiable_observed_only_in_MvsL_channel_replicate_4 +1;
              Apex_values_replicate_4(Do_gaussian_in_channels_position,:)=NaN(1,4);
              
              %Test if observed only in heavy channels
            elseif ~(((sum(isnan(HvsL_Ratio_Apex(1,:))) == 5) | (sum(isnan(HvsL_Ratio_Apex(1,2:4))) >= 2)) & ((sum(isnan(MvsL_Ratio_Apex(1,:))) == 5) | (sum(isnan(MvsL_Ratio_Apex(1,2:4))) >= 2)))
              %Write out what channel was observed
              Do_gaussian_in_channels_change_replicate_4{Gaussian_counter1,internal_counter5}= 'Unquantifiable_observed_only_in HvsL channel';
              Do_gaussian_in_channels_change_for_figure_replicate_4(Gaussian_counter1,internal_counter5)=0;
              log_2_fold_comparsion_of_gaussian_replicate_4(Gaussian_counter1,internal_counter5)= NaN;
              Unquantifiable_observed_only_in_HvsL_channel_replicate_4=Unquantifiable_observed_only_in_HvsL_channel_replicate_4+1;
              Apex_values_replicate_4(Do_gaussian_in_channels_position,:)=NaN(1,4);
              
              %Test if not observed in both channels
            elseif (((sum(isnan(MvsL_Ratio_Apex(1,:))) == 5 | sum(isnan(MvsL_Ratio_Apex(1,2:4))) >= 2)) & ((sum(isnan(HvsL_Ratio_Apex(1,:))) == 5 | sum(isnan(HvsL_Ratio_Apex(1,2:4))) >= 2)))
              %Write out what channel was observed
              Do_gaussian_in_channels_change_replicate_4{Gaussian_counter1,internal_counter5}= 'Unquantifiable_Not_observed_in_replicate';
              Do_gaussian_in_channels_change_for_figure_replicate_4(Gaussian_counter1,internal_counter5)=0;
              log_2_fold_comparsion_of_gaussian_replicate_4(Gaussian_counter1,internal_counter5)= NaN;
              Unquantifiable_Not_observed_in_replicate_replicate_4=Unquantifiable_Not_observed_in_replicate_replicate_4+1;
              Apex_values_replicate_4(Do_gaussian_in_channels_position,:)=NaN(1,4);
            end
            
          end
          
          %Test if observed in both light channels- enabling quantification
        elseif ~(sum(isnan(MvsL_Ratio_Apex(1,:))) == 5 | sum(isnan(MvsL_Ratio_Apex(1,2:4))) >= 2 | sum(isnan(HvsL_Ratio_Apex(1,:))) == 5 | sum(isnan(HvsL_Ratio_Apex(1,2:4))) >= 2)
          
          %Write values to table to output
          if replicate_counter == 1
            
            %Save values for each channel (use of MWW test)
            MvsL_Ratio_replicate_1(Gaussian_counter1,internal_counter5)={MvsL_Ratio_Apex};
            HvsL_Ratio_replicate_1(Gaussian_counter1,internal_counter5)={HvsL_Ratio_Apex};
            HvsM_Ratio_replicate_1(Gaussian_counter1,internal_counter5)={HvsM_Ratio_Apex};
            
            %save values
            MvsL_Ratio_Apex_prefiltered=MvsL_Ratio_Apex;
            HvsL_Ratio_Apex_prefiltered=HvsL_Ratio_Apex;
            
            %Remove nan values
            HvsM_Ratio_Apex(isnan(HvsM_Ratio_Apex(1,:)))=[];
            HvsL_Ratio_Apex(isnan(HvsL_Ratio_Apex(1,:)))=[];
            MvsL_Ratio_Apex(isnan(MvsL_Ratio_Apex(1,:)))=[];
            
            Protein_location_replicate1(Gaussian_counter1,internal_counter5)=location_Protein_in_raw_textdata(replicate_counter);
            %write out apex value for replicate 1
            Apex_values_replicate_1(Do_gaussian_in_channels_position,1)=mean(MvsL_Ratio_Apex_prefiltered(1,:)./HvsL_Ratio_Apex_prefiltered(1,:));
            Apex_values_replicate_1(Do_gaussian_in_channels_position,2)=mean(MvsL_Ratio_Apex(1,:)); %save the mean of MvsL ratios
            Apex_values_replicate_1(Do_gaussian_in_channels_position,3)=mean(HvsL_Ratio_Apex(1,:)); %save the mean of HvsL ratios
            
            %test if the Apex values are high enough to determine chnages
            if ~(Apex_values_replicate_1(Do_gaussian_in_channels_position,2) < 0.2 | Apex_values_replicate_1(Do_gaussian_in_channels_position,3) <0.2 |...
                isnan(Apex_values_replicate_1(Do_gaussian_in_channels_position,2)) | isnan(Apex_values_replicate_1(Do_gaussian_in_channels_position,3)))
              Apex_values_replicate_1(Do_gaussian_in_channels_position,4)=mean(HvsM_Ratio_Apex(1,:)); %save the mean of HvsM ratios
            elseif (Apex_values_replicate_1(Do_gaussian_in_channels_position,2) < 0.2 | Apex_values_replicate_1(Do_gaussian_in_channels_position,3) <0.2 |....
                isnan(Apex_values_replicate_1(Do_gaussian_in_channels_position,2)) | isnan(Apex_values_replicate_1(Do_gaussian_in_channels_position,3)))
              Apex_values_replicate_1(Do_gaussian_in_channels_position,4)=NaN(1);
            end
            
            %Test if apex changes
            if (Apex_values_replicate_1(Do_gaussian_in_channels_position,4) <= 2 && Apex_values_replicate_1(Do_gaussian_in_channels_position,4) >= 0.5)
              Increase_at_apex_in_stimulation_replicate_1(Gaussian_counter1,internal_counter5)= 0;
              Decrease_at_apex_in_stimulation_replicate_1(Gaussian_counter1,internal_counter5)= 0;
              Do_gaussian_in_channels_change_replicate_1{Gaussian_counter1,internal_counter5}= 'No change';
              Do_gaussian_in_channels_change_for_figure_replicate_1(Gaussian_counter1,internal_counter5)=5;
              log_2_fold_comparsion_of_gaussian_replicate_1(Gaussian_counter1,internal_counter5)= log2(Apex_values_replicate_1(Do_gaussian_in_channels_position,4));
              No_change_counter_replicate_1=No_change_counter_replicate_1+1;
            elseif (Apex_values_replicate_1(Do_gaussian_in_channels_position,4) >= 2)
              Increase_at_apex_in_stimulation_replicate_1(Gaussian_counter1,internal_counter5)= 1;
              Decrease_at_apex_in_stimulation_replicate_1(Gaussian_counter1,internal_counter5)= 0;
              Do_gaussian_in_channels_change_replicate_1{Gaussian_counter1,internal_counter5}='Increase';
              Do_gaussian_in_channels_change_for_figure_replicate_1(Gaussian_counter1,internal_counter5)=10;
              log_2_fold_comparsion_of_gaussian_replicate_1(Gaussian_counter1,internal_counter5)= log2(Apex_values_replicate_1(Do_gaussian_in_channels_position,4));
              Increase_counter_replicate_1=Increase_counter_replicate_1+1;
            elseif (Apex_values_replicate_1(Do_gaussian_in_channels_position,4) <= 0.5)
              Decrease_at_apex_in_stimulation_replicate_1(Gaussian_counter1,internal_counter5)= 1;
              Increase_at_apex_in_stimulation_replicate_1(Gaussian_counter1,internal_counter5)= 0;
              Do_gaussian_in_channels_change_replicate_1{Gaussian_counter1,internal_counter5}='Decrease';
              Do_gaussian_in_channels_change_for_figure_replicate_1(Gaussian_counter1,internal_counter5)=1;
              log_2_fold_comparsion_of_gaussian_replicate_1(Gaussian_counter1,internal_counter5)= log2(Apex_values_replicate_1(Do_gaussian_in_channels_position,4));
              Decrease_counter_replicate_1=Decrease_counter_replicate_1+1;
            elseif isnan(Apex_values_replicate_1(Do_gaussian_in_channels_position,4))
              Increase_at_apex_in_stimulation_replicate_1(Gaussian_counter1,internal_counter5)= 0;
              Decrease_at_apex_in_stimulation_replicate_1(Gaussian_counter1,internal_counter5)= 0;
              Do_gaussian_in_channels_change_replicate_1{Gaussian_counter1,internal_counter5}= 'Unquantifiable_filtered';
              Do_gaussian_in_channels_change_for_figure_replicate_1(Gaussian_counter1,internal_counter5)=5;
              log_2_fold_comparsion_of_gaussian_replicate_1(Gaussian_counter1,internal_counter5)= log2(Apex_values_replicate_1(Do_gaussian_in_channels_position,4));
              No_change_counter_replicate_1=No_change_counter_replicate_1+1;
            end
            
          elseif replicate_counter == 2
            
            %Save values for each channel (use of MWW test)
            MvsL_Ratio_replicate_2(Gaussian_counter1,internal_counter5)={MvsL_Ratio_Apex};
            HvsL_Ratio_replicate_2(Gaussian_counter1,internal_counter5)={HvsL_Ratio_Apex};
            HvsM_Ratio_replicate_2(Gaussian_counter1,internal_counter5)={HvsM_Ratio_Apex};
            
            %save values
            MvsL_Ratio_Apex_prefiltered=MvsL_Ratio_Apex;
            HvsL_Ratio_Apex_prefiltered=HvsL_Ratio_Apex;
            
            %Remove nan values
            HvsM_Ratio_Apex(isnan(HvsM_Ratio_Apex(1,:)))=[];
            HvsL_Ratio_Apex(isnan(HvsL_Ratio_Apex(1,:)))=[];
            MvsL_Ratio_Apex(isnan(MvsL_Ratio_Apex(1,:)))=[];
            
            Protein_location_replicate2(Gaussian_counter1,internal_counter5)=location_Protein_in_raw_textdata(replicate_counter);
            %write out apex value for replicate 2
            Apex_values_replicate_2(Do_gaussian_in_channels_position,1)=mean(MvsL_Ratio_Apex_prefiltered(1,:)./HvsL_Ratio_Apex_prefiltered(1,:));
            Apex_values_replicate_2(Do_gaussian_in_channels_position,2)=mean(MvsL_Ratio_Apex(1,:)); %save the mean of MvsL ratios
            Apex_values_replicate_2(Do_gaussian_in_channels_position,3)=mean(HvsL_Ratio_Apex(1,:)); %save the mean of HvsL ratios
            
            %test if the Apex values are high enough to determine chnages
            if ~(Apex_values_replicate_2(Do_gaussian_in_channels_position,2) < 0.2 | Apex_values_replicate_2(Do_gaussian_in_channels_position,3) <0.2 |...
                isnan(Apex_values_replicate_2(Do_gaussian_in_channels_position,2)) | isnan(Apex_values_replicate_2(Do_gaussian_in_channels_position,3)))
              Apex_values_replicate_2(Do_gaussian_in_channels_position,4)=mean(HvsM_Ratio_Apex(1,:)); %save the mean of HvsM ratios
            elseif (Apex_values_replicate_2(Do_gaussian_in_channels_position,2) < 0.2 | Apex_values_replicate_2(Do_gaussian_in_channels_position,3) <0.2 |....
                isnan(Apex_values_replicate_2(Do_gaussian_in_channels_position,2)) | isnan(Apex_values_replicate_2(Do_gaussian_in_channels_position,3)))
              Apex_values_replicate_2(Do_gaussian_in_channels_position,4)=NaN(1);
            end
            
            %Test if apex changes
            if (Apex_values_replicate_2(Do_gaussian_in_channels_position,4) <= 2 && Apex_values_replicate_2(Do_gaussian_in_channels_position,4) >= 0.5)
              Increase_at_apex_in_stimulation_replicate_2(Gaussian_counter1,internal_counter5)= 0;
              Decrease_at_apex_in_stimulation_replicate_2(Gaussian_counter1,internal_counter5)= 0;
              Do_gaussian_in_channels_change_replicate_2{Gaussian_counter1,internal_counter5}= 'No change';
              Do_gaussian_in_channels_change_for_figure_replicate_2(Gaussian_counter1,internal_counter5)=5;
              log_2_fold_comparsion_of_gaussian_replicate_2(Gaussian_counter1,internal_counter5)= log2(Apex_values_replicate_2(Do_gaussian_in_channels_position,4));
              No_change_counter_replicate_2=No_change_counter_replicate_2+1;
            elseif (Apex_values_replicate_2(Do_gaussian_in_channels_position,4) >= 2)
              Increase_at_apex_in_stimulation_replicate_2(Gaussian_counter1,internal_counter5)= 1;
              Decrease_at_apex_in_stimulation_replicate_2(Gaussian_counter1,internal_counter5)= 0;
              Do_gaussian_in_channels_change_replicate_2{Gaussian_counter1,internal_counter5}='Increase';
              Do_gaussian_in_channels_change_for_figure_replicate_2(Gaussian_counter1,internal_counter5)=10;
              log_2_fold_comparsion_of_gaussian_replicate_2(Gaussian_counter1,internal_counter5)= log2(Apex_values_replicate_2(Do_gaussian_in_channels_position,4));
              Increase_counter_replicate_2=Increase_counter_replicate_2+1;
            elseif (Apex_values_replicate_2(Do_gaussian_in_channels_position,4) <= 0.5)
              Decrease_at_apex_in_stimulation_replicate_2(Gaussian_counter1,internal_counter5)= 1;
              Increase_at_apex_in_stimulation_replicate_2(Gaussian_counter1,internal_counter5)= 0;
              Do_gaussian_in_channels_change_replicate_2{Gaussian_counter1,internal_counter5}='Decrease';
              Do_gaussian_in_channels_change_for_figure_replicate_2(Gaussian_counter1,internal_counter5)=1;
              log_2_fold_comparsion_of_gaussian_replicate_2(Gaussian_counter1,internal_counter5)= log2(Apex_values_replicate_2(Do_gaussian_in_channels_position,4));
              Decrease_counter_replicate_2=Decrease_counter_replicate_2+1;
            elseif isnan(Apex_values_replicate_2(Do_gaussian_in_channels_position,4))
              Increase_at_apex_in_stimulation_replicate_2(Gaussian_counter1,internal_counter5)= 0;
              Decrease_at_apex_in_stimulation_replicate_2(Gaussian_counter1,internal_counter5)= 0;
              Do_gaussian_in_channels_change_replicate_2{Gaussian_counter1,internal_counter5}= 'Unquantifiable_filtered';
              Do_gaussian_in_channels_change_for_figure_replicate_2(Gaussian_counter1,internal_counter5)=5;
              log_2_fold_comparsion_of_gaussian_replicate_2(Gaussian_counter1,internal_counter5)= log2(Apex_values_replicate_2(Do_gaussian_in_channels_position,4));
              No_change_counter_replicate_2=No_change_counter_replicate_2+1;
            end
            
          elseif replicate_counter == 3
            
            %Save values for each channel (use of MWW test)
            MvsL_Ratio_replicate_3(Gaussian_counter1,internal_counter5)={MvsL_Ratio_Apex};
            HvsL_Ratio_replicate_3(Gaussian_counter1,internal_counter5)={HvsL_Ratio_Apex};
            HvsM_Ratio_replicate_3(Gaussian_counter1,internal_counter5)={HvsM_Ratio_Apex};
            
            %save values
            MvsL_Ratio_Apex_prefiltered=MvsL_Ratio_Apex;
            HvsL_Ratio_Apex_prefiltered=HvsL_Ratio_Apex;
            
            %Remove nan values
            HvsM_Ratio_Apex(isnan(HvsM_Ratio_Apex(1,:)))=[];
            HvsL_Ratio_Apex(isnan(HvsL_Ratio_Apex(1,:)))=[];
            MvsL_Ratio_Apex(isnan(MvsL_Ratio_Apex(1,:)))=[];
            
            Protein_location_replicate3(Gaussian_counter1,internal_counter5)=location_Protein_in_raw_textdata(replicate_counter);
            %write out apex value for replicate 3
            Apex_values_replicate_3(Do_gaussian_in_channels_position,1)=mean(MvsL_Ratio_Apex_prefiltered(1,:)./HvsL_Ratio_Apex_prefiltered(1,:));
            Apex_values_replicate_3(Do_gaussian_in_channels_position,2)=mean(MvsL_Ratio_Apex(1,:)); %save the mean of MvsL ratios
            Apex_values_replicate_3(Do_gaussian_in_channels_position,3)=mean(HvsL_Ratio_Apex(1,:)); %save the mean of HvsL ratios
            
            %test if the Apex values are high enough to determine chnages
            if ~(Apex_values_replicate_3(Do_gaussian_in_channels_position,2) < 0.2 || Apex_values_replicate_3(Do_gaussian_in_channels_position,3) <0.2 |...
                isnan(Apex_values_replicate_3(Do_gaussian_in_channels_position,2)) || isnan(Apex_values_replicate_3(Do_gaussian_in_channels_position,3)))
              Apex_values_replicate_3(Do_gaussian_in_channels_position,4)=mean(HvsM_Ratio_Apex(1,:)); %save the mean of HvsM ratios
            elseif (Apex_values_replicate_3(Do_gaussian_in_channels_position,2) < 0.2 || Apex_values_replicate_3(Do_gaussian_in_channels_position,3) <0.2 |....
                isnan(Apex_values_replicate_3(Do_gaussian_in_channels_position,2)) || isnan(Apex_values_replicate_3(Do_gaussian_in_channels_position,3)))
              Apex_values_replicate_3(Do_gaussian_in_channels_position,4)=NaN(1);
            end
            
            %Test if apex changes
            if (Apex_values_replicate_3(Do_gaussian_in_channels_position,4) <= 2 && Apex_values_replicate_3(Do_gaussian_in_channels_position,4) >= 0.5)
              Increase_at_apex_in_stimulation_replicate_3(Gaussian_counter1,internal_counter5)= 0;
              Decrease_at_apex_in_stimulation_replicate_3(Gaussian_counter1,internal_counter5)= 0;
              Do_gaussian_in_channels_change_replicate_3{Gaussian_counter1,internal_counter5}= 'No change';
              Do_gaussian_in_channels_change_for_figure_replicate_3(Gaussian_counter1,internal_counter5)=5;
              log_2_fold_comparsion_of_gaussian_replicate_3(Gaussian_counter1,internal_counter5)= log2(Apex_values_replicate_3(Do_gaussian_in_channels_position,4));
              No_change_counter_replicate_3=No_change_counter_replicate_3+1;
            elseif (Apex_values_replicate_3(Do_gaussian_in_channels_position,4) >= 2)
              Increase_at_apex_in_stimulation_replicate_3(Gaussian_counter1,internal_counter5)= 1;
              Decrease_at_apex_in_stimulation_replicate_3(Gaussian_counter1,internal_counter5)= 0;
              Do_gaussian_in_channels_change_replicate_3{Gaussian_counter1,internal_counter5}='Increase';
              Do_gaussian_in_channels_change_for_figure_replicate_3(Gaussian_counter1,internal_counter5)=10;
              log_2_fold_comparsion_of_gaussian_replicate_3(Gaussian_counter1,internal_counter5)= log2(Apex_values_replicate_3(Do_gaussian_in_channels_position,4));
              Increase_counter_replicate_3=Increase_counter_replicate_3+1;
            elseif (Apex_values_replicate_3(Do_gaussian_in_channels_position,4) <= 0.5)
              Decrease_at_apex_in_stimulation_replicate_3(Gaussian_counter1,internal_counter5)= 1;
              Increase_at_apex_in_stimulation_replicate_3(Gaussian_counter1,internal_counter5)= 0;
              Do_gaussian_in_channels_change_replicate_3{Gaussian_counter1,internal_counter5}='Decrease';
              Do_gaussian_in_channels_change_for_figure_replicate_3(Gaussian_counter1,internal_counter5)=1;
              log_2_fold_comparsion_of_gaussian_replicate_3(Gaussian_counter1,internal_counter5)= log2(Apex_values_replicate_3(Do_gaussian_in_channels_position,4));
              Decrease_counter_replicate_3=Decrease_counter_replicate_3+1;
            elseif isnan(Apex_values_replicate_3(Do_gaussian_in_channels_position,4))
              Increase_at_apex_in_stimulation_replicate_3(Gaussian_counter1,internal_counter5)= 0;
              Decrease_at_apex_in_stimulation_replicate_3(Gaussian_counter1,internal_counter5)= 0;
              Do_gaussian_in_channels_change_replicate_3{Gaussian_counter1,internal_counter5}= 'Unquantifiable_filtered';
              Do_gaussian_in_channels_change_for_figure_replicate_3(Gaussian_counter1,internal_counter5)=5;
              log_2_fold_comparsion_of_gaussian_replicate_3(Gaussian_counter1,internal_counter5)= log2(Apex_values_replicate_3(Do_gaussian_in_channels_position,4));
              No_change_counter_replicate_3=No_change_counter_replicate_3+1;
            end
            
          elseif replicate_counter == 4
            
            %Save values for each channel (use of MWW test)
            MvsL_Ratio_replicate_4(Gaussian_counter1,internal_counter5)={MvsL_Ratio_Apex};
            HvsL_Ratio_replicate_4(Gaussian_counter1,internal_counter5)={HvsL_Ratio_Apex};
            HvsM_Ratio_replicate_4(Gaussian_counter1,internal_counter5)={HvsM_Ratio_Apex};
            
            %save values
            MvsL_Ratio_Apex_prefiltered=MvsL_Ratio_Apex;
            HvsL_Ratio_Apex_prefiltered=HvsL_Ratio_Apex;
            
            %Remove nan values
            HvsM_Ratio_Apex(isnan(HvsM_Ratio_Apex(1,:)))=[];
            HvsL_Ratio_Apex(isnan(HvsL_Ratio_Apex(1,:)))=[];
            MvsL_Ratio_Apex(isnan(MvsL_Ratio_Apex(1,:)))=[];
            
            Protein_location_replicate4(Gaussian_counter1,internal_counter5)=location_Protein_in_raw_textdata(replicate_counter);
            %write out apex value for replicate 1
            Apex_values_replicate_4(Do_gaussian_in_channels_position,1)=mean(MvsL_Ratio_Apex_prefiltered(1,:)./HvsL_Ratio_Apex_prefiltered(1,:));
            Apex_values_replicate_4(Do_gaussian_in_channels_position,2)=mean(MvsL_Ratio_Apex(1,:)); %save the mean of MvsL ratios
            Apex_values_replicate_4(Do_gaussian_in_channels_position,3)=mean(HvsL_Ratio_Apex(1,:)); %save the mean of HvsL ratios
            
            %test if the Apex values are high enough to determine chnages
            if ~(Apex_values_replicate_4(Do_gaussian_in_channels_position,2) < 0.2 | Apex_values_replicate_4(Do_gaussian_in_channels_position,3) <0.2 |...
                isnan(Apex_values_replicate_4(Do_gaussian_in_channels_position,2)) | isnan(Apex_values_replicate_4(Do_gaussian_in_channels_position,3)))
              Apex_values_replicate_4(Do_gaussian_in_channels_position,4)=mean(HvsM_Ratio_Apex(1,:)); %save the mean of HvsM ratios
            elseif (Apex_values_replicate_4(Do_gaussian_in_channels_position,2) < 0.2 | Apex_values_replicate_4(Do_gaussian_in_channels_position,3) <0.2 |....
                isnan(Apex_values_replicate_4(Do_gaussian_in_channels_position,2)) | isnan(Apex_values_replicate_4(Do_gaussian_in_channels_position,3)))
              Apex_values_replicate_4(Do_gaussian_in_channels_position,4)=NaN(1);
            end
            
            %Test if apex changes
            if (Apex_values_replicate_4(Do_gaussian_in_channels_position,4) <= 2 && Apex_values_replicate_4(Do_gaussian_in_channels_position,4) >= 0.5)
              Increase_at_apex_in_stimulation_replicate_4(Gaussian_counter1,internal_counter5)= 0;
              Decrease_at_apex_in_stimulation_replicate_4(Gaussian_counter1,internal_counter5)= 0;
              Do_gaussian_in_channels_change_replicate_4{Gaussian_counter1,internal_counter5}= 'No change';
              Do_gaussian_in_channels_change_for_figure_replicate_4(Gaussian_counter1,internal_counter5)=5;
              log_2_fold_comparsion_of_gaussian_replicate_4(Gaussian_counter1,internal_counter5)= log2(Apex_values_replicate_4(Do_gaussian_in_channels_position,4));
              No_change_counter_replicate_4=No_change_counter_replicate_4+1;
            elseif (Apex_values_replicate_4(Do_gaussian_in_channels_position,4) >= 2)
              Increase_at_apex_in_stimulation_replicate_4(Gaussian_counter1,internal_counter5)= 1;
              Decrease_at_apex_in_stimulation_replicate_4(Gaussian_counter1,internal_counter5)= 0;
              Do_gaussian_in_channels_change_replicate_4{Gaussian_counter1,internal_counter5}='Increase';
              Do_gaussian_in_channels_change_for_figure_replicate_4(Gaussian_counter1,internal_counter5)=10;
              log_2_fold_comparsion_of_gaussian_replicate_4(Gaussian_counter1,internal_counter5)= log2(Apex_values_replicate_4(Do_gaussian_in_channels_position,4));
              Increase_counter_replicate_4=Increase_counter_replicate_4+1;
            elseif (Apex_values_replicate_4(Do_gaussian_in_channels_position,4) <= 0.5)
              Decrease_at_apex_in_stimulation_replicate_4(Gaussian_counter1,internal_counter5)= 1;
              Increase_at_apex_in_stimulation_replicate_4(Gaussian_counter1,internal_counter5)= 0;
              Do_gaussian_in_channels_change_replicate_4{Gaussian_counter1,internal_counter5}='Decrease';
              Do_gaussian_in_channels_change_for_figure_replicate_4(Gaussian_counter1,internal_counter5)=1;
              log_2_fold_comparsion_of_gaussian_replicate_4(Gaussian_counter1,internal_counter5)= log2(Apex_values_replicate_4(Do_gaussian_in_channels_position,4));
              Decrease_counter_replicate_4=Decrease_counter_replicate_4+1;
            elseif isnan(Apex_values_replicate_4(Do_gaussian_in_channels_position,4))
              Increase_at_apex_in_stimulation_replicate_4(Gaussian_counter1,internal_counter5)= 0;
              Decrease_at_apex_in_stimulation_replicate_4(Gaussian_counter1,internal_counter5)= 0;
              Do_gaussian_in_channels_change_replicate_4{Gaussian_counter1,internal_counter5}= 'Unquantifiable_filtered';
              Do_gaussian_in_channels_change_for_figure_replicate_4(Gaussian_counter1,internal_counter5)=5;
              log_2_fold_comparsion_of_gaussian_replicate_4(Gaussian_counter1,internal_counter5)= log2(Apex_values_replicate_4(Do_gaussian_in_channels_position,4));
              No_change_counter_replicate_4=No_change_counter_replicate_4+1;
            end
            
          end
        end
        
      end
    end
    %Add one to the row counter
    Do_gaussian_in_channels_position=Do_gaussian_in_channels_position+1;
  end
end

%remove NaN entries to work out mean
manipulation1_log_2_fold_comparsion_of_gaussian_replicate_1=log_2_fold_comparsion_of_gaussian_replicate_1;
manipulation1_log_2_fold_comparsion_of_gaussian_replicate_1(manipulation1_log_2_fold_comparsion_of_gaussian_replicate_1==0)=[];
manipulation2_log_2_fold_comparsion_of_gaussian_replicate_1=manipulation1_log_2_fold_comparsion_of_gaussian_replicate_1;
manipulation2_log_2_fold_comparsion_of_gaussian_replicate_1(isnan(manipulation2_log_2_fold_comparsion_of_gaussian_replicate_1))=[];

manipulation1_log_2_fold_comparsion_of_gaussian_replicate_2=log_2_fold_comparsion_of_gaussian_replicate_2;
manipulation1_log_2_fold_comparsion_of_gaussian_replicate_2(manipulation1_log_2_fold_comparsion_of_gaussian_replicate_2==0)=[];
manipulation2_log_2_fold_comparsion_of_gaussian_replicate_2=manipulation1_log_2_fold_comparsion_of_gaussian_replicate_2;
manipulation2_log_2_fold_comparsion_of_gaussian_replicate_2(isnan(manipulation2_log_2_fold_comparsion_of_gaussian_replicate_2))=[];

manipulation1_log_2_fold_comparsion_of_gaussian_replicate_3=log_2_fold_comparsion_of_gaussian_replicate_3;
manipulation1_log_2_fold_comparsion_of_gaussian_replicate_3(manipulation1_log_2_fold_comparsion_of_gaussian_replicate_3==0)=[];
manipulation2_log_2_fold_comparsion_of_gaussian_replicate_3=manipulation1_log_2_fold_comparsion_of_gaussian_replicate_3;
manipulation2_log_2_fold_comparsion_of_gaussian_replicate_3(isnan(manipulation2_log_2_fold_comparsion_of_gaussian_replicate_3))=[];

manipulation1_log_2_fold_comparsion_of_gaussian_replicate_4=log_2_fold_comparsion_of_gaussian_replicate_4;
manipulation1_log_2_fold_comparsion_of_gaussian_replicate_4(manipulation1_log_2_fold_comparsion_of_gaussian_replicate_4==0)=[];
manipulation2_log_2_fold_comparsion_of_gaussian_replicate_4=manipulation1_log_2_fold_comparsion_of_gaussian_replicate_4;
manipulation2_log_2_fold_comparsion_of_gaussian_replicate_4(isnan(manipulation2_log_2_fold_comparsion_of_gaussian_replicate_4))=[];

%Normalise log2 comparsion
mean_log_2_fold_comparsion_of_gaussian_replicate_1=mean(manipulation2_log_2_fold_comparsion_of_gaussian_replicate_1);
mean_log_2_fold_comparsion_of_gaussian_replicate_2=mean(manipulation2_log_2_fold_comparsion_of_gaussian_replicate_2);
mean_log_2_fold_comparsion_of_gaussian_replicate_3=mean(manipulation2_log_2_fold_comparsion_of_gaussian_replicate_3);
mean_log_2_fold_comparsion_of_gaussian_replicate_4=mean(manipulation2_log_2_fold_comparsion_of_gaussian_replicate_4);

%create varibles
normaliased_log2_comparsion_replicate_1=NaN(Dimension_of_master_gaussian_list(1),Dimension_of_master_gaussian_list(2));
normaliased_log2_comparsion_replicate_2=NaN(Dimension_of_master_gaussian_list(1),Dimension_of_master_gaussian_list(2));
normaliased_log2_comparsion_replicate_3=NaN(Dimension_of_master_gaussian_list(1),Dimension_of_master_gaussian_list(2));
normaliased_log2_comparsion_replicate_4=NaN(Dimension_of_master_gaussian_list(1),Dimension_of_master_gaussian_list(2));
log2_comparsion_replicate_1=NaN(Dimension_of_master_gaussian_list(1),Dimension_of_master_gaussian_list(2));
log2_comparsion_replicate_2=NaN(Dimension_of_master_gaussian_list(1),Dimension_of_master_gaussian_list(2));
log2_comparsion_replicate_3=NaN(Dimension_of_master_gaussian_list(1),Dimension_of_master_gaussian_list(2));
log2_comparsion_replicate_4=NaN(Dimension_of_master_gaussian_list(1),Dimension_of_master_gaussian_list(2));

%Copy Log2 values into same format as Height, Center and width
%replicate 1
for writeout_counter1= 1:Dimension_of_master_gaussian_list(1)
  for writeout_counter2= 1:nnz(log_2_fold_comparsion_of_gaussian_replicate_1(writeout_counter1,:))
    normaliased_log2_comparsion_replicate_1(writeout_counter1,writeout_counter2)=log_2_fold_comparsion_of_gaussian_replicate_1(writeout_counter1,writeout_counter2)-mean_log_2_fold_comparsion_of_gaussian_replicate_1;
    log2_comparsion_replicate_1(writeout_counter1,writeout_counter2)=log_2_fold_comparsion_of_gaussian_replicate_1(writeout_counter1,writeout_counter2);
  end
end

%replicate 2
for writeout_counter1= 1:Dimension_of_master_gaussian_list(1)
  for writeout_counter2= 1:nnz(log_2_fold_comparsion_of_gaussian_replicate_2(writeout_counter1,:))
    normaliased_log2_comparsion_replicate_2(writeout_counter1,writeout_counter2)=log_2_fold_comparsion_of_gaussian_replicate_2(writeout_counter1,writeout_counter2)-mean_log_2_fold_comparsion_of_gaussian_replicate_2;
    log2_comparsion_replicate_2(writeout_counter1,writeout_counter2)=log_2_fold_comparsion_of_gaussian_replicate_2(writeout_counter1,writeout_counter2);
  end
end

%replicate 3
for writeout_counter1= 1:Dimension_of_master_gaussian_list(1)
  for writeout_counter2= 1:nnz(log_2_fold_comparsion_of_gaussian_replicate_3(writeout_counter1,:))
    normaliased_log2_comparsion_replicate_3(writeout_counter1,writeout_counter2)=log_2_fold_comparsion_of_gaussian_replicate_3(writeout_counter1,writeout_counter2)-mean_log_2_fold_comparsion_of_gaussian_replicate_3;
    log2_comparsion_replicate_3(writeout_counter1,writeout_counter2)=log_2_fold_comparsion_of_gaussian_replicate_3(writeout_counter1,writeout_counter2);
  end
end

%replicate 4
for writeout_counter1= 1:Dimension_of_master_gaussian_list(1)
  for writeout_counter2= 1:nnz(log_2_fold_comparsion_of_gaussian_replicate_4(writeout_counter1,:))
    normaliased_log2_comparsion_replicate_4(writeout_counter1,writeout_counter2)=log_2_fold_comparsion_of_gaussian_replicate_4(writeout_counter1,writeout_counter2)-mean_log_2_fold_comparsion_of_gaussian_replicate_4;
    log2_comparsion_replicate_4(writeout_counter1,writeout_counter2)=log_2_fold_comparsion_of_gaussian_replicate_4(writeout_counter1,writeout_counter2);
  end
end


%Determine consistency of changes across replicates
%Determine if there is no change

No_change_across_replicates=zeros(Dimension_of_master_gaussian_list(1),Dimension_of_master_gaussian_list(2));

for writeout_counter1= 1:Dimension_of_master_gaussian_list(1)
  for writeout_counter2= 1:nnz(Finalised_Master_Gaussian_list.Center(writeout_counter1,:))
    No_changes_matrix_across_replicate(4)= (normaliased_log2_comparsion_replicate_4(writeout_counter1,writeout_counter2) >-1 & normaliased_log2_comparsion_replicate_4(writeout_counter1,writeout_counter2) <1);
    No_changes_matrix_across_replicate(3)= (normaliased_log2_comparsion_replicate_3(writeout_counter1,writeout_counter2) >-1 & normaliased_log2_comparsion_replicate_3(writeout_counter1,writeout_counter2) <1);
    No_changes_matrix_across_replicate(2)= (normaliased_log2_comparsion_replicate_2(writeout_counter1,writeout_counter2) >-1 & normaliased_log2_comparsion_replicate_2(writeout_counter1,writeout_counter2) <1);
    No_changes_matrix_across_replicate(1)= (normaliased_log2_comparsion_replicate_1(writeout_counter1,writeout_counter2) >-1 & normaliased_log2_comparsion_replicate_1(writeout_counter1,writeout_counter2) <1);
    if sum(No_changes_matrix_across_replicate) >= 2
      No_change_across_replicates(writeout_counter1,writeout_counter2) =1;
    end
  end
end

Increase_change_across_replicates=zeros(Dimension_of_master_gaussian_list(1),Dimension_of_master_gaussian_list(2));

%Determine if there is an increase change
for writeout_counter1= 1:Dimension_of_master_gaussian_list(1)
  for writeout_counter2= 1:nnz(Finalised_Master_Gaussian_list.Center(writeout_counter1,:))
    Increase_changes_matrix_across_replicate(4)= (normaliased_log2_comparsion_replicate_4(writeout_counter1,writeout_counter2) >=1);
    Increase_changes_matrix_across_replicate(3)= (normaliased_log2_comparsion_replicate_3(writeout_counter1,writeout_counter2) >=1);
    Increase_changes_matrix_across_replicate(2)= (normaliased_log2_comparsion_replicate_2(writeout_counter1,writeout_counter2) >=1);
    Increase_changes_matrix_across_replicate(1)= (normaliased_log2_comparsion_replicate_1(writeout_counter1,writeout_counter2) >=1);
    if sum(Increase_changes_matrix_across_replicate) >= 2
      Increase_change_across_replicates(writeout_counter1,writeout_counter2) =1;
    end
  end
end

Decrease_change_across_replicates=zeros(Dimension_of_master_gaussian_list(1),Dimension_of_master_gaussian_list(2));

%Determine if there is an decrease change
for writeout_counter1= 1:Dimension_of_master_gaussian_list(1)
  for writeout_counter2= 1:nnz(Finalised_Master_Gaussian_list.Center(writeout_counter1,:))
    Decrease_changes_matrix_across_replicate(4)= (normaliased_log2_comparsion_replicate_4(writeout_counter1,writeout_counter2) <=-1);
    Decrease_changes_matrix_across_replicate(3)= (normaliased_log2_comparsion_replicate_3(writeout_counter1,writeout_counter2) <=-1);
    Decrease_changes_matrix_across_replicate(2)= (normaliased_log2_comparsion_replicate_2(writeout_counter1,writeout_counter2) <=-1);
    Decrease_changes_matrix_across_replicate(1)= (normaliased_log2_comparsion_replicate_1(writeout_counter1,writeout_counter2) <=-1);
    if sum(Decrease_changes_matrix_across_replicate) >= 2
      Decrease_change_across_replicates(writeout_counter1,writeout_counter2) =1;
    end
  end
end

Inconsistent_across_replicates=zeros(Dimension_of_master_gaussian_list(1),Dimension_of_master_gaussian_list(2));

%Determine if not detected
for writeout_counter1= 1:Dimension_of_master_gaussian_list(1)
  for writeout_counter2= 1:nnz(Finalised_Master_Gaussian_list.Center(writeout_counter1,:))
    Decrease_changes_matrix_across_replicate(4)= (normaliased_log2_comparsion_replicate_4(writeout_counter1,writeout_counter2) <=-1);
    Decrease_changes_matrix_across_replicate(3)= (normaliased_log2_comparsion_replicate_3(writeout_counter1,writeout_counter2) <=-1);
    Decrease_changes_matrix_across_replicate(2)= (normaliased_log2_comparsion_replicate_2(writeout_counter1,writeout_counter2) <=-1);
    Decrease_changes_matrix_across_replicate(1)= (normaliased_log2_comparsion_replicate_1(writeout_counter1,writeout_counter2) <=-1);
    Increase_changes_matrix_across_replicate(4)= (normaliased_log2_comparsion_replicate_4(writeout_counter1,writeout_counter2) >=1);
    Increase_changes_matrix_across_replicate(3)= (normaliased_log2_comparsion_replicate_3(writeout_counter1,writeout_counter2) >=1);
    Increase_changes_matrix_across_replicate(2)= (normaliased_log2_comparsion_replicate_2(writeout_counter1,writeout_counter2) >=1);
    Increase_changes_matrix_across_replicate(1)= (normaliased_log2_comparsion_replicate_1(writeout_counter1,writeout_counter2) >=1);
    No_changes_matrix_across_replicate(4)= (normaliased_log2_comparsion_replicate_4(writeout_counter1,writeout_counter2) >-1 & normaliased_log2_comparsion_replicate_4(writeout_counter1,writeout_counter2) <1);
    No_changes_matrix_across_replicate(3)= (normaliased_log2_comparsion_replicate_3(writeout_counter1,writeout_counter2) >-1 & normaliased_log2_comparsion_replicate_3(writeout_counter1,writeout_counter2) <1);
    No_changes_matrix_across_replicate(2)= (normaliased_log2_comparsion_replicate_2(writeout_counter1,writeout_counter2) >-1 & normaliased_log2_comparsion_replicate_2(writeout_counter1,writeout_counter2) <1);
    No_changes_matrix_across_replicate(1)= (normaliased_log2_comparsion_replicate_1(writeout_counter1,writeout_counter2) >-1 & normaliased_log2_comparsion_replicate_1(writeout_counter1,writeout_counter2) <1);
    if ~(sum(Decrease_changes_matrix_across_replicate) >= 2 | sum(Increase_changes_matrix_across_replicate) >= 2 | sum(No_changes_matrix_across_replicate) >= 2)
      Inconsistent_across_replicates(writeout_counter1,writeout_counter2) =1;
    end
  end
end

tt = toc;
fprintf('  ...  %.2f seconds\n',tt)



%% 9. Use mwwtest and ttest to determine if complexes change
fprintf('    9. Use mwwtest and ttest to determine if complexes change (BETWEEN REPLICATES)')


%Count the number of unique gaussian detected
Unique_final_gaussian_counter=0;
for writeout_counter1= 1:Dimension_of_master_gaussian_list(1)
  for writeout_counter2= 1:nnz(Finalised_Master_Gaussian_list.Center(writeout_counter1,:))
    Unique_final_gaussian_counter=1+Unique_final_gaussian_counter;
  end
end

%Create variables to copy log2 values into
Average_fold_change=zeros(Dimension_of_master_gaussian_list(1),Dimension_of_master_gaussian_list(2));
Normalised_average_fold_change=zeros(Dimension_of_master_gaussian_list(1),Dimension_of_master_gaussian_list(2));
Stdev_fold_change=zeros(Dimension_of_master_gaussian_list(1),Dimension_of_master_gaussian_list(2));
Average_log2_value=zeros(1,4);
Average_log2_value_across_replicates=zeros(1,4);
log2_comparsion_replicate_1_for_figure=zeros(Unique_final_gaussian_counter,1);
log2_comparsion_replicate_2_for_figure=zeros(Unique_final_gaussian_counter,1);
log2_comparsion_replicate_3_for_figure=zeros(Unique_final_gaussian_counter,1);
log2_comparsion_replicate_4_for_figure=zeros(Unique_final_gaussian_counter,1);

%Write T_test data in form to use in figures
T_test_stats.d_for_figure=nan(Total_number_of_unique_gaussians_master_list,1);
T_test_stats.h_for_figure=nan(Total_number_of_unique_gaussians_master_list,1);
T_test_stats.p_for_figure=nan(Total_number_of_unique_gaussians_master_list,1);
T_test_stats.ci1_for_figure=nan(Total_number_of_unique_gaussians_master_list,1);
T_test_stats.ci2_for_figure=nan(Total_number_of_unique_gaussians_master_list,1);
T_test_stats.tstat_for_figure=nan(Total_number_of_unique_gaussians_master_list,1);
T_test_stats.df_for_figure=nan(Total_number_of_unique_gaussians_master_list,1);
T_test_stats.sd_for_figure=nan(Total_number_of_unique_gaussians_master_list,1);
T_test_stats.num_obserations_for_figure=nan(Total_number_of_unique_gaussians_master_list,1);
T_test_stats.Adjusted_pvalue_for_figure=cell(Total_number_of_unique_gaussians_master_list,1);

T_test_stats.d=nan(Dimension_of_master_gaussian_list(1),Dimension_of_master_gaussian_list(2));
T_test_stats.h=nan(Dimension_of_master_gaussian_list(1),Dimension_of_master_gaussian_list(2));
T_test_stats.p=nan(Dimension_of_master_gaussian_list(1),Dimension_of_master_gaussian_list(2));
T_test_stats.ci1=nan(Dimension_of_master_gaussian_list(1),Dimension_of_master_gaussian_list(2));
T_test_stats.ci2=nan(Dimension_of_master_gaussian_list(1),Dimension_of_master_gaussian_list(2));
T_test_stats.tstat=nan(Dimension_of_master_gaussian_list(1),Dimension_of_master_gaussian_list(2));
T_test_stats.df=nan(Dimension_of_master_gaussian_list(1),Dimension_of_master_gaussian_list(2));
T_test_stats.sd=nan(Dimension_of_master_gaussian_list(1),Dimension_of_master_gaussian_list(2));
T_test_stats.num_obserations=nan(Dimension_of_master_gaussian_list(1),Dimension_of_master_gaussian_list(2));
T_test_stats.Adjusted_pvalue=cell(Dimension_of_master_gaussian_list(1),Dimension_of_master_gaussian_list(2));

MWW_test_stats.T=nan(Dimension_of_master_gaussian_list(1),Dimension_of_master_gaussian_list(2));
MWW_test_stats.U=nan(Dimension_of_master_gaussian_list(1),Dimension_of_master_gaussian_list(2));
MWW_test_stats.mean=nan(Dimension_of_master_gaussian_list(1),Dimension_of_master_gaussian_list(2));
MWW_test_stats.std=nan(Dimension_of_master_gaussian_list(1),Dimension_of_master_gaussian_list(2));
MWW_test_stats.z=nan(Dimension_of_master_gaussian_list(1),Dimension_of_master_gaussian_list(2));
MWW_test_stats.p=nan(Dimension_of_master_gaussian_list(1),Dimension_of_master_gaussian_list(2));
MWW_test_stats.method=cell(Dimension_of_master_gaussian_list(1),Dimension_of_master_gaussian_list(2));

MWW_test_stats.T_for_figure=cell(Total_number_of_unique_gaussians_master_list,1);
MWW_test_stats.U_for_figure=cell(Total_number_of_unique_gaussians_master_list,1);
MWW_test_stats.mean_for_figure=cell(Total_number_of_unique_gaussians_master_list,1);
MWW_test_stats.std_for_figure=cell(Total_number_of_unique_gaussians_master_list,1);
MWW_test_stats.z_for_figure=cell(Total_number_of_unique_gaussians_master_list,1);
MWW_test_stats.p_for_figure=cell(Total_number_of_unique_gaussians_master_list,1);

%Determine Average fold change and signficances across replicates
fold_counter1=1;
for writeout_counter1= 1:Dimension_of_master_gaussian_list(1) % looping over proteins
  for writeout_counter2= 1:nnz(Finalised_Master_Gaussian_list.Center(writeout_counter1,:)) % looping over... what? gaussians in a complex?
    %Non-normlaised values
    Average_log2_value(1,1)=  log2_comparsion_replicate_1(writeout_counter1,writeout_counter2);
    Average_log2_value(1,2)=  log2_comparsion_replicate_2(writeout_counter1,writeout_counter2);
    Average_log2_value(1,3)=  log2_comparsion_replicate_3(writeout_counter1,writeout_counter2);
    Average_log2_value(1,4)=  log2_comparsion_replicate_4(writeout_counter1,writeout_counter2);
    
    Average_log2_value_across_replicates(fold_counter1,1)=  log2_comparsion_replicate_1(writeout_counter1,writeout_counter2);
    Average_log2_value_across_replicates(fold_counter1,2)=  log2_comparsion_replicate_2(writeout_counter1,writeout_counter2);
    Average_log2_value_across_replicates(fold_counter1,3)=  log2_comparsion_replicate_3(writeout_counter1,writeout_counter2);
    Average_log2_value_across_replicates(fold_counter1,4)=  log2_comparsion_replicate_4(writeout_counter1,writeout_counter2);
    
    %normlaised values
    normalised_average_log2_value(1,1)=normaliased_log2_comparsion_replicate_1(writeout_counter1,writeout_counter2);
    normalised_average_log2_value(1,2)=normaliased_log2_comparsion_replicate_2(writeout_counter1,writeout_counter2);
    normalised_average_log2_value(1,3)=normaliased_log2_comparsion_replicate_3(writeout_counter1,writeout_counter2);
    normalised_average_log2_value(1,4)=normaliased_log2_comparsion_replicate_4(writeout_counter1,writeout_counter2);
    
    %Copy log2 to n*1 matrix for
    log2_comparsion_replicate_1_for_figure(fold_counter1)=log2_comparsion_replicate_1(writeout_counter1,writeout_counter2);
    log2_comparsion_replicate_2_for_figure(fold_counter1)=log2_comparsion_replicate_2(writeout_counter1,writeout_counter2);
    log2_comparsion_replicate_3_for_figure(fold_counter1)=log2_comparsion_replicate_3(writeout_counter1,writeout_counter2);
    log2_comparsion_replicate_4_for_figure(fold_counter1)=log2_comparsion_replicate_4(writeout_counter1,writeout_counter2);
    
    %remove NaN values
    average_values_to_keep=~isnan(Average_log2_value(:));
    normalised_average_values_to_keep=~isnan(normalised_average_log2_value(:));
    
    Average_log2_value=Average_log2_value(average_values_to_keep==1);
    normalised_average_log2_value=normalised_average_log2_value(normalised_average_values_to_keep==1);
    
    if ~isempty(Average_log2_value)
      %determine the average fold change
      Average_fold_change(writeout_counter1,writeout_counter2) = mean(Average_log2_value(:));
      Stdev_fold_change(writeout_counter1,writeout_counter2) =std(Average_log2_value(:));
      
      %Convert MvsL and HvsL measurements to a form to do mann-whitney
      %wilcoxon testing
      MvsL_measurements_all_replicates=[MvsL_Ratio_replicate_1{writeout_counter1,writeout_counter2},MvsL_Ratio_replicate_2{writeout_counter1,writeout_counter2},...
        MvsL_Ratio_replicate_3{writeout_counter1,writeout_counter2},MvsL_Ratio_replicate_4{writeout_counter1,writeout_counter2}];
      
      HvsL_measurements_all_replicates=[HvsL_Ratio_replicate_1{writeout_counter1,writeout_counter2},HvsL_Ratio_replicate_2{writeout_counter1,writeout_counter2},...
        HvsL_Ratio_replicate_3{writeout_counter1,writeout_counter2},HvsL_Ratio_replicate_4{writeout_counter1,writeout_counter2}];
      
      %remove nan and zeros
      MvsL_measurements_all_replicates(isnan(MvsL_measurements_all_replicates(1,:)))=[];
      HvsL_measurements_all_replicates(isnan(HvsL_measurements_all_replicates(1,:)))=[];
      MvsL_measurements_all_replicates((MvsL_measurements_all_replicates(1,:))==0)=[];
      HvsL_measurements_all_replicates((HvsL_measurements_all_replicates(1,:))==0)=[];
      
      STATS=mwwtest(MvsL_measurements_all_replicates(:),HvsL_measurements_all_replicates(:));
      
      %Values to array to use later in 1*n array for figure/Tables
      try
        MWW_test_stats.T_for_figure(fold_counter1,1)={STATS.T};
        MWW_test_stats.U_for_figure(fold_counter1,1)={STATS.U};
        MWW_test_stats.mean_for_figure(fold_counter1,1)={STATS.mean};
        MWW_test_stats.std_for_figure(fold_counter1,1)={STATS.std};
        MWW_test_stats.z_for_figure(fold_counter1,1)={STATS.z};
        MWW_test_stats.p_for_figure(fold_counter1,1)={STATS.p};
      catch
        MWW_test_stats.p_for_figure(fold_counter1,1)={STATS.p};
        MWW_test_stats.U_for_figure(fold_counter1,1)={STATS.U};
        MWW_test_stats.T_for_figure(fold_counter1,1)={STATS.T};
        MWW_test_stats.mean_for_figure(fold_counter1,1)={NaN};
        MWW_test_stats.std_for_figure(fold_counter1,1)={NaN};
        MWW_test_stats.z_for_figure(fold_counter1,1)={NaN};
      end
      %in n*m array
      try
        MWW_test_stats.method(writeout_counter1,writeout_counter2)={STATS.method};
        MWW_test_stats.T(writeout_counter1,writeout_counter2)=STATS.T;
        MWW_test_stats.U(writeout_counter1,writeout_counter2)=STATS.U;
        MWW_test_stats.mean(writeout_counter1,writeout_counter2)=STATS.mean;
        MWW_test_stats.std(writeout_counter1,writeout_counter2)=STATS.std;
        MWW_test_stats.z(writeout_counter1,writeout_counter2)=STATS.z;
        MWW_test_stats.p(writeout_counter1,writeout_counter2)=STATS.p;
      catch
        MWW_test_stats.method(writeout_counter1,writeout_counter2)={STATS.method};
        MWW_test_stats.p(writeout_counter1,writeout_counter2)=STATS.p;
        MWW_test_stats.U(writeout_counter1,writeout_counter2)=STATS.U;
        MWW_test_stats.T(writeout_counter1,writeout_counter2)=STATS.T;
        MWW_test_stats.mean(writeout_counter1,writeout_counter2)=nan;
        MWW_test_stats.std(writeout_counter1,writeout_counter2)=nan;
        MWW_test_stats.z(writeout_counter1,writeout_counter2)=nan;
      end
      
      % Do T test
      [ttest_h,ttest_p,ttest_ci,ttest_stats]=ttest2(MvsL_measurements_all_replicates(:),HvsL_measurements_all_replicates(:));
      %Copy T test values to array to use later
      %in 1*n array for figure/Tables
      T_test_stats.h_for_figure(fold_counter1,1)=ttest_h;
      T_test_stats.p_for_figure(fold_counter1,1)=ttest_p;
      T_test_stats.ci1_for_figure(fold_counter1,1)=ttest_ci(1);
      T_test_stats.ci2_for_figure(fold_counter1,1)=ttest_ci(2);
      T_test_stats.tstat_for_figure(fold_counter1,1)=ttest_stats.tstat(1);
      T_test_stats.df_for_figure(fold_counter1,1)=ttest_stats.df(1);
      T_test_stats.sd_for_figure(fold_counter1,1)=ttest_stats.sd(1);
      %in n*m array
      T_test_stats.h(writeout_counter1,writeout_counter2)=ttest_h;
      T_test_stats.p(writeout_counter1,writeout_counter2)=ttest_p;
      T_test_stats.ci1(writeout_counter1,writeout_counter2)=ttest_ci(1);
      T_test_stats.ci2(writeout_counter1,writeout_counter2)=ttest_ci(2);
      T_test_stats.tstat(writeout_counter1,writeout_counter2)=ttest_stats.tstat(1);
      T_test_stats.df(writeout_counter1,writeout_counter2)=ttest_stats.df(1);
      T_test_stats.sd(writeout_counter1,writeout_counter2)=ttest_stats.sd(1);
      
      %Number of observation of the gaussian
      T_test_stats.num_obserations(writeout_counter1,writeout_counter2)=length(Average_log2_value);
      T_test_stats.num_obserations_for_figure(fold_counter1,1)=length(Average_log2_value);
      
      %Multiple hyopthesis correction, using Bonferroni correction based on the
      %number of replicates in with the Gaussian was observed, if signficant
      %add 1, if not 0
      if ttest_p<0.05/Total_number_of_unique_gaussians_master_list
        T_test_stats.Adjusted_pvalue_for_figure(fold_counter1,1)= {'+'};
        T_test_stats.Adjusted_pvalue(writeout_counter1,writeout_counter2)= {'+'};
      else
        T_test_stats.Adjusted_pvalue_for_figure(fold_counter1,1)= {' '};
        T_test_stats.Adjusted_pvalue(writeout_counter1,writeout_counter2)= {' '};
      end
      
      if MWW_test_stats.p(writeout_counter1,writeout_counter2)<0.05/Total_number_of_unique_gaussians_master_list
        MWW_test_stats.Adjusted_pvalue_for_figure(fold_counter1,1)= {'+'};
        MWW_test_stats.Adjusted_pvalue(writeout_counter1,writeout_counter2)= {'+'};
      else
        MWW_test_stats.Adjusted_pvalue_for_figure(fold_counter1,1)= {' '};
        MWW_test_stats.Adjusted_pvalue(writeout_counter1,writeout_counter2)= {' '};
      end
      
      Normalised_average_fold_change(writeout_counter1,writeout_counter2)=mean(normalised_average_log2_value(:));
      
    elseif isempty(Average_log2_value)
      %determine the average fold change
      Average_fold_change(writeout_counter1,writeout_counter2)= 0;
    end
    fold_counter1=fold_counter1+1;
  end
end


%Write comparsion output to combined structure array
Finalised_Master_Gaussian_list.Do_gaussian_in_channels_change_replicate_1=Do_gaussian_in_channels_change_replicate_1;
Finalised_Master_Gaussian_list.Do_gaussian_in_channels_change_replicate_2=Do_gaussian_in_channels_change_replicate_2;
Finalised_Master_Gaussian_list.Do_gaussian_in_channels_change_replicate_3=Do_gaussian_in_channels_change_replicate_3;
Finalised_Master_Gaussian_list.Do_gaussian_in_channels_change_replicate_4=Do_gaussian_in_channels_change_replicate_4;

Finalised_Master_Gaussian_list.normaliased_log2_comparsion_replicate_1= normaliased_log2_comparsion_replicate_1;
Finalised_Master_Gaussian_list.normaliased_log2_comparsion_replicate_2= normaliased_log2_comparsion_replicate_2;
Finalised_Master_Gaussian_list.normaliased_log2_comparsion_replicate_3= normaliased_log2_comparsion_replicate_3;
Finalised_Master_Gaussian_list.normaliased_log2_comparsion_replicate_4= normaliased_log2_comparsion_replicate_4;

Finalised_Master_Gaussian_list.log2_comparsion_replicate_1=log2_comparsion_replicate_1;
Finalised_Master_Gaussian_list.log2_comparsion_replicate_2=log2_comparsion_replicate_2;
Finalised_Master_Gaussian_list.log2_comparsion_replicate_3=log2_comparsion_replicate_3;
Finalised_Master_Gaussian_list.log2_comparsion_replicate_4=log2_comparsion_replicate_4;

Finalised_Master_Gaussian_list.Protein_location_replicate1=Protein_location_replicate1;
Finalised_Master_Gaussian_list.Protein_location_replicate2=Protein_location_replicate2;
Finalised_Master_Gaussian_list.Protein_location_replicate3=Protein_location_replicate3;
Finalised_Master_Gaussian_list.Protein_location_replicate4=Protein_location_replicate4;

Finalised_Master_Gaussian_list.No_change_across_replicates=No_change_across_replicates;
Finalised_Master_Gaussian_list.Increase_change_across_replicates=Increase_change_across_replicates;
Finalised_Master_Gaussian_list.Decrease_change_across_replicates=Decrease_change_across_replicates;
Finalised_Master_Gaussian_list.Inconsistent_across_replicates=Inconsistent_across_replicates;

Finalised_Master_Gaussian_list.Average_fold_change=Average_fold_change;
Finalised_Master_Gaussian_list.Normalised_average_fold_change=Normalised_average_fold_change;
Finalised_Master_Gaussian_list.Stdev_fold_change=Stdev_fold_change;

%Add ttest output
Finalised_Master_Gaussian_list.Ttest_p_values=T_test_stats.p;
Finalised_Master_Gaussian_list.degree_freedom=T_test_stats.df;
Finalised_Master_Gaussian_list.Number_observation=T_test_stats.num_obserations;
Finalised_Master_Gaussian_list.Satifies_Adjusted_pvalue=T_test_stats.Adjusted_pvalue;

%add MWW test output
Finalised_Master_Gaussian_list.MWW_p_values=MWW_test_stats.p;
Finalised_Master_Gaussian_list.MWW_U=MWW_test_stats.U;
Finalised_Master_Gaussian_list.MWW_Satifies_Adjusted_pvalue=MWW_test_stats.Adjusted_pvalue;

%Create varibles to write to
Decrease_change_across_protein_replicates_persaus_input=cell(Dimension_of_master_gaussian_list(1),1);
Decrease_change_across_replicates_persaus_input=cell(Dimension_of_master_gaussian_list(1),Dimension_of_master_gaussian_list(2));
Increase_change_across_protein_replicates_persaus_input=cell(Dimension_of_master_gaussian_list(1),1);
Increase_change_across_replicates_persaus_input=cell(Dimension_of_master_gaussian_list(1),Dimension_of_master_gaussian_list(2));
Changes__across_protein_replicates_persaus_input=cell(Dimension_of_master_gaussian_list(1),3);

%Convert Increase/Decrease changes to a format usible by Perseus
for writeout_counter1= 1:Dimension_of_master_gaussian_list(1)
  for writeout_counter2= 1:nnz(log_2_fold_comparsion_of_gaussian_replicate_3(writeout_counter1,:))
    %test if positions is equal to 1
    if Finalised_Master_Gaussian_list.Increase_change_across_replicates(writeout_counter1,writeout_counter2) == 1
      Increase_change_across_replicates_persaus_input(writeout_counter1,writeout_counter2)= {'+'};
      Increase_change_across_protein_replicates_persaus_input(writeout_counter1)= {'+'};
    else
      Increase_change_across_replicates_persaus_input(writeout_counter1,writeout_counter2)= {' '};
    end
    %test if positions is equal to 1
    if Finalised_Master_Gaussian_list.Decrease_change_across_replicates(writeout_counter1,writeout_counter2) == 1
      Decrease_change_across_replicates_persaus_input(writeout_counter1,writeout_counter2)= {'+'};
      Decrease_change_across_protein_replicates_persaus_input(writeout_counter1)= {'+'};
    else
      Decrease_change_across_replicates_persaus_input(writeout_counter1,writeout_counter2)= {' '};
    end
  end
  
  %Test if chnages observed at the protein level
  if (isequal(Decrease_change_across_protein_replicates_persaus_input(writeout_counter1),{'+'}) &...
      isequal(Increase_change_across_protein_replicates_persaus_input(writeout_counter1),{'+'}))
    Changes__across_protein_replicates_persaus_input(writeout_counter1,1)= {'+'};
    Changes__across_protein_replicates_persaus_input(writeout_counter1,2)= {' '};
    Changes__across_protein_replicates_persaus_input(writeout_counter1,3)= {' '};
  elseif (isequal(Increase_change_across_protein_replicates_persaus_input(writeout_counter1),{'+'}))
    Changes__across_protein_replicates_persaus_input(writeout_counter1,1)= {' '};
    Changes__across_protein_replicates_persaus_input(writeout_counter1,2)= {'+'};
    Changes__across_protein_replicates_persaus_input(writeout_counter1,3)= {' '};
  elseif (isequal(Decrease_change_across_protein_replicates_persaus_input(writeout_counter1),{'+'}))
    Changes__across_protein_replicates_persaus_input(writeout_counter1,1)= {' '};
    Changes__across_protein_replicates_persaus_input(writeout_counter1,2)= {' '};
    Changes__across_protein_replicates_persaus_input(writeout_counter1,3)= {'+'};
  elseif ~(isequal(Decrease_change_across_protein_replicates_persaus_input(writeout_counter1),{'+'})|...
      isequal(Increase_change_across_protein_replicates_persaus_input(writeout_counter1),{'+'}))
    Changes__across_protein_replicates_persaus_input(writeout_counter1,1)= {' '};
    Changes__across_protein_replicates_persaus_input(writeout_counter1,2)= {' '};
    Changes__across_protein_replicates_persaus_input(writeout_counter1,3)= {' '};
  end
  
end

%Determine if changes were detected by multiple hypothesis testing in
%protein group
Singficant_change_ttest_persaus_protein=cell(Dimension_of_master_gaussian_list(1),1);
Singficant_change_mww_persaus_protein=cell(Dimension_of_master_gaussian_list(1),1);
Singficant_change_ttest_persaus_protein(cellfun(@isempty, Singficant_change_ttest_persaus_protein))={' '};
Singficant_change_mww_persaus_protein(cellfun(@isempty, Singficant_change_mww_persaus_protein))={' '};

for writeout_counter1= 1:Dimension_of_master_gaussian_list(1)
  for writeout_counter2= 1:nnz(log_2_fold_comparsion_of_gaussian_replicate_3(writeout_counter1,:))
    %Check MWW test
    if isequal(Finalised_Master_Gaussian_list.MWW_Satifies_Adjusted_pvalue(writeout_counter1,writeout_counter2),{'+'})
      Singficant_change_mww_persaus_protein(writeout_counter1)= {'+'};
    end
    %Check ttest
    if isequal(Finalised_Master_Gaussian_list.Satifies_Adjusted_pvalue(writeout_counter1,writeout_counter2),{'+'})
      Singficant_change_ttest_persaus_protein(writeout_counter1)= {'+'};
    end
  end
end


%Copy to Final Gaussian list
Finalised_Master_Gaussian_list.Increase_change_across_replicates_persaus_input=Increase_change_across_replicates_persaus_input;
Finalised_Master_Gaussian_list.Decrease_change_across_replicates_persaus_input=Decrease_change_across_replicates_persaus_input;
Finalised_Master_Gaussian_list.Increase_and_Decrease_protein_persaus_input=Changes__across_protein_replicates_persaus_input(:,1);
Finalised_Master_Gaussian_list.Increase_protein_persaus_input=Changes__across_protein_replicates_persaus_input(:,2);
Finalised_Master_Gaussian_list.Decrease_protein_persaus_input=Changes__across_protein_replicates_persaus_input(:,3);

tt = toc;
fprintf('  ...  %.2f seconds\n',tt)



%% 10. Determine the coverage (number of datapoint observed) for reach protein of each replicate
fprintf('    10. Determine the coverage (number of datapoint observed) (BETWEEN REPLICATES)')

%Determine size of values to caliculate percentage
size_of_num_val_MvsL_for_figures=size(num_val_MvsL_for_figures);
size_of_num_val_HvsL_for_figures=size(num_val_HvsL_for_figures);

%generate array to write values to
Coverage_array_MvsL=zeros(Dimension_of_master_gaussian_list(1),4);
Coverage_array_HvsL=zeros(Dimension_of_master_gaussian_list(1),4);
Coverage_array_summed_value_MvsL=zeros(Dimension_of_master_gaussian_list(1),4);
Coverage_array_summed_value_HvsL=zeros(Dimension_of_master_gaussian_list(1),4);
Coverage_array_summed_value_HvsM=zeros(Dimension_of_master_gaussian_list(1),4);

%MvsL values
for Gaussian_counter1= 1:Dimension_of_master_gaussian_list(1)
  
  %locate protein to compare
  [location_Protein_in_raw_textdata]=ind2sub(size(txt_MvsL(:,2)),...
    strmatch(Finalised_Master_Gaussian_list.Protein_name{Gaussian_counter1},txt_MvsL(:,2),'exact'));
  
  %count how many replicates protein was observed in
  number_of_replicates_observed=length(location_Protein_in_raw_textdata);
  
  for replicate_counter = 1:number_of_replicates_observed
    %Detemine replicate observed in
    replicate_temp_value=num_val_MvsL_for_figures((location_Protein_in_raw_textdata(replicate_counter)-1),1);
    %Determine number of non na values
    number_of_observation=sum(~isnan(num_val_MvsL_for_figures((location_Protein_in_raw_textdata(replicate_counter)-1),2:end)));
    
    %Write value to array
    Coverage_array_MvsL(Gaussian_counter1,replicate_temp_value)= number_of_observation/(size_of_num_val_MvsL_for_figures(2)-1);
  end
end

%MvsL values to sum values
for Gaussian_counter1= 1:Dimension_of_master_gaussian_list(1)
  
  %locate protein to compare
  [location_Protein_in_raw_textdata]=ind2sub(size(txt_MvsL(:,2)),...
    strmatch(Finalised_Master_Gaussian_list.Protein_name{Gaussian_counter1},txt_MvsL(:,2),'exact'));
  
  %count how many replicates protein was observed in
  number_of_replicates_observed=length(location_Protein_in_raw_textdata);
  
  for replicate_counter = 1:number_of_replicates_observed
    %Detemine replicate observed in
    replicate_temp_value=num_val_MvsL_for_figures((location_Protein_in_raw_textdata(replicate_counter)-1),1);
    %Determine number of non na values
    values_MvsL_for_figures=num_val_MvsL_for_figures((location_Protein_in_raw_textdata(replicate_counter)-1),2:end);
    number_of_observation=isnan(values_MvsL_for_figures);
    
    %remove missing values
    values_MvsL_for_figures(number_of_observation)=[];
    
    %Write value to array
    Coverage_array_summed_value_MvsL(Gaussian_counter1,replicate_temp_value)= sum(values_MvsL_for_figures)/(fraction_to_plot*1/Diltuion_factor_master_mix);
  end
end

%HvsL values
for Gaussian_counter1= 1:Dimension_of_master_gaussian_list(1)
  
  %locate protein to compare
  [location_Protein_in_raw_textdata]=ind2sub(size(txt_HvsL(:,2)),...
    strmatch(Finalised_Master_Gaussian_list.Protein_name{Gaussian_counter1},txt_HvsL(:,2),'exact'));
  
  %count how many replicates protein was observed in
  number_of_replicates_observed=length(location_Protein_in_raw_textdata);
  
  for replicate_counter = 1:number_of_replicates_observed
    %Detemine replicate observed in
    replicate_temp_value=num_val_HvsL_for_figures((location_Protein_in_raw_textdata(replicate_counter)-1),1);
    %Determine number of non na values
    number_of_observation=sum(~isnan(num_val_HvsL_for_figures((location_Protein_in_raw_textdata(replicate_counter)-1),2:end)));
    
    %Write value to array
    Coverage_array_HvsL(Gaussian_counter1,replicate_temp_value)= number_of_observation/(size_of_num_val_HvsL_for_figures(2)-1);
  end
end

%HvsL values to sum values
for Gaussian_counter1= 1:Dimension_of_master_gaussian_list(1)
  
  %locate protein to compare
  [location_Protein_in_raw_textdata]=ind2sub(size(txt_HvsL(:,2)),...
    strmatch(Finalised_Master_Gaussian_list.Protein_name{Gaussian_counter1},txt_HvsL(:,2),'exact'));
  
  %count how many replicates protein was observed in
  number_of_replicates_observed=length(location_Protein_in_raw_textdata);
  
  for replicate_counter = 1:number_of_replicates_observed
    %Detemine replicate observed in
    replicate_temp_value=num_val_HvsL_for_figures((location_Protein_in_raw_textdata(replicate_counter)-1),1);
    %Determine number of non na values
    values_HvsL_for_figures=num_val_HvsL_for_figures((location_Protein_in_raw_textdata(replicate_counter)-1),2:end);
    number_of_observation=isnan(values_HvsL_for_figures);
    
    %remove missing values
    values_HvsL_for_figures(number_of_observation)=[];
    
    %Write value to array
    Coverage_array_summed_value_HvsL(Gaussian_counter1,replicate_temp_value)= sum(values_HvsL_for_figures)/(fraction_to_plot*1/Diltuion_factor_master_mix);
  end
end


%Determine the global ratio of HvsM
for Gaussian_counter1= 1:Dimension_of_master_gaussian_list(1)
  
  %locate protein to compare
  [location_Protein_in_raw_textdata]=ind2sub(size(txt_HvsL(:,2)),...
    strmatch(Finalised_Master_Gaussian_list.Protein_name{Gaussian_counter1},txt_HvsL(:,2),'exact'));
  
  %count how many replicates protein was observed in
  number_of_replicates_observed=length(location_Protein_in_raw_textdata);
  
  for replicate_counter = 1:number_of_replicates_observed
    %Detemine replicate observed in
    replicate_temp_value=num_val_HvsM((location_Protein_in_raw_textdata(replicate_counter)-1),1);
    %Determine number of non na values
    num_val_HvsM_temp=num_val_HvsM((location_Protein_in_raw_textdata(replicate_counter)-1),2:end);
    number_of_observation=isnan(num_val_HvsM_temp);
    
    %remove missing values
    num_val_HvsM_temp(number_of_observation)=[];
    
    %Write value to array
    Coverage_array_summed_value_HvsM(Gaussian_counter1,replicate_temp_value)= sum(num_val_HvsM_temp)/(fraction_to_plot*1/Diltuion_factor_master_mix);
  end
end


%Write values to Master list
Finalised_Master_Gaussian_list.Summed_Gaussian_areaMvsL=Coverage_array_summed_value_MvsL;
Finalised_Master_Gaussian_list.Summed_Gaussian_areaHvsL=Coverage_array_summed_value_HvsL;
Finalised_Master_Gaussian_list.Ratio_Gaussian_area=Coverage_array_summed_value_HvsM;
Finalised_Master_Gaussian_list.SEC_coverageMvsL=Coverage_array_MvsL;
Finalised_Master_Gaussian_list.SEC_coverageHvsL=Coverage_array_HvsL;


%Determine the global trend of guassians (as observed within atleast two
%replicates)
Finalised_Master_Gaussian_list.global_change_trends=cell(Dimension_of_master_gaussian_list(1),1);
global_no_change=0;
global_increase=0;
global_decrease=0;
global_increase_decrease=0;
global_inconsistent_across_replicates=0;

for writeout_counter1= 1:Dimension_of_master_gaussian_list(1)
  %Reset value
  global_trend_signficant=0;
  %determine how many Gaussians to assess in the row
  global_trend1=nnz(Finalised_Master_Gaussian_list.No_change_across_replicates(writeout_counter1,:));
  global_trend2=nnz(Finalised_Master_Gaussian_list.Increase_change_across_replicates(writeout_counter1,:));
  global_trend3=nnz(Finalised_Master_Gaussian_list.Decrease_change_across_replicates(writeout_counter1,:));
  global_trend4=nnz(Finalised_Master_Gaussian_list.Inconsistent_across_replicates(writeout_counter1,:));
  global_trend5=nnz(Finalised_Master_Gaussian_list.Center(writeout_counter1,:));
  global_trend_signficant=0;
  
  significance_Bon_variable =zeros(1,1);
  significance_increase_variable =zeros(1,1);
  significance_decrease_variable =zeros(1,1);
  
  significance_temp_variable=zeros(global_trend5,1);
  for Signficance_loop=1:global_trend5
    significance_Bon_variable(Signficance_loop,1)=isequal(Finalised_Master_Gaussian_list.MWW_Satifies_Adjusted_pvalue(writeout_counter1, Signficance_loop),{'+'});
    significance_increase_variable(Signficance_loop,1)=Finalised_Master_Gaussian_list.Average_fold_change(writeout_counter1,Signficance_loop) >=0;
    significance_decrease_variable(Signficance_loop,1)=Finalised_Master_Gaussian_list.Average_fold_change(writeout_counter1,Signficance_loop) <=0;
  end
  
  if sum(significance_Bon_variable)==global_trend5
    global_trend_signficant=1;
  end
  
  summed_number_increase=sum(significance_increase_variable);
  summed_number_decrease=sum(significance_decrease_variable);
  
  if global_trend_signficant == 0
    Finalised_Master_Gaussian_list.global_change_trends(writeout_counter1,1)= {'No Change'};
    global_no_change=1+global_no_change;
    
  elseif  global_trend_signficant ==1 & summed_number_increase == global_trend5
    Finalised_Master_Gaussian_list.global_change_trends(writeout_counter1,1)= {'Increase'};
    global_increase=1+global_increase;
    
  elseif global_trend_signficant ==1 & summed_number_decrease == global_trend5
    Finalised_Master_Gaussian_list.global_change_trends(writeout_counter1,1)= {'Decrease'};
    global_decrease=1+global_decrease;
    
  elseif global_trend_signficant ==1 & ((summed_number_increase + summed_number_decrease) == global_trend5)
    Finalised_Master_Gaussian_list.global_change_trends(writeout_counter1,1)= {'+/-'};
    global_increase_decrease=1+ global_increase_decrease;
    
  elseif ~((summed_number_increase + summed_number_decrease) == global_trend5) | global_trend_signficant == 0
    Finalised_Master_Gaussian_list.global_change_trends(writeout_counter1,1)= {'No Change'}; %If a change is not considered signficant it is denoted as not change
    global_inconsistent_across_replicates=1 + global_inconsistent_across_replicates;
  end
  
end

% Count how many data points quantifed across all four biological
% replicates
replicate_1_counter=0;
replicate_2_counter=0;
replicate_3_counter=0;
replicate_4_counter=0;

for writeout_counter1= 1:Dimension_of_master_gaussian_list(1)
  for writeout_counter2= 1:nnz(Finalised_Master_Gaussian_list.Center(writeout_counter1,:))
    if ~isnan(Finalised_Master_Gaussian_list.normaliased_log2_comparsion_replicate_1(writeout_counter1,writeout_counter2))
      replicate_1_counter=replicate_1_counter+1;
    end
    if ~isnan(Finalised_Master_Gaussian_list.normaliased_log2_comparsion_replicate_2(writeout_counter1,writeout_counter2))
      replicate_2_counter=replicate_2_counter+1;
    end
    if ~isnan(Finalised_Master_Gaussian_list.normaliased_log2_comparsion_replicate_3(writeout_counter1,writeout_counter2))
      replicate_3_counter=replicate_3_counter+1;
    end
    if ~isnan(Finalised_Master_Gaussian_list.normaliased_log2_comparsion_replicate_4(writeout_counter1,writeout_counter2))
      replicate_4_counter=replicate_4_counter+1;
    end
    
  end
end

tt = toc;
fprintf('  ...  %.2f seconds\n',tt)



%% 11. Write out tables
fprintf('    11. Write out tables')

% #output
% write out table of master guassians list
fid_master_list1= fopen([datadir1 'Master_guassian_list_fold_comparsion.csv'],'wt'); % create the output file with the header infromation
fprintf (fid_master_list1,'%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,\n',... %header for OutputGaus output
  'Protein name', 'Unique_identifier (of best gaussian)','Replicate (of best gaussian)',...
  'Channel (of best gaussian)', 'Guassian_index_number (of best gaussian)',...
  'Center (of best gaussian)', 'Height (of best gaussian)', 'Width (of best gaussian)',...
  'SSE (of best gaussian)', 'adjrsquare (of best gaussian)','Gaussian Area (of best gaussian)',...
  'Complex Size (of best gaussian)', 'Change observed replicate 1',...
  'Change observed replicate 2','Change observed replicate 3',...
  'Change observed replicate 4','Std Dev of fold change','Normalised fold change',...
  'Normalised fold change replicate 1','Normalised fold change replicate 2',...
  'Normalised fold change replicate 3','Normalised fold change replicate 4',...
  'Average fold change','Fold change replicate 1', 'Fold change replicate 2',...
  'Fold change replicate 3', 'Fold change replicate 4', 'p-value (ttest)',...
  'Below Adjusted_pvalue(bonferroni correction) (ttest)','p-value (MWW)',...
  'Below Adjusted_pvalue(bonferroni correction) (MWW)'); %Write Header
for writeout_counter1= 1:Dimension_of_master_gaussian_list(1)
  for writeout_counter2= 1:nnz(Finalised_Master_Gaussian_list.Center(writeout_counter1,:))
    fprintf(fid_master_list1,'%s,%s,%6.3f,%s,%6.3f,%6.3f,%6.3f,%6.3f,%6.3f,%6.3f,%6.3f,%6.3f,%s,%s,%s,%s,%6.3f,%6.3f,%6.3f,%6.3f,%6.3f,%6.3f,%6.3f,%6.3f,%6.3f,%6.3f,%6.3f,%6.9f,%s,%6.9f,%s\n',...
      Finalised_Master_Gaussian_list.Protein_name{writeout_counter1},...
      Finalised_Master_Gaussian_list.Unique_identifier{writeout_counter1,writeout_counter2},...
      Finalised_Master_Gaussian_list.Replicate(writeout_counter1,writeout_counter2),...
      Finalised_Master_Gaussian_list.Channel{writeout_counter1,writeout_counter2},...
      Finalised_Master_Gaussian_list.Guassian_index_number(writeout_counter1,writeout_counter2),...
      Finalised_Master_Gaussian_list.Center(writeout_counter1,writeout_counter2),...
      Finalised_Master_Gaussian_list.Height(writeout_counter1,writeout_counter2),...
      Finalised_Master_Gaussian_list.Width(writeout_counter1,writeout_counter2),...
      Finalised_Master_Gaussian_list.SSE(writeout_counter1,writeout_counter2),...
      Finalised_Master_Gaussian_list.adjrsquare(writeout_counter1,writeout_counter2),...
      Finalised_Master_Gaussian_list.Gaussian_area(writeout_counter1,writeout_counter2),...
      Finalised_Master_Gaussian_list.Complex_size(writeout_counter1,writeout_counter2),...
      Finalised_Master_Gaussian_list.Do_gaussian_in_channels_change_replicate_1{writeout_counter1,writeout_counter2},...
      Finalised_Master_Gaussian_list.Do_gaussian_in_channels_change_replicate_2{writeout_counter1,writeout_counter2},...
      Finalised_Master_Gaussian_list.Do_gaussian_in_channels_change_replicate_3{writeout_counter1,writeout_counter2},...
      Finalised_Master_Gaussian_list.Do_gaussian_in_channels_change_replicate_4{writeout_counter1,writeout_counter2},...
      Finalised_Master_Gaussian_list.Stdev_fold_change(writeout_counter1,writeout_counter2),...
      Finalised_Master_Gaussian_list.Normalised_average_fold_change(writeout_counter1,writeout_counter2),...
      Finalised_Master_Gaussian_list.normaliased_log2_comparsion_replicate_1(writeout_counter1,writeout_counter2),...
      Finalised_Master_Gaussian_list.normaliased_log2_comparsion_replicate_2(writeout_counter1,writeout_counter2),...
      Finalised_Master_Gaussian_list.normaliased_log2_comparsion_replicate_3(writeout_counter1,writeout_counter2),...
      Finalised_Master_Gaussian_list.normaliased_log2_comparsion_replicate_4(writeout_counter1,writeout_counter2),...
      Finalised_Master_Gaussian_list.Average_fold_change(writeout_counter1,writeout_counter2),...
      Finalised_Master_Gaussian_list.log2_comparsion_replicate_1(writeout_counter1,writeout_counter2),...
      Finalised_Master_Gaussian_list.log2_comparsion_replicate_2(writeout_counter1,writeout_counter2),...
      Finalised_Master_Gaussian_list.log2_comparsion_replicate_3(writeout_counter1,writeout_counter2),...
      Finalised_Master_Gaussian_list.log2_comparsion_replicate_4(writeout_counter1,writeout_counter2),...
      Finalised_Master_Gaussian_list.Ttest_p_values(writeout_counter1,writeout_counter2),...
      Finalised_Master_Gaussian_list.Satifies_Adjusted_pvalue{writeout_counter1,writeout_counter2},...
      Finalised_Master_Gaussian_list.MWW_p_values(writeout_counter1,writeout_counter2),...
      Finalised_Master_Gaussian_list.MWW_Satifies_Adjusted_pvalue{writeout_counter1,writeout_counter2});
  end
end
fclose(fid_master_list1);


%Write out a table of all T-test information
fid_ttest= fopen([datadir1 'Student_ttests_results.csv'],'wt'); % create the output file with the header infromation
fprintf (fid_ttest,'%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,\n',... %header for OutputGaus outpufid_master_list1= fopen('Master_guassian_list_fold_comparsion.csv','wt'); % create the output file with the header infromation
  'Protein name', 'Unique_identifier (of best gaussian)','Replicate (of best gaussian)',...
  'Channel (of best gaussian)', 'Guassian_index_number (of best gaussian)',...
  'Center (of best gaussian)', 'Height (of best gaussian)', 'Width (of best gaussian)',...
  'SSE (of best gaussian)', 'adjrsquare (of best gaussian)','Gaussian Area (of best gaussian)',...
  'Complex Size (of best gaussian)', 'Fold change replicate 1', 'Fold change replicate 2',...
  'Fold change replicate 3', 'Fold change replicate 4', 'p-value (ttest)','Std Dev of fold change',...
  'Number of observation', 'Degrees of freedom','Below Adjusted_pvalue '); %Write Header
for writeout_counter1= 1:Dimension_of_master_gaussian_list(1)
  for writeout_counter2= 1:nnz(Finalised_Master_Gaussian_list.Center(writeout_counter1,:))
    fprintf(fid_ttest,'%s,%s,%6.3f,%s,%6.3f,%6.3f,%6.3f,%6.3f,%6.3f,%6.3f,%6.3f,%6.3f,%6.3f,%6.3f,%6.3f,%6.3f,%6.3f,%6.3f,%6.3f,%6.3f,%s\n',...
      Finalised_Master_Gaussian_list.Protein_name{writeout_counter1},...
      Finalised_Master_Gaussian_list.Unique_identifier{writeout_counter1,writeout_counter2},...
      Finalised_Master_Gaussian_list.Replicate(writeout_counter1,writeout_counter2),...
      Finalised_Master_Gaussian_list.Channel{writeout_counter1,writeout_counter2},...
      Finalised_Master_Gaussian_list.Guassian_index_number(writeout_counter1,writeout_counter2),...
      Finalised_Master_Gaussian_list.Center(writeout_counter1,writeout_counter2),...
      Finalised_Master_Gaussian_list.Height(writeout_counter1,writeout_counter2),...
      Finalised_Master_Gaussian_list.Width(writeout_counter1,writeout_counter2),...
      Finalised_Master_Gaussian_list.SSE(writeout_counter1,writeout_counter2),...
      Finalised_Master_Gaussian_list.adjrsquare(writeout_counter1,writeout_counter2),...
      Finalised_Master_Gaussian_list.Gaussian_area(writeout_counter1,writeout_counter2),...
      Finalised_Master_Gaussian_list.Complex_size(writeout_counter1,writeout_counter2),...
      Finalised_Master_Gaussian_list.log2_comparsion_replicate_1(writeout_counter1,writeout_counter2),...
      Finalised_Master_Gaussian_list.log2_comparsion_replicate_2(writeout_counter1,writeout_counter2),...
      Finalised_Master_Gaussian_list.log2_comparsion_replicate_3(writeout_counter1,writeout_counter2),...
      Finalised_Master_Gaussian_list.log2_comparsion_replicate_4(writeout_counter1,writeout_counter2),...
      Finalised_Master_Gaussian_list.Ttest_p_values(writeout_counter1,writeout_counter2),...
      Finalised_Master_Gaussian_list.Stdev_fold_change(writeout_counter1,writeout_counter2),...
      Finalised_Master_Gaussian_list.Number_observation(writeout_counter1,writeout_counter2),...
      Finalised_Master_Gaussian_list.degree_freedom(writeout_counter1,writeout_counter2),...
      Finalised_Master_Gaussian_list.Satifies_Adjusted_pvalue{writeout_counter1,writeout_counter2});
  end
end
fclose(fid_ttest);

fid_summary_master_protein_comparsion1= fopen([datadir1 'Summary_changes_Protein_across_replicate.csv'],'wt'); % create the output file with the header infromation
fprintf (fid_summary_master_protein_comparsion1,'%s,%s,%s,%s,%s\n',... %header for OutputGaus output
  'Increase(in atleast two replicates)', 'Decrease (in atleast two replicates)', 'No change (in atleast two replicates)',...
  '+/- (in atleast two replicates)', 'Inconsistent across replicates'); %Write Header
fprintf (fid_summary_master_protein_comparsion1,'%6.4f,%6.4f,%6.4f,%6.4f,%6.4f\n',...
  global_increase, global_decrease, global_no_change,...
  global_increase_decrease, global_inconsistent_across_replicates);
fclose(fid_summary_master_protein_comparsion1);

fid_summary_master_comparsion1= fopen([datadir1 'Summary_changes_based_master_gaussian_list.csv'],'wt'); % create the output file with the header infromation
fprintf (fid_summary_master_comparsion1,'%s,%s,%s,%s\n',... %header for OutputGaus output
  'Increase(in atleast two replicates without observation of no changes)', 'Decrease (in atleast two replicates without observation of no changes)',...
  'No change (in atleast two replicates)',...
  'Inconsistent across replicates'); %Write Header
fprintf (fid_summary_master_comparsion1,'%6.4f,%6.4f,%6.4f,%6.4f,%6.4f,%6.4f,%6.4f,%6.4f,%6.4f,%6.4f,%6.4f,%6.4f\n',...
  sum(sum(Increase_change_across_replicates & ~(Decrease_change_across_replicates | No_change_across_replicates | Inconsistent_across_replicates))),...
  sum(sum(Decrease_change_across_replicates & ~(Increase_change_across_replicates | No_change_across_replicates | Inconsistent_across_replicates))),...
  sum(sum(No_change_across_replicates)),...
  sum(sum(Inconsistent_across_replicates & ~(Increase_change_across_replicates | Decrease_change_across_replicates | No_change_across_replicates))));
fclose(fid_summary_master_comparsion1);

fid_summary_master_comparsion1= fopen([datadir1 'Summary_changes_based_master_gaussian_list_across_replicate_1.csv'],'wt'); % create the output file with the header infromation
fprintf (fid_summary_master_comparsion1,'%s,%s,%s,%s,%s,%s\n',... %header for OutputGaus output
  'Increase in replicate 1', 'Decrease in replicate 1', 'No change in replicate 1',...
  'Unquantifiable observed only in MvsL channel replicate 1',...
  'Unquantifiable observed only in HvsL channel replicate 1',...
  'Unquantifiable not observed in replicate 1'); %Write Header
fprintf (fid_summary_master_comparsion1,'%6.4f,%6.4f,%6.4f,%6.4f,%6.4f,%6.4f\n',...
  Increase_counter_replicate_1, Decrease_counter_replicate_1,...
  No_change_counter_replicate_1,Unquantifiable_observed_only_in_MvsL_channel_replicate_1,...
  Unquantifiable_observed_only_in_HvsL_channel_replicate_1,...
  Unquantifiable_Not_observed_in_replicate_replicate_1);
fclose(fid_summary_master_comparsion1);

fid_summary_master_comparsion2= fopen([datadir1 'Summary_changes_based_master_gaussian_list_across_replicate_2.csv'],'wt'); % create the output file with the header infromation
fprintf (fid_summary_master_comparsion2,'%s,%s,%s,%s,%s,%s\n',... %header for OutputGaus output
  'Increase in replicate 2', 'Decrease in replicate 2', 'No change in replicate 2',...
  'Unquantifiable observed only in MvsL channel replicate 2',...
  'Unquantifiable observed only in HvsL channel replicate 2',...
  'Unquantifiable not observed in replicate 2'); %Write Header
fprintf (fid_summary_master_comparsion2,'%6.4f,%6.4f,%6.4f,%6.4f,%6.4f,%6.4f\n',...
  Increase_counter_replicate_2,Decrease_counter_replicate_2,...
  No_change_counter_replicate_2,Unquantifiable_observed_only_in_MvsL_channel_replicate_2,...
  Unquantifiable_observed_only_in_HvsL_channel_replicate_2,...
  Unquantifiable_Not_observed_in_replicate_replicate_2);
fclose(fid_summary_master_comparsion2);

fid_summary_master_comparsion3= fopen([datadir1 'Summary_changes_based_master_gaussian_list_across_replicate_3.csv'],'wt'); % create the output file with the header infromation
fprintf (fid_summary_master_comparsion3,'%s,%s,%s,%s,%s,%s\n',... %header for OutputGaus output
  'Increase in replicate 3', 'Decrease in replicate 3', 'No change in replicate 3',...
  'Unquantifiable observed only in MvsL channel replicate 3',...
  'Unquantifiable observed only in HvsL channel replicate 3',...
  'Unquantifiable not observed in replicate 3'); %Write Header
fprintf (fid_summary_master_comparsion3,'%6.4f,%6.4f,%6.4f,%6.4f,%6.4f,%6.4f\n',...
  Increase_counter_replicate_3,Decrease_counter_replicate_3,...
  No_change_counter_replicate_3,Unquantifiable_observed_only_in_MvsL_channel_replicate_3,...
  Unquantifiable_observed_only_in_HvsL_channel_replicate_3,...
  Unquantifiable_Not_observed_in_replicate_replicate_3);
fclose(fid_summary_master_comparsion3);

fid_summary_master_comparsion4= fopen([datadir1 'Summary_changes_based_master_gaussian_list_across_replicate_4.csv'],'wt'); % create the output file with the header infromation
fprintf (fid_summary_master_comparsion4,'%s,%s,%s,%s,%s,%s\n',... %header for OutputGaus output
  'Increase in replicate 4', 'Decrease in replicate 4', 'No change in replicate 3',...
  'Unquantifiable observed only in MvsL channel replicate 4',...
  'Unquantifiable observed only in HvsL channel replicate 4',...
  'Unquantifiable not observed in replicate 4'); %Write Header
fprintf (fid_summary_master_comparsion4,'%6.4f,%6.4f,%6.4f,%6.4f,%6.4f,%6.4f\n',...
  Increase_counter_replicate_4,Decrease_counter_replicate_4,...
  No_change_counter_replicate_4,Unquantifiable_observed_only_in_MvsL_channel_replicate_4,...
  Unquantifiable_observed_only_in_HvsL_channel_replicate_4,...
  Unquantifiable_Not_observed_in_replicate_replicate_4);
fclose(fid_summary_master_comparsion4);

fid_summary_master_coverage= fopen([datadir1 'Summary_observed_protein_PCP_SILAC_coverage_within_experiments.csv'],'wt'); % create the output file with the header infromation
fprintf (fid_summary_master_coverage,'%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,\n',... %header for OutputGaus output
  'Replicate number','Isotopologue Channel',...
  'Average PCP_SILAC coverage','Average coverage of the SEC', 'Median PCP_SILAC coverage',...
  'Number of proteins with a coverage less then 0.5 ',...
  'Number of proteins with a coverage less then 0.7',...
  'Number of proteins with a coverage less then 0.85',...
  'Number of proteins with a coverage less then 1',...
  'Number of proteins with a coverage less then 1.5',...
  'Number of proteins with a coverage less then 2',...
  'Number of proteins with a coverage greater then 2',...
  'Minimum coverage observed',...
  'Maximum coverage observed'); %Write Header
%Replicate 1 information
%MvsL
fprintf (fid_summary_master_coverage,'%s,%s,%6.4f,%6.4f,%6.4f,%6.4f,%6.4f,%6.4f,%6.4f,%6.4f,%6.4f,%6.4f,%6.4f,%6.4f,\n',...
  'Replicate 1','MvsL',...
  mean(Finalised_Master_Gaussian_list.Summed_Gaussian_areaMvsL(:,1)),...
  mean(Finalised_Master_Gaussian_list.SEC_coverageMvsL(:,1)),...
  median(Finalised_Master_Gaussian_list.Summed_Gaussian_areaMvsL(:,1)),...
  sum((Finalised_Master_Gaussian_list.Summed_Gaussian_areaMvsL(:,1)<=0.5)&(Finalised_Master_Gaussian_list.Summed_Gaussian_areaMvsL(:,1)>=0)),...
  sum((Finalised_Master_Gaussian_list.Summed_Gaussian_areaMvsL(:,1)<=0.7)&(Finalised_Master_Gaussian_list.Summed_Gaussian_areaMvsL(:,1)>=0.5)),...
  sum((Finalised_Master_Gaussian_list.Summed_Gaussian_areaMvsL(:,1)<=0.85)&(Finalised_Master_Gaussian_list.Summed_Gaussian_areaMvsL(:,1)>=0.7)),...
  sum((Finalised_Master_Gaussian_list.Summed_Gaussian_areaMvsL(:,1)<=1)&(Finalised_Master_Gaussian_list.Summed_Gaussian_areaMvsL(:,1)>=0.85)),...
  sum((Finalised_Master_Gaussian_list.Summed_Gaussian_areaMvsL(:,1)<=1.5)&(Finalised_Master_Gaussian_list.Summed_Gaussian_areaMvsL(:,1)>=1)),...
  sum((Finalised_Master_Gaussian_list.Summed_Gaussian_areaMvsL(:,1)<=2)&(Finalised_Master_Gaussian_list.Summed_Gaussian_areaMvsL(:,1)>=1.5)),...
  sum((Finalised_Master_Gaussian_list.Summed_Gaussian_areaMvsL(:,1)>=2)),...
  min(Finalised_Master_Gaussian_list.Summed_Gaussian_areaMvsL(:,1)),...
  max(Finalised_Master_Gaussian_list.Summed_Gaussian_areaMvsL(:,1)));
%HvsL
fprintf (fid_summary_master_coverage,'%s,%s,%6.4f,%6.4f,%6.4f,%6.4f,%6.4f,%6.4f,%6.4f,%6.4f,%6.4f,%6.4f,%6.4f,%6.4f,\n',...
  'Replicate 1','HvsL',...
  mean(Finalised_Master_Gaussian_list.Summed_Gaussian_areaHvsL(:,1)),...
  mean(Finalised_Master_Gaussian_list.SEC_coverageHvsL(:,1)),...
  median(Finalised_Master_Gaussian_list.Summed_Gaussian_areaHvsL(:,1)),...
  sum((Finalised_Master_Gaussian_list.Summed_Gaussian_areaHvsL(:,1)<=0.5)&(Finalised_Master_Gaussian_list.Summed_Gaussian_areaHvsL(:,1)>=0)),...
  sum((Finalised_Master_Gaussian_list.Summed_Gaussian_areaHvsL(:,1)<=0.7)&(Finalised_Master_Gaussian_list.Summed_Gaussian_areaHvsL(:,1)>=0.5)),...
  sum((Finalised_Master_Gaussian_list.Summed_Gaussian_areaHvsL(:,1)<=0.85)&(Finalised_Master_Gaussian_list.Summed_Gaussian_areaHvsL(:,1)>=0.7)),...
  sum((Finalised_Master_Gaussian_list.Summed_Gaussian_areaHvsL(:,1)<=1)&(Finalised_Master_Gaussian_list.Summed_Gaussian_areaHvsL(:,1)>=0.85)),...
  sum((Finalised_Master_Gaussian_list.Summed_Gaussian_areaHvsL(:,1)<=1.5)&(Finalised_Master_Gaussian_list.Summed_Gaussian_areaHvsL(:,1)>=1)),...
  sum((Finalised_Master_Gaussian_list.Summed_Gaussian_areaHvsL(:,1)<=2)&(Finalised_Master_Gaussian_list.Summed_Gaussian_areaHvsL(:,1)>=1.5)),...
  sum((Finalised_Master_Gaussian_list.Summed_Gaussian_areaHvsL(:,1)>=2)),...
  min(Finalised_Master_Gaussian_list.Summed_Gaussian_areaHvsL(:,1)),...
  max(Finalised_Master_Gaussian_list.Summed_Gaussian_areaHvsL(:,1)));

%Replicate 2 information
%MvsL
fprintf (fid_summary_master_coverage,'%s,%s,%6.4f,%6.4f,%6.4f,%6.4f,%6.4f,%6.4f,%6.4f,%6.4f,%6.4f,%6.4f,%6.4f,%6.4f,\n',...
  'Replicate 2','MvsL',...
  mean(Finalised_Master_Gaussian_list.Summed_Gaussian_areaMvsL(:,2)),...
  mean(Finalised_Master_Gaussian_list.SEC_coverageMvsL(:,2)),...
  median(Finalised_Master_Gaussian_list.Summed_Gaussian_areaMvsL(:,2)),...
  sum((Finalised_Master_Gaussian_list.Summed_Gaussian_areaMvsL(:,2)<=0.5)&(Finalised_Master_Gaussian_list.Summed_Gaussian_areaMvsL(:,2)>=0)),...
  sum((Finalised_Master_Gaussian_list.Summed_Gaussian_areaMvsL(:,2)<=0.7)&(Finalised_Master_Gaussian_list.Summed_Gaussian_areaMvsL(:,2)>=0.5)),...
  sum((Finalised_Master_Gaussian_list.Summed_Gaussian_areaMvsL(:,2)<=0.85)&(Finalised_Master_Gaussian_list.Summed_Gaussian_areaMvsL(:,2)>=0.7)),...
  sum((Finalised_Master_Gaussian_list.Summed_Gaussian_areaMvsL(:,2)<=1)&(Finalised_Master_Gaussian_list.Summed_Gaussian_areaMvsL(:,2)>=0.85)),...
  sum((Finalised_Master_Gaussian_list.Summed_Gaussian_areaMvsL(:,2)<=1.5)&(Finalised_Master_Gaussian_list.Summed_Gaussian_areaMvsL(:,2)>=1)),...
  sum((Finalised_Master_Gaussian_list.Summed_Gaussian_areaMvsL(:,2)<=2)&(Finalised_Master_Gaussian_list.Summed_Gaussian_areaMvsL(:,2)>=1.5)),...
  sum((Finalised_Master_Gaussian_list.Summed_Gaussian_areaMvsL(:,2)>=2)),...
  min(Finalised_Master_Gaussian_list.Summed_Gaussian_areaMvsL(:,2)),...
  max(Finalised_Master_Gaussian_list.Summed_Gaussian_areaMvsL(:,2)));
%HvsL
fprintf (fid_summary_master_coverage,'%s,%s,%6.4f,%6.4f,%6.4f,%6.4f,%6.4f,%6.4f,%6.4f,%6.4f,%6.4f,%6.4f,%6.4f,%6.4f,\n',...
  'Replicate 2','HvsL',...
  mean(Finalised_Master_Gaussian_list.Summed_Gaussian_areaHvsL(:,2)),...
  mean(Finalised_Master_Gaussian_list.SEC_coverageHvsL(:,2)),...
  median(Finalised_Master_Gaussian_list.Summed_Gaussian_areaHvsL(:,2)),...
  sum((Finalised_Master_Gaussian_list.Summed_Gaussian_areaHvsL(:,2)<=0.5)&(Finalised_Master_Gaussian_list.Summed_Gaussian_areaHvsL(:,2)>=0)),...
  sum((Finalised_Master_Gaussian_list.Summed_Gaussian_areaHvsL(:,2)<=0.7)&(Finalised_Master_Gaussian_list.Summed_Gaussian_areaHvsL(:,2)>=0.5)),...
  sum((Finalised_Master_Gaussian_list.Summed_Gaussian_areaHvsL(:,2)<=0.85)&(Finalised_Master_Gaussian_list.Summed_Gaussian_areaHvsL(:,2)>=0.7)),...
  sum((Finalised_Master_Gaussian_list.Summed_Gaussian_areaHvsL(:,2)<=1)&(Finalised_Master_Gaussian_list.Summed_Gaussian_areaHvsL(:,2)>=0.85)),...
  sum((Finalised_Master_Gaussian_list.Summed_Gaussian_areaHvsL(:,2)<=1.5)&(Finalised_Master_Gaussian_list.Summed_Gaussian_areaHvsL(:,2)>=1)),...
  sum((Finalised_Master_Gaussian_list.Summed_Gaussian_areaHvsL(:,2)<=2)&(Finalised_Master_Gaussian_list.Summed_Gaussian_areaHvsL(:,2)>=1.5)),...
  sum((Finalised_Master_Gaussian_list.Summed_Gaussian_areaHvsL(:,2)>=2)),...
  min(Finalised_Master_Gaussian_list.Summed_Gaussian_areaHvsL(:,2)),...
  max(Finalised_Master_Gaussian_list.Summed_Gaussian_areaHvsL(:,2)));


%Replicate 3 information
%MvsL
fprintf (fid_summary_master_coverage,'%s,%s,%6.4f,%6.4f,%6.4f,%6.4f,%6.4f,%6.4f,%6.4f,%6.4f,%6.4f,%6.4f,%6.4f,%6.4f,\n',...
  'Replicate 3','MvsL',...
  mean(Finalised_Master_Gaussian_list.Summed_Gaussian_areaMvsL(:,3)),...
  mean(Finalised_Master_Gaussian_list.SEC_coverageMvsL(:,3)),...
  median(Finalised_Master_Gaussian_list.Summed_Gaussian_areaMvsL(:,3)),...
  sum((Finalised_Master_Gaussian_list.Summed_Gaussian_areaMvsL(:,3)<=0.5)&(Finalised_Master_Gaussian_list.Summed_Gaussian_areaMvsL(:,3)>=0)),...
  sum((Finalised_Master_Gaussian_list.Summed_Gaussian_areaMvsL(:,3)<=0.7)&(Finalised_Master_Gaussian_list.Summed_Gaussian_areaMvsL(:,3)>=0.5)),...
  sum((Finalised_Master_Gaussian_list.Summed_Gaussian_areaMvsL(:,3)<=0.85)&(Finalised_Master_Gaussian_list.Summed_Gaussian_areaMvsL(:,3)>=0.7)),...
  sum((Finalised_Master_Gaussian_list.Summed_Gaussian_areaMvsL(:,3)<=1)&(Finalised_Master_Gaussian_list.Summed_Gaussian_areaMvsL(:,3)>=0.85)),...
  sum((Finalised_Master_Gaussian_list.Summed_Gaussian_areaMvsL(:,3)<=1.5)&(Finalised_Master_Gaussian_list.Summed_Gaussian_areaMvsL(:,3)>=1)),...
  sum((Finalised_Master_Gaussian_list.Summed_Gaussian_areaMvsL(:,3)<=2)&(Finalised_Master_Gaussian_list.Summed_Gaussian_areaMvsL(:,3)>=1.5)),...
  sum((Finalised_Master_Gaussian_list.Summed_Gaussian_areaMvsL(:,3)>=2)),...
  min(Finalised_Master_Gaussian_list.Summed_Gaussian_areaMvsL(:,3)),...
  max(Finalised_Master_Gaussian_list.Summed_Gaussian_areaMvsL(:,3)));
%HvsL
fprintf (fid_summary_master_coverage,'%s,%s,%6.4f,%6.4f,%6.4f,%6.4f,%6.4f,%6.4f,%6.4f,%6.4f,%6.4f,%6.4f,%6.4f,%6.4f,\n',...
  'Replicate 3','HvsL',...
  mean(Finalised_Master_Gaussian_list.Summed_Gaussian_areaHvsL(:,3)),...
  mean(Finalised_Master_Gaussian_list.SEC_coverageHvsL(:,3)),...
  median(Finalised_Master_Gaussian_list.Summed_Gaussian_areaHvsL(:,3)),...
  sum((Finalised_Master_Gaussian_list.Summed_Gaussian_areaHvsL(:,3)<=0.5)&(Finalised_Master_Gaussian_list.Summed_Gaussian_areaHvsL(:,3)>=0)),...
  sum((Finalised_Master_Gaussian_list.Summed_Gaussian_areaHvsL(:,3)<=0.7)&(Finalised_Master_Gaussian_list.Summed_Gaussian_areaHvsL(:,3)>=0.5)),...
  sum((Finalised_Master_Gaussian_list.Summed_Gaussian_areaHvsL(:,3)<=0.85)&(Finalised_Master_Gaussian_list.Summed_Gaussian_areaHvsL(:,3)>=0.7)),...
  sum((Finalised_Master_Gaussian_list.Summed_Gaussian_areaHvsL(:,3)<=1)&(Finalised_Master_Gaussian_list.Summed_Gaussian_areaHvsL(:,3)>=0.85)),...
  sum((Finalised_Master_Gaussian_list.Summed_Gaussian_areaHvsL(:,3)<=1.5)&(Finalised_Master_Gaussian_list.Summed_Gaussian_areaHvsL(:,3)>=1)),...
  sum((Finalised_Master_Gaussian_list.Summed_Gaussian_areaHvsL(:,3)<=2)&(Finalised_Master_Gaussian_list.Summed_Gaussian_areaHvsL(:,3)>=1.5)),...
  sum((Finalised_Master_Gaussian_list.Summed_Gaussian_areaHvsL(:,3)>=2)),...
  min(Finalised_Master_Gaussian_list.Summed_Gaussian_areaHvsL(:,3)),...
  max(Finalised_Master_Gaussian_list.Summed_Gaussian_areaHvsL(:,3)));

%Replicate 4 information
%MvsL
fprintf (fid_summary_master_coverage,'%s,%s,%6.4f,%6.4f,%6.4f,%6.4f,%6.4f,%6.4f,%6.4f,%6.4f,%6.4f,%6.4f,%6.4f,%6.4f,\n',...
  'Replicate 4','MvsL',...
  mean(Finalised_Master_Gaussian_list.Summed_Gaussian_areaMvsL(:,4)),...
  mean(Finalised_Master_Gaussian_list.SEC_coverageMvsL(:,4)),...
  median(Finalised_Master_Gaussian_list.Summed_Gaussian_areaMvsL(:,4)),...
  sum((Finalised_Master_Gaussian_list.Summed_Gaussian_areaMvsL(:,4)<=0.5)&(Finalised_Master_Gaussian_list.Summed_Gaussian_areaMvsL(:,4)>=0)),...
  sum((Finalised_Master_Gaussian_list.Summed_Gaussian_areaMvsL(:,4)<=0.7)&(Finalised_Master_Gaussian_list.Summed_Gaussian_areaMvsL(:,4)>=0.5)),...
  sum((Finalised_Master_Gaussian_list.Summed_Gaussian_areaMvsL(:,4)<=0.85)&(Finalised_Master_Gaussian_list.Summed_Gaussian_areaMvsL(:,4)>=0.7)),...
  sum((Finalised_Master_Gaussian_list.Summed_Gaussian_areaMvsL(:,4)<=1)&(Finalised_Master_Gaussian_list.Summed_Gaussian_areaMvsL(:,4)>=0.85)),...
  sum((Finalised_Master_Gaussian_list.Summed_Gaussian_areaMvsL(:,4)<=1.5)&(Finalised_Master_Gaussian_list.Summed_Gaussian_areaMvsL(:,4)>=1)),...
  sum((Finalised_Master_Gaussian_list.Summed_Gaussian_areaMvsL(:,4)<=2)&(Finalised_Master_Gaussian_list.Summed_Gaussian_areaMvsL(:,4)>=1.5)),...
  sum((Finalised_Master_Gaussian_list.Summed_Gaussian_areaMvsL(:,4)>=2)),...
  min(Finalised_Master_Gaussian_list.Summed_Gaussian_areaMvsL(:,4)),...
  max(Finalised_Master_Gaussian_list.Summed_Gaussian_areaMvsL(:,4)));
%HvsL
fprintf (fid_summary_master_coverage,'%s,%s,%6.4f,%6.4f,%6.4f,%6.4f,%6.4f,%6.4f,%6.4f,%6.4f,%6.4f,%6.4f,%6.4f,%6.4f,\n',...
  'Replicate 4','HvsL',...
  mean(Finalised_Master_Gaussian_list.Summed_Gaussian_areaHvsL(:,4)),...
  mean(Finalised_Master_Gaussian_list.SEC_coverageHvsL(:,4)),...
  median(Finalised_Master_Gaussian_list.Summed_Gaussian_areaHvsL(:,4)),...
  sum((Finalised_Master_Gaussian_list.Summed_Gaussian_areaHvsL(:,4)<=0.5)&(Finalised_Master_Gaussian_list.Summed_Gaussian_areaHvsL(:,4)>=0)),...
  sum((Finalised_Master_Gaussian_list.Summed_Gaussian_areaHvsL(:,4)<=0.7)&(Finalised_Master_Gaussian_list.Summed_Gaussian_areaHvsL(:,4)>=0.5)),...
  sum((Finalised_Master_Gaussian_list.Summed_Gaussian_areaHvsL(:,4)<=0.85)&(Finalised_Master_Gaussian_list.Summed_Gaussian_areaHvsL(:,4)>=0.7)),...
  sum((Finalised_Master_Gaussian_list.Summed_Gaussian_areaHvsL(:,4)<=1)&(Finalised_Master_Gaussian_list.Summed_Gaussian_areaHvsL(:,4)>=0.85)),...
  sum((Finalised_Master_Gaussian_list.Summed_Gaussian_areaHvsL(:,4)<=1.5)&(Finalised_Master_Gaussian_list.Summed_Gaussian_areaHvsL(:,4)>=1)),...
  sum((Finalised_Master_Gaussian_list.Summed_Gaussian_areaHvsL(:,4)<=2)&(Finalised_Master_Gaussian_list.Summed_Gaussian_areaHvsL(:,4)>=1.5)),...
  sum((Finalised_Master_Gaussian_list.Summed_Gaussian_areaHvsL(:,4)>=2)),...
  min(Finalised_Master_Gaussian_list.Summed_Gaussian_areaHvsL(:,4)),...
  max(Finalised_Master_Gaussian_list.Summed_Gaussian_areaHvsL(:,4)));
%Close file
fclose(fid_summary_master_coverage);

fid_summary_protein_coverage= fopen([datadir1 'Summary_PCP_SILAC_coverage_Individual_proteins_within_experiments.csv'],'wt'); % create the output file with the header infromation
fprintf (fid_summary_protein_coverage,'%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,\n',... %header for OutputGaus output
  'Protein name','Summed isotopologue area MvsL replicate 1',...
  'Summed isotopologue area MvsL replicate 2',...
  'Summed isotopologue area MvsL replicate 3',...
  'Summed isotopologue area MvsL replicate 4',...
  'Summed isotopologue area HvsL replicate 1',...
  'Summed isotopologue area HvsL replicate 2',...
  'Summed isotopologue area HvsL replicate 3',...
  'Summed isotopologue area HvsL replicate 4',...
  'Summed isotopologue HvsM ratio replicate 1',...
  'Summed isotopologue HvsM ratio replicate 2',...
  'Summed isotopologue HvsM ratio replicate 3',...
  'Summed isotopologue HvsM ratio replicate 4',...
  'Number of Gaussians');
for writeout_counter1= 1:Dimension_of_master_gaussian_list(1)
  fprintf (fid_summary_protein_coverage,'%s,%6.4f,%6.4f,%6.4f,%6.4f,%6.4f,%6.4f,%6.4f,%6.4f,%6.4f,%6.4f,%6.4f,%6.4f,%6f,\n',...
    Finalised_Master_Gaussian_list.Protein_name{writeout_counter1},...
    Finalised_Master_Gaussian_list.Summed_Gaussian_areaMvsL(writeout_counter1,1),...
    Finalised_Master_Gaussian_list.Summed_Gaussian_areaMvsL(writeout_counter1,2),...
    Finalised_Master_Gaussian_list.Summed_Gaussian_areaMvsL(writeout_counter1,3),...
    Finalised_Master_Gaussian_list.Summed_Gaussian_areaMvsL(writeout_counter1,4),...
    Finalised_Master_Gaussian_list.Summed_Gaussian_areaHvsL(writeout_counter1,1),...
    Finalised_Master_Gaussian_list.Summed_Gaussian_areaHvsL(writeout_counter1,2),...
    Finalised_Master_Gaussian_list.Summed_Gaussian_areaHvsL(writeout_counter1,3),...
    Finalised_Master_Gaussian_list.Summed_Gaussian_areaHvsL(writeout_counter1,4),...
    Finalised_Master_Gaussian_list.Ratio_Gaussian_area(writeout_counter1,1),...
    Finalised_Master_Gaussian_list.Ratio_Gaussian_area(writeout_counter1,2),...
    Finalised_Master_Gaussian_list.Ratio_Gaussian_area(writeout_counter1,3),...
    Finalised_Master_Gaussian_list.Ratio_Gaussian_area(writeout_counter1,4),...
    nnz(Finalised_Master_Gaussian_list.Center(writeout_counter1,:)));
end
fclose(fid_summary_protein_coverage);


fid_Quantation_numbers= fopen([datadir1 'Quantation_values_identified.csv'],'wt'); % create the output file with the header infromation
fprintf (fid_Quantation_numbers,'%s,%s,%s,%s,%s\n',... %header for OutputGaus output
  'Total number of values across replicates','Replicate 1',...
  'Replicate 2','Replicate 3', 'Replicate 4'); %Write Header
fprintf (fid_Quantation_numbers,'%6.4f,%6.4f,%6.4f,%6.4f,%6.4f\n',...
  nnz(Finalised_Master_Gaussian_list.Center(:,:)),...
  replicate_1_counter,replicate_2_counter,...
  replicate_3_counter,replicate_4_counter);
fclose(fid_Quantation_numbers);

tt = toc;
fprintf('  ...  %.2f seconds\n',tt)



%% 12. Write out txt for Go enrichment
fprintf('    12. Write out txt for Go enrichment')

%Aim: Write out file compatible with Perseus to rapidly allow the go
%enrichment
%Guassian level
fid_Go_enrichment= fopen([datadir1 'Perseus_enrichment_Gaussian_level_file.csv'],'wt'); % create the output file with the header infromation
fprintf (fid_Go_enrichment,'%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,\n',... %header for OutputGaus output
  'Center','Protein name','Normalised Average fold change', 'Std dev fold change',...
  'Normalised fold change replicate 1','Normalised fold change replicate 2',...
  'Normalised fold change replicate 3','Normalised fold change replicate 3',...
  'Increase (>= two replicates)', 'Decrease (>= two replicates)',...
  'p-value (t-test)','Below Adjusted_pvalue (bonferroni correction) (ttest)',...
  'p-value (MWW)','Below Adjusted_pvalue(bonferroni correction) (MWW)'); %Write Header
for writeout_counter1= 1:Dimension_of_master_gaussian_list(1)
  for writeout_counter2= 1:nnz(Finalised_Master_Gaussian_list.Center(writeout_counter1,:))
    fprintf(fid_Go_enrichment,'%6.3f,%s,%6.3f,%6.3f,%6.3f,%6.3f,%6.3f,%6.3f,%s,%s,%6.9f,%s,%6.9f,%s,\n',...
      Finalised_Master_Gaussian_list.Center(writeout_counter1,writeout_counter2),...
      Finalised_Master_Gaussian_list.Protein_name{writeout_counter1},...
      Finalised_Master_Gaussian_list.Normalised_average_fold_change(writeout_counter1,writeout_counter2),...
      Finalised_Master_Gaussian_list.Stdev_fold_change(writeout_counter1,writeout_counter2),...
      Finalised_Master_Gaussian_list.normaliased_log2_comparsion_replicate_1(writeout_counter1,writeout_counter2),...
      Finalised_Master_Gaussian_list.normaliased_log2_comparsion_replicate_2(writeout_counter1,writeout_counter2),...
      Finalised_Master_Gaussian_list.normaliased_log2_comparsion_replicate_3(writeout_counter1,writeout_counter2),...
      Finalised_Master_Gaussian_list.normaliased_log2_comparsion_replicate_4(writeout_counter1,writeout_counter2),...
      Finalised_Master_Gaussian_list.Increase_change_across_replicates_persaus_input{writeout_counter1,writeout_counter2},...
      Finalised_Master_Gaussian_list.Decrease_change_across_replicates_persaus_input{writeout_counter1,writeout_counter2},...
      Finalised_Master_Gaussian_list.Ttest_p_values(writeout_counter1,writeout_counter2),...
      Finalised_Master_Gaussian_list.Satifies_Adjusted_pvalue{writeout_counter1,writeout_counter2},...
      Finalised_Master_Gaussian_list.MWW_p_values(writeout_counter1,writeout_counter2),...
      Finalised_Master_Gaussian_list.MWW_Satifies_Adjusted_pvalue{writeout_counter1,writeout_counter2});
  end
end
fclose(fid_Go_enrichment);


%protein level
fid_Go_enrichment= fopen([datadir1 'Perseus_enrichment_Protein_level_file.csv'],'wt'); % create the output file with the header infromation
fprintf (fid_Go_enrichment,'%s,%s,%s,%s,%s,%s,\n',... %header for OutputGaus output
  'Protein name','Increase AND Decreases (>= two replicates)',...
  'Increase (>= two replicates)', 'Decrease (>= two replicates)',...
  'Gaussian fitted with change below adjusted_pvalue(bonferroni correction) (ttest)',...
  'Gaussian fitted with change below adjusted_pvalue(bonferroni correction) (MWW)');
for writeout_counter1= 1:Dimension_of_master_gaussian_list(1)
  fprintf(fid_Go_enrichment,'%s,%s,%s,%s,%s,%s,\n',...
    Finalised_Master_Gaussian_list.Protein_name{writeout_counter1},...
    Finalised_Master_Gaussian_list.Increase_and_Decrease_protein_persaus_input{writeout_counter1},...
    Finalised_Master_Gaussian_list.Increase_protein_persaus_input{writeout_counter1},...
    Finalised_Master_Gaussian_list.Decrease_protein_persaus_input{writeout_counter1},...
    Singficant_change_ttest_persaus_protein{writeout_counter1},...
    Singficant_change_mww_persaus_protein{writeout_counter1});
end
fclose(fid_Go_enrichment);

tt = toc;
tt = round(tt/10)*10;
fprintf('  ...  %.2f seconds\n',tt)



%% 13. Make figures
fprintf('    12. Make figures')


%Figure showing changes in response to treatment

%Graph data as log2 scatter
f5=figure;
log2_sorted=sort(normaliased_log_2_fold_comparsion_of_gaussian(:));
f5=subplot(2,1,1)
P4 =scatter([1:Total_number_of_unique_gaussians],log2_sorted(:),8,'fill');
title('Log2 changes between replicate: BN-PAGE Mitochondira Fas treated','FontSize', 10);
ylabel('Normalised Log2 ratio','FontSize', 10);
xlabel('Gaussians identified','FontSize', 10);
%Set limits based on observed values
%ylim([log2_sorted(1)*1.05,log2_sorted(Total_number_of_unique_gaussians)*1.05]);
%Set limits based on user defined cut off
ylim([-6,6]);
xlim([-100,Total_number_of_unique_gaussians+100]);
hold on
hold on %add 1 dotted line
PP1= plot([-100 Total_number_of_unique_gaussians+100],[1 1],':');
hold on %add 1 dotted line
PP1= plot([-100 Total_number_of_unique_gaussians+100],[-1 -1],':');

%plot distribution of protein chnages across the fractions
%define colours
subplot(2,1,2);
f6_figure=bar(1:fraction_number(2), [Hist_array(:,1) Hist_array(:,2) Hist_array(:,3)],0.6, 'stack');
for k=1:3
  set(f6_figure(k),'facecolor', myC(k,:), 'EdgeColor', 'k' )
end
legend('No change', 'Increase', 'Decrease', 'FontSize',8, 'Location', 'Best');
xlim([-2,fraction_to_plot+2]);
title('Changes in Gaussians observed across fractions');
xlabel('Fractions');
ylabel('# of Guassian centers');
saveas(f5, [figdir1 'All_gaussians_changes_observed.png']);



for figure_counter=1:replicate_num
  
  %Determine Threshold to use for determine signficant change, use 1.96*stdev
  for writeout_counter2= 1:Dimension_of_master_gaussian_list(1)
    percentage_derivation(writeout_counter2)= log2(Finalised_Master_Gaussian_list.Summed_Gaussian_areaMvsL(writeout_counter2,figure_counter)/...
      Finalised_Master_Gaussian_list.Summed_Gaussian_areaHvsL(writeout_counter2,figure_counter));
    if isinf(percentage_derivation(writeout_counter2))==1
      percentage_derivation(writeout_counter2)=NaN;
    end
  end
  
  %Remove nan values
  percentage_derivation(isnan(percentage_derivation))=[];
  
  %calculate stdev
  std_of_area_measurment=std(percentage_derivation);
  if isnan(std_of_area_measurment)
    threshold_of_change=max(abs(percentage_derivation));
  else
    threshold_of_change=1.96*std_of_area_measurment;
  end
  
  %Figure showing coverage within PCP-SEC-SILAC experiments
  %%Count how many proteins change across the observed SEC fraction
  bin_number_hist_MvsL= 0.05:0.05:ceil(max(max(Finalised_Master_Gaussian_list.Summed_Gaussian_areaMvsL)));
  bin_number_hist_MvsL_plus1= 0.05:0.05:(ceil(max(max(Finalised_Master_Gaussian_list.Summed_Gaussian_areaMvsL)))+0.05);
  bin_number_hist2_MvsL=length(bin_number_hist_MvsL);
  
  Hist_array3A=zeros(Dimension_of_master_gaussian_list(1),(bin_number_hist2_MvsL+1));
  Hist_array3B=zeros(Dimension_of_master_gaussian_list(1),(bin_number_hist2_MvsL+1));
  Hist_array3C=zeros(Dimension_of_master_gaussian_list(1),(bin_number_hist2_MvsL+1));
  Hist_array3D=zeros(Dimension_of_master_gaussian_list(1),(bin_number_hist2_MvsL+1));
  
  %Count the number of gaussian detected vs the area
  for hist2_counter1= 1:(bin_number_hist2_MvsL)+1
    for writeout_counter2= 1:Dimension_of_master_gaussian_list(1)
      %determine bin values
      Coverage_Bin_number=ceil(Finalised_Master_Gaussian_list.Summed_Gaussian_areaMvsL(writeout_counter2,figure_counter)/0.05);
      %determine if the fold change is within the bin value
      if  hist2_counter1-1==Coverage_Bin_number
        HvsL_values=Finalised_Master_Gaussian_list.Summed_Gaussian_areaHvsL(writeout_counter2,figure_counter);
        MvsL_values=Finalised_Master_Gaussian_list.Summed_Gaussian_areaMvsL(writeout_counter2,figure_counter);
        %determine bin of the total observed isotopologue 'volumn' in MvsL
        Hist_array3A(writeout_counter2,hist2_counter1)=...
          Hist_array3A(writeout_counter2,hist2_counter1) + 1;
        %determine if the 'volumn' is different from HvsL
        if log2(MvsL_values/HvsL_values)<=threshold_of_change &...
            log2(MvsL_values/HvsL_values)>=-threshold_of_change
          Hist_array3B(writeout_counter2,hist2_counter1)=...
            Hist_array3B(writeout_counter2,hist2_counter1) + 1;
        elseif log2(MvsL_values/HvsL_values)>threshold_of_change
          Hist_array3C(writeout_counter2,hist2_counter1)=...
            Hist_array3C(writeout_counter2,hist2_counter1) + 1;
        elseif log2(MvsL_values/HvsL_values)<-threshold_of_change
          Hist_array3D(writeout_counter2,hist2_counter1)=...
            Hist_array3D(writeout_counter2,hist2_counter1) + 1;
        end
      end
    end
  end
  
  %Figure showing coverage within PCP-SEC-SILAC experiments
  %%Count how many proteins change across the observed SEC fraction
  bin_number_hist_HvsL= 0.05:0.05:ceil(max(max(Finalised_Master_Gaussian_list.Summed_Gaussian_areaHvsL)));
  bin_number_hist_HvsL_plus1= 0.05:0.05:(ceil(max(max(Finalised_Master_Gaussian_list.Summed_Gaussian_areaHvsL)))+0.05);
  bin_number_hist2_HvsL=length(bin_number_hist_HvsL);
  
  Hist_array4A=zeros(Dimension_of_master_gaussian_list(1),(bin_number_hist2_HvsL+1));
  Hist_array4B=zeros(Dimension_of_master_gaussian_list(1),(bin_number_hist2_HvsL+1));
  Hist_array4C=zeros(Dimension_of_master_gaussian_list(1),(bin_number_hist2_HvsL+1));
  Hist_array4D=zeros(Dimension_of_master_gaussian_list(1),(bin_number_hist2_HvsL+1));
  
  %Count the number of gaussian detected vs the area
  for hist2_counter1= 1:(bin_number_hist2_HvsL)+1
    for writeout_counter2= 1:Dimension_of_master_gaussian_list(1)
      %determine bin values
      Coverage_Bin_number=ceil(Finalised_Master_Gaussian_list.Summed_Gaussian_areaHvsL(writeout_counter2,figure_counter)/0.05);
      %determine if the fold change is within the bin value
      if  hist2_counter1-1==Coverage_Bin_number
        HvsL_values=Finalised_Master_Gaussian_list.Summed_Gaussian_areaHvsL(writeout_counter2,figure_counter);
        MvsL_values=Finalised_Master_Gaussian_list.Summed_Gaussian_areaMvsL(writeout_counter2,figure_counter);
        %determine bin of the total observed isotopologue 'volumn' in MvsL
        Hist_array4A(writeout_counter2,hist2_counter1)=...
          Hist_array4A(writeout_counter2,hist2_counter1) + 1;
        %determine if the 'volumn' is different from HvsL
        if log2(HvsL_values/MvsL_values)<=threshold_of_change &...
            log2(HvsL_values/MvsL_values)>=-threshold_of_change
          Hist_array4B(writeout_counter2,hist2_counter1)=...
            Hist_array4B(writeout_counter2,hist2_counter1) + 1;
        elseif log2(HvsL_values/MvsL_values)>threshold_of_change
          Hist_array4C(writeout_counter2,hist2_counter1)=...
            Hist_array4C(writeout_counter2,hist2_counter1) + 1;
        elseif log2(HvsL_values/MvsL_values)<-threshold_of_change
          Hist_array4D(writeout_counter2,hist2_counter1)=...
            Hist_array4D(writeout_counter2,hist2_counter1) + 1;
        end
      end
    end
  end
  
  %add value to a single array to grpah
  Isotopologue_areaA(:,1)=sum(Hist_array3A)';
  Isotopologue_areaA(:,2)=sum(Hist_array3B)';
  Isotopologue_areaA(:,3)=sum(Hist_array3C)';
  Isotopologue_areaA(:,4)=sum(Hist_array3D)';
  
  Isotopologue_areaB(:,1)=sum(Hist_array4A)';
  Isotopologue_areaB(:,2)=sum(Hist_array4B)';
  Isotopologue_areaB(:,3)=sum(Hist_array4C)';
  Isotopologue_areaB(:,4)=sum(Hist_array4D)';
  
  try
    %Find bin which contains the max
    [fitted_valuesMvsL fitted_statsMvsL] = fit(bin_number_hist_MvsL_plus1.',Isotopologue_areaA(:,1),'gauss1');
  catch
    fitted_valuesMvsL.a1=NaN;
    fitted_valuesMvsL.b1=NaN;
    fitted_valuesMvsL.c1=NaN;
    fitted_statsMvsL.rsquare=NaN;
  end
  
  %Generate Gaussian for MvsL and HvsL
  for test_value=1:(bin_number_hist2_MvsL+1)
    %MvsL
    gaussian_fit_MvsL_bin(test_value)=(fitted_valuesMvsL.a1*exp(-((0.05*(test_value)- fitted_valuesMvsL.b1)...
      /fitted_valuesMvsL.c1).^2));
  end
  
  try
    %Find bin which contains the max
    [fitted_valuesHvsL fitted_statsHvsL] = fit(bin_number_hist_HvsL_plus1.',Isotopologue_areaB(:,1),'gauss3');
  catch
    fitted_valuesHvsL.a1=NaN;
    fitted_valuesHvsL.a2=NaN;
    fitted_valuesHvsL.a3=NaN;
    fitted_valuesHvsL.b1=NaN;
    fitted_valuesHvsL.b2=NaN;
    fitted_valuesHvsL.b3=NaN;
    fitted_valuesHvsL.c1=NaN;
    fitted_valuesHvsL.c2=NaN;
    fitted_valuesHvsL.c3=NaN;
    fitted_statsHvsL.rsquare=NaN;
  end
  
  for test_value=1:(bin_number_hist2_HvsL+1)
    %HvsL
    gaussian_fit_HvsL_bin1(test_value)=(fitted_valuesHvsL.a1*exp(-((0.05*(test_value)- fitted_valuesHvsL.b1)...
      /fitted_valuesHvsL.c1).^2));
    gaussian_fit_HvsL_bin2(test_value)=(fitted_valuesHvsL.a2*exp(-((0.05*(test_value)- fitted_valuesHvsL.b2)...
      /fitted_valuesHvsL.c2).^2));
    gaussian_fit_HvsL_bin3(test_value)=(fitted_valuesHvsL.a3*exp(-((0.05*(test_value)- fitted_valuesHvsL.b3)...
      /fitted_valuesHvsL.c3).^2));
  end
  
  
  %Figure showing coverage within PCP-SEC-SILAC experiments
  %%Count how many proteins change across the observed SEC fraction
  bin_number_hist_SEC=0.05:0.05:ceil(max(max(Finalised_Master_Gaussian_list.SEC_coverageMvsL)));
  bin_number_hist_SEC_plus1= 0.05:0.05:(ceil(max(max(Finalised_Master_Gaussian_list.SEC_coverageMvsL)))+0.05);
  bin_number_hist2_SEC=length(bin_number_hist_SEC);
  
  Hist_array5A_SEC=zeros(Dimension_of_master_gaussian_list(1),(bin_number_hist2_SEC+1));
  Hist_array5B_SEC=zeros(Dimension_of_master_gaussian_list(1),(bin_number_hist2_SEC+1));
  Hist_array5C_SEC=zeros(Dimension_of_master_gaussian_list(1),(bin_number_hist2_SEC+1));
  Hist_array5D_SEC=zeros(Dimension_of_master_gaussian_list(1),(bin_number_hist2_SEC+1));
  
  %Count the number of gaussian detected vs the area
  for hist2_counter1= 1:(bin_number_hist2_SEC)+1
    for writeout_counter2= 1:Dimension_of_master_gaussian_list(1)
      %determine bin values
      Coverage_Bin_number=ceil(Finalised_Master_Gaussian_list.SEC_coverageMvsL(writeout_counter2,figure_counter)/0.05);
      %determine if the fold change is within the bin value
      if  hist2_counter1-1==Coverage_Bin_number
        HvsL_values=Finalised_Master_Gaussian_list.Summed_Gaussian_areaHvsL(writeout_counter2,figure_counter);
        MvsL_values=Finalised_Master_Gaussian_list.Summed_Gaussian_areaMvsL(writeout_counter2,figure_counter);
        %determine bin of the total observed isotopologue 'volumn' in MvsL
        Hist_array5A_SEC(writeout_counter2,hist2_counter1)=...
          Hist_array5A_SEC(writeout_counter2,hist2_counter1)+ 1;
        %determine if the 'volumn' is different from HvsL
        if log2(MvsL_values/HvsL_values)<=threshold_of_change &...
            log2(MvsL_values/HvsL_values)>=-threshold_of_change
          Hist_array5B_SEC(writeout_counter2,hist2_counter1)=...
            Hist_array5B_SEC(writeout_counter2,hist2_counter1) + 1;
        elseif log2(MvsL_values/HvsL_values)>threshold_of_change
          Hist_array5C_SEC(writeout_counter2,hist2_counter1)=...
            Hist_array5C_SEC(writeout_counter2,hist2_counter1) + 1;
        elseif log2(MvsL_values/HvsL_values)<-threshold_of_change
          Hist_array5D_SEC(writeout_counter2,hist2_counter1)=...
            Hist_array5D_SEC(writeout_counter2,hist2_counter1) + 1;
        end
      end
    end
  end
  
  Hist_array6A_SEC=zeros(Dimension_of_master_gaussian_list(1),(bin_number_hist2_SEC+1));
  Hist_array6B_SEC=zeros(Dimension_of_master_gaussian_list(1),(bin_number_hist2_SEC+1));
  Hist_array6C_SEC=zeros(Dimension_of_master_gaussian_list(1),(bin_number_hist2_SEC+1));
  Hist_array6D_SEC=zeros(Dimension_of_master_gaussian_list(1),(bin_number_hist2_SEC+1));
  
  %Count the number of gaussian detected vs the area
  for hist2_counter1= 1:(bin_number_hist2_SEC)+1
    for writeout_counter2= 1:Dimension_of_master_gaussian_list(1)
      %determine bin values
      Coverage_Bin_number=ceil(Finalised_Master_Gaussian_list.SEC_coverageHvsL(writeout_counter2,figure_counter)/0.05);
      %determine if the fold change is within the bin value
      if  hist2_counter1-1==Coverage_Bin_number
        HvsL_values=Finalised_Master_Gaussian_list.Summed_Gaussian_areaHvsL(writeout_counter2,figure_counter);
        MvsL_values=Finalised_Master_Gaussian_list.Summed_Gaussian_areaMvsL(writeout_counter2,figure_counter);
        %determine bin of the total observed isotopologue 'volumn' in MvsL
        Hist_array6A_SEC(writeout_counter2,hist2_counter1)=...
          Hist_array6A_SEC(writeout_counter2,hist2_counter1)+ 1;
        %determine if the 'volumn' is different from HvsL
        if log2(HvsL_values/MvsL_values)<=threshold_of_change &...
            log2(HvsL_values/MvsL_values)>=-threshold_of_change
          Hist_array6B_SEC(writeout_counter2,hist2_counter1)=...
            Hist_array6B_SEC(writeout_counter2,hist2_counter1) + 1;
        elseif  log2(HvsL_values/MvsL_values)>threshold_of_change
          Hist_array6C_SEC(writeout_counter2,hist2_counter1)=...
            Hist_array6C_SEC(writeout_counter2,hist2_counter1) + 1;
        elseif log2(HvsL_values/MvsL_values)<-threshold_of_change
          Hist_array6D_SEC(writeout_counter2,hist2_counter1)=...
            Hist_array6D_SEC(writeout_counter2,hist2_counter1) + 1;
        end
      end
    end
  end
  
  %add value to a single array to grpah
  Coverage_areaA(:,1)=sum(Hist_array5A_SEC)';
  Coverage_areaA(:,2)=sum(Hist_array5B_SEC)';
  Coverage_areaA(:,3)=sum(Hist_array5C_SEC)';
  Coverage_areaA(:,4)=sum(Hist_array5D_SEC)';
  
  Coverage_areaB(:,1)=sum(Hist_array6A_SEC)';
  Coverage_areaB(:,2)=sum(Hist_array6B_SEC)';
  Coverage_areaB(:,3)=sum(Hist_array6C_SEC)';
  Coverage_areaB(:,4)=sum(Hist_array6D_SEC)';
  
  
  %Compare coverage to observed guassian within replicates
  %MvsL
  Number_Gaussian_MvsL=zeros(Dimension_of_master_gaussian_list(1),1);
  
  for writeout_counter2= 1:Dimension_of_master_gaussian_list(1)
    %Find position of protein in MvsL array
    [position_inMvsL]=ind2sub(length(MvsL_Gaussians.Protein_name),...
      strmatch(Finalised_Master_Gaussian_list.Protein_name(writeout_counter2), MvsL_Gaussians.Protein_name, 'exact'));
    %determine the length of position_inMvsL
    number_observation=length(position_inMvsL);
    
    %Find how many observations were recorded within each isotopologue channel
    for counter=1:number_observation
      if MvsL_Gaussians.Replicate(position_inMvsL(counter))==figure_counter
        temp_value=MvsL_Gaussians.Unique_identifier(position_inMvsL(counter),:);
        Number_Gaussian_MvsL(writeout_counter2)=sum(~cellfun('isempty',temp_value));
      end
    end
  end
  
  %HvsL
  Number_Gaussian_HvsL=zeros(Dimension_of_master_gaussian_list(1),1);
  
  for writeout_counter2= 1:Dimension_of_master_gaussian_list(1)
    %Find position of protein in MvsL array
    [position_inHvsL]=ind2sub(length(HvsL_Gaussians.Protein_name),...
      strmatch(Finalised_Master_Gaussian_list.Protein_name(writeout_counter2), HvsL_Gaussians.Protein_name, 'exact'));
    %determine the length of position_inMvsL
    number_observation=length(position_inHvsL);
    
    %Find how many observations were recorded within each isotopologue channel
    for counter=1:number_observation
      if HvsL_Gaussians.Replicate(position_inHvsL(counter))==figure_counter
        temp_value=HvsL_Gaussians.Unique_identifier(position_inHvsL(counter),:);
        Number_Gaussian_HvsL(writeout_counter2)=sum(~cellfun('isempty',temp_value));
      end
    end
  end
  
  %Count how many of each Guassian were observed for each bin
  Hist_array7A=zeros(bin_number_hist2_MvsL,6);
  Hist_array7B=zeros(bin_number_hist2_SEC,6);
  Hist_array7C=zeros(bin_number_hist2_HvsL,6);
  Hist_array7D=zeros(bin_number_hist2_SEC,6);
  
  %Count the number of gaussian detected vs the area based on isotopologue area
  for hist2_counter1= 1:(bin_number_hist2_MvsL)+1
    for writeout_counter2= 1:Dimension_of_master_gaussian_list(1)
      %determine bin values
      Coverage_Bin_number=ceil(Finalised_Master_Gaussian_list.Summed_Gaussian_areaMvsL(writeout_counter2,figure_counter)/0.05);
      %determine if the fold change is within the bin value
      if  hist2_counter1-1==Coverage_Bin_number
        %Create matrix couting how Gaussian were detected
        Hist_array7A(hist2_counter1, (Number_Gaussian_MvsL(writeout_counter2)+1))=...
          Hist_array7A(hist2_counter1, (Number_Gaussian_MvsL(writeout_counter2)+1))+1;
      end
    end
  end
  
  %Count the number of gaussian detected vs the area based on isotopologue area
  for hist2_counter1= 1:(bin_number_hist2_HvsL)+1
    for writeout_counter2= 1:Dimension_of_master_gaussian_list(1)
      %determine bin values
      Coverage_Bin_number=ceil(Finalised_Master_Gaussian_list.Summed_Gaussian_areaHvsL(writeout_counter2,figure_counter)/0.05);
      %determine if the fold change is within the bin value
      if  hist2_counter1-1==Coverage_Bin_number
        %Create matrix couting how Gaussian were detected
        Hist_array7C(hist2_counter1, (Number_Gaussian_HvsL(writeout_counter2)+1))=...
          Hist_array7C(hist2_counter1, (Number_Gaussian_HvsL(writeout_counter2)+1))+1;
      end
    end
  end
  
  %Count the number of gaussian detected vs the area based on SEC coverage
  for hist2_counter1= 1:(bin_number_hist2_SEC)+1
    for writeout_counter2= 1:Dimension_of_master_gaussian_list(1)
      %determine bin values
      Coverage_Bin_number=ceil(Finalised_Master_Gaussian_list.SEC_coverageMvsL(writeout_counter2,figure_counter)/0.05);
      %determine if the fold change is within the bin value
      if  hist2_counter1-1==Coverage_Bin_number
        %Create matrix couting how Gaussian were detected
        Hist_array7B(hist2_counter1, (Number_Gaussian_MvsL(writeout_counter2)+1))=...
          Hist_array7B(hist2_counter1, (Number_Gaussian_MvsL(writeout_counter2)+1))+1;
      end
    end
  end
  
  %Count the number of gaussian detected vs the area based on SEC coverage
  for hist2_counter1= 1:(bin_number_hist2_SEC)+1
    for writeout_counter2= 1:Dimension_of_master_gaussian_list(1)
      %determine bin values
      Coverage_Bin_number=ceil(Finalised_Master_Gaussian_list.SEC_coverageHvsL(writeout_counter2,figure_counter)/0.05);
      %determine if the fold change is within the bin value
      if  hist2_counter1-1==Coverage_Bin_number
        %Create matrix couting how Gaussian were detected
        Hist_array7D(hist2_counter1, (Number_Gaussian_HvsL(writeout_counter2)+1))=...
          Hist_array7D(hist2_counter1, (Number_Gaussian_HvsL(writeout_counter2)+1))+1;
      end
    end
  end
  
  %write out figures of coverage
  Text_for_figureMvsL1=strcat('Apex: ', mat2str(round(fitted_valuesMvsL.b1*100)/100));
  Text_for_figureMvsL2=strcat('r-square: ', mat2str(round(fitted_statsMvsL.rsquare*100)/100));
  Text_for_figureHvsL1=strcat('Apex1: ', mat2str(round(fitted_valuesHvsL.b1*100)/100));
  Text_for_figureHvsL2=strcat('Apex2: ', mat2str(round(fitted_valuesHvsL.b2*100)/100));
  Text_for_figureHvsL3=strcat('Apex3: ', mat2str(round(fitted_valuesHvsL.b3*100)/100));
  Text_for_figureHvsL4=strcat('r-square: ', mat2str(round(fitted_statsHvsL.rsquare*100)/100));
  
  %Write out histogram of complex coverage within sample
  figure
  f7=subplot(2,1,1);
  f7_1_figure=bar(bin_number_hist_MvsL_plus1,Isotopologue_areaA(:,1), 0.6, 'stack');
  set(f7_1_figure(1),'facecolor', myC(1,:), 'EdgeColor', 'k' );
  hold on
  plot(bin_number_hist_MvsL_plus1.',gaussian_fit_MvsL_bin,'linewidth', 2, 'Color', colour_to_use(3,:));
  try
    xlim([0,max(Finalised_Master_Gaussian_list.Summed_Gaussian_areaHvsL(:,figure_counter))*1.1])
  catch
  end
  try
    ylim([0,max(Isotopologue_areaA(:,1))*1.1]);
  catch
  end
  title('Total observed protein PCP-SILAC coverage within the Lys4Arg6 sample','FontSize', 12);
  ylabel('Number of proteins','FontSize', 8);
  xlabel('Coverage (expected area under the curve vs observed)','FontSize', 8);
  try
    text(max(Finalised_Master_Gaussian_list.Summed_Gaussian_areaHvsL(:,figure_counter))*0.9,max(Isotopologue_areaA(:,1)), Text_for_figureMvsL1, 'FontSize', 8);
    text(max(Finalised_Master_Gaussian_list.Summed_Gaussian_areaHvsL(:,figure_counter))*0.9,max(Isotopologue_areaA(:,1))*0.90, Text_for_figureMvsL2, 'FontSize', 8);
  catch
  end
  
  f7=subplot(2,1,2);
  f7_2_figure=bar(bin_number_hist_HvsL_plus1,Isotopologue_areaB(:,1) , 0.6, 'stack');
  set(f7_2_figure(1),'facecolor', myC(1,:), 'EdgeColor', 'k' );
  hold on
  fit_fig1=plot(bin_number_hist_HvsL_plus1.',gaussian_fit_HvsL_bin1,'linewidth', 2, 'Color', colour_to_use(6,:));
  hold on
  fit_fig2=plot(bin_number_hist_HvsL_plus1.',gaussian_fit_HvsL_bin2,'linewidth', 2,'Color', colour_to_use(3,:));
  hold on
  fit_fig3=plot(bin_number_hist_HvsL_plus1.',gaussian_fit_HvsL_bin3,'linewidth', 2,'Color', colour_to_use(5,:));
  try
    xlim([0,max(Finalised_Master_Gaussian_list.Summed_Gaussian_areaHvsL(:,figure_counter))*1.1])
  catch
  end
  try
    ylim([0,max(Isotopologue_areaB(:,1))*1.1]);
  catch
  end
  title('Total observed protein PCP-SILAC coverage within the Lys8Arg10 samples','FontSize', 12);
  ylabel('Number of proteins','FontSize', 8);
  xlabel('Coverage (expected area under the curve vs observed)','FontSize', 8);
  try
    text(max(Finalised_Master_Gaussian_list.Summed_Gaussian_areaHvsL(:,figure_counter))*0.9,max(Isotopologue_areaB(:,1)), Text_for_figureHvsL1, 'FontSize', 8);
    text(max(Finalised_Master_Gaussian_list.Summed_Gaussian_areaHvsL(:,figure_counter))*0.9,max(Isotopologue_areaB(:,1))*0.90, Text_for_figureHvsL2, 'FontSize', 8);
    text(max(Finalised_Master_Gaussian_list.Summed_Gaussian_areaHvsL(:,figure_counter))*0.9,max(Isotopologue_areaB(:,1))*0.80, Text_for_figureHvsL3, 'FontSize', 8);
    text(max(Finalised_Master_Gaussian_list.Summed_Gaussian_areaHvsL(:,figure_counter))*0.9,max(Isotopologue_areaB(:,1))*0.70, Text_for_figureHvsL4, 'FontSize', 8);
  catch
  end
  
  % save Histogram
  Save_name_plot=strcat(figdir1,'Isotopologue_coverage_Histogram_replicate',mat2str(figure_counter),'.png');
  saveas(f7, Save_name_plot);
  
  figure
  f8=subplot(2,1,1);
  f8_1_figure=bar(bin_number_hist_SEC_plus1,Coverage_areaA(:,1), 0.6, 'stack');
  set(f8_1_figure(1),'facecolor', myC(1,:), 'EdgeColor', 'k' );
  xlim([0,1]);
  ylim([0,max(Coverage_areaA(:,1))*1.1]);
  title('Total observed protein coverage across SEC fractions within the Lys4Arg6 sample','FontSize', 12);
  ylabel('Number of proteins','FontSize', 8);
  xlabel('Coverage (% of fracations analyzed)','FontSize', 8);
  % Convert y-axis values to percentage values by multiplication
  a=[cellstr(num2str(get(gca,'xtick')'*100))];
  % Create a vector of '%' signs
  pct = char(ones(size(a,1),1)*'%');
  % Append the '%' signs after the percentage values
  new_xticks = [char(a),pct];
  % 'Reflect the changes on the plot
  set(gca,'xticklabel',new_xticks);
  
  f8=subplot(2,1,2);
  f8_2_figure=bar(bin_number_hist_SEC_plus1,Coverage_areaB(:,1), 0.6, 'stack');
  set(f8_2_figure(1),'facecolor', myC(1,:), 'EdgeColor', 'k' );
  xlim([0,1]);
  ylim([0,max(Coverage_areaA(:,1))*1.1]);
  title('Total observed protein coverage across SEC fractions within the Lys8Arg10 sample','FontSize', 12);
  ylabel('Number of proteins','FontSize', 8);
  xlabel('Coverage (% of fracations analyzed)','FontSize', 8);
  % Convert y-axis values to percentage values by multiplication
  a=[cellstr(num2str(get(gca,'xtick')'*100))];
  % Create a vector of '%' signs
  pct = char(ones(size(a,1),1)*'%');
  % Append the '%' signs after the percentage values
  new_xticks = [char(a),pct];
  % 'Reflect the changes on the plot
  set(gca,'xticklabel',new_xticks);
  
  % save Histogram
  Save_name_plot=strcat(figdir1, 'SEC_fraction_coverage_Histogram_replicate',mat2str(figure_counter),'.png');
  saveas(f8, Save_name_plot);
  
  %Create figure of observed changes across Isotopologue_coverage
  figure
  f9=subplot(2,1,1);
  f9_1_figure=bar(bin_number_hist_MvsL_plus1,[Isotopologue_areaA(:,4),Isotopologue_areaA(:,2),Isotopologue_areaA(:,3)] , 0.6, 'stack');
  order_of_columns=[1 3 6];
  for k=1:3
    set(f9_1_figure(k),'facecolor', colour_to_use(order_of_columns(k),:), 'EdgeColor', 'k' );
  end
  legend('Decrease','No change','Increase');
  xlim([0,max(Finalised_Master_Gaussian_list.Summed_Gaussian_areaHvsL(:,figure_counter))*1.1])
  title('Changes (95%), PCP-SILAC coverage within experiments (Lys4Arg6 sample)','FontSize', 12);
  ylabel('Number of proteins','FontSize', 8);
  xlabel('Coverage (expected area under the curve vs observed)','FontSize', 8);
  
  f9=subplot(2,1,2)
  f9_2_figure=bar(bin_number_hist_SEC_plus1,[Coverage_areaA(:,4),Coverage_areaA(:,2),Coverage_areaA(:,3)], 0.6, 'stack');
  for k=1:3
    set(f9_2_figure(k),'facecolor', colour_to_use(order_of_columns(k),:), 'EdgeColor', 'k' );
  end
  xlim([0,1]);
  ylim([0,max(Coverage_areaA(:,1))*1.1]);
  title('Changes (95%), coverage across SEC fractions within experiments (Lys4Arg6 sample)','FontSize', 12);
  ylabel('Number of proteins','FontSize', 8);
  xlabel('Coverage (% of fracations analyzed)','FontSize', 8);
  % Convert y-axis values to percentage values by multiplication
  a=[cellstr(num2str(get(gca,'xtick')'*100))];
  % Create a vector of '%' signs
  pct = char(ones(size(a,1),1)*'%');
  % Append the '%' signs after the percentage values
  new_xticks = [char(a),pct];
  % 'Reflect the changes on the plot
  set(gca,'xticklabel',new_xticks);
  
  % save Histogram
  Save_name_plot=strcat(figdir1, 'Changes_Histogram_Lys4Arg6_replicate',mat2str(figure_counter),'.png');
  saveas(f9, Save_name_plot);
  
  %Create figure of observed changes across Isotopologue_coverage
  figure
  f9B=subplot(2,1,1);
  f9B_1_figure=bar(bin_number_hist_HvsL_plus1,[Isotopologue_areaB(:,4),Isotopologue_areaB(:,2),Isotopologue_areaB(:,3)] , 0.6, 'stack');
  order_of_columns=[1 3 6];
  for k=1:3
    set(f9B_1_figure(k),'facecolor', colour_to_use(order_of_columns(k),:), 'EdgeColor', 'k' );
  end
  legend('Decrease','No change','Increase');
  xlim([0,max(Finalised_Master_Gaussian_list.Summed_Gaussian_areaHvsL(:,figure_counter))*1.1])
  title('Changes (95%), PCP-SILAC coverage within experiments (Lys8Arg10 sample)','FontSize', 12);
  ylabel('Number of proteins','FontSize', 8);
  xlabel('Coverage (expected area under the curve vs observed)','FontSize', 8);
  
  f9B=subplot(2,1,2);
  f9B_2_figure=bar(bin_number_hist_SEC_plus1,[Coverage_areaB(:,4),Coverage_areaB(:,2),Coverage_areaB(:,3)], 0.6, 'stack');
  for k=1:3
    set(f9B_2_figure(k),'facecolor', colour_to_use(order_of_columns(k),:), 'EdgeColor', 'k' );
  end
  xlim([0,1]);
  ylim([0,max(Coverage_areaB(:,1))*1.1]);
  title('Changes (95%), coverage across SEC fractions within experiments (Lys8Arg10 sample)','FontSize', 12);
  ylabel('Number of proteins','FontSize', 8);
  xlabel('Coverage (% of fracations analyzed)','FontSize', 8);
  % Convert y-axis values to percentage values by multiplication
  a=[cellstr(num2str(get(gca,'xtick')'*100))];
  % Create a vector of '%' signs
  pct = char(ones(size(a,1),1)*'%');
  % Append the '%' signs after the percentage values
  new_xticks = [char(a),pct];
  % 'Reflect the changes on the plot
  set(gca,'xticklabel',new_xticks);
  
  % save Histogram
  Save_name_plot=strcat(figdir1, 'Changes_Histogram_Lys8Arg10_replicate ',mat2str(figure_counter),'.png');
  saveas(f9B, Save_name_plot);
  
  figure
  f10=subplot(2,1,1)
  f10_1_figure=bar(bin_number_hist_MvsL,[Hist_array7A(:,1),Hist_array7A(:,2),Hist_array7A(:,3),...
    Hist_array7A(:,4),Hist_array7A(:,5),Hist_array7A(:,6)], 0.6, 'stack');
  order_of_columns=[1 3 6];
  for k=1:6
    set(f10_1_figure(k),'facecolor', colour_to_use((k),:), 'EdgeColor', 'k' );
  end
  legend('zero','one','two','three','four','five','FontSize',8, 'Location', 'EastOutside');
  xlim([0,max(Finalised_Master_Gaussian_list.Summed_Gaussian_areaHvsL(:,figure_counter))*1.1])
  title('Number of fitted Gaussian compared to PCP-SILAC coverage','FontSize', 12);
  ylabel('Number of proteins','FontSize', 8);
  xlabel('Coverage (expected area under the curve vs observed)','FontSize', 8);
  
  f10=subplot(2,1,2);
  f10_2_figure=bar(bin_number_hist_SEC,[Hist_array7B(:,1),Hist_array7B(:,2),Hist_array7B(:,3),...
    Hist_array7B(:,4),Hist_array7B(:,5),Hist_array7B(:,6)], 0.6, 'stack');
  for k=1:6
    set(f10_2_figure(k),'facecolor', colour_to_use((k),:), 'EdgeColor', 'k' );
  end
  xlim([0,1]);
  ylim([0,max(Coverage_areaA(:,1))*1.1]);
  legend('zero','one','two','three','four','five','FontSize',8, 'Location', 'EastOutside');
  title('Number of fitted Gaussian compared to SEC fractions coverage','FontSize', 12);
  ylabel('Number of proteins','FontSize', 8);
  xlabel('Coverage (% of fracations analyzed)','FontSize', 8);
  % Convert y-axis values to percentage values by multiplication
  a=[cellstr(num2str(get(gca,'xtick')'*100))];
  % Create a vector of '%' signs
  pct = char(ones(size(a,1),1)*'%');
  % Append the '%' signs after the percentage values
  new_xticks = [char(a),pct];
  % 'Reflect the changes on the plot
  set(gca,'xticklabel',new_xticks);
  
  % save Histogram
  Save_name_plot=strcat(figdir1, 'Gaussians_Histogram_Lys4Arg6_replicate ',mat2str(figure_counter),'.png');
  saveas(f10, Save_name_plot);
  
  figure
  f10B=subplot(2,1,1);
  f10B_1_figure=bar(bin_number_hist_HvsL,[Hist_array7C(:,1),Hist_array7C(:,2),Hist_array7C(:,3),...
    Hist_array7C(:,4),Hist_array7C(:,5),Hist_array7C(:,6)], 0.6, 'stack');
  order_of_columns=[1 3 6];
  for k=1:6
    set(f10B_1_figure(k),'facecolor', colour_to_use((k),:), 'EdgeColor', 'k' );
  end
  legend('zero','one','two','three','four','five','FontSize',8, 'Location', 'EastOutside');
  xlim([0,max(Finalised_Master_Gaussian_list.Summed_Gaussian_areaHvsL(:,figure_counter))*1.1])
  title('Number of fitted Gaussian compared to PCP-SILAC coverage','FontSize', 12);
  ylabel('Number of proteins','FontSize', 8);
  xlabel('Coverage (expected area under the curve vs observed)','FontSize', 8);
  
  f10B=subplot(2,1,2);
  f10B_2_figure=bar(bin_number_hist_SEC,[Hist_array7D(:,1),Hist_array7D(:,2),Hist_array7D(:,3),...
    Hist_array7D(:,4),Hist_array7D(:,5),Hist_array7D(:,6)], 0.6, 'stack');
  for k=1:6
    set(f10B_2_figure(k),'facecolor', colour_to_use((k),:), 'EdgeColor', 'k' );
  end
  xlim([0,1]);
  ylim([0,max(Coverage_areaB(:,1))*1.1]);
  legend('zero','one','two','three','four','five','FontSize',8, 'Location', 'EastOutside');
  title('Number of fitted Gaussian compared to SEC fractions coverage','FontSize', 12);
  ylabel('Number of proteins','FontSize', 8);
  xlabel('Coverage (% of fracations analyzed)','FontSize', 8);
  % Convert y-axis values to percentage values by multiplication
  a=[cellstr(num2str(get(gca,'xtick')'*100))];
  % Create a vector of '%' signs
  pct = char(ones(size(a,1),1)*'%');
  % Append the '%' signs after the percentage values
  new_xticks = [char(a),pct];
  % 'Reflect the changes on the plot
  set(gca,'xticklabel',new_xticks);
  
  % save Histogram
  Save_name_plot=strcat(figdir1, 'Gaussians_Histogram_Lys8Arg10_replicate ',mat2str(figure_counter),'.png');
  saveas(f10B, Save_name_plot);
  
  figure
  f11=subplot(2,1,1);
  f11_1_figure=bar(bin_number_hist_HvsL_plus1,Isotopologue_areaB(:,1), 0.6, 'stack');
  set(f11_1_figure(1),'facecolor', myC(1,:), 'EdgeColor', 'k' );
  hold on
  fit_fig1=plot(bin_number_hist_HvsL_plus1.',gaussian_fit_HvsL_bin1,'linewidth', 2, 'Color', colour_to_use(6,:));
  hold on
  fit_fig2=plot(bin_number_hist_HvsL_plus1.',gaussian_fit_HvsL_bin2,'linewidth', 2,'Color', colour_to_use(3,:));
  hold on
  fit_fig3=plot(bin_number_hist_HvsL_plus1.',gaussian_fit_HvsL_bin3,'linewidth', 2,'Color', colour_to_use(5,:));
  try
    xlim([0,max(Finalised_Master_Gaussian_list.Summed_Gaussian_areaHvsL(:,figure_counter))*1.1])
  catch
  end
  try
    ylim([0,max(Isotopologue_areaB(:,1))*1.1]);
  catch
  end
  title('Total observed protein PCP-SILAC coverage within the Lys8Arg10 samples','FontSize', 12);
  ylabel('Number of proteins','FontSize', 8);
  xlabel('Coverage (expected area under the curve vs observed)','FontSize', 8);
  try
    text(max(Finalised_Master_Gaussian_list.Summed_Gaussian_areaHvsL(:,figure_counter))*0.9,max(Isotopologue_areaB(:,1)), Text_for_figureHvsL1, 'FontSize', 8);
    text(max(Finalised_Master_Gaussian_list.Summed_Gaussian_areaHvsL(:,figure_counter))*0.9,max(Isotopologue_areaB(:,1))*0.90, Text_for_figureHvsL2, 'FontSize', 8);
    text(max(Finalised_Master_Gaussian_list.Summed_Gaussian_areaHvsL(:,figure_counter))*0.9,max(Isotopologue_areaB(:,1))*0.80, Text_for_figureHvsL3, 'FontSize', 8);
    text(max(Finalised_Master_Gaussian_list.Summed_Gaussian_areaHvsL(:,figure_counter))*0.9,max(Isotopologue_areaB(:,1))*0.70, Text_for_figureHvsL4, 'FontSize', 8);
  catch
  end
  
  f11=subplot(2,1,2)
  f11_1_figure=bar(bin_number_hist_HvsL_plus1,[Isotopologue_areaB(:,4),Isotopologue_areaB(:,2),Isotopologue_areaB(:,3)] , 0.6, 'stack');
  order_of_columns=[1 3 6];
  for k=1:3
    set(f11_1_figure(k),'facecolor', colour_to_use(order_of_columns(k),:), 'EdgeColor', 'k' );
  end
  legend('Decrease','No change','Increase');
  xlim([0,max(Finalised_Master_Gaussian_list.Summed_Gaussian_areaHvsL(:,figure_counter))*1.1])
  title('Changes (95%), PCP-SILAC coverage within experiments (Lys8Arg10 sample)','FontSize', 12);
  ylabel('Number of proteins','FontSize', 8);
  xlabel('Coverage (expected area under the curve vs observed)','FontSize', 8);
  
  %save Histogram
  Save_name_plot=strcat(figdir1,'Comparsion_Histogram_Lys8Arg10_replicate ',mat2str(figure_counter),'.png');
  saveas(f11, Save_name_plot);
end

%%Count how many guassian change in each fraction
Hist_array4_MvsL=zeros(fraction_to_plot,3);
Hist_array4_HvsL=zeros(fraction_to_plot,3);

for fraction_counter1= 1:Dimension_of_master_gaussian_list(1)
  for fraction_counter2= 1:nnz(Finalised_Master_Gaussian_list.Center(fraction_counter1,:))
    %Find Center
    Center_of_gaussain= floor(Finalised_Master_Gaussian_list.Center(fraction_counter1,fraction_counter2));
    %Determine if Center is greater then zero
    if  Center_of_gaussain>0
      
      if Finalised_Master_Gaussian_list.Average_fold_change(fraction_counter1,fraction_counter2) >0 &&...
          Center_of_gaussain <=fraction_to_plot && isequal(Finalised_Master_Gaussian_list.Satifies_Adjusted_pvalue(fraction_counter1,fraction_counter2),{'+'})
        
        Hist_array4_MvsL(Center_of_gaussain,1)=Hist_array4_MvsL(Center_of_gaussain,1)+1;
        Hist_array4_MvsL(Center_of_gaussain,2)=Hist_array4_MvsL(Center_of_gaussain,2)+1;
      elseif Finalised_Master_Gaussian_list.Average_fold_change(fraction_counter1,fraction_counter2) <0 &&...
          Center_of_gaussain <=fraction_to_plot && isequal(Finalised_Master_Gaussian_list.Satifies_Adjusted_pvalue(fraction_counter1,fraction_counter2),{'+'})
        
        Hist_array4_MvsL(Center_of_gaussain,1)=Hist_array4_MvsL(Center_of_gaussain,1)+1;
        Hist_array4_MvsL(Center_of_gaussain,3)=Hist_array4_MvsL(Center_of_gaussain,3)+1;
      elseif Center_of_gaussain <=fraction_to_plot
        
        Hist_array4_MvsL(Center_of_gaussain,1)=Hist_array4_MvsL(Center_of_gaussain,1)+1;
      end
    end
  end
end

%Create array to plot
plotting_counter=1;

for writeout_counter1= 1:Dimension_of_master_gaussian_list(1)
  for writeout_counter2= 1:nnz(Finalised_Master_Gaussian_list.Center(writeout_counter1,:))
    data_scatter_plot(plotting_counter,1)=Finalised_Master_Gaussian_list.Average_fold_change(writeout_counter1,writeout_counter2);
    data_scatter_plot(plotting_counter,2)=Finalised_Master_Gaussian_list.Stdev_fold_change(writeout_counter1,writeout_counter2);
    if isequal(Finalised_Master_Gaussian_list.Satifies_Adjusted_pvalue(writeout_counter1,writeout_counter2),{'+'})
      data_scatter_plot(plotting_counter,3)=1;
    else
      data_scatter_plot(plotting_counter,3)=0;
    end
    plotting_counter=plotting_counter+1;
  end
end

%Sort arrya for plot and divided
[unique_log2_sorted1,unique_log2_sorted2]=sort(data_scatter_plot(:,1));
data_scatter_plot=[data_scatter_plot(unique_log2_sorted2,1) data_scatter_plot(unique_log2_sorted2,2) data_scatter_plot(unique_log2_sorted2,3)];

%Create array
length_values_plot=length(data_scatter_plot);
Multiple_measurement_Adjusted_pvalue_corrected=zeros(length_values_plot,2);

%Work out mean stdev
work_out_mean=data_scatter_plot(:,2);
work_out_mean(work_out_mean==0)=[];
mean_stdev=mean(work_out_mean);

%Divide dataset to colour code individual changes
Gaussian_ploted=0;
Gaussian_not_ploted=1;
for counter= 1:length(data_scatter_plot)
  if data_scatter_plot(counter,3) == 1 && (data_scatter_plot(counter,1) >0.2 | data_scatter_plot(counter,1) <-0.2 )
    Multiple_measurement_Adjusted_pvalue_corrected(counter,1)=data_scatter_plot(counter,1);
    Multiple_measurement_Adjusted_pvalue_corrected(counter,2)=data_scatter_plot(counter,2);
    Gaussian_ploted=1+Gaussian_ploted;
  elseif ~(data_scatter_plot(counter,3) == 1 && (data_scatter_plot(counter,1) >0.2 | data_scatter_plot(counter,1) <-0.2 ))
    Multiple_measurement_Adjusted_pvalue_corrected(counter,1)=NaN;
    Multiple_measurement_Adjusted_pvalue_corrected(counter,2)=NaN;
    Gaussian_not_ploted=1+Gaussian_not_ploted;
  end
end

%Text for figure
Text_for_figure1=strcat('Mean std dev: ', mat2str((round(mean_stdev*1000))/1000));
Text_for_figure2=strcat('Number of Guassian below bonferroni adjusted p-value: ', mat2str(Gaussian_ploted));

%Graph data as log2 scatter
figure;
f_all=subplot(2,1,1);
P4A =scatter([1:length_values_plot],data_scatter_plot(:,1),5,'fill', 'markerfacecolor',colour_to_use(1,:));
hold on
P4B =scatter([1:length_values_plot],Multiple_measurement_Adjusted_pvalue_corrected(:,1), 5,'fill', 'markerfacecolor',colour_to_use(6,:));
title('Log2 changes in Gaussians: Cyotplasmic- Fas treated','FontSize', 12);
ylabel('Log2 ratio','FontSize', 10);
xlabel('Gaussians identified','FontSize', 10);
ylim([data_scatter_plot(1,1)*1.05,data_scatter_plot(end,1)*1.05]);
xlim([-100,length_values_plot+100]);
hold on %add 1 dotted line
PP1= plot([-100,length_values_plot+100],[1 1],':');
hold on %add 1 dotted line
PP1= plot([-100,length_values_plot+100],[-1 -1],':');
text(100, 3.5, Text_for_figure1, 'FontSize', 8);
text(100, 3.0, Text_for_figure2, 'FontSize', 8);

%plot distribution of protein chnages across the fractions
%define colours
subplot(2,1,2);
f_all_figure=bar(1:fraction_to_plot, [Hist_array4_MvsL(:,1) Hist_array4_MvsL(:,2) Hist_array4_MvsL(:,3)],0.6, 'stack');
for k=1:3
  set(f_all_figure(k),'facecolor', myC(k,:), 'EdgeColor', 'k' )
end
legend('No change', 'Increase', 'Decrease', 'FontSize',8, 'Location', 'Best');
xlim([-2,fraction_to_plot+2]);
title('Changes in Gaussians observed across fractions');
xlabel('Fractions');
ylabel('# of Guassian centers');
saveas(f_all, [figdir1 'Unique_gaussian_changes_observed.png']);

%Create figures for all proteins showing overlay of
%Create a file to write pdf name out to
List_of_pdf=cell(Dimension_of_master_gaussian_list(1)*2,2);
List_of_pdf_counter=1;


%Determine which Gaussians change by comparing fitted data with raw data
for Gaussian_counter1= 1:Dimension_of_master_gaussian_list(1)
  
  %Protein being plotted
  Protein_to_plot=Finalised_Master_Gaussian_list.Protein_name{Gaussian_counter1};
  
  %Replicate 1
  %MvsL Data to plot for replicate
  MvsL_to_plot_replicate1=num_val_MvsL_for_figures((Finalised_Master_Gaussian_list.Protein_location_replicate1(Gaussian_counter1,1)-1),(1+position_fraction1:fraction_number(2)-5));
  for NaN_removal_counter=1:fraction_to_plot
    if isnan(MvsL_to_plot_replicate1(NaN_removal_counter))
      MvsL_to_plot_replicate1(NaN_removal_counter)= 0;
    end
  end
  %Smooth data to use in figure for replicate 1
  MvsL_to_plot_replicate1=smooth(MvsL_to_plot_replicate1,2);
  Max_MvsL_replicate1=max(MvsL_to_plot_replicate1);
  
  %HvsL Data to plot for replicate 1
  HvsL_to_plot_replicate1=num_val_HvsL_for_figures((Finalised_Master_Gaussian_list.Protein_location_replicate1(Gaussian_counter1,1)-1),(1+position_fraction1:fraction_number(2)-5));
  for NaN_removal_counter=1:fraction_to_plot
    if isnan(HvsL_to_plot_replicate1(NaN_removal_counter))
      HvsL_to_plot_replicate1(NaN_removal_counter)= 0;
    end
  end
  %Smooth data to use in figure for replicate 1
  HvsL_to_plot_replicate1=smooth(HvsL_to_plot_replicate1,2);
  Max_HvsL_replicate1=max(HvsL_to_plot_replicate1);
  
  %Maximum of both mvsl and HvsL replicate 1
  real_max_replicate1=max([Max_MvsL_replicate1 Max_HvsL_replicate1]);
  
  if replicate_num >= 2
    %Replicate 2
    %MvsL Data to plot for replicate 2
    MvsL_to_plot_replicate2=num_val_MvsL_for_figures((Finalised_Master_Gaussian_list.Protein_location_replicate2(Gaussian_counter1,1)-1),(1+position_fraction1:fraction_number(2)-5));
    for NaN_removal_counter=1:fraction_to_plot
      if isnan(MvsL_to_plot_replicate2(NaN_removal_counter))
        MvsL_to_plot_replicate2(NaN_removal_counter)= 0;
      end
    end
    %Smooth data to use in figure for replicate 2
    MvsL_to_plot_replicate2=smooth(MvsL_to_plot_replicate2,2);
    Max_MvsL_replicate2=max(MvsL_to_plot_replicate2);
    
    %HvsL Data to plot for replicate 2
    HvsL_to_plot_replicate2=num_val_HvsL_for_figures((Finalised_Master_Gaussian_list.Protein_location_replicate2(Gaussian_counter1,1)-1),(1+position_fraction1:fraction_number(2)-5));
    for NaN_removal_counter=1:fraction_to_plot
      if isnan(HvsL_to_plot_replicate2(NaN_removal_counter))
        HvsL_to_plot_replicate2(NaN_removal_counter)= 0;
      end
    end
    %Smooth data to use in figure for replicate 2
    HvsL_to_plot_replicate2=smooth(HvsL_to_plot_replicate2,2);
    Max_HvsL_replicate2=max(HvsL_to_plot_replicate2);
    
    %Maximum of both mvsl and HvsL replicate 2
    real_max_replicate2=max([Max_MvsL_replicate2 Max_HvsL_replicate2]);
  else
    real_max_replicate2=0;
  end
  
  if replicate_num >= 3
    %Replicate 3
    %MvsL Data to plot for replicate 3
    MvsL_to_plot_replicate3=num_val_MvsL_for_figures((Finalised_Master_Gaussian_list.Protein_location_replicate3(Gaussian_counter1,1)-1),(1+position_fraction1:fraction_number(2)-5));
    for NaN_removal_counter=1:fraction_to_plot
      if isnan(MvsL_to_plot_replicate3(NaN_removal_counter))
        MvsL_to_plot_replicate3(NaN_removal_counter)= 0;
      end
    end
    %Smooth data to use in figure for replicate 3
    MvsL_to_plot_replicate3=smooth(MvsL_to_plot_replicate3,2);
    Max_MvsL_replicate3=max(MvsL_to_plot_replicate3);
    
    %HvsL Data to plot for replicate 3
    HvsL_to_plot_replicate3=num_val_HvsL_for_figures((Finalised_Master_Gaussian_list.Protein_location_replicate3(Gaussian_counter1,1)-1),(1+position_fraction1:fraction_number(2)-5));
    for NaN_removal_counter=1:fraction_to_plot
      if isnan(HvsL_to_plot_replicate3(NaN_removal_counter))
        HvsL_to_plot_replicate3(NaN_removal_counter)= 0;
      end
    end
    
    %Smooth data to use in figure for replicate 3
    HvsL_to_plot_replicate3=smooth(HvsL_to_plot_replicate3,2);
    Max_HvsL_replicate3=max(HvsL_to_plot_replicate3);
    
    %Maximum of both mvsl and HvsL replicate 3
    real_max_replicate3=max([Max_MvsL_replicate3 Max_HvsL_replicate3]);
  else
    real_max_replicate3=0;
  end
  
  if replicate_num >=4
    %Replicate 4
    %MvsL Data to plot for replicate 4
    MvsL_to_plot_replicate4=num_val_MvsL_for_figures((Finalised_Master_Gaussian_list.Protein_location_replicate4(Gaussian_counter1,1)-1),(1+position_fraction1:fraction_number(2)-5));
    for NaN_removal_counter=1:fraction_to_plot
      if isnan(MvsL_to_plot_replicate4(NaN_removal_counter))
        MvsL_to_plot_replicate4(NaN_removal_counter)= 0;
      end
    end
    %Smooth data to use in figure for replicate 4
    MvsL_to_plot_replicate4=smooth(MvsL_to_plot_replicate4,2);
    Max_MvsL_replicate4=max(MvsL_to_plot_replicate4);
    
    
    %HvsL Data to plot for replicate 4
    HvsL_to_plot_replicate4=num_val_HvsL_for_figures((Finalised_Master_Gaussian_list.Protein_location_replicate4(Gaussian_counter1,1)-1),(1+position_fraction1:fraction_number(2)-5));
    for NaN_removal_counter=1:fraction_to_plot
      if isnan(HvsL_to_plot_replicate4(NaN_removal_counter))
        HvsL_to_plot_replicate4(NaN_removal_counter)= 0;
      end
    end
    
    %Smooth data to use in figure for replicate 4
    HvsL_to_plot_replicate4=smooth(HvsL_to_plot_replicate4,2);
    Max_HvsL_replicate4=max(HvsL_to_plot_replicate4);
    
    %Maximum of both mvsl and HvsL replicate 4
    real_max_replicate4=max([Max_MvsL_replicate4 Max_HvsL_replicate4]);
  else
    real_max_replicate4=0;
  end
  
  %plot Master gaussian map
  f8 = subplot((replicate_num+1),1,1);
  %Count the number of gaussians detected
  number_of_gaussian_to_plot= nnz(Finalised_Master_Gaussian_list.Center(Gaussian_counter1,:));
  
  for hold_on_counter=1:number_of_gaussian_to_plot
    hold on  %Graph Gaus
    Center_output=Finalised_Master_Gaussian_list.Center(Gaussian_counter1, hold_on_counter);
    Height_output=Finalised_Master_Gaussian_list.Height(Gaussian_counter1, hold_on_counter);
    Width_output=Finalised_Master_Gaussian_list.Width(Gaussian_counter1, hold_on_counter);
    Fitted_gaus=1:0.1:fraction_to_plot;
    number_of_data_points=length(Fitted_gaus);
    gaus_colour_set=colour_to_use(hold_on_counter,:);
    for fill_in_counter=1:number_of_data_points
      hold on  %Graph Gaus fill area
      y1 =  Height_output*exp(-((Fitted_gaus(fill_in_counter)-Center_output)/Width_output).^2);
      patch([Fitted_gaus(fill_in_counter) Fitted_gaus(fill_in_counter)], [0 y1], 'w','EdgeColor',gaus_colour_set(:),'EdgeAlpha',0.2,'LineWidth',2);
    end
    y1 =  Height_output*exp(-((Fitted_gaus-Center_output)/Width_output).^2);
    P1 = plot(Fitted_gaus,y1);
    set(P1,'Color','black','LineWidth',1);
    xlim([0,fraction_to_plot])
  end
  title_name_plot=strcat('Gaussian curves identified across replicates of :',Protein_to_plot);
  title(title_name_plot,'FontSize', 12);
  ylabel('Isotopologue ratio','FontSize', 6);
  xlabel('Fractions','FontSize', 6);
  
  
  %test if detected in replicate
  if ~(real_max_replicate1==0)
    %Plot replicate 1 raw data
    hold on %Graph points identified in medium channel
    f8 = subplot((replicate_num+1),1,2);
    P2A =plot((1:fraction_to_plot),MvsL_to_plot_replicate1/real_max_replicate1);
    set(P2A,'Color', [0.5 0.5 0.5],'LineWidth',0.5);
    hold on
    P2B =plot((1:fraction_to_plot),MvsL_to_plot_replicate1/real_max_replicate1,'s');
    set(P2B,'MarkerFaceColor',colour_to_use(2,:),'MarkerSize',4);
    hold on  %Graph points identified in heavy channel
    P3A =plot((1:fraction_to_plot),HvsL_to_plot_replicate1/real_max_replicate1);
    set(P3A,'Color', [0.5 0.5 0.5],'LineWidth',0.5);
    hold on
    P3B =plot((1:fraction_to_plot),HvsL_to_plot_replicate1/real_max_replicate1,'d');
    set(P3B,'MarkerFaceColor',colour_to_use(5,:),'MarkerSize',5);
    xlim([0,fraction_to_plot])
    ylim([0,(real_max_replicate1/real_max_replicate1)*1.05]); % Set to adjust to maximum
    ylabel('Isotopologue ratio','FontSize', 6);
    xlabel('Fractions','FontSize', 6);
    hold on %add 0.2 dotted line gaussian under which qunation is not considered
    PP1B= plot([0 fraction_to_plot],[0.2/real_max_replicate1 0.2/real_max_replicate1],':');
    set(PP1B,'Color','red','LineWidth',1);
    hold on %legend
    legend([P2B P3B],{'MvsL', 'HvsL'}, 'Location', 'Best');
    legend boxoff;
  end
  
  %test if detected in replicate
  if ~(real_max_replicate2==0)
    %Plot replicate 2 raw data
    hold on %Graph points identified in medium channel
    f8 = subplot((replicate_num+1),1,3);
    P4A =plot((1:fraction_to_plot),MvsL_to_plot_replicate2/real_max_replicate2);
    set(P4A,'Color', [0.5 0.5 0.5],'LineWidth',0.5);
    hold on
    P4B =plot((1:fraction_to_plot),MvsL_to_plot_replicate2/real_max_replicate2,'s');
    set(P4B,'MarkerFaceColor',colour_to_use(2,:),'MarkerSize',4);
    hold on  %Graph points identified in heavy channel
    P5A =plot((1:fraction_to_plot),HvsL_to_plot_replicate2/real_max_replicate2);
    set(P5A,'Color', [0.5 0.5 0.5],'LineWidth',0.5);
    hold on
    P5B =plot((1:fraction_to_plot),HvsL_to_plot_replicate2/real_max_replicate2,'d');
    set(P5B,'MarkerFaceColor',colour_to_use(5,:),'MarkerSize',5);
    xlim([0,fraction_to_plot])
    ylim([0,(real_max_replicate2/real_max_replicate2)*1.05]);% Set to adjust to maximum
    ylabel('Isotopologue realtive %','FontSize', 6);
    xlabel('Fractions','FontSize', 6);
    hold on %add 0.2 dotted line gaussian under which qunation is not considered
    PP1B= plot([0 fraction_to_plot],[0.2/real_max_replicate2 0.2/real_max_replicate2],':');
    set(PP1B,'Color','red','LineWidth',1);
    hold on %legend
    legend([P4B P5B],{'MvsL', 'HvsL'}, 'Location', 'Best');
    legend boxoff;
  end
  
  if ~(real_max_replicate3==0)
    %Plot replicate 3 raw data
    hold on %Graph points identified in medium channel
    f8 = subplot((replicate_num+1),1,4);
    P6A =plot((1:fraction_to_plot),MvsL_to_plot_replicate3/real_max_replicate3);
    set(P6A,'Color', [0.5 0.5 0.5],'LineWidth',0.5);
    hold on
    P6B =plot((1:fraction_to_plot),MvsL_to_plot_replicate3/real_max_replicate3,'s');
    set(P6B,'MarkerFaceColor',colour_to_use(2,:),'MarkerSize',4);
    hold on  %Graph points identified in heavy channel
    P7A =plot((1:fraction_to_plot),HvsL_to_plot_replicate3/real_max_replicate3);
    set(P7A,'Color', [0.5 0.5 0.5],'LineWidth',0.5);
    hold on
    P7B =plot((1:fraction_to_plot),HvsL_to_plot_replicate3/real_max_replicate3,'d');
    set(P7B,'MarkerFaceColor',colour_to_use(5,:),'MarkerSize',5);
    xlim([0,fraction_to_plot])
    ylim([0,(real_max_replicate3/real_max_replicate3)*1.05]); % Set to adjust to maximum
    ylabel('Isotopologue realtive %','FontSize', 6);
    xlabel('Fractions','FontSize', 6);
    hold on %add 0.2 dotted line gaussian under which qunation is not considered
    PP1B= plot([0 fraction_to_plot],[0.2/real_max_replicate3 0.2/real_max_replicate3],':');
    set(PP1B,'Color','red','LineWidth',1);
    hold on %legend
    legend([P6B P7B],{'MvsL', 'HvsL'}, 'Location', 'Best');
    legend boxoff;
  end
  
  if ~(real_max_replicate4)==0
    %Plot replicate 4 raw data
    hold on %Graph points identified in medium channel
    f8 = subplot((replicate_num+1),1,5);
    P8A =plot((1:fraction_to_plot),MvsL_to_plot_replicate4/real_max_replicate4);
    set(P8A,'Color', [0.5 0.5 0.5],'LineWidth',0.5);
    hold on
    P8B =plot((1:fraction_to_plot),MvsL_to_plot_replicate4/real_max_replicate4,'s');
    set(P8B,'MarkerFaceColor',colour_to_use(2,:),'MarkerSize',4);
    hold on  %Graph points identified in heavy channel
    P9A =plot((1:fraction_to_plot),HvsL_to_plot_replicate4/real_max_replicate4);
    set(P9A,'Color', [0.5 0.5 0.5],'LineWidth',0.5);
    hold on
    P9B =plot((1:fraction_to_plot),HvsL_to_plot_replicate4/real_max_replicate4,'d');
    set(P9B,'MarkerFaceColor',colour_to_use(5,:),'MarkerSize',5);
    xlim([0,fraction_to_plot])
    ylim([0,(real_max_replicate4/real_max_replicate4)*1.05]); % Set to adjust to maximum
    ylabel('Isotopologue relative %','FontSize', 6);
    xlabel('Fractions','FontSize', 6);
    hold on %add 0.2 dotted line gaussian under which qunation is not considered
    PP1B= plot([0 fraction_to_plot1],[0.2/real_max_replicate4 0.2/real_max_replicate4],':');
    set(PP1B,'Color','red','LineWidth',1);
    hold on %legend
    legend([P8B P9B],{'MvsL', 'HvsL'}, 'Location', 'Best');
    legend boxoff;
  end
  
  Save_name_replicates=strcat(figdir2, mat2str(Gaussian_counter1),'_1_PCP_SEC_Profiles of_',Protein_to_plot,'.png');
  List_of_pdf{List_of_pdf_counter,1}=Save_name_replicates;
  List_of_pdf{List_of_pdf_counter,2}=List_of_pdf_counter;
  List_of_pdf_counter=List_of_pdf_counter+1;
  saveas(f8, Save_name_replicates);
  %print('-dpdf', '-r600', Save_name_replicates);
  close 'all';
  
  % Determine if Gaussians should be plotted
  
  %Count the number of gaussians detected
  number_of_gaussian_to_plot= nnz(Finalised_Master_Gaussian_list.Center(Gaussian_counter1,:));
  
  %create counter
  hold_on_counter1=1;
  
  %rest varibles
  Center_test=[];
  Height_test=[];
  Width_test=[];
  
  %determine if Center is less then fraction_to_plot
  for hold_on_counter=1:number_of_gaussian_to_plot
    if ~(Finalised_Master_Gaussian_list.Center(Gaussian_counter1, hold_on_counter)>fraction_to_plot-2) % minus two add for comsetics
      Center_test(hold_on_counter1)=Finalised_Master_Gaussian_list.Center(Gaussian_counter1, hold_on_counter);
      Height_test(hold_on_counter1)=Finalised_Master_Gaussian_list.Height(Gaussian_counter1, hold_on_counter);
      Width_test(hold_on_counter1)=Finalised_Master_Gaussian_list.Width(Gaussian_counter1, hold_on_counter);
      hold_on_counter1=hold_on_counter1+1;
    end
  end
  
  number_of_gaussian_to_plot=length(Center_test);
  
  if ~isempty(Center_test)
    %Plot Quantiation_of proteins
    f9 = subplot(2,1,1);
    for hold_on_counter=1:number_of_gaussian_to_plot
      hold on  %Graph Gaus
      Center_output=Center_test(hold_on_counter);
      Height_output=Height_test(hold_on_counter);
      Width_output=Width_test(hold_on_counter);
      Fitted_gaus=1:0.1:fraction_number(2);
      number_of_data_points=length(Fitted_gaus);
      gaus_colour_set=colour_to_use(hold_on_counter,:);
      for fill_in_counter=1:number_of_data_points
        hold on  %Graph Gaus fill area
        y1 =  Height_output*exp(-((Fitted_gaus(fill_in_counter)-Center_output)/Width_output).^2);
        patch([Fitted_gaus(fill_in_counter) Fitted_gaus(fill_in_counter)], [0 y1], 'w','EdgeColor',gaus_colour_set(:),'EdgeAlpha',0.2,'LineWidth',2);
      end
      y1 =  Height_output*exp(-((Fitted_gaus-Center_output)/Width_output).^2);
      P1 = plot(Fitted_gaus,y1);
      set(P1,'Color','black','LineWidth',1);
      xlim([0,fraction_to_plot]);
    end
    title_name_plot=strcat('Gaussian curves identified across replicates of :',Protein_to_plot);
    title(title_name_plot,'FontSize', 12);
    ylabel('Isotopologue ratio','FontSize', 10);
    xlabel('Fractions','FontSize', 10);
    
    
    f9 = subplot(2,1,2);
    %Format log2 values into ascending order to plot
    Values_for_bar_graph=zeros(number_of_gaussian_to_plot*replicate_num,1);
    Names_for_bar_graph=cell(number_of_gaussian_to_plot*replicate_num,1);
    bar_position=[1:replicate_num:(number_of_gaussian_to_plot*replicate_num)];
    
    Center_to_plot =Center_test(:);
    Center_counter1=1;
    Gaussian_number_counter1=1;
    
    %Copy value to matrix to manipulate
    Replicate1_normalised_raw_data=Finalised_Master_Gaussian_list.normaliased_log2_comparsion_replicate_1(Gaussian_counter1,1:number_of_gaussian_to_plot);
    Replicate2_normalised_raw_data=Finalised_Master_Gaussian_list.normaliased_log2_comparsion_replicate_2(Gaussian_counter1,1:number_of_gaussian_to_plot);
    Replicate3_normalised_raw_data=Finalised_Master_Gaussian_list.normaliased_log2_comparsion_replicate_3(Gaussian_counter1,1:number_of_gaussian_to_plot);
    Replicate4_normalised_raw_data=Finalised_Master_Gaussian_list.normaliased_log2_comparsion_replicate_4(Gaussian_counter1,1:number_of_gaussian_to_plot);
    
    while ~isempty(Center_to_plot)
      %Find minimum center to plot
      [index_minimum,index_minimum]=min(Center_to_plot);
      %replicate 1
      Values_for_bar_graph(bar_position(Center_counter1),1)=Replicate1_normalised_raw_data(index_minimum);
      Bar_bin_name_replicate_1=strcat('G_',mat2str(Gaussian_number_counter1),'_R_1');
      Names_for_bar_graph{bar_position(Center_counter1),1}=Bar_bin_name_replicate_1;
      
      if replicate_num >=2
        %replicate 2
        Values_for_bar_graph(bar_position(Center_counter1)+1,1)=Replicate2_normalised_raw_data(index_minimum);
        Bar_bin_name_replicate_2=strcat('G_',mat2str(Gaussian_number_counter1),'_R_2');
        Names_for_bar_graph{bar_position(Center_counter1)+1,1}=Bar_bin_name_replicate_2;
      end
      
      if replicate_num >=3
        %replicate 3
        Values_for_bar_graph(bar_position(Center_counter1)+2,1)=Replicate3_normalised_raw_data(index_minimum);
        Bar_bin_name_replicate_3=strcat('G_',mat2str(Gaussian_number_counter1),'_R_3');
        Names_for_bar_graph{bar_position(Center_counter1)+2,1}=Bar_bin_name_replicate_3;
      end
      
      if replicate_num ==4
        %replicate 4
        Values_for_bar_graph(bar_position(Center_counter1)+3,1)=Replicate4_normalised_raw_data(index_minimum);
        Bar_bin_name_replicate_4=strcat('G_',mat2str(Gaussian_number_counter1),'_R_4');
        Names_for_bar_graph{bar_position(Center_counter1)+3,1}=Bar_bin_name_replicate_4;
      end
      
      %Remove values of data already plotted
      Replicate1_normalised_raw_data(index_minimum)=[];
      Replicate2_normalised_raw_data(index_minimum)=[];
      Replicate3_normalised_raw_data(index_minimum)=[];
      Replicate4_normalised_raw_data(index_minimum)=[];
      Center_to_plot(index_minimum)=[];
      
      Center_counter1=Center_counter1+1;
      Gaussian_number_counter1=Gaussian_number_counter1+1;
    end
    
    
    P8= bar(Values_for_bar_graph);
    title_name_bar=strcat('Log2 ratio of treated to untreated at gaussian apex of :',Protein_to_plot);
    hold on %add 1 dotted line
    PP1= plot([0 (length(Values_for_bar_graph)+1)],[1 1],':');
    set(PP1,'Color','black','LineWidth',1);
    hold on %add 1 dotted line
    PP2= plot([0 (length(Values_for_bar_graph)+1)],[-1 -1],':');
    set(PP2,'Color','black','LineWidth',1);
    title(title_name_bar,'FontSize', 12);
    ylabel('Log2 ratio','FontSize', 10);
    xlabel('Fractions','FontSize', 10);
    set(gca, 'XTickLabel',Names_for_bar_graph, 'XTick',1:numel(Names_for_bar_graph),'FontSize', 10/(number_of_gaussian_to_plot*0.5));
    xlim([0 (length(Values_for_bar_graph)+1)]);
    
    %Save image
    Save_name_plot=strcat(figdir2, mat2str(Gaussian_counter1),'_2_PCP_SEC_Profiles of of_',Protein_to_plot,'.png');
    List_of_pdf{List_of_pdf_counter,1}=Save_name_plot;
    List_of_pdf{List_of_pdf_counter,2}=List_of_pdf_counter;
    List_of_pdf_counter=List_of_pdf_counter+1;
    saveas(f9, Save_name_plot);
    %print('-dpdf', '-r600', Save_name_plot);
    close 'all';
  end
  
end


% 4. Make pie chart figure



%Create array to store which replicate each protein
Proteins_observed_in_replicates=zeros(number_of_unique_protein_with_gaussians, 8);

for count_shared_guassians3=1:number_of_unique_protein_with_gaussians
  [internal_location_of_protein_of_interest_MvsL]=ind2sub(Number_All_guassians_identified_nonredunent_list(1), strmatch(Unique_protein_names(count_shared_guassians3), MvsL_Gaussians.Protein_name, 'exact'));
  [internal_location_of_protein_of_interest_HvsL]=ind2sub(Number_All_guassians_identified_nonredunent_list(1), strmatch(Unique_protein_names(count_shared_guassians3), HvsL_Gaussians.Protein_name, 'exact'));
  %Record which repliates the protein of interest was seen in MvsL replicates
  Number_of_times_seen_MvsL=length(internal_location_of_protein_of_interest_MvsL);
  for Number_of_times_seen_counter1= 1:Number_of_times_seen_MvsL
    if MvsL_Gaussians.Replicate(internal_location_of_protein_of_interest_MvsL(Number_of_times_seen_counter1)) == 1;
      Proteins_observed_in_replicates(count_shared_guassians3,1)=1;
    elseif MvsL_Gaussians.Replicate(internal_location_of_protein_of_interest_MvsL(Number_of_times_seen_counter1)) == 2;
      Proteins_observed_in_replicates(count_shared_guassians3,2)=1;
    elseif MvsL_Gaussians.Replicate(internal_location_of_protein_of_interest_MvsL(Number_of_times_seen_counter1)) == 3;
      Proteins_observed_in_replicates(count_shared_guassians3,3)=1;
    elseif MvsL_Gaussians.Replicate(internal_location_of_protein_of_interest_MvsL(Number_of_times_seen_counter1)) == 4;
      Proteins_observed_in_replicates(count_shared_guassians3,4)=1;
    elseif MvsL_Gaussians.Replicate(internal_location_of_protein_of_interest_MvsL(Number_of_times_seen_counter1)) == 0;
      Proteins_observed_in_replicates(count_shared_guassians3,4)=0;
      Proteins_observed_in_replicates(count_shared_guassians3,3)=0;
      Proteins_observed_in_replicates(count_shared_guassians3,2)=0;
      Proteins_observed_in_replicates(count_shared_guassians3,1)=0;
    end
  end
  
  %Record which repliates the protein of interest was seen in HvsL replicates
  Number_of_times_seen_HvsL=length(internal_location_of_protein_of_interest_HvsL);
  for Number_of_times_seen_counter2= 1:Number_of_times_seen_HvsL
    if HvsL_Gaussians.Replicate(internal_location_of_protein_of_interest_HvsL(Number_of_times_seen_counter2)) == 1;
      Proteins_observed_in_replicates(count_shared_guassians3,5)=1;
    elseif HvsL_Gaussians.Replicate(internal_location_of_protein_of_interest_HvsL(Number_of_times_seen_counter2)) == 2;
      Proteins_observed_in_replicates(count_shared_guassians3,6)=1;
    elseif HvsL_Gaussians.Replicate(internal_location_of_protein_of_interest_HvsL(Number_of_times_seen_counter2)) == 3;
      Proteins_observed_in_replicates(count_shared_guassians3,7)=1;
    elseif HvsL_Gaussians.Replicate(internal_location_of_protein_of_interest_HvsL(Number_of_times_seen_counter2)) == 4;
      Proteins_observed_in_replicates(count_shared_guassians3,8)=1;
    elseif HvsL_Gaussians.Replicate(internal_location_of_protein_of_interest_HvsL(Number_of_times_seen_counter2)) == 0;
      Proteins_observed_in_replicates(count_shared_guassians3,5)=0;
      Proteins_observed_in_replicates(count_shared_guassians3,6)=0;
      Proteins_observed_in_replicates(count_shared_guassians3,7)=0;
      Proteins_observed_in_replicates(count_shared_guassians3,8)=0;
    end
  end
  
end

%Write Protein information to Protein structure array
Protein_information.Protein_name=Unique_protein_names;
Protein_information.Observed_MvsL_replicate1= Proteins_observed_in_replicates(:,1);
Protein_information.Observed_MvsL_replicate2= Proteins_observed_in_replicates(:,2);
Protein_information.Observed_MvsL_replicate3= Proteins_observed_in_replicates(:,3);
Protein_information.Observed_MvsL_replicate4= Proteins_observed_in_replicates(:,4);
Protein_information.Observed_HvsL_replicate1= Proteins_observed_in_replicates(:,5);
Protein_information.Observed_HvsL_replicate2= Proteins_observed_in_replicates(:,6);
Protein_information.Observed_HvsL_replicate3= Proteins_observed_in_replicates(:,7);
Protein_information.Observed_HvsL_replicate4= Proteins_observed_in_replicates(:,8);

%Number of times observed in specific channels
Protein_information.Total_number_channels_observed= (sum(Proteins_observed_in_replicates(:,:)'))';
Protein_information.Observed_in_MvsL= (sum(Proteins_observed_in_replicates(:,1:4)')>0)';
Protein_information.Observed_in_HvsL= (sum(Proteins_observed_in_replicates(:,5:8)')>0)';
Protein_information.Number_Observed_in_MvsL= (sum(Proteins_observed_in_replicates(:,1:4)'))';
Protein_information.Number_Observed_in_HvsL= (sum(Proteins_observed_in_replicates(:,5:8)'))';
Protein_information.Observed_only_HvsL= ((Protein_information.Observed_in_HvsL>0) & (Protein_information.Observed_in_MvsL==0));
Protein_information.Observed_only_MvsL= ((Protein_information.Observed_in_MvsL>0) & (Protein_information.Observed_in_HvsL==0));
Protein_information.Observed_both= ((Protein_information.Observed_in_MvsL>0) & (Protein_information.Observed_in_HvsL>0));

%Number of time observed in replicates
Protein_information.Replicate1=Protein_information.Observed_MvsL_replicate1 & Protein_information.Observed_HvsL_replicate1;
Protein_information.Replicate2=Protein_information.Observed_MvsL_replicate2 & Protein_information.Observed_HvsL_replicate2;
Protein_information.Replicate3=Protein_information.Observed_MvsL_replicate3 & Protein_information.Observed_HvsL_replicate3;
Protein_information.Replicate4=Protein_information.Observed_MvsL_replicate4 & Protein_information.Observed_HvsL_replicate4;

%Number of times observed between replicates
Protein_information.Observed_in_zero_HvsL= (Protein_information.Number_Observed_in_HvsL==0);
Protein_information.Observed_in_zero_MvsL= (Protein_information.Number_Observed_in_MvsL==0);
Protein_information.Observed_in_one_HvsL= (Protein_information.Number_Observed_in_HvsL==1);
Protein_information.Observed_in_one_MvsL= (Protein_information.Number_Observed_in_MvsL==1);
Protein_information.Observed_in_two_HvsL= (Protein_information.Number_Observed_in_HvsL==2);
Protein_information.Observed_in_two_MvsL= (Protein_information.Number_Observed_in_MvsL==2);
Protein_information.Observed_in_three_HvsL= (Protein_information.Number_Observed_in_HvsL==3);
Protein_information.Observed_in_three_MvsL= (Protein_information.Number_Observed_in_MvsL==3);
Protein_information.Observed_in_four_HvsL= (Protein_information.Number_Observed_in_HvsL==4);
Protein_information.Observed_in_four_MvsL= (Protein_information.Number_Observed_in_MvsL==4);

%Number of times observed between replicates and channels
%HvsL
Protein_information.Observed_in_one_HvsL_both_channels= Protein_information.Observed_in_one_HvsL & Protein_information.Observed_both;
Protein_information.Observed_in_two_HvsL_both_channels= Protein_information.Observed_in_two_HvsL & Protein_information.Observed_both;
Protein_information.Observed_in_three_HvsL_both_channels= Protein_information.Observed_in_three_HvsL & Protein_information.Observed_both;
Protein_information.Observed_in_four_HvsL_both_channels= Protein_information.Observed_in_four_HvsL & Protein_information.Observed_both;
%MvsL
Protein_information.Observed_in_one_MvsL_both_channels= Protein_information.Observed_in_one_MvsL & Protein_information.Observed_both;
Protein_information.Observed_in_two_MvsL_both_channels= Protein_information.Observed_in_two_MvsL & Protein_information.Observed_both;
Protein_information.Observed_in_three_MvsL_both_channels= Protein_information.Observed_in_three_MvsL & Protein_information.Observed_both;
Protein_information.Observed_in_four_MvsL_both_channels= Protein_information.Observed_in_four_MvsL & Protein_information.Observed_both;
%HvsL
Protein_information.Observed_in_zero_MvsL_HvsL_channels= Protein_information.Observed_in_zero_MvsL & Protein_information.Observed_only_HvsL;
Protein_information.Observed_in_one_HvsL_HvsL_channels= Protein_information.Observed_in_one_HvsL & Protein_information.Observed_only_HvsL;
Protein_information.Observed_in_two_HvsL_HvsL_channels= Protein_information.Observed_in_two_HvsL & Protein_information.Observed_only_HvsL;
Protein_information.Observed_in_three_HvsL_HvsL_channels= Protein_information.Observed_in_three_HvsL & Protein_information.Observed_only_HvsL;
Protein_information.Observed_in_four_HvsL_HvsL_channels= Protein_information.Observed_in_four_HvsL & Protein_information.Observed_only_HvsL;
%MvsL
Protein_information.Observed_in_zero_HvsL_MvsL_channels= Protein_information.Observed_in_zero_HvsL & Protein_information.Observed_only_MvsL;
Protein_information.Observed_in_one_MvsL_MvsL_channels= Protein_information.Observed_in_one_MvsL & Protein_information.Observed_only_MvsL;
Protein_information.Observed_in_two_MvsL_MvsL_channels= Protein_information.Observed_in_two_MvsL & Protein_information.Observed_only_MvsL;
Protein_information.Observed_in_three_MvsL_MvsL_channels= Protein_information.Observed_in_three_MvsL & Protein_information.Observed_only_MvsL;
Protein_information.Observed_in_four_MvsL_MvsL_channels= Protein_information.Observed_in_four_MvsL & Protein_information.Observed_only_MvsL;

%Determine overlap of proteins identifications between replicates

%Compare how many proteins were observed in each channel (use for Figure 1)
pie_figure1(1,1)=nnz(Protein_information.Total_number_channels_observed(Protein_information.Total_number_channels_observed==1));
pie_figure1(1,2)=nnz(Protein_information.Total_number_channels_observed(Protein_information.Total_number_channels_observed==2));
pie_figure1(1,3)=nnz(Protein_information.Total_number_channels_observed(Protein_information.Total_number_channels_observed==3));
pie_figure1(1,4)=nnz(Protein_information.Total_number_channels_observed(Protein_information.Total_number_channels_observed==4));
pie_figure1(1,5)=nnz(Protein_information.Total_number_channels_observed(Protein_information.Total_number_channels_observed==5));
pie_figure1(1,6)=nnz(Protein_information.Total_number_channels_observed(Protein_information.Total_number_channels_observed==6));
pie_figure1(1,7)=nnz(Protein_information.Total_number_channels_observed(Protein_information.Total_number_channels_observed==7));
pie_figure1(1,8)=nnz(Protein_information.Total_number_channels_observed(Protein_information.Total_number_channels_observed==8));

%Comparing what channels the proteins were observed in across all
%replicates
pie_figure2(1,1)=nnz(Protein_information.Observed_both);
pie_figure2(1,2)=nnz(Protein_information.Observed_only_HvsL);
pie_figure2(1,3)=nnz(Protein_information.Observed_only_MvsL);

%Compare what the number of times observed between replicates (Figure 3)
pie_figure3(1,1)=nnz(Protein_information.Observed_in_zero_MvsL);
pie_figure3(1,3)=nnz(Protein_information.Observed_in_one_MvsL);
pie_figure3(1,5)=nnz(Protein_information.Observed_in_two_MvsL);
pie_figure3(1,7)=nnz(Protein_information.Observed_in_three_MvsL);
pie_figure3(1,9)=nnz(Protein_information.Observed_in_four_MvsL);
pie_figure3(1,2)=nnz(Protein_information.Observed_in_zero_HvsL);
pie_figure3(1,4)=nnz(Protein_information.Observed_in_one_HvsL);
pie_figure3(1,6)=nnz(Protein_information.Observed_in_two_HvsL);
pie_figure3(1,8)=nnz(Protein_information.Observed_in_three_HvsL);
pie_figure3(1,10)=nnz(Protein_information.Observed_in_four_HvsL);

%Compare the Number of times observed between replicates and channels
%(Figure 4)
pie_figure4(1,1)=nnz(Protein_information.Observed_in_zero_HvsL_MvsL_channels);
pie_figure4(1,3)=nnz(Protein_information.Observed_in_one_HvsL_HvsL_channels);
pie_figure4(1,5)=nnz(Protein_information.Observed_in_two_HvsL_HvsL_channels);
pie_figure4(1,2)=nnz(Protein_information.Observed_in_zero_MvsL_HvsL_channels);
pie_figure4(1,4)=nnz(Protein_information.Observed_in_one_MvsL_MvsL_channels);
pie_figure4(1,6)=nnz(Protein_information.Observed_in_two_MvsL_MvsL_channels);
pie_figure4(1,7)=nnz(Protein_information.Observed_in_three_HvsL_HvsL_channels);
pie_figure4(1,8)=nnz(Protein_information.Observed_in_three_MvsL_MvsL_channels);
pie_figure4(1,9)=nnz(Protein_information.Observed_in_four_HvsL_HvsL_channels);
pie_figure4(1,10)=nnz(Protein_information.Observed_in_four_MvsL_MvsL_channels);

%(Figure 5)
pie_figure5(1,1)=nnz(Protein_information.Observed_in_one_HvsL_both_channels);
pie_figure5(1,3)=nnz(Protein_information.Observed_in_two_HvsL_both_channels);
pie_figure5(1,5)=nnz(Protein_information.Observed_in_three_HvsL_both_channels);
pie_figure5(1,7)=nnz(Protein_information.Observed_in_four_HvsL_both_channels);
pie_figure5(1,2)=nnz(Protein_information.Observed_in_one_MvsL_both_channels);
pie_figure5(1,4)=nnz(Protein_information.Observed_in_two_MvsL_both_channels);
pie_figure5(1,6)=nnz(Protein_information.Observed_in_three_MvsL_both_channels);
pie_figure5(1,8)=nnz(Protein_information.Observed_in_four_MvsL_both_channels);


%Compare what the number of times proteins observed between replicates (Figure 3)
pie_figure6(1,1)=nnz(Protein_information.Replicate1);
pie_figure6(1,2)=nnz(Protein_information.Replicate2);
pie_figure6(1,3)=nnz(Protein_information.Replicate3);
pie_figure6(1,4)=nnz(Protein_information.Replicate4);

% #output
%write out list of which proteins were observed (as gaussians in each replicate)
fid_Proteins_in_each_rep = fopen([datadir1 'Protein_gaussian_observed_in_each_replicate.csv'],'w');
fprintf (fid_Proteins_in_each_rep,'%s,%s,%s,%s,%s,%s,%s,%s,%s,\n',...
  'Protein_name','Replicate 1 MvsL channel','Replicate 2 MvsL channel','Replicate 3 MvsL channel','Replicate 4 MvsL channel','Replicate 1 HvsL channel','Replicate 2 HvsL channel','Replicate 3 HvsL channel','Replicate 4 HvsL channel'); %Write Header
for write_out_Proteins_in_rep=1:number_of_unique_protein_with_gaussians(1)
  fprintf(fid_Proteins_in_each_rep,'%s,', Unique_protein_names{write_out_Proteins_in_rep});
  fprintf(fid_Proteins_in_each_rep,'%6.4g,', Proteins_observed_in_replicates(write_out_Proteins_in_rep,:));
  fprintf(fid_Proteins_in_each_rep,'\n');
end
fclose(fid_Proteins_in_each_rep);

%write out analysis of observation
fid_Proteins_in_each_rep2 = fopen([datadir1 'Protein_observed_in_each_replicate.csv'],'w');
fprintf (fid_Proteins_in_each_rep2,'%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,\n',...
  'Protein_Observed_in_both_channels','Protein_Observed_only_in_MvsL_channel','Protein_Observed_only_in_HvsL_channel',...
  'Protein_Observed_in_MvsL_once_channels','Protein_Observed_in_MvsL_twice_channel','Protein_Observedin_MvsL_three_times_channel','Protein_Observedin_MvsL_four_times_channel',...
  'Protein_Observed_in_HvsL_once_channels','Protein_Observed_in_HvsL_twice_channel','Protein_Observedin_HvsL_three_times_channel','Protein_Observedin_HvsL_four_times_channel'); %Write Header
fprintf(fid_Proteins_in_each_rep2,'%6.4f,', pie_figure2(1,1));
fprintf(fid_Proteins_in_each_rep2,'%6.4f,', pie_figure2(1,3));
fprintf(fid_Proteins_in_each_rep2,'%6.4f,', pie_figure2(1,2));
fprintf(fid_Proteins_in_each_rep2,'%6.4f,', pie_figure3(1,3));
fprintf(fid_Proteins_in_each_rep2,'%6.4f,', pie_figure3(1,5));
fprintf(fid_Proteins_in_each_rep2,'%6.4f,', pie_figure3(1,7));
fprintf(fid_Proteins_in_each_rep2,'%6.4f,', pie_figure3(1,9));
fprintf(fid_Proteins_in_each_rep2,'%6.4f,', pie_figure3(1,4));
fprintf(fid_Proteins_in_each_rep2,'%6.4f,', pie_figure3(1,6));
fprintf(fid_Proteins_in_each_rep2,'%6.4f,', pie_figure3(1,8));
fprintf(fid_Proteins_in_each_rep2,'%6.4f,', pie_figure3(1,10));
fclose(fid_Proteins_in_each_rep2);

%write out analysis of observation
fid_Proteins_in_each_rep3 = fopen([datadir1 'Protein_observed_in_each_replicate_and_channels.csv'],'w');
fprintf (fid_Proteins_in_each_rep3,'%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,\n',...
  'Protein_Observed_only_in_MvsL_channel','Protein_Observed_only_in_HvsL_channel','Protein_Observed_only_in_HvsL_channel_once',...
  'Protein_Observed_only_in_HvsL_channel_twice','Protein_Observed_only_in_HvsL_channel_three_times','Protein_Observed_only_in_HvsL_channel_four_times',...
  'Protein_Observed_only_in_MvsL_channel_once','Protein_Observed_only_in_MvsL_channesl_twice','Protein_Observed_only_in_MvsL_channel_three_times',...
  'Protein_Observed_only_in_MvsL_channel_four_times','Protein_Observed_Both_channesl_once_HvsL','Protein_Observed_Both_channesl_twice_HvsL',...
  'Protein_Observed_Both_channesl_three_times_HvsL','Protein_Observed_Both_channesl_four_times_HvsL','Protein_Observed_Both_channesl_once_MvsL',...
  'Protein_Observed_Both_channesl_twice_MvsL','Protein_Observed_Both_channesl_three_times_MvsL','Protein_Observed_Both_channesl_four_times_MvsL'); %Write Header
fprintf(fid_Proteins_in_each_rep3,'%6.4f,', pie_figure4(1,1));
fprintf(fid_Proteins_in_each_rep3,'%6.4f,', pie_figure4(1,2));
fprintf(fid_Proteins_in_each_rep3,'%6.4f,', pie_figure4(1,3));
fprintf(fid_Proteins_in_each_rep3,'%6.4f,', pie_figure4(1,5));
fprintf(fid_Proteins_in_each_rep3,'%6.4f,', pie_figure4(1,7));
fprintf(fid_Proteins_in_each_rep3,'%6.4f,', pie_figure4(1,9));
fprintf(fid_Proteins_in_each_rep3,'%6.4f,', pie_figure4(1,4));
fprintf(fid_Proteins_in_each_rep3,'%6.4f,', pie_figure4(1,6));
fprintf(fid_Proteins_in_each_rep3,'%6.4f,', pie_figure4(1,8));
fprintf(fid_Proteins_in_each_rep3,'%6.4f,', pie_figure4(1,10));
fprintf(fid_Proteins_in_each_rep3,'%6.4f,', pie_figure5(1,1));
fprintf(fid_Proteins_in_each_rep3,'%6.4f,', pie_figure5(1,3));
fprintf(fid_Proteins_in_each_rep3,'%6.4f,', pie_figure5(1,5));
fprintf(fid_Proteins_in_each_rep3,'%6.4f,', pie_figure5(1,7));
fprintf(fid_Proteins_in_each_rep3,'%6.4f,', pie_figure5(1,2));
fprintf(fid_Proteins_in_each_rep3,'%6.4f,', pie_figure5(1,4));
fprintf(fid_Proteins_in_each_rep3,'%6.4f,', pie_figure5(1,6));
fprintf(fid_Proteins_in_each_rep3,'%6.4f,', pie_figure5(1,8));
fclose(fid_Proteins_in_each_rep3);


%Figure parameters, Figure 1
%Note: Option are present to also output a visulation on the overlap
%between replicates and isotopologue channels but this has found to be
%largely uninformative
f1 = figure;
P(1)=subplot(1,1,1);
Figure1= pie(pie_figure1);
hp = findobj(P(1), 'Type', 'patch');
figure_pos = get(P(1),'position');

%create varible for legend
combinedstrings=cell(length(hp),1);

for colour_count=1:length(hp)
  set(hp(colour_count), 'FaceColor', colour_to_use(colour_count,:)); % set colour according to defined palate
  combinedstrings{colour_count,1}=mat2str(colour_count);
end

set(P(1), 'position',[figure_pos(1)*1.6 figure_pos(2)*1.0 figure_pos(3:4)]);

%Create label
hText = findobj(Figure1,'Type','text'); % text handles
percentValues = get(hText,'String'); % percent values

%Position text
textPositions_cell = get(hText,'Position'); % cell array
textPositions = cell2mat(textPositions_cell); % numeric array
textPositions(:,1) = textPositions(:,1)*1.2; % add offset
%set(hText,'Position',num2cell(textPositions,[3,2]));
for jj = 1:size(hText,1)
  set(hText(jj),'Position',textPositions(jj,:));
end

%# add legend
Legend_1_position = legend(P(1), combinedstrings);
Legend_pos = get(Legend_1_position,'position');
set(Legend_1_position, 'position',[0.05 0.25 0.15 0.5]);

saveas(f1, [figdir1 'Protein_analysis.png']);

tt = toc;
tt = round(tt/10)*10;
fprintf('  ...  %.2f seconds\n',tt)