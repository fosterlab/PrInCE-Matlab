%%%%%%%%%%%%%%% Logic:
% Since replicates can be run on different columns, etc., it's necessary to align replicates. Each
% replicate will be offset by a constant amount, so it's sufficient to find offset between a few
% landmarks
% i. Read chromatogram data (MaxQuant verbose output) and Gaussian fitting data.
% ii. Impute missing chromatogram values.
% iii. Find all proteins with a single Gaussian in all replicates.
% iv. Choose the best replicate, i.e. the one with the most overlap of single-Gaussian proteins with
%     the other two. Align to this replicate.
% v. For every other replicate, fit the line: y_mu = m * x_mu + b, where y_mu is the centers of the
%     best replicate, x_mu is the centers of the other replicate.
% vi. Scale the other replicates, x, by the above function. Use linear interpolation.
%
% NB1 - Best replicate doesn't seem to be chosen using SSE of fits.
% NB2 - x and y likely need to be rescaled to reflect fraction spacing (e.g. 20 in 10 minutes, etc.)
%     - That is, the spacing b/w fractions is not constant. That seems major, right?


%%%%%%%%%%%%%%% Small changes:
% 1. Make all file references and directories absolute. Remove 'cd'. Make this directory-independent!
% 2. Replace all '\' with '/'
% 3. Replace all 'if test==1' with 'if test'
% 4. 'for replicate_to_align_against= y_maximum' does not need a for loop.


%%%%%%%%%%%%%%% Big changes:
% 1. Impute multiple missing values (more than just one, as is currently being done)
% 2. Don't write csv files (these seem useful for debugging, but otherwise aren't needed).
% 3. Include an option not to plot (separate plotting script?).
% 4. Why does Gauss_Build clean the chromatograms more thoroughly? Here we just impute single
%    values. Should we do more?


%%%%%%%%%%%%%%% Questions:
% 1. Why do the chromatograms have 5 extra points on either side?
%
%

disp('Alignment.m')


%% 0. Initialize
tic
fprintf('\n    0. Initialize')

% User defined variables
User_alignment_window1 = 20; %for initial alignment with proteins fitted to a single Gaussian
User_alignment_window2 = 8; %for Second round of alignment
Experimental_channels = {'MvsL' 'HvsL'};
Number_of_experimental_channels = length(Experimental_channels);  % Defines the number of experiments to be compared
Nreplicates = 3; % number of replicates
Alignment_to_user= 'MvsL'; %Define the experimental channel to use for alignment
User_defined_zero_value = 0.2; %lowest value to be shown in Adjusted Chromatograms


% Define folders, i.e. define where everything lives.
maindir = '/Users/Mercy/Academics/Foster/NickCodeData/GregPCP-SILAC/'; % where everything lives
codedir = [maindir 'Code/']; % where this script lives
funcdir = [maindir 'Code/Functions/']; % where small pieces of code live
datadir = [maindir 'Data/']; % where data files live
datadir1 = [maindir 'Data/Alignment/']; % where data files live
datadir2 = [maindir 'Data/GaussBuild/']; % where data files live
figdir1 = [maindir 'Figures/Alignment/']; % where figures live
% Make folders if necessary
if ~exist(codedir, 'dir'); mkdir(codedir); end
if ~exist(funcdir, 'dir'); mkdir(funcdir); end
if ~exist(datadir, 'dir'); mkdir(datadir); end
if ~exist(datadir1, 'dir'); mkdir(datadir1); end
if ~exist(datadir2, 'dir'); error('\nData files from Gauss_Build.m not found\n'); end
if ~exist(figdir1, 'dir'); mkdir(figdir1); end


% List all input files. These contain data that will be read by this script.
clear InputFile
InputFile{1} = [datadir 'Combined_replicates_2014_04_22_contaminates_removed_for_MvsL_scripts.xlsx'];
InputFile{2} = [datadir 'Combined_replicates_2014_04_22_contaminates_removed_for_HvsL_scripts.xlsx'];
InputFile{3} = [datadir 'Combined_replicates_2014_04_22_contaminates_removed_for_HvsM_scripts.xlsx'];
InputFile{4} = [datadir 'SEC_alignment.xlsx'];

% for now, read these in from Nick's data
% in the future read them in from my data
if 1 %nick's output
  pw = '/Users/Mercy/Academics/Foster/NickCodeData/2_Alignment processing/MvsL/';
  InputFile{5} = '/Users/Mercy/Academics/Foster/NickCodeData/2_Alignment processing/MvsL/MvsL_Combined_OutputGaus.csv'; % From Gauss_build
  InputFile{6} = '/Users/Mercy/Academics/Foster/NickCodeData/2_Alignment processing/MvsL/MvsL_Summary_Gausians_for_individual_proteins.csv';   % ''
  InputFile{7} = '/Users/Mercy/Academics/Foster/NickCodeData/2_Alignment processing/HvsL/HvsL_Combined_OutputGaus.csv';  % ''
  InputFile{8} = '/Users/Mercy/Academics/Foster/NickCodeData/2_Alignment processing/HvsL/HvsL_Summary_Gausians_for_individual_proteins.csv'; % ''
  GaussInputFile = cell(Number_of_experimental_channels, Nreplicates);
  GassSumInputFile = cell(Number_of_experimental_channels, Nreplicates);
  for ei=1:Number_of_experimental_channels
    tmp = Experimental_channels{ei};
    for replicates= 1:Nreplicates
      GaussInputFile{ei,replicates} = ['/Users/Mercy/Academics/Foster/NickCodeData/2_Alignment processing/' tmp '_alignment/Processed Gaussian/' tmp '_Combined_OutputGaus_rep' num2str(replicates) '.csv'];
      GassSumInputFile{ei,replicates} = ['/Users/Mercy/Academics/Foster/NickCodeData/2_Alignment processing/' tmp '_alignment/Processed Gaussian/' tmp '_Summary_Gausians_for_individual_proteins_rep' num2str(replicates) '.csv'];
    end
  end
  
else %my output
  InputFile{5} = [datadir 'MvsL_Combined_OutputGaus.csv'];                        % From Gauss_build
  InputFile{6} = [datadir 'MvsL_Summary_Gausians_for_individual_proteins.csv'];   % ''
  InputFile{7} = [datadir 'HvsL_Combined_OutputGaus.csv'];                        % ''
  InputFile{8} = [datadir 'HvsL_Summary_Gausians_for_individual_proteins.csv'];   % ''
  GaussInputFile = cell(Number_of_experimental_channels, Nreplicates);
  GassSumInputFile = cell(Number_of_experimental_channels, Nreplicates);
  for ei=1:Number_of_experimental_channels
    for replicates= 1:Nreplicates
      GaussInputFile{ei,replicates} = [datadir2 Experimental_channels{ei} '_Combined_OutputGaus_rep' num2str(replicates) '.csv'];
      GassSumInputFile{ei,replicates} = [datadir2 Experimental_channels{ei} '_Summary_Gausians_for_individual_proteins_rep' num2str(replicates) '.csv'];
    end
  end
  
end
% *vsL_Combined_OutputGaus.csv: data on all the fitted Gaussians
% *vsL_Summary_Gausians_for_individual_proteins.csv: how many Gaussians were fitted per protein


tt = toc;
fprintf('  ...  %.2f seconds\n',tt)



%% 1. Read input
tic
fprintf('    1. Read input')

% Import MaxQuant data files
[num_val_MvsL,txt_MvsL] = xlsread(InputFile{1}); %Import file raw Maxqaunt output
[num_val_HvsL,txt_HvsL] = xlsread(InputFile{2}); %Import file raw Maxqaunt output
[num_val_HvsM,txt_HvsM] = xlsread(InputFile{3}); %Import file raw Maxqaunt output
[SEC_size_alignment]=xlsread(InputFile{4});

% Import Gaussian fits
f1 = fopen(InputFile{5});
Processed_data1{1} = textscan(f1, '%s', 'Delimiter',','); %#
fclose(f1);
f2 = fopen(InputFile{6});
Processed_data2{1} = textscan(f2, '%s', 'Delimiter',','); %#
fclose(f2);
f1 = fopen(InputFile{7});
Processed_data1{2} = textscan(f1, '%s', 'Delimiter',','); %#
fclose(f1);
f2 = fopen(InputFile{8});
Processed_data2{2} = textscan(f2, '%s', 'Delimiter',','); %#
fclose(f2);

% Import Gauss fits for each replicate
%   Gaus_import: mx6, where m is the number of proteins with a fitted Gaussian
%   Summary_gausian_infomration: nx6, where n is the unique protein number (1-3217)
Gaus_import = cell(Number_of_experimental_channels, Nreplicates);
Summary_gausian_infomration = cell(Number_of_experimental_channels, Nreplicates);
for ci = 1:Number_of_experimental_channels
  for replicates= 1:Nreplicates
    Gaus_import{ci,replicates} = importdata(GaussInputFile{ci,replicates});
    Summary_gausian_infomration{ci,replicates} = importdata(GassSumInputFile{ci,replicates});
  end
end

% Calibration
SEC_fit=polyfit(SEC_size_alignment(1,:),SEC_size_alignment(2,:),1);

%Number of fractions
fraction_number=size(num_val_MvsL);
fraction_number(2) = fraction_number(2)-1;

% number of fractions,equal to # fractions + 10
Chromatograms_axis = fraction_number(2)+10;

% Clean chromatograms
for ri = 1:size(num_val_MvsL,1) % loop over proteins
  num_val_MvsL(ri,:) = cleanChromatogram(num_val_MvsL(ri,:),[1 3]);
  num_val_HvsL(ri,:) = cleanChromatogram(num_val_HvsL(ri,:),[1 3]);
  num_val_HvsM(ri,:) = cleanChromatogram(num_val_HvsM(ri,:),[1 3]);
end

tt = toc;
fprintf('  ...  %.2f seconds\n',tt)



%% 2. Find the best replicates to align to
% i) Find names of proteins with a single Gaussian
% ii) Find the overlap in single Gaussians b/w replicates
% iii) Choose the replicate with maximum overlap

replicate_to_align_against = nan(Number_of_experimental_channels,1);
summerised_names_G1 = cell(Number_of_experimental_channels,Nreplicates); % proteins with 1 Gaussian
for ci = 1:Number_of_experimental_channels
  
  % i) Find names of proteins with a single Gaussian in each replicate
  %summerised_protein_number_G1 = cell(Nreplicates,1);
  for rr= 1:Nreplicates
    Ngauss = size(Summary_gausian_infomration{ci,rr},1);
    
    Inotsingle = find(Summary_gausian_infomration{ci,rr}.data(:,2) ~= 1) + 1;
    %summerised_protein_number_G1{rr} = Summary_gausian_infomration{ci,rr}.textdata(Isingle+1,1);
    %summerised_protein_number_G1{rr} = cellfun(@str2num,summerised_protein_number_G1{rr});
    %summerised_names_G1{ei,rr} = Summary_gausian_infomration{ci,rr}.textdata(Isingle+1,2);
    summerised_names_G1{ci,rr} = Summary_gausian_infomration{ci,rr}.textdata(:,2);
    summerised_names_G1{ci,rr}(Inotsingle) = {'-1'};
  end
  
  % ii) Find the overlap in single Gaussians b/w replicates
  Nintersect = zeros(Nreplicates,Nreplicates);
  for rr1 = 1:Nreplicates
    for rr2 = 1:Nreplicates
      Nintersect(rr1,rr2) = length(intersect(summerised_names_G1{ci,rr1},summerised_names_G1{ci,rr2}));
    end
  end
  Nintersect(eye(Nreplicates)==1) = 0; % remove diagonal, i.e. self-overlap
  Nintersect = sum(Nintersect); % collapse to one dimension
  
  [~,replicate_to_align_against(ci)] = max(Nintersect); % ding ding ding!
end



%% 3. Calculate best fit lines for adjustment
% i) Find the overlapping proteins
% ii) Find the Centers of these proteins
% iii) Fit a line

pfit = nan(Number_of_experimental_channels,Nreplicates,2);

for ci = 1:Number_of_experimental_channels
  align_rep = replicate_to_align_against(ci);
  for rr = 1:Nreplicates
    
    rep_to_align = rr;
    
    % i) Find the overlapping proteins
    overlap = intersect(summerised_names_G1{ci,align_rep},summerised_names_G1{ci,rep_to_align});
    overlap([1 2]) = [];
    
    % ii) Find their centers
    Ia = find(ismember(Gaus_import{ci,align_rep}.textdata(:,4),overlap));
    Ib = find(ismember(Gaus_import{ci,rep_to_align}.textdata(:,4),overlap));
    Ca = Gaus_import{ci,align_rep}.data(Ia-1,2); % align to this replicate, x
    Cb = Gaus_import{ci,rep_to_align}.data(Ib-1,2);% align this replicate, y
    
    % quality control...
    if length(overlap)~=length(Ia) || length(overlap)~=length(Ib)
      error('Error: Alignment.m: uhhhh')
    end
    
    % iii) Fit a line
    I = abs(Ca - Cb)<User_alignment_window1;
    pfit(ci,rr,:) = robustfit(Ca(I), Cb(I));
    
  end
end




%% 4. Using fitted curves, adjust replicate data
% i) Gaussian fits: shift Center
% ii) Chromatograms: shift data points

Adjusted_Gaus_import = Gaus_import;

for ci = 1:Number_of_experimental_channels
  align_rep = replicate_to_align_against(ci);
  for rr = 1:Nreplicates
    
    b = pfit(ci,rr,1); % intercept
    m = pfit(ci,rr,2); % slope
    
    % i) Gaussian fits: shift Center
    adjustedCenters = (Gaus_import{ci,rr}.data(:,2) - b) / m;
    Adjusted_Gaus_import{ci,rr}.data(:,2) = adjustedCenters;    
  end  
end

% ii) Chromatograms: shift data points

% MvsL
x = 1:55;
adjusted_raw_data{1} = nan(size(num_val_MvsL));
for ri=1:size(num_val_MvsL)
  y = num_val_MvsL(ri,2:end);
  
  rr = num_val_MvsL(ri,1);
  b = pfit(ci,rr,1); % intercept
  m = pfit(ci,rr,2); % slope
  
  x2 = (x-b)/m;
  
  adjusted_raw_data{1}(ri,2:end) = interp1(x,y,x2);
end



%% 5. Write output
%    fid9_Name = strcat('Adjusted_',Experimental_channel,'_Raw_data_maxquant_rep',mat2str(alignment_counter),'.csv');
%    fid9B_Name = strcat('Adjusted_',Experimental_channel,'_Raw_for_ROC_analysis_rep',mat2str(alignment_counter),'.csv');
%    fid9B_Name = strcat('Adjusted_HvsM_Raw_data_maxquant_rep',mat2str(alignment_counter),'.csv');
%    fid6_Name = strcat('Adjusted_Chromatograms_vobose_rep',mat2str(alignment_counter),'_','.csv');
%    fid7_Name = strcat('Adjusted_Combined_OutputGaus_rep',mat2str(alignment_counter),'.csv');
%  fid10_Name = strcat('Adjusted_',Experimental_channel,'_Combined_OutputGaus.csv');
%  fid11_Name = strcat('Adjusted_',Experimental_channel,'_Raw_data_maxquant.csv');




%% 6. Make figures



%% Is each G1 name really a G1?

clear A1 A2 N1 N2

listofnames = Summary_gausian_infomration{ci,rr}.textdata(:,2);
kk=0;
overlap = intersect(summerised_names_G1{1,2},summerised_names_G1{1,3});
for ri = 2:length(overlap)
  
  protname = overlap{ri};
  if strmatch('-1',protname);continue;end
  kk = kk+1;

  % how many G?
  I = strmatch(protname,listofnames,'exact');
  A1(kk,:) = Summary_gausian_infomration{1,2}.data(I-1,:);
  A2(kk,:) = Summary_gausian_infomration{1,3}.data(I-1,:);
  
  % how many times does this name show up in Gaus_import?
  N1(kk) = length(strmatch(protname,Gaus_import{1,2}.textdata(:,4),'exact'));
  N2(kk) = length(strmatch(protname,Gaus_import{1,3}.textdata(:,4),'exact'));
  
end

