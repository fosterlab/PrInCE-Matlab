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


%%%%%%%%%%%%%%% BUG BUG BUG BUG UBGSS BUGGGGGSSS!!!!
% 1. My Alignment numbers are off. Nick's are closer to the raw data (mine are quite far away).
% 2. I'm adding too many leading nans (two too many?).


%%%%%%%%%%%%%%% Questions:
% 1. Why do the chromatograms have 5 extra points on either side?
%
%

disp('Alignment.m')


%% 0. Initialize
tic
fprintf('\n    0. Initialize')


% Load user settings
maindir = user.maindir;
Experimental_channels = user.silacratios;
User_alignment_window1 = user.User_alignment_window1;
clear InputFile
for ii = 1:length(user.MQfiles)
  InputFile{1} = user.MQfiles{1};
  InputFile{2} = user.MQfiles{2};
end
InputFile{4} = user.calfile;



%User_alignment_window2 = 8; %for Second round of alignment
Number_of_experimental_channels = length(Experimental_channels);  % Defines the number of experiments to be compared
%Alignment_to_user= 'MvsL'; %Define the experimental channel to use for alignment
%User_defined_zero_value = 0.2; %lowest value to be shown in Adjusted Chromatograms


% Define folders, i.e. define where everything lives.
%maindir = '/Users/Mercy/Academics/Foster/NickCodeData/GregPCP-SILAC/'; % where everything lives
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
% InputFile{1} = [datadir 'Combined_replicates_2014_04_22_contaminates_removed_for_MvsL_scripts.xlsx'];
% InputFile{2} = [datadir 'Combined_replicates_2014_04_22_contaminates_removed_for_HvsL_scripts.xlsx'];
% InputFile{3} = [datadir 'Combined_replicates_2014_04_22_contaminates_removed_for_HvsM_scripts.xlsx'];
% InputFile{4} = [datadir 'SEC_alignment.xlsx'];
flag3 = 1;
try ls(InputFile{3});
catch 
  flag3 = 0;
end


% for now, read files in from Nick's data
% in the future read them in from my data
if user.nickflag %nick's output
  pw = '/Users/Mercy/Academics/Foster/NickCodeData/2_Alignment processing/MvsL/';
  Nreplicates = 3;
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
  dd = dir([datadir2 '*Summary_Gausians_for_individual_proteins_rep*csv']);
  Nreplicates = length(dd) / length(Experimental_channels);
  InputFile{5} = [datadir2 'MvsL_Combined_OutputGaus.csv'];                        % From Gauss_build
  InputFile{6} = [datadir2 'MvsL_Summary_Gausians_for_individual_proteins.csv'];   % ''
  InputFile{7} = [datadir2 'HvsL_Combined_OutputGaus.csv'];                        % ''
  InputFile{8} = [datadir2 'HvsL_Summary_Gausians_for_individual_proteins.csv'];   % ''
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
if flag3; [num_val_HvsM,txt_HvsM] = xlsread(InputFile{3});end %Import file raw Maxqaunt output
[SEC_size_alignment] = xlsread(InputFile{4});

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

% make protgausI
protgausI = cell(Number_of_experimental_channels, Nreplicates);
for ci = 1:Number_of_experimental_channels
  for rr = 1:Nreplicates
    protgausI{ci,rr} = nan(size(Gaus_import{ci,rr}.textdata,1)-1);
    for gi = 1:size(Gaus_import{ci,rr}.textdata,1)-1
      protgausI{ci,rr}(gi) = str2num(Gaus_import{ci,rr}.textdata{gi+1,3});
    end
  end
end

% Calibration
SEC_fit=polyfit(SEC_size_alignment(1,:),SEC_size_alignment(2,:),1);

%Number of fractions
%fraction_number=size(num_val_MvsL);
%fraction_number(2) = fraction_number(2)-1;
[~, fraction_number]=size(num_val_MvsL);
fraction_number = fraction_number-1;

% number of fractions,equal to # fractions + 10
Chromatograms_axis = fraction_number+10;

% replicates
replicates =  num_val_MvsL(:,1);

% Clean chromatograms
cleandata{1} = nan(size(num_val_MvsL,1),size(num_val_MvsL,2)+10);
cleandata{2} = nan(size(num_val_HvsL,1),size(num_val_HvsL,2)+10);
if flag3; cleandata{3} = nan(size(num_val_HvsM,1),size(num_val_HvsM,2)+10);end
cleandata{1}(:,1) = num_val_MvsL(:,1);
cleandata{2}(:,2) = num_val_HvsL(:,1);
if flag3; cleandata{3}(:,3) = num_val_HvsM(:,1);end
for ri = 1:size(num_val_MvsL,1) % loop over proteins
  cleandata{1}(ri,2:end) = cleanChromatogram(num_val_MvsL(ri,2:end),[1 3 5 7]);
  cleandata{2}(ri,2:end) = cleanChromatogram(num_val_HvsL(ri,2:end),[1 3 5 7]);
  if flag3; cleandata{3}(ri,2:end) = cleanChromatogram(num_val_HvsM(ri,2:end),[1 3 5 7]);end
end

tt = toc;
fprintf('  ...  %.2f seconds\n',tt)



%% 2. Find the best replicates to align to
% i) Find names of proteins with a single Gaussian
% ii) Find the overlap in single Gaussians b/w replicates
% iii) Choose the replicate with maximum overlap
tic
fprintf('    2. Find the best replicates to align to')

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

tt = toc;
fprintf('  ...  %.2f seconds\n',tt)



%% 3. Calculate best fit lines for adjustment
% i) Find the overlapping proteins
% ii) Find the Centers of these proteins
% iii) Fit a line
tic
fprintf('    3. Calculate best fit lines for adjustment')

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
      disp('Error: Alignment.m: uhhhh')
    end
    
    % iii) Fit a line
    I = abs(Ca - Cb)<User_alignment_window1;
    pfit(ci,rr,:) = robustfit(Ca(I), Cb(I));
    
  end
end

tt = toc;
fprintf('  ...  %.2f seconds\n',tt)



%% 4. Using fitted curves, adjust replicate data
% i) Gaussian fits: shift Center
% ii) Chromatograms: shift data points
tic
fprintf('    4. Using fitted curves, adjust replicate data')

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
x = -4:fraction_number+5;
adjusted_raw_data{1} = nan(size(cleandata{1}));
adjusted_raw_data{2} = nan(size(cleandata{2}));
for ri=1:size(num_val_MvsL)
  % MvsL
  y = cleandata{1}(ri,2:end);
  rr = num_val_MvsL(ri,1);
  b = pfit(ci,rr,1); % intercept
  m = pfit(ci,rr,2); % slope
  x2 = (x-b)/m;
  y2 = interp1(x,y,x2);
  y2(y2<.2) = nan;
  adjusted_raw_data{1}(ri,2:end) = y2;
  
  % HvsL
  y = cleandata{2}(ri,2:end);
  rr = num_val_HvsL(ri,1);
  b = pfit(ci,rr,1); % intercept
  m = pfit(ci,rr,2); % slope
  x2 = (x-b)/m;
  y2 = interp1(x,y,x2);
  y2(y2<.2) = nan;
  adjusted_raw_data{2}(ri,2:end) = y2;
end

tt = toc;
fprintf('  ...  %.2f seconds\n',tt)



%% 5. Make some summary statistics
% Delta_Center
% Delta_Height
% Delta_Width
% EuDist
% Divergence

tic
fprintf('    5. Make some summary statistics')

ci = 1;

Delta_center = cell(Number_of_experimental_channels,Nreplicates,Nreplicates);
Delta_height = cell(Number_of_experimental_channels,Nreplicates,Nreplicates);
Delta_width = cell(Number_of_experimental_channels,Nreplicates,Nreplicates);
EuDis = cell(Number_of_experimental_channels,Nreplicates,Nreplicates);

x = 1:fraction_number+10;

for ci = 1:Number_of_experimental_channels,
  for rr1 = 1:Nreplicates
    for rr2 = 1:Nreplicates
      
      % i) Find the overlapping proteins
      overlap = intersect(summerised_names_G1{ci,rr1},summerised_names_G1{ci,rr2});
      overlap([1 2]) = [];
      Ia = find(ismember(Gaus_import{ci,rr1}.textdata(:,4),overlap));
      Ib = find(ismember(Gaus_import{ci,rr2}.textdata(:,4),overlap));
      
      % ii) Calculate Gaussian curves
      G1 = zeros(fraction_number+10,length(overlap));
      G2 = zeros(fraction_number+10,length(overlap));
      for ri = 1:length(overlap)
        c1 = Gaus_import{ci,rr1}.data(Ia(ri)-1,1:3);
        c2 = Gaus_import{ci,rr2}.data(Ib(ri)-1,1:3);
        G1(:,ri) = c1(1)*exp( -(x-(c1(2)+5)).^2 /c1(3).^2 /2);
        G2(:,ri) = c2(1)*exp( -(x-(c2(2)+5)).^2 /c2(3).^2 /2);
      end
      
      % ii3) Calculate statistics
      Delta_center{ci,rr1,rr2} = abs(Gaus_import{ci,rr1}.data(Ia-1,2) - Gaus_import{ci,rr2}.data(Ib-1,2));
      Delta_height{ci,rr1,rr2} = abs(Gaus_import{ci,rr1}.data(Ia-1,1) - Gaus_import{ci,rr2}.data(Ib-1,1));
      Delta_width{ci,rr1,rr2} = abs(Gaus_import{ci,rr1}.data(Ia-1,3) - Gaus_import{ci,rr2}.data(Ib-1,3));
      EuDis{ci,rr1,rr2} = sqrt(sum((G1 - G2) .^ 2));
      
      % quality control...
      if length(overlap)~=length(Ia) || length(overlap)~=length(Ib)
        disp('Error: Alignment.m: uhhhh 2')
      end
      
    end
  end
end

tt = toc;
fprintf('  ...  %.2f seconds\n',tt)



%% 6. Write output
%    fid9_Name = strcat('Adjusted_',Experimental_channel,'_Raw_data_maxquant_rep',mat2str(alignment_counter),'.csv');
%    fid9B_Name = strcat('Adjusted_',Experimental_channel,'_Raw_for_ROC_analysis_rep',mat2str(alignment_counter),'.csv');
%    fid9B_Name = strcat('Adjusted_HvsM_Raw_data_maxquant_rep',mat2str(alignment_counter),'.csv');
%    fid6_Name = strcat('Adjusted_Chromatograms_vobose_rep',mat2str(alignment_counter),'_','.csv');
%    fid7_Name = strcat('Adjusted_Combined_OutputGaus_rep',mat2str(alignment_counter),'.csv');
%  fid10_Name = strcat('Adjusted_',Experimental_channel,'_Combined_OutputGaus.csv');
%  fid11_Name = strcat('Adjusted_',Experimental_channel,'_Raw_data_maxquant.csv');
tic
fprintf('    6. Write output')

if flag3
  protnames = {txt_MvsL, txt_HvsL, txt_HvsM};
else
  protnames = {txt_MvsL, txt_HvsL};
end

writeOutput_alignment

tt = toc;
fprintf('  ...  %.2f seconds\n',tt)



%% 7. Make figures

tic
fprintf('    7. Make figures')

makeFigures_alignment

tt = toc;
fprintf('  ...  %.2f seconds\n',tt)


