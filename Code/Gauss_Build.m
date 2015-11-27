%%%%%%%%%%%%%%% Instructions for running Gauss_Build.m:
% - Rename maindir for your computer. (The command 'pwd' will give you the full path to a folder.)
% - Add maindir/Code and maindir/Code/Functions to Matlab path.
%
%%%%%%%%%%%%%%% Logic:
% 0. Initalize
% 1. Read input data (MaxQuant output)
% 2. Clean the chromatograms
% 3. Fit 1-5 Gaussians on each cleaned chromatogram
% 4. Write output
%
%%%%%%%%%%%%%%% Custom functions called by Gauss_build.m:
% - cleanChromatogram.m
% - choosemodel_holdout.m
% - gaussfitICs.m
% - fitgaussmodel.m
%
%%%%%%%%%%%%%%% To do:
% - Include the parfoor loops
% - Right now the logic is a bit roundabout:
%       First select a model and then fit that model. This is dangerous because the final fit
%       could be bad! Better to keep the parameters of the best fit from model selection.


%% 0. Initialize
% - Set home folder (homedir), folder for subfunctions (funcdir), folder containing data files
%   (datdir), folder for figures (figdir), and finally the folder containing everything (maindir).
% - List all input files.
% - List all output files as either DebugOutputFiles or MainOutputFiles.
% - Define some variables.

% Define folders, i.e. define where everything lives.
maindir = '/Users/Mercy/Academics/Foster/NickCodeData/GregPCP-SILAC/'; % where everything lives
codedir = [maindir 'Code/']; % where this script lives
funcdir = [maindir 'Code/Functions/']; % where small pieces of code live
datadir = [maindir 'DataFiles/']; % where data files live
datadir1 = [datadir 'Output_Chromatograms/'];
datadir2 = [datadir 'Output_Chromatograms_filtered_out/'];
datadir3 = [datadir 'OutputGaus/'];
datadir4 = [datadir 'OutputGaus_filtered_out/'];
figdir = [maindir 'Figures/']; % where figures live
% Make folders if necessary
if ~exist(datadir1, 'dir'); mkdir(datadir1); end
if ~exist(datadir2, 'dir'); mkdir(datadir2); end
if ~exist(datadir3, 'dir'); mkdir(datadir3); end
if ~exist(datadir4, 'dir'); mkdir(datadir4); end


% List all input files. These contain data that will be read by this script.
InputFile{1} = [datadir 'Combined_replicates_2014_04_22_contaminates_removed_for_MvsL_scripts.xlsx'];
InputFile{2} = [datadir 'Combined_replicates_2014_04_22_contaminates_removed_for_HvsL_scripts.xlsx'];
InputFile{3} = [datadir 'SEC_alignment.xlsx'];

% List output files.
% MainOutputFiles contain data that will be used by a downstream program.
% DebugOutputFiles are not used by another program.
MainOutputFile{1} = [datadir 'HvsL_Combined_OutputGaus.csv'];
MainOutputFile{2} = [datadir 'MvsL_Combined_OutputGaus.csv'];
MainOutputFile{3} = [datadir 'HvsL_Summary_Gausians_for_individual_proteins.csv'];
MainOutputFile{4} = [datadir 'MvsL_Summary_Gausians_for_individual_proteins.csv'];
DebugOutputFile{1} = [datadir 'Combined_Chromatograms.csv'];
DebugOutputFile{2} = [datadir 'Combined_Chromatograms_filtered_out.csv'];
DebugOutputFile{3} = [datadir 'Combined_OutputGaus_filtered_out.csv'];
DebugOutputFile{4} = [datadir 'Summary_Gausians_identifed.csv'];
DebugOutputFile{5} = [datadir 'Summary_Proteins_with_Gausians.csv'];
DebugOutputFile{6} = [datadir 'Proteins_not_fitted_to_gaussian_Hvsl.csv'];
DebugOutputFile{7} = [datadir 'Proteins_not_fitted_to_gaussian_Mvsl.csv'];

% Define some variables.
MaxIter = 500;                          % Number of iterations in holdout analysis
experimental_channels = {'MvsL' 'HvsL'};
Nchannels = length(experimental_channels);



%% 1. Read input data

% Import data from all input files
[num_val_MvsL,txt_MvsL] = xlsread(InputFile{1}); %Import file MvsL
[num_val_HvsL,txt_HvsL] = xlsread(InputFile{2}); %Import file HvsL
SEC_size_alignment = xlsread(InputFile{3});

% Remove first column, this is the replicate number
replicate = num_val_MvsL(:,1);
num_val_MvsL = num_val_MvsL(:,2:end);
num_val_HvsL = num_val_HvsL(:,2:end);

% How many proteins and fractions are there?
[Nproteins, Nfractions]=size(num_val_MvsL);

% Put raw data into a single variable for easy access
rawdata{1} = num_val_MvsL;
rawdata{2} = num_val_HvsL;



%% 2. Clean the chromatograms.

% store all clean chromatograms in cleandata
cleandata = cell(size(rawdata));

% a couple of housekeeping variables
tmp1 = cell(size(rawdata));
tmp2 = cell(size(rawdata));

for ci = 1:Nchannels % loop over channels
  cleandata{ci} = zeros(size(rawdata{ci},1),size(rawdata{ci},2)+10);
  for ri = 1:Nproteins % loop over proteins
    
    raw_chromatogram = rawdata{ci}(ri,1:Nfractions);
    [clean_chromatogram, x1, x2] = cleanChromatogram(raw_chromatogram);
    
    % store clean_chromatogram in cleandata
    cleandata{ci}(ri,:) = clean_chromatogram;
    
    % store housekeeping variables
    tmp1{ci}(ri,:) = x1;
    tmp2{ci}(ri,:) = x2;
  end
end



%% 3. Fit 1-5 Gaussians on each cleaned chromatogram

% pre-allocate some variables
Coef = cell(Nchannels, Nproteins);
SSE = nan(size(Coef));
adjrsquare = nan(size(Coef));
fit_flag = nan(size(Coef));
Try_Fit = zeros(size(Coef));

for ci = 1:Nchannels % loop over channels
  for ri = 1:Nproteins % loop over proteins
    
    % get a single clean chromatogram
    clean_chromatogram = cleandata{ci}(ri,:);
    
    % don't fit Gaussians if there are less than 5 good data points in the chromatogram
    if sum(clean_chromatogram > 0.05)<5
      continue
    end
    Try_Fit(ci,ri)=1;
    
    % choose the best model to fit
    model = choosemodel_holdout(clean_chromatogram,MaxIter);
    disp(['    fit protein number ' num2str(ri) ' with ' num2str(model.Ngauss) ' Gaussians...'])
    
    % fit that model
    [Coef{ci,ri},SSE(ci,ri),adjrsquare(ci,ri),fit_flag(ci,ri)] = fitgaussmodel(clean_chromatogram,model);

  end
end



%% 4. Write output

writeOutput_gaussbuild(datadir,MainOutputFile,DebugOutputFile,...
  Coef,SSE,adjrsquare,Try_Fit,...
  txt_MvsL,txt_HvsL,replicate,SEC_size_alignment,experimental_channels,...
  cleandata,rawdata,tmp1,tmp2);

