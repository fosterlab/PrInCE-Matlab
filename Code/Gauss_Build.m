%%%%%%%%%%%%%%% Logic:
% 0. Initalize
% 1. Read input data (MaxQuant output)
% 2. Clean the chromatograms
% 3. Fit 1-5 Gaussians on each cleaned chromatogram
% 4. Write output

%%%%%%%%%%%%%%% Instructions:
% - Rename maindir for your computer. (Try using the command 'pwd' to find the full path to a folder.)
% - Add maindir/Code and maindir/Code/Functions to Matlab path.

%%%%%%%%%%%%%%% To do:
% - Include the parfoor loops



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
figdir = [maindir 'Figures/']; % where figures live

% List all input files. These contain data that will be read by this script.
InputFile{1} = [datadir 'Combined_replicates_2014_04_22_contaminates_removed_for_MvsL_scripts.xlsx'];
InputFile{2} = [datadir 'Combined_replicates_2014_04_22_contaminates_removed_for_HvsL_scripts.xlsx'];
InputFile{3} = [datadir 'SEC_alignment.xlsx'];

% List all output files.
% MainOutputFiles contain data that will be used by a downstream program.
% DebugOutputFiles are not used by another program.
MainOutputFile{1} = [datadir 'HvsL_Combined_OutputGaus.csv'];
MainOutputFile{2} = [datadir 'HvsL_Summary_Gausians_for_individual_proteins.csv'];
DebugOutputFile{1} = [datadir 'Combined_Chromatograms.csv'];
DebugOutputFile{2} = [datadir 'Combined_Chromatograms_filtered_out.csv'];
DebugOutputFile{3} = [datadir 'Combined_OutputGaus_filtered_out.csv'];
DebugOutputFile{4} = [datadir 'Summary_Gausians_identifed.csv'];
DebugOutputFile{5} = [datadir 'Summary_Proteins_with_Gausians.csv'];
DebugOutputFile{6} = [datadir 'Proteins_not_fitted_to_gaussian.csv'];

% Define some variables.
MaxIter = 500;                          % Number of iterations in crossvalidation
firstfrac= 2;                       % The first fraction to be analyzed
%lastfrac = 55;                      % The last fraction to be analyzed
experimental_channels = {'MvsL' 'HvsL'};
Nchannels = length(experimental_channels);
%Protein_number=1:Proteins;



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

for ci = 1:Nchannels % loop over channels
  cleandata{ci} = zeros(size(rawdata{ci},1),size(rawdata{ci},2)+10);
  for ri = 1:Nproteins % loop over proteins
    
    raw_chromatogram = rawdata{ci}(ri,1:Nfractions);
    clean_chromatogram = cleanChromatogram(raw_chromatogram);
    
    % store clean_chromatogram in cleandata
    cleandata{ci}(ri,:) = clean_chromatogram;
    
  end
end



%% 3. Fit the Gaussians on the clean chromatograms

for ci = 1:Nchannels % loop over channels
  Try_Fit=zeros(1,Nproteins); % housekeeping variable
  for ri = 1:Nproteins % loop over proteins
    
    clean_chromatogram = cleandata{ci}(ri,:);
    
    % Don't fit Gaussians if there are less than 5 good data points
    if sum(clean_chromatogram > 0.05)<5
      continue
    end
    Try_Fit(1,ri)=1;
    
    ICs = gaussfitICs(clean_chromatogram); % 1x5 cell
    Output_SSE = holdoutSSE(MaxIter,clean_chromatogram);
    
  end
end




