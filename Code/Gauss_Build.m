%%%%%%%%%%%%%%% Logic:
% 0. Initalize
% 1. Read input data (MaxQuant output)
% 2. Clean the chromatograms
% 3. On each profile, fit 1-5 Gaussians
% 4. Write output

%%%%%%%%%%%%%%% Instructions:
% - Rename maindir for your computer. (Try using the command 'pwd' to find the full path to a folder.)
% - Add maindir/Code and maindir/Code/Functions to Matlab path.


%% 0. Initialize
% - Set home folder (homedir), folder for subfunctions (funcdir), folder containing data files
%   (datdir), folder for figures (figdir), and finally the folder containing everything (maindir).
% - List all input files.
% - List all output files as either DebugOutputFiles or MainOutputFiles.
% - Define some variables.

% Set folders, i.e. define where everything lives.
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
% Ite = 500;                          % Number of iterations in crossvalidation
firstfrac= 2;                       % The first fraction to be analyzed
%lastfrac = 55;                      % The last fraction to be analyzed
experimental_channels = ['MvsL' 'HvsL'];
Nchannels = length(experimental_channels);
%Protein_number=1:Proteins;



%% 1. Read input data

% Import data from all input files
[num_val_MvsL,txt_MvsL] = xlsread(InputFile{1}); %Import file MvsL
[num_val_HvsL,txt_HvsL] = xlsread(InputFile{2}); %Import file HvsL
[SEC_size_alignment]=xlsread(InputFile{3});
% How many proteins and fractions are there?
[Nproteins, Nfractions]=size(num_val_MvsL);



%% 2. Clean the chromatograms.


for ci = 1:Nchannels % loop over channels
  
  % choose rawdata, i.e. data for the current channel
  if ci==1
    rawdata = num_val_MvsL;
  elseif ci==2
    rawdata = num_val_HvsL;
  end
  
  for ri = 1:Nproteins % loop over proteins
    
    raw_chromatogram = rawdata(ri,firstfrac:Nfractions);
    pause
    clean_chromatogram = cleanChromatogram(raw_chromatogram);
    
  end
end









