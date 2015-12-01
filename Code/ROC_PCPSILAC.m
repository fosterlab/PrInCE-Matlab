%%%%%%%%%%%%%%% Instructions for running Gauss_Build.m:
%
%
%%%%%%%%%%%%%%% Logic:
% 0. Initalize
%
%%%%%%%%%%%%%%% Custom functions called by Gauss_build.m:
% -
%
%%%%%%%%%%%%%%% To do:
% - Get all the input/output files sorted out.



%% 0. Initialize

% user defined settings
number_of_replicates=3;
number_of_channels=2;

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
InputFile{3} = [datadir 'Input/Corum_2012_human.xlsx'];
InputFile{4} = [datadir 'Input/'];

%define Raw SILAC ratios data, This is the output from the alignment script
MvsL_filename_Raw_rep1=[tmpdir 'MvsL_alignment\Realignment\Adjusted_MvsL_Raw_for_ROC_analysis_rep1.csv'];
MvsL_filename_Raw_rep2=[tmpdir 'MvsL_alignment\Realignment\Adjusted_MvsL_Raw_for_ROC_analysis_rep2.csv'];
MvsL_filename_Raw_rep3=[tmpdir 'MvsL_alignment\Realignment\Adjusted_MvsL_Raw_for_ROC_analysis_rep3.csv'];
HvsL_filename_Raw_rep1=[tmpdir 'HvsL_alignment\Realignment\Adjusted_HvsL_Raw_for_ROC_analysis_rep1.csv'];
HvsL_filename_Raw_rep2=[tmpdir 'HvsL_alignment\Realignment\Adjusted_HvsL_Raw_for_ROC_analysis_rep2.csv'];
HvsL_filename_Raw_rep3=[tmpdir 'HvsL_alignment\Realignment\Adjusted_HvsL_Raw_for_ROC_analysis_rep3.csv'];
MvsL_filename_gaus_rep1=[tmpdir 'MvsL_alignment\Realignment\Adjusted_Combined_OutputGaus_rep1.csv'];
MvsL_filename_gaus_rep2=[tmpdir 'MvsL_alignment\Realignment\Adjusted_Combined_OutputGaus_rep2.csv'];
MvsL_filename_gaus_rep3=[tmpdir 'MvsL_alignment\Realignment\Adjusted_Combined_OutputGaus_rep3.csv'];
HvsL_filename_gaus_rep1=[tmpdir 'HvsL_alignment\Realignment\Adjusted_Combined_OutputGaus_rep1.csv'];
HvsL_filename_gaus_rep2=[tmpdir 'HvsL_alignment\Realignment\Adjusted_Combined_OutputGaus_rep2.csv'];
HvsL_filename_gaus_rep3=[tmpdir 'HvsL_alignment\Realignment\Adjusted_Combined_OutputGaus_rep3.csv'];
List_of_Raw_filename={MvsL_filename_Raw_rep1,MvsL_filename_Raw_rep2,MvsL_filename_Raw_rep3,...
  HvsL_filename_Raw_rep1,HvsL_filename_Raw_rep2,HvsL_filename_Raw_rep3};
List_of_Gaus_filename={MvsL_filename_gaus_rep1,MvsL_filename_gaus_rep2,MvsL_filename_gaus_rep3,...
  HvsL_filename_gaus_rep1,HvsL_filename_gaus_rep2,HvsL_filename_gaus_rep3};



%% 1. Read input data

[~, Protein_IDs] = xlsread(InputFile{1});
fid=fopen(InputFile{2}, 'rt');    %Input corum data base.
Corum_Import= textscan(fid,'%s\t', 'Delimiter',',');
fclose(fid);
[~, ~, Corum_complexes] = xlsread('Corum_2012_human.xlsx','Sheet1','A2:L1847');


No=length(Corum_Import{1})/2;
Corum_Protein_names = (reshape(Corum_Import{1,1},2,No)');
Unique_Corum = unique(Corum_Protein_names);


Corum_complexes(cellfun(@(x) ~isempty(x) && isnumeric(x) && isnan(x),Corum_complexes)) = {''};
cellVectors = Corum_complexes(:,[2,3,4,5,7,9,10,11,12]);
Corum_complexes = Corum_complexes(:,[1,6,8]);
% Replace non-numeric cells with NaN
R = cellfun(@(x) ~isnumeric(x) && ~islogical(x),Corum_complexes); % Find non-numeric cells
Corum_complexes(R) = {NaN}; % Replace non-numeric cells
Corum_complexes_data = reshape([Corum_complexes{:}],size(Corum_complexes)); % Create output variable


Master_Gaussian_list=importdata('Master_guassian_list.csv',',');


% Allocate imported array to column variable names
Corum.Complexid = Corum_complexes_data(:,1);
Corum.Complexname = cellVectors(:,1);
Corum.Synonyms = cellVectors(:,2);
Corum.organism = cellVectors(:,3);
Corum.subunitsUniProtIDs = cellVectors(:,4);
Corum.subunitsEntrezIDs = Corum_complexes_data(:,2);
Corum.proteincomplexpurificationmethod = cellVectors(:,5);
Corum.PubMedid = Corum_complexes_data(:,3);
Corum.FunCatcategories = cellVectors(:,6);
Corum.functionalcomment = cellVectors(:,7);
Corum.diseasecomment = cellVectors(:,8);
Corum.subunitcomment = cellVectors(:,9);

%Split UniprotID to individual entries
length_of_corum=length(Corum.subunitsUniProtIDs);
Corum.subunitsUniProtIDs_split=cell(length_of_corum,1);
Split_Complex=regexp(Corum.subunitsUniProtIDs,',','split');
for Complex_no=1:length_of_corum
  %using reg expression split entries in corum
  number_proteins_in_complex=length(Split_Complex{Complex_no});
  if number_proteins_in_complex > 1
    Complex_to_write_out=Split_Complex{Complex_no};
    for write_out_counter1=1:number_proteins_in_complex
      Corum.subunitsUniProtIDs_split{Complex_no,write_out_counter1}=Complex_to_write_out(1,write_out_counter1);
    end
  elseif number_proteins_in_complex == 1
    Corum.subunitsUniProtIDs_split{Complex_no,1}=Split_Complex{Complex_no};
  end
end

%Create file to save replicate information to
Save_replicate_folder_name=cell((number_of_replicates*number_of_channels),1);

