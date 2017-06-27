function user_new = readExperimentFile(user)
% READEXPERIMENTFILE Reads experimental_design.rtf.
%
%   READEXPERIMENTFILE(user) reads the experimental design and plotting
%   parameters in experimental_design.rtf. These include the number of
%   conditions, replicates, and the number of plots generated. user is a
%   structure used by master script PRINCE.
%
%   See also PRINCE


%%%%%%% Hard coded parameters
user.userwindow = 2;
user.separateByReplicate = 0;
user.separateByChannel = 1;
user.fdr = 0.01;
user.User_alignment_window1 = 20;
user.nickflag = 0;
tmp = dir([user.maindir 'Input/*ondition*.csv']);
if isempty(tmp)
  error('No data files in Input folder. Data files must be named in this format: conditionX.csv, where X is the condition number.')
end
user.MQfiles = cell(size(tmp));
user.silacratios = cell(size(tmp));
for ii = 1:length(tmp)
  user.MQfiles{ii} = tmp(ii).name;
  if ~strfind(tmp(ii).name,'csv')
    error('Data files (conditionX.csv) must in csv format. Please Save As csv.')
  end
  ss = strrep(tmp(ii).name,'.csv','');
  user.silacratios{ii} = ss;
end
if length(user.MQfiles)<=2
  user.comparisonpairs = user.silacratios;
else
  user.comparisonpairs = user.silacratios(1:2);
  warning('Fold changes will be calculated between condition1 and condition2...')
end


%%%%%%% Read parameters
code_strings = {'precision' ...
  'Major protein groups filename' ...
  'Reference database filename' ...
  'Treatment condition' ...
  'No-treatment condition' ...
  'Number of fractions' ...
  'Number of replicates' ...
  'Skip Alignment' ...
  'Skip FoldChanges' ...
  'Skip Interactions' ...
  'Skip Complexes' ...
  'Plots of Gaussians' ...
  'Plots of fold changes'};

fn = [user.maindir 'experimental_design.rtf'];
fid = fopen(fn);

while ~feof(fid)
  
  t1 = fgetl(fid);
  
  Ifield = zeros(size(code_strings));
  for ii = 1:length(Ifield)
    %I(ii) = ismember(code_strings{ii},t1);
    tmp = strfind(t1, code_strings{ii});
    if ~isempty(tmp)
      Ifield(ii) = tmp;
    end
  end
  
  Ifield = find(Ifield>0);
  
  if sum(Ifield)==0
    continue;
  elseif length(Ifield)==1
    usable_string = strsplit(t1,':');
    usable_string = usable_string{2};
    usable_string = strsplit(usable_string,'\');
    usable_string = usable_string{1};
    
    if Ifield==1 % desired precision
      I = ismember(usable_string,['1' '2' '3' '4' '5' '6' '7' '8' '9' '0' '.']);
      desiredPrecision = str2double(usable_string(I));
      if desiredPrecision>0 & desiredPrecision<=1
        user.desiredPrecision = desiredPrecision;
      elseif desiredPrecision>=1 & desiredPrecision<100
        user.desiredPrecision = desiredPrecision/100;
      else
        error('Desired precision of interaction list must be a percentage or a number between 0 and 1.')
      end
      
    elseif Ifield==2 % 'Major protein groups filename'
      user.majorproteingroupsfile = process_field(usable_string,'');
      
    elseif Ifield==3 % 'Reference database filename'
      user.corumfile = process_field(usable_string,[user.maindir 'Input/allComplexes.txt']);
      
    elseif Ifield==4 % 'Treatment condition'
      user.treatmentcondition = process_field(usable_string,'condition1');
      
    elseif Ifield==5 % 'No-treatment condition'
      user.notreatmentcondition = process_field(usable_string,'condition2');
      if length(user.MQfiles)==1
        user.notreatmentcondition = '';
      end
      
    elseif Ifield==6 % 'Number of fractions'
      user.Nfraction = process_field(usable_string,'-1');
      user.Nfraction = str2double(user.Nfraction);
      
    elseif Ifield==7 % 'Number of replicates'
      user.Nreplicate = process_field(usable_string,'-1');
      user.Nreplicate = str2double(user.Nreplicate);
      
    elseif Ifield==8 % 'Skip GaussBuild'
      user.skipalignment = process_field(usable_string,'0');
      user.skipalignment = str2double(user.skipalignment);
      
    elseif Ifield==9 % 'Skip FoldChanges'
      user.skipcomparison = process_field(usable_string,'0');
      user.skipcomparison = str2double(user.skipcomparison);
    
    elseif Ifield==10 % 'Skip FoldChanges'
      user.skipinteractions = process_field(usable_string,'0');
      user.skipinteractions = str2double(user.skipinteractions);
      
    elseif Ifield==11 % 'Skip FoldChanges'
      user.skipcomplexes = process_field(usable_string,'0');
      user.skipcomplexes = str2double(user.skipcomplexes);
      
    elseif Ifield==12 % 'Plots of Gaussians'
      user.fastgaussbuild = process_field(usable_string,'0');
      user.fastgaussbuild = str2double(user.fastgaussbuild);
    
    elseif Ifield==13 % 'Plots of fold changes'
      user.fastcomparison = process_field(usable_string,'0');
      user.fastcomparison = str2double(user.fastcomparison);
      
    end
    
  else
    error('Bad experimental_design.rtf file.')
  end
  
end

user_new = user;
end


function ss = process_field(usable_string,default)
ss = strrep(usable_string,'\','');
ss = strrep(ss,'}','');
ss = strrep(ss,'{','');
ss = strtrim(ss);
if isempty(ss)
  ss = default;
end
if isequal(lower(ss),'no')
  ss = '0';
elseif isequal(lower(ss),'yes')
  ss = '1';
end
end
