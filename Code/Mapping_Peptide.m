%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                               Script to map coverage of peptide within each protein group
%                                               Created by Nichollas Scott,
%                                                    Foster lab,UBC, 2015
%Description: Code to generate peptide and protein profile as well as
%determine uneveness and diversity from fasta file.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all

%Define number of replicates
replicate_number=3;

%Experiment boundries
Define_experiment_fractions_start=[1,56,111];
Define_experiment_fractions_end=[55,110,165];

%colour palatte
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic

%% Import Contamiation file
Mouse_fasta_file=fastaread('UP_human_canonical_isoforms_TrEMBL_Oct2013.fasta');
Mouse_fasta_header=cell(length(Mouse_fasta_file),1);
%Convert header into form which can be used
for counter=1:length(Mouse_fasta_file)
  seperated_values=regexp(Mouse_fasta_file(counter).Header,'\|+','split');
  Mouse_fasta_header{counter}=seperated_values{2};
end

%% Read in peptide results
peptides_results='modificationSpecificPeptides.txt';
inFile=fopen(peptides_results);
headerLine=fgetl(inFile);

%Read in header
tempCellArray=textscan(headerLine,'%s','DELIMITER','\t','BUFSIZE',100000);

%look for Calibrated retention time column
%Convert header to upper case, this is to ensure matching for caps
fileCols=upper(tempCellArray{1,1});
columnNames = upper({'Sequence','Modifications','Mass','Proteins','Calibrated retention time'});

%look up if position of column file
idx = ismember(fileCols,columnNames);

%Determine which columns to read in, if no the correct column assign %*s
format='';
for ci=1:length(fileCols)
  if idx(ci)==0
    format=strcat(format,'%*s ');
  elseif idx(ci)==1
    format=strcat(format,'%s ');
  end;
end;

% read rest of matched_features file containing columns of interest
Peptide_content=textscan(inFile,format,'DELIMITER','\t','BUFSIZE',1000000);

%% unpack data to alterative varibles
Peptide_content=[Peptide_content{:,:}];

%% Import MvsL from peptide results
expression_ratio_MvsL = upper(strcat('Ratio ','M/L'));
expression_MvsL=strcat(expression_ratio_MvsL,'(\d+)');

%ignore whitespaces
column_name_temp = strrep(upper(tempCellArray{1,1}), ' ', '');

%Find all MvsL ratios
MvsL_ratios=regexp(column_name_temp,expression_MvsL);

%find the columns which contain the isotopologue values of interest
idx_MvsL = cellfun(@(x) isequal(x, 1), MvsL_ratios);

%Determine which columns to read in, if no the correct column assign %*s
format_MvsL='';
for ci=1:length(idx_MvsL)
  if idx_MvsL(ci)==0
    format_MvsL=strcat(format_MvsL,'%*s ');
  elseif idx_MvsL(ci)==1
    format_MvsL=strcat(format_MvsL,'%s ');
  end;
end;

%re-open matched feature text file
inFile=fopen(peptides_results);

%read rest of matched_features file containing columns of interest
MvsL_peptide_content=textscan(inFile,format_MvsL,'DELIMITER','\t','treatAsEmpty','NULL', 'EmptyValue',NaN,'BUFSIZE',1000000);

%Re-order columns, create list of headers
header_names=tempCellArray{1,1};
header_names(idx_MvsL==0) = [];

%Extract numbers within headers
header_names_numbers=regexp(header_names,'(\d+)?(\d+)?(\d+)','match');
header_numbers=str2double([header_names_numbers{:}]);

corrected_ordered = zeros(1, length(header_numbers));  % Pre-allocate positions
for ii = 1:length(header_numbers)
  [value, index] = min(header_numbers);
  corrected_ordered(ii) = index;
  header_numbers(index) = Inf;
end

%Re-order column to be used downstream
MvsL_peptide_content=MvsL_peptide_content(1,corrected_ordered);

%Save
save('MvsL_peptide_content','MvsL_peptide_content');

clearvars MvsL_peptide_content

%% Import HvsL from peptide results
expression_ratio_HvsL = upper(strcat('Ratio ','H/L'));
expression_HvsL=strcat(expression_ratio_HvsL,'(\d+)');

%Find all MvsL ratios
HvsL_ratios=regexp(column_name_temp,expression_HvsL);

%find the columns which contain the isotopologue values of interest
idx_HvsL = cellfun(@(x) isequal(x, 1), HvsL_ratios);

%Determine which columns to read in, if no the correct column assign %*s
format_HvsL='';
for ci=1:length(idx_HvsL)
  if idx_HvsL(ci)==0
    format_HvsL=strcat(format_HvsL,'%*s ');
  elseif idx_HvsL(ci)==1
    format_HvsL=strcat(format_HvsL,'%s ');
  end;
end;

%re-open matched feature text file
inFile=fopen(peptides_results);

%read rest of matched_features file containing columns of interest
HvsL_peptide_content=textscan(inFile,format_HvsL,'DELIMITER','\t','treatAsEmpty','NULL', 'EmptyValue',NaN,'BUFSIZE',1000000);

%Re-order columns, create list of headers
header_names=tempCellArray{1,1};
header_names(idx_HvsL==0) = [];

%Extract numbers within headers
header_names_numbers=regexp(header_names,'(\d+)?(\d+)?(\d+)','match');
header_numbers=str2double([header_names_numbers{:}]);

corrected_ordered = zeros(1, length(header_numbers));  % Pre-allocate positions
for ii = 1:length(header_numbers)
  [value, index] = min(header_numbers);
  corrected_ordered(ii) = index;
  header_numbers(index) = Inf;
end

%Re-order column to be used downstream
HvsL_peptide_content=HvsL_peptide_content(1,corrected_ordered);

%Save
save('HvsL_peptide_content', 'HvsL_peptide_content');

clearvars HvsL_peptide_content

%% Import HvsM from peptide results
expression_ratio_HvsM = upper(strcat('Ratio ','H/M'));
expression_HvsM=strcat(expression_ratio_HvsM,'(\d+)');

%Find all MvsL ratios
HvsM_ratios=regexp(column_name_temp,expression_HvsM);

%find the columns which contain the isotopologue values of interest
idx_HvsM = cellfun(@(x) isequal(x, 1), HvsM_ratios);

%Determine which columns to read in, if no the correct column assign %*s
format_HvsM='';
for ci=1:length(idx_HvsM)
  if idx_HvsM(ci)==0
    format_HvsM=strcat(format_HvsM,'%*s ');
  elseif idx_HvsM(ci)==1
    format_HvsM=strcat(format_HvsM,'%s ');
  end;
end;

%re-open matched feature text file
inFile=fopen(peptides_results);

%read rest of matched_features file containing columns of interest
HvsM_peptide_content=textscan(inFile,format_HvsM,'DELIMITER','\t','treatAsEmpty','NULL', 'EmptyValue',NaN,'BUFSIZE',1000000);

%Re-order columns, create list of headers
header_names=tempCellArray{1,1};
header_names(idx_HvsM==0) = [];

%Extract numbers within headers
header_names_numbers=regexp(header_names,'(\d+)?(\d+)?(\d+)','match');
header_numbers=str2double([header_names_numbers{:}]);

corrected_ordered = zeros(1, length(header_numbers));  % Pre-allocate positions
for ii = 1:length(header_numbers)
  [value, index] = min(header_numbers);
  corrected_ordered(ii) = index;
  header_numbers(index) = Inf;
end

%Re-order column to be used downstream
HvsM_peptide_content=HvsM_peptide_content(1,corrected_ordered);

%Save
save('HvsM_peptide_content', 'HvsM_peptide_content');

clearvars HvsM_peptide_content

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Read in peptide results
Protein_results='proteinGroups.txt';
inFile=fopen(Protein_results);
headerLine=fgetl(inFile);

%Read in header
tempCellArray=textscan(headerLine,'%s','DELIMITER','\t','BUFSIZE',100000);

%look for Calibrated retention time column
%Convert header to upper case, this is to ensure matching for caps
fileCols=upper(tempCellArray{1,1});
columnNames = upper({'Majority protein IDs', 'Peptides'});

%look up if position of column file
idx = ismember(fileCols,columnNames);

%Determine which columns to read in, if no the correct column assign %*s
format='';
for ci=1:length(fileCols)
  if idx(ci)==0
    format=strcat(format,'%*s ');
  elseif idx(ci)==1
    format=strcat(format,'%s ');
  end;
end;

% read rest of matched_features file containing columns of interest
Protein_content=textscan(inFile,format,'DELIMITER','\t','BUFSIZE',1000000);


%% Import MvsL from peptide results
expression_ratio_MvsL = upper(strcat('Ratio ','M/L'));
expression_MvsL=strcat(expression_ratio_MvsL,'(\d+)');

%ignore whitespaces
column_name_temp = strrep(upper(tempCellArray{1,1}), ' ', '');

%Find all MvsL ratios
MvsL_ratios=regexp(column_name_temp,expression_MvsL);

%find the columns which contain the isotopologue values of interest
idx_MvsL = cellfun(@(x) isequal(x, 1), MvsL_ratios);

%Determine which columns to read in, if no the correct column assign %*s
format_MvsL='';
for ci=1:length(idx_MvsL)
  if idx_MvsL(ci)==0
    format_MvsL=strcat(format_MvsL,'%*s ');
  elseif idx_MvsL(ci)==1
    format_MvsL=strcat(format_MvsL,'%s ');
  end;
end;

%re-open matched feature text file
inFile=fopen(Protein_results);

%read rest of matched_features file containing columns of interest
MvsL_Protein_content=textscan(inFile,format_MvsL,'DELIMITER','\t','treatAsEmpty','NULL', 'EmptyValue',NaN,'BUFSIZE',1000000);

%Re-order columns, create list of headers
header_names=tempCellArray{1,1};
header_names(idx_MvsL==0) = [];

%Extract numbers within headers
header_names_numbers=regexp(header_names,'(\d+)?(\d+)?(\d+)','match');
header_numbers=str2double([header_names_numbers{:}]);

corrected_ordered = zeros(1, length(header_numbers));  % Pre-allocate positions
for ii = 1:length(header_numbers)
  [value, index] = min(header_numbers);
  corrected_ordered(ii) = index;
  header_numbers(index) = Inf;
end

%Re-order column to be used downstream
MvsL_Protein_content=MvsL_Protein_content(1,corrected_ordered);


%% Import HvsL from peptide results
expression_ratio_HvsL = upper(strcat('Ratio ','H/L'));
expression_HvsL=strcat(expression_ratio_HvsL,'(\d+)');

%Find all MvsL ratios
HvsL_ratios=regexp(column_name_temp,expression_HvsL);

%find the columns which contain the isotopologue values of interest
idx_HvsL = cellfun(@(x) isequal(x, 1), HvsL_ratios);

%Determine which columns to read in, if no the correct column assign %*s
format_HvsL='';
for ci=1:length(idx_HvsL)
  if idx_HvsL(ci)==0
    format_HvsL=strcat(format_HvsL,'%*s ');
  elseif idx_HvsL(ci)==1
    format_HvsL=strcat(format_HvsL,'%s ');
  end;
end;

%re-open matched feature text file
inFile=fopen(Protein_results);

%read rest of matched_features file containing columns of interest
HvsL_Protein_content=textscan(inFile,format_HvsL,'DELIMITER','\t','treatAsEmpty','NULL', 'EmptyValue',NaN,'BUFSIZE',1000000);

%Re-order columns, create list of headers
header_names=tempCellArray{1,1};
header_names(idx_HvsL==0) = [];

%Extract numbers within headers
header_names_numbers=regexp(header_names,'(\d+)?(\d+)?(\d+)','match');
header_numbers=str2double([header_names_numbers{:}]);

corrected_ordered = zeros(1, length(header_numbers));  % Pre-allocate positions
for ii = 1:length(header_numbers)
  [value, index] = min(header_numbers);
  corrected_ordered(ii) = index;
  header_numbers(index) = Inf;
end

%Re-order column to be used downstream
HvsL_Protein_content=HvsL_Protein_content(1,corrected_ordered);



%% Import HvsL from peptide results
expression_ratio_HvsM = upper(strcat('Ratio ','H/L'));
expression_HvsM=strcat(expression_ratio_HvsM,'(\d+)');

%Find all MvsL ratios
HvsM_ratios=regexp(column_name_temp,expression_HvsM);

%find the columns which contain the isotopologue values of interest
idx_HvsM = cellfun(@(x) isequal(x, 1), HvsM_ratios);

%Determine which columns to read in, if no the correct column assign %*s
format_HvsM='';
for ci=1:length(idx_HvsM)
  if idx_HvsM(ci)==0
    format_HvsM=strcat(format_HvsM,'%*s ');
  elseif idx_HvsM(ci)==1
    format_HvsM=strcat(format_HvsM,'%s ');
  end;
end;

%re-open matched feature text file
inFile=fopen(Protein_results);

%read rest of matched_features file containing columns of interest
HvsM_Protein_content=textscan(inFile,format_HvsM,'DELIMITER','\t','treatAsEmpty','NULL', 'EmptyValue',NaN,'BUFSIZE',1000000);

%Re-order columns, create list of headers
header_names=tempCellArray{1,1};
header_names(idx_HvsM==0) = [];

%Extract numbers within headers
header_names_numbers=regexp(header_names,'(\d+)?(\d+)?(\d+)','match');
header_numbers=str2double([header_names_numbers{:}]);

corrected_ordered = zeros(1, length(header_numbers));  % Pre-allocate positions
for ii = 1:length(header_numbers)
  [value, index] = min(header_numbers);
  corrected_ordered(ii) = index;
  header_numbers(index) = Inf;
end

%Re-order column to be used downstream
HvsM_Protein_content=HvsM_Protein_content(1,corrected_ordered);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% unpack data to alterative varibles
Protein_content=[Protein_content{:,:}];

%% Import Majority Protein ID list
[~, Protein_IDs] = xlsread('Major_protein_groups.xlsx');

%determine all the peptides which belong to a protein group
Dimension_of_Protein_IDs=size(Protein_IDs);

% For each Major protein group look up every peptide, create a text data
%file off all the read outs and make a figure

mkdir('ProteinGroups text_output');
mkdir('ProteinGroups image_output');


%Load in require isotopologue channel
for isotope_channel=3 %1:3
  clearvars Protein_isotope_ratio_values
  clearvars Peptide_isotope_ratio_values
  
  %Import Peptide data, this file may be very large so code is designed to
  %read the files in as few times as possible
  
  if isotope_channel==1
    
    %load in protein ratio (unested form)
    load 'MvsL_peptide_content' MvsL_peptide_content
    Peptide_isotope_ratio_values=[MvsL_peptide_content{:,:}];
    
    %delete varible to free up memory
    %clearvars MvsL_peptide_content
    
    isotope_name_value={'MvsL'};
  elseif isotope_channel==2
    
    %load in protein ratio (unested form)
    load 'HvsL_peptide_content' HvsL_peptide_content
    
    Peptide_isotope_ratio_values=[HvsL_peptide_content{:,:}];
    
    %delete varible to free up memory
    %clearvars HvsL_peptide_content
    
    isotope_name_value={'HvsL'};
  elseif isotope_channel==3
    %load in protein ratio (unested form)
    load 'HvsM_peptide_content' HvsM_peptide_content
    load 'MvsL_peptide_content' MvsL_peptide_content
    load 'HvsL_peptide_content' HvsL_peptide_content
    
    Peptide_isotope_ratio_values=[HvsM_peptide_content{:,:}];
    Peptide_isotope_ratio_values_HvsL=[HvsL_peptide_content{:,:}];
    Peptide_isotope_ratio_values_MvsL=[MvsL_peptide_content{:,:}];
    %delete varible to free up memory
    clearvars HvsM_peptide_content
    clearvars HvsL_peptide_content
    clearvars MvsL_peptide_content
    isotope_name_value={'HvsM'};
  end
  %Read out so users can know which isotopology channels are being processed
  strcat('now processing_',isotope_name_value,'_Profiles')
  
  %Create table of all the sequences used for mapping, this is for the users
  %records
  fileName = strcat('Summary_of_',isotope_name_value,'_protein_identified_from_uniprot.csv');
  fid_Summary_uniprot = fopen(fileName{:},'a+');
  fprintf (fid_Summary_uniprot,'%s,%s,%s,%s\n',...
    'Major Protein group', 'Alterative Protein name 1','Alterative Protein name 2','Sequence');
  
  %Create summary table to
  fileName = strcat('Peptide_Summary_',isotope_name_value,'_table.csv');
  fid_summary = fopen(fileName{:},'a+');
  fprintf (fid_summary,'%s,%s,%s,%s,%s,%s,%s,%s,%s,\n',...
    'Protein','Isotopic channel','Replicate', 'Cophenetic correlation (raw values from maxquant)',...
    'Cophenetic correlation (fixed values)','Contains peptide which differ based on stdev (raw values from maxquant)',...
    'Contains peptide which differ based on stdev (fixed values)', 'Shannon diversity index(H)','Evenness ');  %Write Header
  
  
  %Determine which peptide belong to which protein
  for Protein_ID_counter1 =2:Dimension_of_Protein_IDs(1)
    
    
    %create a logic array to determine which peptides belong to this protein group
    Found_peptide=zeros(length(Peptide_content),1);
    
    for Protein_ID_counter2 =1:Dimension_of_Protein_IDs(2)-sum(cellfun('isempty',Protein_IDs(Protein_ID_counter1,:)))
      
      %find all peptides belonging to protein group
      
      Find_peptide=cellfun(@(IDX) ~isempty(IDX), strfind(Peptide_content(:,4), Protein_IDs{Protein_ID_counter1,Protein_ID_counter2}));
      
      %Combine two logic list
      Found_peptide=(Find_peptide | Found_peptide);
    end
    
    %find the position of all peptide from this protein group
    idx = find(Found_peptide);
    
    number_idx=length(idx);
    
    %Find protein in protein
    Find_protein=strfind(Protein_content(:,1), Protein_IDs{Protein_ID_counter1,Protein_ID_counter2});
    
    %Find all non empty cells
    Find_protein1=find(~cellfun(@isempty,Find_protein));
    
    %Save protein group information
    Protein_group.Proteins_assigned_in_Protein_table_assigned_in_peptide_table=Protein_content(Find_protein1,1);
    
    %% import protein data
    if isotope_channel==1
      %load in protein ratio (unested form)
      Protein_isotope_ratio_values=[MvsL_Protein_content{:,:}];
      
      %write protein ratio
      for ci2=1:length(Protein_isotope_ratio_values(Find_protein1(1),:))
        protein_ratio(1,ci2)=str2double(Protein_isotope_ratio_values{Find_protein1(1)+1,ci2});
      end
      
      %delete varible to free up memory
      clearvars Protein_isotope_ratio_values
    elseif isotope_channel==2
      %load in protein ratio (unested form)
      Protein_isotope_ratio_values=[HvsL_Protein_content{:,:}];
      
      %write protein ratio
      for ci2=1:length(Protein_isotope_ratio_values(Find_protein1(1),:))
        protein_ratio(1,ci2)=str2double(Protein_isotope_ratio_values{Find_protein1(1)+1,ci2});
      end
      
      %delete varible to free up memory
      clearvars Protein_isotope_ratio_values
      
    elseif isotope_channel==3
      %load in protein ratio (unested form)
      Protein_isotope_ratio_values=[HvsM_Protein_content{:,:}];
      
      %write protein ratio
      for ci2=1:length(Protein_isotope_ratio_values(Find_protein1(1),:))
        protein_ratio(1,ci2)=str2double(Protein_isotope_ratio_values{Find_protein1(1)+1,ci2});
      end
      
      %delete varible to free up memory
      clearvars Protein_isotope_ratio_values
      
      %load in protein ratio (unested form)
      Protein_isotope_ratio_values=[MvsL_Protein_content{:,:}];
      
      %write protein ratio
      for ci2=1:length(Protein_isotope_ratio_values(Find_protein1(1),:))
        protein_ratio_MvsL(1,ci2)=str2double(Protein_isotope_ratio_values{Find_protein1(1)+1,ci2});
      end
      
      %delete varible to free up memory
      clearvars Protein_isotope_ratio_values
      
      %load in protein ratio (unested form)
      Protein_isotope_ratio_values=[HvsL_Protein_content{:,:}];
      
      %write protein ratio
      for ci2=1:length(Protein_isotope_ratio_values(Find_protein1(1),:))
        protein_ratio_HvsL(1,ci2)=str2double(Protein_isotope_ratio_values{Find_protein1(1)+1,ci2});
      end
      
      %delete varible to free up memory
      clearvars Protein_isotope_ratio_values
      
    end
    
    
    for replicate_counter=1:replicate_number
      
      %Copy peptide information to struture to use for figure
      Protein_group.Sequence=cell(number_idx,1);
      Protein_group.Proteins_assigned_in_peptide_table_assigned_in_peptide_table=cell(number_idx,1);
      Protein_group.mass=cell(number_idx,1);
      Protein_group.calibrated_time=cell(number_idx,1);
      
      %Create array to store ratio information
      Protein_group.protein_ratio_raw=zeros(1,length(protein_ratio));
      Protein_group.protein_ratio_fixed=zeros(1,length(protein_ratio));
      Protein_group.peptide_ratio_raw=zeros(1,length(protein_ratio));
      Protein_group.peptide_ratio_fixed=zeros(1,length(protein_ratio));
      
      
      %Save protein information before fixing values
      Protein_group.protein_ratio_raw=protein_ratio;
      
      %Save major protein group
      Protein_group.Major_protein_group_name=Protein_IDs(Protein_ID_counter1,1);
      
      %% Fix missing values and extreme values for protein
      %Define variable
      Count =numel(protein_ratio(1,:));
      AK2= isnan(protein_ratio(1,:));
      
      %find NaN values flanked by numbers
      for i=2:(Count-1)
        if 1==isnan(protein_ratio(1,i))
          Test1=AK2((i-1):(i+1));
          if  1==isequal(~Test1,[1 0 1])
            Number=protein_ratio(1,(i-1):(i+1));
            Number2 = Number(~Test1);
            miss = mean(Number2); %average if  data point is missing
          else
            miss=NaN;
          end
          protein_ratio(1,i) = miss;
        else
        end
      end
      
      %Fix extreme value
      for i=2:(Count-1)
        Test2=AK2((i-1):(i+1));
        if  1==isequal(~Test2,[1 1 1])
          
          %Create value correspond to a value in front and behide
          Number2 =protein_ratio(1,(i-1):(i+1));
          
          if (protein_ratio(1,i)/sum(Number2))>0.95
            %If the value is extreme, filter by avaering the
            %flanking numbers
            protein_ratio(1,i) = (protein_ratio(1,i-1)+protein_ratio(1,i+1))/2;
            
          end
        end
      end
      
      %Unmodified values for plotting
      protein_ratio_for_figure=protein_ratio;
      
      %Remove values below 0.2
      protein_ratio(protein_ratio<0.2) = NaN;
      
      %% Consecutive numbers, if less then 5 consecutive number removes
      %% chromogram and replace with 0.05
      Multi_NaN = NaN(1,10);
      Sig=[Multi_NaN protein_ratio Multi_NaN];
      AKK= ~(isnan(Sig));
      Count2= numel(Sig);
      
      Sig2=zeros(1,Count2); %Allocate variable
      for i= 11:(Count2-10)
        if 5== sum(AKK((i-(4)):i));
          Sig2(1,i) = Sig(i);
        elseif 5== sum(AKK((i-(3)):(i+1)));
          Sig2(1,i) = Sig(i);
        elseif 5== sum(AKK((i-(2)):(i+2)));
          Sig2(1,i) = Sig(i);
        elseif 5== sum(AKK((i-(1)):(i+3)));
          Sig2(1,i) = Sig(i);
        elseif 5== sum(AKK((i):(i+4)));
          Sig2(1,i) = Sig(i);
        else Sig2(1,i) = 0.2;
        end
        
      end
      Out= Sig2(1,((10+1):(Count2-10)));
      
      % Smooth data, change vary the number of data points used for smoothing
      protein_ratio= (smooth(Out,3))';
      
      protein_ratio(protein_ratio<0.3) = NaN;
      Protein_group.protein_ratio_fixed=protein_ratio;
      
      %% If looking at the HvsM channel, also process the MvsL and Hvsl to fix
      %missing values BUT this is only so the profiles can be mapped
      if isotope_channel==3
        
        %Fix missing values and extreme values for MvsL protein values
        %Define variable
        Count =numel(protein_ratio_MvsL(1,:));
        AK2= isnan(protein_ratio_MvsL(1,:));
        
        %find NaN values flanked by numbers
        for i=2:(Count-1)
          if 1==isnan(protein_ratio_MvsL(1,i))
            Test1=AK2((i-1):(i+1));
            if  1==isequal(~Test1,[1 0 1])
              Number=protein_ratio_MvsL(1,(i-1):(i+1));
              Number2 = Number(~Test1);
              miss = mean(Number2); %average if  data point is missing
            else
              miss=NaN;
            end
            protein_ratio_MvsL(1,i) = miss;
          else
          end
        end
        
        %Fix extreme value
        for i=2:(Count-1)
          Test2=AK2((i-1):(i+1));
          if  1==isequal(~Test2,[1 1 1])
            
            %Create value correspond to a value in front and behide
            Number2 =protein_ratio_MvsL(1,(i-1):(i+1));
            
            if (protein_ratio_MvsL(1,i)/sum(Number2))>0.95
              %If the value is extreme, filter by avaering the
              %flanking numbers
              protein_ratio_MvsL(1,i) = (protein_ratio_MvsL(1,i-1)+protein_ratio_MvsL(1,i+1))/2;
              
            end
          end
        end
        
        %Remove values below 0.2
        protein_ratio_MvsL(protein_ratio_MvsL<0.2) = NaN;
        
        %% Consecutive numbers, if less then 5 consecutive number removes
        %% chromogram and replace with 0.05
        Multi_NaN = NaN(1,10);
        Sig=[Multi_NaN protein_ratio_MvsL Multi_NaN];
        AKK= ~(isnan(Sig));
        Count2= numel(Sig);
        
        Sig2=zeros(1,Count2); %Allocate variable
        for i= 11:(Count2-10)
          if 5== sum(AKK((i-(4)):i));
            Sig2(1,i) = Sig(i);
          elseif 5== sum(AKK((i-(3)):(i+1)));
            Sig2(1,i) = Sig(i);
          elseif 5== sum(AKK((i-(2)):(i+2)));
            Sig2(1,i) = Sig(i);
          elseif 5== sum(AKK((i-(1)):(i+3)));
            Sig2(1,i) = Sig(i);
          elseif 5== sum(AKK((i):(i+4)));
            Sig2(1,i) = Sig(i);
          else Sig2(1,i) = 0.2;
          end
          
        end
        Out= Sig2(1,((10+1):(Count2-10)));
        
        % Smooth data, change vary the number of data points used for smoothing
        protein_ratio_MvsL= (smooth(Out,3))';
        protein_ratio_MvsL(protein_ratio_MvsL<0.3) = NaN;
        
        
        
        %Fix missing values and extreme values for HvsL protein values
        %Define variable
        Count =numel(protein_ratio_HvsL(1,:));
        AK2= isnan(protein_ratio_HvsL(1,:));
        
        %find NaN values flanked by numbers
        for i=2:(Count-1)
          if 1==isnan(protein_ratio_HvsL(1,i))
            Test1=AK2((i-1):(i+1));
            if  1==isequal(~Test1,[1 0 1])
              Number=protein_ratio_HvsL(1,(i-1):(i+1));
              Number2 = Number(~Test1);
              miss = mean(Number2); %average if  data point is missing
            else
              miss=NaN;
            end
            protein_ratio_HvsL(1,i) = miss;
          else
          end
        end
        
        %Fix extreme value
        for i=2:(Count-1)
          Test2=AK2((i-1):(i+1));
          if  1==isequal(~Test2,[1 1 1])
            
            %Create value correspond to a value in front and behide
            Number2 =protein_ratio_HvsL(1,(i-1):(i+1));
            
            if (protein_ratio_HvsL(1,i)/sum(Number2))>0.95
              %If the value is extreme, filter by avaering the
              %flanking numbers
              protein_ratio_HvsL(1,i) = (protein_ratio_HvsL(1,i-1)+protein_ratio_HvsL(1,i+1))/2;
              
            end
          end
        end
        
        %Remove values below 0.2
        protein_ratio_HvsL(protein_ratio_HvsL<0.2) = NaN;
        
        %% Consecutive numbers, if less then 5 consecutive number removes
        %% chromogram and replace with 0.05
        Multi_NaN = NaN(1,10);
        Sig=[Multi_NaN protein_ratio_HvsL Multi_NaN];
        AKK= ~(isnan(Sig));
        Count2= numel(Sig);
        
        Sig2=zeros(1,Count2); %Allocate variable
        for i= 11:(Count2-10)
          if 5== sum(AKK((i-(4)):i));
            Sig2(1,i) = Sig(i);
          elseif 5== sum(AKK((i-(3)):(i+1)));
            Sig2(1,i) = Sig(i);
          elseif 5== sum(AKK((i-(2)):(i+2)));
            Sig2(1,i) = Sig(i);
          elseif 5== sum(AKK((i-(1)):(i+3)));
            Sig2(1,i) = Sig(i);
          elseif 5== sum(AKK((i):(i+4)));
            Sig2(1,i) = Sig(i);
          else Sig2(1,i) = 0.2;
          end
          
        end
        Out= Sig2(1,((10+1):(Count2-10)));
        
        % Smooth data, change vary the number of data points used for smoothing
        protein_ratio_HvsL= (smooth(Out,3))';
        protein_ratio_HvsL(protein_ratio_HvsL<0.3) = NaN;
        
      end
      
      %Define peptide array to write to
      peptide_ratio=NaN(number_idx,length(protein_ratio));
      
      %loop to extract ratio all peptides
      for ci=1:number_idx
        Protein_group.Sequence(ci,1)=Peptide_content(idx(ci),1);
        Protein_group.Proteins_assigned_in_peptide_table(ci,1)=Peptide_content(idx(ci),4);
        Protein_group.mass(ci,1)=Peptide_content(idx(ci),3);
        Protein_group.calibrated_time(ci,1)=Peptide_content(idx(ci),5);
        Protein_group.modification(ci,1)=Peptide_content(idx(ci),2);
        
        %write ratio for each peptide
        for ci3=1:length(protein_ratio);
          peptide_ratio(ci,ci3)=str2double(Peptide_isotope_ratio_values{idx(ci)+1,ci3});
        end
        
        if isotope_channel==3
          for ci3=1:length(protein_ratio);
            peptide_ratio_HvsL(ci,ci3)=str2double(Peptide_isotope_ratio_values_HvsL{idx(ci)+1,ci3});
            peptide_ratio_MvsL(ci,ci3)=str2double(Peptide_isotope_ratio_values_MvsL{idx(ci)+1,ci3});
          end
        end
        
        
        % Fix missing values and extreme values for peptide
        %Define variable
        Count =numel(peptide_ratio(ci,:));
        AK2= isnan(peptide_ratio(ci,:));
        Output_mean=zeros(1,Count); %Allocate Output_mean
        
        %find NaN values flanked by numbers
        for i=2:(Count-1)
          if 1==isnan(peptide_ratio(ci,i))
            Test1=AK2((i-1):(i+1));
            if  1==isequal(~Test1,[1 0 1])
              Number=peptide_ratio(ci,(i-1):(i+1));
              Number2 = Number(~Test1);
              miss = mean(Number2); %average if  data point is missing
            else
              miss=NaN;
            end
            peptide_ratio(ci,i) = miss;
          else
          end
        end
        
        %Fix extreme value
        for i=2:(Count-1)
          Test2=AK2((i-1):(i+1));
          if  1==isequal(~Test2,[1 1 1])
            
            %Create value correspond to a value in front and behide
            Number2 =peptide_ratio(ci,(i-1):(i+1));
            
            if (peptide_ratio(ci,i)/sum(Number2))>0.95
              %If the value is extreme, filter by avaering the
              %flanking numbers
              peptide_ratio(ci,i) = (peptide_ratio(ci,i-1)+peptide_ratio(ci,i+1))/2;
            end
          end
        end
      end
      
      %Copy raw data to strcuture
      Protein_group.peptide_ratio_raw=peptide_ratio;
      Protein_group.peptide_ratio_MvsL=peptide_ratio_MvsL;
      Protein_group.peptide_ratio_HvsL=peptide_ratio_HvsL;
      
      %Unmodified values for plotting
      peptide_ratio_for_figure=peptide_ratio;
      
      
      %Remove values below 0.2
      peptide_ratio(peptide_ratio<0.1) = NaN;
      Protein_group.peptide_ratio_fixed=peptide_ratio;
      
      
      %% create Figures
      [Find_FASTA_location]=ind2sub(size(Mouse_fasta_header), strmatch(Protein_IDs{Protein_ID_counter1,1}, Mouse_fasta_header, 'exact'));
      %Copy sequence
      sequence=Mouse_fasta_file(Find_FASTA_location).Sequence;
      
      try
        fprintf(fid_Summary_uniprot, '%s,%s,%s,\n',...
          Protein_group.Major_protein_group_name{1}, Protein_IDs{Protein_ID_counter1,1},sequence);
      catch
        fprintf(fid_Summary_uniprot, '%s,%s,%s,%s\n',...
          Protein_group.Major_protein_group_name{1}, Protein_IDs{Protein_ID_counter1,1});
      end
      
      %Find start position of peptide in protein
      Protein_group.Start_position=cell(length(Protein_group.Sequence),1);
      Protein_group.End_position=cell(length(Protein_group.Sequence),1);
      Protein_group.Semi_tryptic=cell(length(Protein_group.Sequence),3);
      
      % (For each peptide?)
      for i=1:length(Protein_group.Sequence)
        %Determine position of peptide in protein
        Protein_group.Start_position{i} = strfind(sequence, Protein_group.Sequence{i});
        Protein_group.End_position{i}= (Protein_group.Start_position{i}+length(Protein_group.Sequence{i})-1);
        
        %test if the proceeding and last amino acid
        if Protein_group.Start_position{i}>1 & Protein_group.End_position{i}<length(sequence)
          
          if ((sequence(Protein_group.Start_position{i}-1)=='R'|  sequence(Protein_group.Start_position{i}-1)=='K')...
              & (sequence(Protein_group.End_position{i})=='R'|  sequence(Protein_group.End_position{i})=='K'))
            Protein_group.Semi_tryptic(i,1)={'No'};
          elseif ~((sequence(Protein_group.Start_position{i}-1)=='R'|  sequence(Protein_group.Start_position{i}-1)=='K')...
              & (sequence(Protein_group.End_position{i})=='R'|  sequence(Protein_group.End_position{i})=='K'))
            Protein_group.Semi_tryptic(i,1)={'Yes'};
          end
        else Protein_group.Semi_tryptic(i,1)={'NaN'};
          
        end
        
        %test if the proceeding amino acid is tryptic
        if Protein_group.Start_position{i}>1
          
          if (sequence(Protein_group.Start_position{i}-1)=='R'|  sequence(Protein_group.Start_position{i}-1)=='K')
            Protein_group.Semi_tryptic(i,2)={'No'};
          elseif ~(sequence(Protein_group.Start_position{i}-1)=='R'|  sequence(Protein_group.Start_position{i}-1)=='K')
            Protein_group.Semi_tryptic(i,2)={'Yes'};
          end
        else Protein_group.Semi_tryptic(i,2)={'N-terminal'};
        end
        
        %test if the proceeding amino acid is tryptic
        if Protein_group.End_position{i}<length(sequence)
          
          if (sequence(Protein_group.End_position{i})=='R'|  sequence(Protein_group.End_position{i})=='K')
            Protein_group.Semi_tryptic(i,3)={'No'};
          elseif ~(sequence(Protein_group.End_position{i})=='R'|  sequence(Protein_group.End_position{i})=='K')
            Protein_group.Semi_tryptic(i,3)={'Yes'};
          end
        else Protein_group.Semi_tryptic(i,3)={'C-terminal'};
        end
      end
      
      
      %Determine if the observed peptides has a different profile
      Protein_group.mean_peptide_fix=zeros(1,length(protein_ratio));
      Protein_group.std_peptide_fix=zeros(1,length(protein_ratio));
      Protein_group.Number_peptide_pre_fraction_fix=zeros(1,length(protein_ratio));
      
      %determine the average and std of all peptides observed from fixed data
      
      for ci4=1:length(protein_ratio)
        
        %create temp values
        temp_values=Protein_group.peptide_ratio_fixed(:,ci4);
        
        %remove NaN values
        remove_values=isnan(temp_values);
        temp_values(remove_values)=[];
        
        Protein_group.Number_peptide_pre_fraction_fix(ci4)=numel(temp_values);
        Protein_group.mean_peptide_fix(ci4)=mean(temp_values);
        Protein_group.std_peptide_fix(ci4)=std(temp_values);
        
      end
      
      %Cover all values equal to zero to NaN
      Protein_group.std_peptide_fix(Protein_group.std_peptide_fix==0)=NaN;
      Protein_group.Number_peptide_pre_fraction_fix(Protein_group.Number_peptide_pre_fraction_fix==0)=NaN;
      
      %Determine if the observed peptides has a different profile
      Protein_group.mean_peptide=zeros(1,length(protein_ratio));
      Protein_group.std_peptide=zeros(1,length(protein_ratio));
      Protein_group.Number_peptide_pre_fraction=zeros(1,length(protein_ratio));
      
      
      %determine the average and std of all peptides observed from raw data
      for ci4=1:length(protein_ratio)
        
        %create temp values
        temp_values=Protein_group.peptide_ratio_raw(:,ci4);
        
        %remove NaN values
        remove_values=isnan(temp_values);
        temp_values(remove_values)=[];
        
        Protein_group.Number_peptide_pre_fraction(ci4)=numel(temp_values);
        Protein_group.mean_peptide(ci4)=mean(temp_values);
        Protein_group.std_peptide(ci4)=std(temp_values);
        
      end
      
      %Cover all values equal to zero to NaN
      Protein_group.std_peptide(Protein_group.std_peptide==0)=NaN;
      Protein_group.Number_peptide_pre_fraction(Protein_group.Number_peptide_pre_fraction==0)=NaN;
      
      %Compare measurement to the mean and std, raw data
      Protein_group.Differ_mean=NaN(number_idx,length(protein_ratio));
      Protein_group.Differ_mean_std=zeros(number_idx,length(protein_ratio));
      for ci5=1:number_idx
        for ci6=1:length(protein_ratio)
          Protein_group.Differ_mean(ci5,ci6)=abs(Protein_group.peptide_ratio_raw(ci5,ci6)-Protein_group.mean_peptide(1,ci6));
          if Protein_group.Differ_mean(ci5,ci6)> 1.96*Protein_group.std_peptide(1,ci6)
            Protein_group.Differ_mean_std(ci5,ci6)=1;
          end
        end
      end
      
      %Compare measurement to the mean and std, fixed data
      Protein_group.Differ_mean_fix=NaN(number_idx,length(protein_ratio));
      Protein_group.Differ_mean_std_fix=zeros(number_idx,length(protein_ratio));
      for ci5=1:number_idx
        for ci6=1:length(protein_ratio)
          Protein_group.Differ_mean_fix(ci5,ci6)=abs(Protein_group.peptide_ratio_fixed(ci5,ci6)-Protein_group.mean_peptide_fix(1,ci6));
          if Protein_group.Differ_mean_fix(ci5,ci6)> 1.96*Protein_group.std_peptide_fix(1,ci6)
            Protein_group.Differ_mean_std_fix(ci5,ci6)=1;
          end
        end
      end
      
      
      %Check if three of more values in the peptide measurements are greater then
      %the std raw data
      Protein_group.Peptide_different_to_mean=cell(number_idx,1);
      Protein_group.Global_Peptide_different_to_mean=cell(1,1);
      for ci5=1:number_idx
        for ci6=3:(Define_experiment_fractions_end(replicate_counter)-(Define_experiment_fractions_start(replicate_counter)+2))
          
          %check all values and write + if peptide measurements are different
          if 3==sum(Protein_group.Differ_mean_std(ci5,((Define_experiment_fractions_start(replicate_counter)+ci6-2):(Define_experiment_fractions_start(replicate_counter)+ci6))));
            Protein_group.Peptide_different_to_mean(ci5,1)={'+'};
            Protein_group.Global_Peptide_different_to_mean(1,1)={'+'};
          elseif 3==sum(Protein_group.Differ_mean_std(ci5,((Define_experiment_fractions_start(replicate_counter)+ci6-1):(Define_experiment_fractions_start(replicate_counter)+ci6+1))));
            Protein_group.Peptide_different_to_mean(ci5,1)={'+'};
            Protein_group.Global_Peptide_different_to_mean(1,1)={'+'};
          elseif 3==sum(Protein_group.Differ_mean_std(ci5,((Define_experiment_fractions_start(replicate_counter)+ci6):(Define_experiment_fractions_start(replicate_counter)+ci6+2))));
            Protein_group.Peptide_different_to_mean(ci5,1)={'+'};
            Protein_group.Global_Peptide_different_to_mean(1,1)={'+'};
          end
          
        end
      end
      
      %Check if three of more values in the peptide measurements are greater then
      %the std fixed
      Protein_group.Peptide_different_to_mean_fixed=cell(number_idx,1);
      Protein_group.Global_Peptide_different_to_mean_fix=cell(1,1);
      for ci5=1:number_idx
        for ci6=3:(Define_experiment_fractions_end(replicate_counter)-(Define_experiment_fractions_start(replicate_counter)+2))
          
          %check all values and write + if peptide measurements are different
          if 3==sum(Protein_group.Differ_mean_std_fix(ci5,((Define_experiment_fractions_start(replicate_counter)+ci6-2):(Define_experiment_fractions_start(replicate_counter)+ci6))));
            Protein_group.Peptide_different_to_mean_fix(ci5,1)={'+'};
            Protein_group.Global_Peptide_different_to_mean_fix(1,1)={'+'};
          elseif 3==sum(Protein_group.Differ_mean_std_fix(ci5,((Define_experiment_fractions_start(replicate_counter)+ci6-1):(Define_experiment_fractions_start(replicate_counter)+ci6+1))));
            Protein_group.Peptide_different_to_mean_fix(ci5,1)={'+'};
            Protein_group.Global_Peptide_different_to_mean_fix(1,1)={'+'};
          elseif 3==sum(Protein_group.Differ_mean_std_fix(ci5,((Define_experiment_fractions_start(replicate_counter)+ci6):(Define_experiment_fractions_start(replicate_counter)+ci6+2))));
            Protein_group.Peptide_different_to_mean_fix(ci5,1)={'+'};
            Protein_group.Global_Peptide_different_to_mean_fix(1,1)={'+'};
          end
          
        end
      end
      
      %Clustering approach to assess if peptides are different
      Temp_value=Protein_group.peptide_ratio_raw(Define_experiment_fractions_start(replicate_counter):Define_experiment_fractions_end(replicate_counter));
      Temp_value(isnan(Temp_value))=0
      try
        cluster_analysis=cophenet(linkage(squareform(pdist(Temp_value)),'complete','jaccard'),squareform(pdist(Temp_value)));
        
      catch
        cluster_analysis=NaN
      end
      
      %Clustering approach to assess if peptides are different
      Temp_value=Protein_group.peptide_ratio_fixed(Define_experiment_fractions_start(replicate_counter):Define_experiment_fractions_end(replicate_counter));
      Temp_value(isnan(Temp_value))=0
      
      try
        cluster_analysis_fix=cophenet(linkage(squareform(pdist(Temp_value)),'complete','jaccard'),squareform(pdist(Temp_value)));
      catch
        cluster_analysis_fix=NaN
      end
      
      
      %Clear data to remove inf values
      if isinf(cluster_analysis)==1
        cluster_analysis=NaN;
      end
      
      %Clear data to remove inf values
      if isinf(cluster_analysis_fix)==1
        cluster_analysis_fix=NaN;
      end
      
      %Use a modified Shannon Wiener Diversity Index to determine evenest of
      %peptides within each protein group
      
      %Formula:
      %H = -SUM[(pi) * ln(pi)]
      %E=H/Hmax
      
      %create arrya to store values
      H_value=(Define_experiment_fractions_end(replicate_counter)-(Define_experiment_fractions_start(replicate_counter)-1));
      H_max=(Define_experiment_fractions_end(replicate_counter)-(Define_experiment_fractions_start(replicate_counter)-1));
      Combined_H_value=NaN((Define_experiment_fractions_end(replicate_counter)-(Define_experiment_fractions_start(replicate_counter)-1)),1);
      Combined_Evenness_values=NaN((Define_experiment_fractions_end(replicate_counter)-(Define_experiment_fractions_start(replicate_counter)-1)),1);
      
      Protein_group.Shannon_wiener_evenest=zeros(1,1);
      Protein_group.Shannon_wiener_diversity=zeros(1,1);
      
      for ci6=Define_experiment_fractions_start(replicate_counter):Define_experiment_fractions_end(replicate_counter)
        %determine value to test
        Temp_value=Protein_group.peptide_ratio_raw(:,ci6);
        %remove NaN values
        Temp_value(isnan(Temp_value))=[];
        %create internal counter to write to
        counter=ci6-(Define_experiment_fractions_start(replicate_counter)-1);
        %store H value
        H_value = Temp_value .* log(Temp_value);
        %store H_max value
        H_max= max(Protein_group.peptide_ratio_raw(:,ci6));
        %caliculate the even-ness of observed peptides
        Combined_Evenness_values(counter)=-sum(H_value)/H_max;
        if ~isempty(H_value)
          Combined_H_value(counter)=mean(H_value);
        end
      end
      %remove NaN values
      Combined_Evenness_values(isnan(Combined_Evenness_values))=[];
      Combined_H_value(isnan(Combined_H_value))=[];
      Protein_group.Shannon_wiener_evenest(1,1)=abs(mean(Combined_Evenness_values));
      Protein_group.Shannon_wiener_diversity(1,1)=abs(mean(Combined_H_value));
      
      
      %For proteins with high unevenness extract peptide information and
      %generate new profile
      
      if isotope_channel==3
        %Determine figure dimension
        Peptide_figure_counter=size(Protein_group.peptide_ratio_fixed);
        %Define rations to consider
        Experiment1_fraction_num=Define_experiment_fractions_start(replicate_counter):Define_experiment_fractions_end(replicate_counter);
        
        %Copy raw data to a form which we can determine the mean of protein
        %fragments which change
        Protein_group.peptide_ratio_raw2=Protein_group.peptide_ratio_raw;
        Protein_group.peptide_ratio_raw2(Protein_group.peptide_ratio_raw2>0.5 & Protein_group.peptide_ratio_raw2<2)=NaN;
        
        %Create array to write to if protein is uneveness with high diversity
        Protein_group.protein_ratio_fixed_subdivided=nan(1,length(protein_ratio));
        Protein_group.peptide_ratio_fixed_subdivided_HvsM=nan(Peptide_figure_counter(1),Peptide_figure_counter(2));
        Protein_group.peptide_ratio_fixed_subdivided_MvsL=nan(Peptide_figure_counter(1),Peptide_figure_counter(2));
        Protein_group.peptide_ratio_fixed_subdivided_HvsL=nan(Peptide_figure_counter(1),Peptide_figure_counter(2));
        
        if abs(mean(Combined_Evenness_values))> 1 & abs(mean(Combined_H_value)) > 1
          
          for ii=1:Peptide_figure_counter(1)
            for iii=Experiment1_fraction_num(1):Experiment1_fraction_num(end)
              %check if peptide has a ratio values
              if ~isnan(Protein_group.peptide_ratio_raw2(ii,iii)) & (nanmean(Protein_group.peptide_ratio_raw(:,iii))>=2 | nanmean(Protein_group.peptide_ratio_raw(:,iii))<=0.5)
                
                Protein_group.peptide_ratio_fixed_subdivided_HvsM(ii,iii)=Protein_group.peptide_ratio_raw(ii,iii);
                Protein_group.peptide_ratio_fixed_subdivided_HvsL(ii,iii)=Protein_group.peptide_ratio_HvsL(ii,iii);
                Protein_group.peptide_ratio_fixed_subdivided_MvsL(ii,iii)=Protein_group.peptide_ratio_MvsL(ii,iii);
                
              end
            end
          end
          
          %Determine if changes observed in the peptides are consistent across
          %peptide measurement within a singel fraction
          for i=1:Peptide_figure_counter(1)
            AK2=isnan(Protein_group.peptide_ratio_fixed_subdivided_HvsM(i,:));
            for iii=(Experiment1_fraction_num(1)+1):(Experiment1_fraction_num(end)-1)
              if 1==isnan(Protein_group.peptide_ratio_fixed_subdivided_HvsM(i,iii))
                Test1=AK2((iii-1):(iii+1));
                if  1==isequal(~Test1,[0 0 1]) | 1==isequal(~Test1,[1 0 1])
                  Protein_group.peptide_ratio_fixed_subdivided_HvsM(i,iii) = NaN;
                  Protein_group.peptide_ratio_fixed_subdivided_HvsL(i,iii) = NaN;
                  Protein_group.peptide_ratio_fixed_subdivided_MvsL(i,iii) = NaN;
                else
                end
              end
            end
          end
          
          
          for ii=1:Peptide_figure_counter(2)
            Protein_group.protein_ratio_fixed_subdivided_HvsM(1,ii)=nanmedian(Protein_group.peptide_ratio_fixed_subdivided_HvsM(:,ii));
            Protein_group.protein_ratio_fixed_subdivided_HvsL(1,ii)=nanmedian(Protein_group.peptide_ratio_fixed_subdivided_HvsL(:,ii));
            Protein_group.protein_ratio_fixed_subdivided_MvsL(1,ii)=nanmedian(Protein_group.peptide_ratio_fixed_subdivided_MvsL(:,ii));
            
          end
          
        end
        
      end
      
      
      if isotope_channel~=3
        
        Figure_isotope=figure;
        subplot1=subplot(1,5,1);
        P1A =plot(protein_ratio_for_figure(Define_experiment_fractions_start(replicate_counter):Define_experiment_fractions_end(replicate_counter)),...
          (Define_experiment_fractions_start(replicate_counter):Define_experiment_fractions_end(replicate_counter)));
        set(P1A,'Color', [0.5 0.5 0.5],'LineWidth',0.5);
        hold on
        P1B =plot(protein_ratio_for_figure(Define_experiment_fractions_start(replicate_counter):Define_experiment_fractions_end(replicate_counter)),...
          (Define_experiment_fractions_start(replicate_counter):Define_experiment_fractions_end(replicate_counter)),'d');
        set(P1B,'MarkerFaceColor',colour_to_use(2,:),'MarkerSize',4);
        hold on %add 0.2 dotted line gaussian under which qunation is not considered
        PP1B= plot([repmat(0.2,...
          ((Define_experiment_fractions_end(replicate_counter)+2)-(Define_experiment_fractions_start(replicate_counter)-3)),1)],...
          (Define_experiment_fractions_start(replicate_counter)-2:Define_experiment_fractions_end(replicate_counter)+2),':');
        set(PP1B,'Color','red','LineWidth',1);
        ylim([Define_experiment_fractions_start(replicate_counter),Define_experiment_fractions_end(replicate_counter)]);
        %Set Ytick to fraction values
        set(gca,'YTick', linspace(Define_experiment_fractions_start(replicate_counter),Define_experiment_fractions_end(replicate_counter),7));
        
        %Set tick values as fraction 1 to 55
        xData=mat2cell(linspace(1,55,7),1);
        set(gca,'YTickLabel',xData);
        
        set(gca,'xdir','reverse')
        
        %Define rations to consider
        Experiment1_fraction_num=Define_experiment_fractions_start(replicate_counter):Define_experiment_fractions_end(replicate_counter);
        
        %Change the Axis of figure
        try
          axis([0,round(max(max(protein_ratio_for_figure(:,Experiment1_fraction_num(1):Experiment1_fraction_num(end)))))*1.2,...
            Define_experiment_fractions_start(replicate_counter)-2, Define_experiment_fractions_end(replicate_counter)+2]);
          
        catch
          axis([0,10,...
            Define_experiment_fractions_start(replicate_counter)-2, Define_experiment_fractions_end(replicate_counter)+2]);
        end
        
        %add figure information
        xlabel('Isotopologue ratio','FontSize', 6);
        ylabel('Fractions','FontSize', 6);
        
        title({'Protein Quantitation'},'FontSize', 8);
        
        % draw map of peptides covering protein
        subplot2=subplot(1,5,2:5);
        axis([-10 length(sequence) Define_experiment_fractions_start(replicate_counter)-2 Define_experiment_fractions_end(replicate_counter)+2]);
        
        %Determine figure dimension
        Peptide_figure_counter=size(Protein_group.peptide_ratio_fixed);
        
        custome_scale = [-3 -2.5 -2.25 -2 -1.8 -1.6 -1.4 -1.2 -1.0 -0.9 -0.8 -0.7 -0.6 -0.5 -0.4 -0.3 -0.2 -0.1 0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1 1.1...
          1.2 1.3 1.4 1.5 1.6 1.7 1.8 1.9 2 2.25 2.5 2.75 3 3.25 3.5 3.75 4 4.5 5];
        custome_scale_number=[1:1:length(custome_scale)];
        custome_scale_for_figure=[0.1 0.25 0.5 1 4 9 16 25 100];
        cmap=jet(length(custome_scale));
        
        %Define maximum values
        Maximum_observed_value=0;
        
        for ii=1:Peptide_figure_counter(1)
          try
            hold on
            for iii=Experiment1_fraction_num(1):Experiment1_fraction_num(end)
              %check if peptide has a ratio values
              if ~isnan(Protein_group.peptide_ratio_fixed(ii,iii))
                
                %Create colour map index
                cm_index=log2(Protein_group.peptide_ratio_fixed(ii,iii));
                
                %find closest number
                tmp = abs(custome_scale-cm_index)
                [idx_for_scale idx_for_scale] = min(tmp) %index of closest value
                closest = custome_scale_number(idx_for_scale) %closest value
                
                %create figure
                rectangle('Position',[Protein_group.Start_position{ii},iii,...
                  (Protein_group.End_position{ii}-Protein_group.Start_position{ii}),0.8],'FaceColor',cmap(closest,:));
                
                %test if value is higher then current maximum value
                if Maximum_observed_value<Protein_group.peptide_ratio_fixed(ii,iii)
                  Maximum_observed_value=Protein_group.peptide_ratio_fixed(ii,iii)
                end
                
              end
            end
            
          catch % if can't find peptide sequence will skip
          end
          %Add semi-tryptic lines
          if isequal(Protein_group.Semi_tryptic(ii,1),'Yes')
            hold on
            % if N-terminal semi-tryptic cleavage
            if  isequal(Protein_group.Semi_tryptic(ii,2),'Yes')
              semitryptic_line= plot([repmat((Protein_group.Start_position{ii}-1), 1,...
                (Define_experiment_fractions_end(replicate_counter)-(Define_experiment_fractions_start(replicate_counter)-1)))],...
                (Define_experiment_fractions_start(replicate_counter):Define_experiment_fractions_end(replicate_counter)),':');
            end
            if  isequal(Protein_group.Semi_tryptic(ii,3),'Yes')
              semitryptic_line= plot([repmat((Protein_group.End_position{ii}+1), 1,...
                (Define_experiment_fractions_end(replicate_counter)-(Define_experiment_fractions_start(replicate_counter)-1)))],...
                (Define_experiment_fractions_start(replicate_counter):Define_experiment_fractions_end(replicate_counter)),':');
            end
          end
          
        end
        
        %Set Ytick to fraction values
        set(gca,'YTick', linspace(Define_experiment_fractions_start(replicate_counter),Define_experiment_fractions_end(replicate_counter),7));
        
        %Set tick values as fraction 1 to 55
        xData=mat2cell(repmat({''},1,7),1);
        set(gca,'YTickLabel',xData{:});
        
        %add figure information
        xlabel('Protein Amino Acid Sequence','FontSize', 6);
        
        title({'Peptide Quantitation'},'FontSize', 8);
        
        %adds the colorbar to your plot
        colour_bar=colorbar
        colormap(cmap);
        initpos = get(colour_bar,'Position');
        set(colour_bar,'YTickLabel',round(linspace(0,Maximum_observed_value,11)*10)/10,'Position',[(initpos(1)+initpos(3)*0.95) (initpos(2)+initpos(4)/4)...
          initpos(3)*0.5 initpos(4)*0.5], 'FontSize', 6);
        ylabel(colour_bar,'Peptide Isotopologue ratio','Rotation',90.0,'FontSize', 8);
        
        %Change size of subplot
        subplot1_positions = get(subplot1, 'position');
        set(subplot1, 'position', [subplot1_positions(1), subplot1_positions(2), subplot1_positions(3)*1.22, subplot1_positions(4)]);
        
        subplot2_positions = get(subplot2, 'position');
        set(subplot2, 'position', [subplot2_positions(1), subplot1_positions(2), subplot2_positions(3)*0.9, subplot2_positions(4)]);
        
        %Add master title of figure
        suptitle(strcat('Protein :',Protein_IDs{Protein_ID_counter1,1},' Biological replicate :',mat2str(replicate_counter),' isotoplogue channel :',isotope_name_value{1}));
        
        %Save imageProtein_IDs(Protein_ID_counter1,:)))
        filename1 = strcat(Protein_IDs{Protein_ID_counter1,1},'_replicate',mat2str(replicate_counter),'_',isotope_name_value{1},'_.pdf');
        print(Figure_isotope(1),'-dpdf','-r600', filename1)
        
        %Close finger to ensure the script does not run out of memory
        clearvars Figure
        close all
        
        movefile(filename1, 'ProteinGroups image_output');
        
      elseif isotope_channel==3
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Creat figure showing MvsL and HvsL with only peptides that changed
        %highlighted
        Figure_comparsion=figure;
        subplot1=subplot(1,5,1);
        %MvsL
        P1A =plot(protein_ratio_MvsL(Define_experiment_fractions_start(replicate_counter):Define_experiment_fractions_end(replicate_counter)),...
          (Define_experiment_fractions_start(replicate_counter):Define_experiment_fractions_end(replicate_counter)));
        set(P1A,'Color', [0.5 0.5 0.5],'LineWidth',0.5);
        hold on
        P1B =plot(protein_ratio_MvsL(Define_experiment_fractions_start(replicate_counter):Define_experiment_fractions_end(replicate_counter)),...
          (Define_experiment_fractions_start(replicate_counter):Define_experiment_fractions_end(replicate_counter)),'d');
        set(P1B,'MarkerFaceColor',colour_to_use(2,:),'MarkerSize',4);
        hold on
        %HvsL
        P2A =plot(protein_ratio_HvsL(Define_experiment_fractions_start(replicate_counter):Define_experiment_fractions_end(replicate_counter)),...
          (Define_experiment_fractions_start(replicate_counter):Define_experiment_fractions_end(replicate_counter)));
        set(P2A,'Color', [0.5 0.5 0.5],'LineWidth',0.5);
        hold on
        P2B =plot(protein_ratio_HvsL(Define_experiment_fractions_start(replicate_counter):Define_experiment_fractions_end(replicate_counter)),...
          (Define_experiment_fractions_start(replicate_counter):Define_experiment_fractions_end(replicate_counter)),'d');
        set(P2B,'MarkerFaceColor',colour_to_use(5,:),'MarkerSize',4);
        hold on
        %HvsL after filtering
        try
          P2C =plot(Protein_group.protein_ratio_fixed_subdivided_HvsL(Define_experiment_fractions_start(replicate_counter):Define_experiment_fractions_end(replicate_counter)),...
            (Define_experiment_fractions_start(replicate_counter):Define_experiment_fractions_end(replicate_counter)));
          set(P2C,'Color', [0.5 0.5 0.5],'LineWidth',0.5);
          hold on
          P2D =plot(Protein_group.protein_ratio_fixed_subdivided_HvsL(Define_experiment_fractions_start(replicate_counter):Define_experiment_fractions_end(replicate_counter)),...
            (Define_experiment_fractions_start(replicate_counter):Define_experiment_fractions_end(replicate_counter)),'d');
          set(P2D,'MarkerFaceColor',colour_to_use(6,:),'MarkerSize',4);
          hold on
        catch
        end
        %MvsL after filtering
        try
          P2E =plot(Protein_group.protein_ratio_fixed_subdivided_MvsL(Define_experiment_fractions_start(replicate_counter):Define_experiment_fractions_end(replicate_counter)),...
            (Define_experiment_fractions_start(replicate_counter):Define_experiment_fractions_end(replicate_counter)));
          set(P2E,'Color', [0.5 0.5 0.5],'LineWidth',0.5);
          hold on
          P2F =plot(Protein_group.protein_ratio_fixed_subdivided_MvsL(Define_experiment_fractions_start(replicate_counter):Define_experiment_fractions_end(replicate_counter)),...
            (Define_experiment_fractions_start(replicate_counter):Define_experiment_fractions_end(replicate_counter)),'d');
          set(P2F,'MarkerFaceColor',colour_to_use(7,:),'MarkerSize',4);
        catch
        end
        hold on
        %add 0.2 dotted line gaussian under which qunation is not considered
        PP1B= plot([repmat(0.2,...
          ((Define_experiment_fractions_end(replicate_counter)+2)-(Define_experiment_fractions_start(replicate_counter)-3)),1)],...
          (Define_experiment_fractions_start(replicate_counter)-2:Define_experiment_fractions_end(replicate_counter)+2),':');
        set(PP1B,'Color','red','LineWidth',1);
        ylim([Define_experiment_fractions_start(replicate_counter),Define_experiment_fractions_end(replicate_counter)]);
        hold on;
        %Set Ytick to fraction values
        set(gca,'YTick', linspace(Define_experiment_fractions_start(replicate_counter),Define_experiment_fractions_end(replicate_counter),7));
        
        %Set tick values as fraction 1 to 55
        xData=mat2cell(linspace(1,55,7),1);
        set(gca,'YTickLabel',xData);
        
        set(gca,'xdir','reverse')
        
        %add figure information
        xlabel('Isotopologue ratio','FontSize', 6);
        ylabel('Fractions','FontSize', 6);
        
        title({'Protein Quantitation'},'FontSize', 8);
        
        %Define rations to consider
        Experiment1_fraction_num=Define_experiment_fractions_start(replicate_counter):Define_experiment_fractions_end(replicate_counter);
        
        %Change the Axis of figure
        try
          axis([0,ceil(max([max(protein_ratio_MvsL) max(protein_ratio_HvsL)])*1.2),...
            Define_experiment_fractions_start(replicate_counter)-2, Define_experiment_fractions_end(replicate_counter)+2]);
          
        catch
          axis([0,10,...
            Define_experiment_fractions_start(replicate_counter)-2, Define_experiment_fractions_end(replicate_counter)+2]);
        end
        
        % draw map of peptides covering protein
        subplot2=subplot(1,5,2:5);
        axis([-10 length(sequence) Define_experiment_fractions_start(replicate_counter)-2 Define_experiment_fractions_end(replicate_counter)+2]);
        
        %Determine figure dimension
        Peptide_figure_counter=size(Protein_group.peptide_ratio_fixed);
        
        
        for ii=1:Peptide_figure_counter(1)
          
          try
            hold on
            for iii=Experiment1_fraction_num(1):Experiment1_fraction_num(end)
              %check if peptide has a ratio values
              if ~isnan(Protein_group.peptide_ratio_fixed(ii,iii))
                
                %Create colour map index
                cm_index=log2(Protein_group.peptide_ratio_fixed(ii,iii));
                
                if cm_index>1
                  %Set values for HvsM over 1
                  custome_scale1 = [1 1.2 1.4 1.6 1.8 2 2.25 2.5 3 3.5 4 4.5 5];
                  custome_scale_number=[1:1:length(custome_scale1)];
                  cmap1=hot(length(custome_scale1));
                  custome_scale_for_figure=[1 5];
                  caxis(custome_scale_for_figure)
                  
                  %find closest number
                  tmp = abs(custome_scale1-cm_index)
                  [idx_for_scale idx_for_scale] = min(tmp) %index of closest value
                  closest_high = custome_scale_number(idx_for_scale) %closest value
                  
                  %create figure
                  rectangle('Position',[Protein_group.Start_position{ii},iii,...
                    (Protein_group.End_position{ii}-Protein_group.Start_position{ii}),0.8],'FaceColor',cmap1(closest_high,:));
                elseif cm_index<-1
                  %Set values for HvsM over 1
                  custome_scale2 = [-1 -1.2 -1.4 -1.6 -1.8 -2 -2.25 -2.5 -3 -3.5 -4 -4.5 -5];
                  custome_scale_number=[1:1:length(custome_scale2)];
                  cmap2=bone(length(custome_scale2));
                  custome_scale_for_figure=[1 5];
                  caxis(custome_scale_for_figure)
                  
                  %find closest number
                  tmp = abs(custome_scale2-cm_index)
                  [idx_for_scale idx_for_scale] = min(tmp) %index of closest value
                  closest_low = custome_scale_number(idx_for_scale) %closest value
                  
                  %create figure
                  rectangle('Position',[Protein_group.Start_position{ii},iii,...
                    (Protein_group.End_position{ii}-Protein_group.Start_position{ii}),0.8],'FaceColor',cmap2(closest_low,:));
                  
                elseif cm_index>=-1 & cm_index<=1
                  %create figure
                  rectangle('Position',[Protein_group.Start_position{ii},iii,...
                    (Protein_group.End_position{ii}-Protein_group.Start_position{ii}),0.8],'FaceColor',[1 1 1]);
                  
                end
              end
            end
            
          catch % if can't find peptide sequence will skip
          end
          
          %Add semi-tryptic lines
          if isequal(Protein_group.Semi_tryptic(ii,1),'Yes')
            hold on
            % if N-terminal semi-tryptic cleavage
            if  isequal(Protein_group.Semi_tryptic(ii,2),'Yes')
              semitryptic_line= plot([repmat((Protein_group.Start_position{ii}-1), 1,...
                (Define_experiment_fractions_end(replicate_counter)-(Define_experiment_fractions_start(replicate_counter)-1)))],...
                (Define_experiment_fractions_start(replicate_counter):Define_experiment_fractions_end(replicate_counter)),':');
            end
            if  isequal(Protein_group.Semi_tryptic(ii,3),'Yes')
              semitryptic_line= plot([repmat((Protein_group.End_position{ii}+1), 1,...
                (Define_experiment_fractions_end(replicate_counter)-(Define_experiment_fractions_start(replicate_counter)-1)))],...
                (Define_experiment_fractions_start(replicate_counter):Define_experiment_fractions_end(replicate_counter)),':');
            end
          end
          
        end
        
        %Set Ytick to fraction values
        set(gca,'YTick', linspace(Define_experiment_fractions_start(replicate_counter),Define_experiment_fractions_end(replicate_counter),7));
        
        %Set tick values as fraction 1 to 55
        xData=mat2cell(repmat({''},1,7),1);
        set(gca,'YTickLabel',xData{:});
        
        %add figure information
        xlabel('Protein Amino Acid Sequence','FontSize', 6);
        
        title({'Peptide Quantitation'},'FontSize', 8);
        
        %adds the colorbar to plot
        hold on
        colour_bar=colorbar
        
        % Build a colormap that consists of three separate
        % colormaps.
        cmap_combined = flipud([cmap1;repmat([1 1 1],5,1);flipud(cmap2)]);
        custome_scale_combined=flipud(horzcat(custome_scale1,repmat([0],1,3),custome_scale2));
        
        %Set postion and range
        initpos = get(colour_bar,'Position');
        set(colour_bar,'YTicklabel',[linspace(min(custome_scale_combined),max(custome_scale_combined),11)],....
          'YTick', [linspace(min(custome_scale_combined),max(custome_scale_combined),11)],....
          'Position',[(initpos(1)+initpos(3)*0.95) (initpos(2))...
          initpos(3)*0.5 initpos(4)], 'FontSize', 6);
        
        %label colourbar
        ylabel(colour_bar,'Peptide Isotopologue ratio (log2)','Rotation',90.0,'FontSize', 8);
        
        %set colour map to use
        colormap(cmap_combined);
        caxis([min(custome_scale_combined) max(custome_scale_combined)])
        
        
        %Change size of subplot
        subplot1_positions = get(subplot1, 'position');
        set(subplot1, 'position', [subplot1_positions(1), subplot1_positions(2), subplot1_positions(3)*1.22, subplot1_positions(4)]);
        
        subplot2_positions = get(subplot2, 'position');
        set(subplot2, 'position', [subplot2_positions(1), subplot1_positions(2), subplot2_positions(3)*0.9, subplot2_positions(4)]);
        
        %Add master title of figure
        suptitle(strcat('Protein :',Protein_IDs{Protein_ID_counter1,1},' Biological replicate :',mat2str(replicate_counter)));
        
        %Save imageProtein_IDs(Protein_ID_counter1,:)))
        filename1 = strcat(Protein_IDs{Protein_ID_counter1,1},'_replicate',mat2str(replicate_counter),'_',isotope_name_value{1},'_.pdf');
        print(Figure_comparsion(1),'-dpdf','-r600', filename1)
        
        %Close finger to ensure the script does not run out of memory
        clearvars Figure
        close all
        
        movefile(filename1, 'ProteinGroups image_output');
        
      end
      
      %% Save text related data
      
      %Experimental channel information
      for column_counter=1:(Define_experiment_fractions_end(replicate_counter)-Define_experiment_fractions_start(replicate_counter)+1)
        Column_name{column_counter}=strcat('Ratio_MvsL_',mat2str(column_counter));
      end
      
      %Create header to write out values
      Column_header1=repmat('%s,', 1, (Define_experiment_fractions_end(replicate_counter)-Define_experiment_fractions_start(replicate_counter)+14));
      Column_header2=['%s,%s,%s,',repmat('%s,', 1, (Define_experiment_fractions_end(replicate_counter)-Define_experiment_fractions_start(replicate_counter)+11))];
      
      %Create file to write all Raw protein and peptide data, raw data
      Protein_group_name1= strcat(Protein_group.Major_protein_group_name(1),'_replicate',mat2str(replicate_counter),'_',isotope_name_value{1},'_Raw_Maxquant_data.csv');
      Protein_group_output1 = fopen(Protein_group_name1{:},'a+');
      
      fprintf (Protein_group_output1,[Column_header1, '\n'],...
        'Sequence','Major Proteingroup','Proteins assignement of peptide',...
        'Modifications','Start position','End position',...
        'Semi_tryptic','N-terminal','C-terminal',...
        'Peptide with atleast three consecative values which differ by 1.96*stdev',...
        'Cophenetic correlation','Mass','Calibrated retention time', Column_name{:});
      
      %Write out protein ratios
      fprintf (Protein_group_output1,[Column_header2, '\n'],...
        'Total Protein',Protein_IDs{Protein_ID_counter1,1},...
        '','','','','','','','','','','', Protein_group.protein_ratio_raw(1,Define_experiment_fractions_start(replicate_counter):Define_experiment_fractions_end(replicate_counter)));
      
      fprintf (Protein_group_output1,[Column_header2, '\n'],...
        'Std of Peptides in each fraction',Protein_IDs{Protein_ID_counter1,1},...
        '','','','','','','','','','','', Protein_group.std_peptide(1,(Define_experiment_fractions_start(replicate_counter):Define_experiment_fractions_end(replicate_counter))));
      
      %The data strcuture of measurement row requires a convoluted write out
      %method
      fprintf (Protein_group_output1,[Column_header2, '\n'],...
        'Measurements',Protein_IDs{Protein_ID_counter1,1},...
        '','','','','','','','','','','');
      
      for counter=1:(Define_experiment_fractions_end(replicate_counter)-(Define_experiment_fractions_start(replicate_counter)-1))
        fprintf (Protein_group_output1, '%s,',...
          num2str(Protein_group.Number_peptide_pre_fraction(1,(Define_experiment_fractions_start(replicate_counter)+(counter-1)))));
      end
      
      %Jump to next line
      fprintf (Protein_group_output1,'\n'),...
        
    
    %Write peptide ratio
    number_of_peptides=length(Protein_group.Sequence);
    
    for counter=1:number_of_peptides
      fprintf (Protein_group_output1,[Column_header2, '\n'],...
        Protein_group.Sequence{counter,1}, Protein_group.Major_protein_group_name{1},...
        Protein_group.Proteins_assigned_in_peptide_table{counter,1},...
        Protein_group.modification{counter,1},num2str(Protein_group.Start_position{counter,1}),...
        num2str(Protein_group.End_position{counter,1}),Protein_group.Semi_tryptic{counter,1},...
        Protein_group.Semi_tryptic{counter,2},Protein_group.Semi_tryptic{counter,3},...
        Protein_group.Peptide_different_to_mean{counter,1},num2str(cluster_analysis(1)),...
        Protein_group.mass{counter,1},Protein_group.calibrated_time{counter,1},...
        Protein_group.peptide_ratio_raw(counter,(Define_experiment_fractions_start(replicate_counter):Define_experiment_fractions_end(replicate_counter))));
    end
    fclose(Protein_group_output1);
    
    
    %Save csv of peptide information for fixed data
    Protein_group_name2= strcat(Protein_IDs(Protein_ID_counter1, 1),'_replicate',mat2str(replicate_counter),'_',isotope_name_value{1},'_Fixed_data.csv');
    Protein_group_output2 = fopen(Protein_group_name2{:},'a+');
    
    fprintf (Protein_group_output2,[Column_header1, '\n'],...
      'Sequence','Major Proteingroup','Proteins assignement of peptide',...
      'Modifications','Start position','End position',...
      'Semi_tryptic','N-terminal','C-terminal',...
      'Peptide with atleast three consecative values which differ by 1.96*stdev',...
      'Cophenetic correlation','Mass','Calibrated retention time', Column_name{:});
    
    %Write out protein ratios
    fprintf (Protein_group_output2,[Column_header2, '\n'],...
      'Total Protein',Protein_IDs{Protein_ID_counter1,1},...
      '','','','','','','','','','','', Protein_group.protein_ratio_fixed(1,Define_experiment_fractions_start(replicate_counter):Define_experiment_fractions_end(replicate_counter)));
    
    fprintf (Protein_group_output2,[Column_header2, '\n'],...
      'Std of Peptides in each fraction',Protein_IDs{Protein_ID_counter1,1},...
      '','','','','','','','','','','', Protein_group.std_peptide_fix(1,(Define_experiment_fractions_start(replicate_counter):Define_experiment_fractions_end(replicate_counter))));
    
    %The data strcuture of measurement row requires a convoluted write out
    %method
    fprintf (Protein_group_output2,[Column_header2, '\n'],...
      'Measurements',Protein_IDs{Protein_ID_counter1,1},...
      '','','','','','','','','','','');
    
    for counter=1:(Define_experiment_fractions_end(replicate_counter)-(Define_experiment_fractions_start(replicate_counter)-1))
      fprintf (Protein_group_output2, '%s,',...
        num2str(Protein_group.Number_peptide_pre_fraction_fix(1,(Define_experiment_fractions_start(replicate_counter)+(counter-1)))));
    end
    
    %Jump to next line
    fprintf (Protein_group_output2,'\n'),...
      
  
  %Write peptide ratio
  number_of_peptides=length(Protein_group.Sequence);
  
  for counter=1:number_of_peptides
    fprintf (Protein_group_output2,[Column_header2, '\n'],...
      Protein_group.Sequence{counter,1}, Protein_group.Major_protein_group_name{1},...
      Protein_group.Proteins_assigned_in_peptide_table{counter,1},...
      Protein_group.modification{counter,1},num2str(Protein_group.Start_position{counter,1}),...
      num2str(Protein_group.End_position{counter,1}),Protein_group.Semi_tryptic{counter,1},...
      Protein_group.Semi_tryptic{counter,2},Protein_group.Semi_tryptic{counter,3},...
      Protein_group.Peptide_different_to_mean_fixed{counter,1},num2str(cluster_analysis_fix(1)),...
      Protein_group.mass{counter,1},Protein_group.calibrated_time{counter,1},...
      Protein_group.peptide_ratio_fixed(counter,(Define_experiment_fractions_start(replicate_counter):Define_experiment_fractions_end(replicate_counter))));
  end
  fclose(Protein_group_output2);
  
  
  
  if isotope_channel==3
    %Write out protein profile for just the fractions in which HvsM is
    %different
    
    Protein_group_name3= strcat(Protein_IDs(Protein_ID_counter1, 1),'_replicate',mat2str(replicate_counter),'_',isotope_name_value{1},'_Differing_HvsM_measuremnt_data.csv');
    Protein_group_output3 = fopen(Protein_group_name3{:},'a+');
    
    fprintf (Protein_group_output3,[Column_header1, '\n'],...
      'Sequence','Major Proteingroup','Proteins assignement of peptide',...
      'Modifications','Start position','End position',...
      'Semi_tryptic','N-terminal','C-terminal',...
      'Peptide with atleast three consecative values which differ by 1.96*stdev',...
      'Cophenetic correlation','Mass','Calibrated retention time', Column_name{:});
    
    %Write out protein ratios
    fprintf (Protein_group_output3,[Column_header2, '\n'],...
      'Total Protein profile for fractions where HvsM is not near zero HvsL',Protein_IDs{Protein_ID_counter1,1},...
      '','','','','','','','','','','', Protein_group.protein_ratio_fixed_subdivided_HvsL(1,Define_experiment_fractions_start(replicate_counter):Define_experiment_fractions_end(replicate_counter)));
    
    fprintf (Protein_group_output3,[Column_header2, '\n'],...
      'Total Protein profile for fractions where HvsM is not near zero MvsL',Protein_IDs{Protein_ID_counter1,1},...
      '','','','','','','','','','','', Protein_group.protein_ratio_fixed_subdivided_MvsL(1,Define_experiment_fractions_start(replicate_counter):Define_experiment_fractions_end(replicate_counter)));
    
    fprintf (Protein_group_output3,[Column_header2, '\n'],...
      'Std of Peptides in each fraction',Protein_IDs{Protein_ID_counter1,1},...
      '','','','','','','','','','','', Protein_group.std_peptide_fix(1,(Define_experiment_fractions_start(replicate_counter):Define_experiment_fractions_end(replicate_counter))));
    
    %The data strcuture of measurement row requires a convoluted write out
    %method
    fprintf (Protein_group_output3,[Column_header2, '\n'],...
      'Measurements',Protein_IDs{Protein_ID_counter1,1},...
      '','','','','','','','','','','');
    
    for counter=1:(Define_experiment_fractions_end(replicate_counter)-(Define_experiment_fractions_start(replicate_counter)-1))
      fprintf (Protein_group_output3, '%s,',...
        num2str(Protein_group.Number_peptide_pre_fraction_fix(1,(Define_experiment_fractions_start(replicate_counter)+(counter-1)))));
    end
    
    %Jump to next line
    fprintf (Protein_group_output3,'\n'),...
      
  
  %Write peptide ratio
  number_of_peptides=length(Protein_group.Sequence);
  
  for counter=1:number_of_peptides
    fprintf (Protein_group_output3,[Column_header2, '\n'],...
      Protein_group.Sequence{counter,1}, Protein_group.Major_protein_group_name{1},...
      Protein_group.Proteins_assigned_in_peptide_table{counter,1},...
      Protein_group.modification{counter,1},num2str(Protein_group.Start_position{counter,1}),...
      num2str(Protein_group.End_position{counter,1}),Protein_group.Semi_tryptic{counter,1},...
      Protein_group.Semi_tryptic{counter,2},Protein_group.Semi_tryptic{counter,3},...
      Protein_group.Peptide_different_to_mean_fixed{counter,1},num2str(cluster_analysis_fix(1)),...
      Protein_group.mass{counter,1},Protein_group.calibrated_time{counter,1},...
      Protein_group.peptide_ratio_raw2(counter,(Define_experiment_fractions_start(replicate_counter):Define_experiment_fractions_end(replicate_counter))));
  end
  fclose(Protein_group_output3);
  
  end
  
  
  %move file to desired folder
  try
    movefile(Protein_group_name1{:}, 'ProteinGroups text_output');
    
  catch
    fclose(Protein_group_output1);
    movefile(Protein_group_name1{:}, 'ProteinGroups text_output');
  end
  
  try
    movefile(Protein_group_name2{:}, 'ProteinGroups text_output');
  catch
    fclose(Protein_group_output2);
    movefile(Protein_group_name2{:}, 'ProteinGroups text_output');
  end
  
  if isotope_channel==3
    try
      movefile(Protein_group_name3{:}, 'ProteinGroups text_output');
    catch
      fclose(Protein_group_output3);
      movefile(Protein_group_name3{:}, 'ProteinGroups text_output');
    end
  end
  
  fprintf(fid_summary, '%s,%s,%6.4f,%6.4f,%6.4f,%s,%s,%6.4f,%6.4f,\n',...
    Protein_group.Major_protein_group_name{1,1}, isotope_name_value{1}, replicate_counter,...
    cluster_analysis(1),cluster_analysis_fix(1),Protein_group.Global_Peptide_different_to_mean{1,1},...
    Protein_group.Global_Peptide_different_to_mean_fix{1,1},Protein_group.Shannon_wiener_evenest(1),...
    Protein_group.Shannon_wiener_diversity);
  
  
    end
    
  end
  
  
end
toc

fclose(fid_summary);
fclose all