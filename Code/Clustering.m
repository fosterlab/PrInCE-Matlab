%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Using Markov clusting to determine the number of unique complexes
%                          Created by Nichollas Scott,
%                              Foster lab,UBC, 2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%clear all
diary([user.maindir 'logfile.txt'])


%% 0. Initialize

%Define Colours to use
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Define fraction number
fraction_number=[2 user.Nfraction];

fn = [user.maindir 'Data/ROC/CombinedResults/'];
dd = dir([fn 'Final_Interactions_list_*_precision*']);
Interaction_file_names = cell(length(dd),1);
for di = 1:length(dd)
  Interaction_file_names{di} = [fn dd(di).name];
end
% Interaction_file_names= {'Final_Interactions_list_70_precision.csv',....
%   'Final_Interactions_list_60_precision.csv',....
%   'Final_Interactions_list_50_precision.csv'};

% Make directories
datadir = [maindir 'Data/Clustering/']; % where data files live
figdir1 = [maindir 'Figures/Clustering/']; % where figures live
%tmpdir = '/Users/Mercy/Academics/Foster/Jenny_PCPSILAC/PCPSILAC_Analysis/Data/Alignment/';
% Make folders if necessary
if ~exist(datadir, 'dir'); mkdir(datadir); end
if ~exist(datadir, 'dir'); mkdir(datadir); end

Interaction_pc_values=user.desiredPrecision;


%% Import master Gaussian list
filename = [user.maindir 'Data/Comparison/Master_guassian_list.csv'];
fileID = fopen(filename,'r');
dataArray = textscan(fileID, '%s%s%f%s%f%f%f%f%f%f%f%f%[^\n\r]', 'Delimiter', ',', 'EmptyValue' ,NaN,'HeaderLines' ,1, 'ReturnOnError', false);

%% Close the text file.
fclose(fileID);

Master_Gaussian.Proteinname = dataArray{:, 1};
Master_Gaussian.Unique_identifierofbestgaussian = dataArray{:, 2};
Master_Gaussian.Replicateofbestgaussian = dataArray{:, 3};
Master_Gaussian.Channelofbestgaussian = dataArray{:, 4};
Master_Gaussian.Guassian_index_numberofbestgaussian = dataArray{:, 5};
Master_Gaussian.Centerofbestgaussian = dataArray{:, 6};
Master_Gaussian.Heightofbestgaussian = dataArray{:, 7};
Master_Gaussian.Widthofbestgaussian = dataArray{:, 8};
Master_Gaussian.SSEofbestgaussian = dataArray{:, 9};
Master_Gaussian.adjrsquareofbestgaussian = dataArray{:, 10};
Master_Gaussian.GaussianAreaofbestgaussian = dataArray{:, 11};
Master_Gaussian.ComplexSizeofbestgaussian = dataArray{:, 12};

number_unique_guassian=length(Master_Gaussian.Proteinname);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% read in annotation data, note as this is so large need to use the headers
% to find out what columns you would like to import

%% Read in peptide results
Annotation_file=[user.maindir 'Data/mainAnnot.homo_sapiens.txt'];
inFile=fopen(Annotation_file);
headerLine=fgetl(inFile);

%Read in header
tempCellArray=textscan(headerLine,'%s','DELIMITER','\t','BUFSIZE',100000);

%look for Calibrated retention time column
%Convert header to upper case, this is to ensure matching for caps
fileCols=upper(tempCellArray{1,1});
columnNames = upper({'UNIPROT','GOMF NAME','GOBP NAME','CORUM'});

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
Annotation_column_content=textscan(inFile,format,'DELIMITER','\t','BUFSIZE',1000000);

%Unpack vaules
Uniprot_annotations=Annotation_column_content{1,1};
GOMF_annotations=Annotation_column_content{1,2};
GOBP_annotations=Annotation_column_content{1,3};
CORUM_annotations=Annotation_column_content{1,4};

%Remove Annotation_column_content to clear up memory
clearvars Annotation_column_content;

Uniprot_annotations_unpacked=regexp(Uniprot_annotations(2:end), ';','split');
for unpack_counter1=1:length(Uniprot_annotations_unpacked)
  %Reset variable
  Temp=cell(0,0);
  %Write out to temp variable
  Temp_var=Uniprot_annotations_unpacked{unpack_counter1};
  %determine length
  Lenght_temp=length(Temp_var);
  for unpack_counter2=1:Lenght_temp
    Uniprot_annotations_unpacked(unpack_counter1,unpack_counter2)=Temp_var(unpack_counter2);
  end
end

%Replace empty cell with NaN
Uniprot_annotations_unpacked(cellfun(@isempty, Uniprot_annotations_unpacked))={'NaN'};

%Create varible
Master_Gaussian.Corum_terms=cell(number_unique_guassian,1);
Master_Gaussian.GOBP_terms=cell(number_unique_guassian,1);
Master_Gaussian.GOMF_terms=cell(number_unique_guassian,1);

%Create smaller list based on all observed proteins
for protein_counter=1:number_unique_guassian
  %look up times an interaction observed
  [col row]=ind2sub(size(Uniprot_annotations_unpacked), strmatch(Master_Gaussian.Proteinname{protein_counter}, Uniprot_annotations_unpacked, 'exact'));
  
  %Copy 'GOMF NAME','GOBP NAME','CORUM'
  
  %Write out GOMF annotation
  if ~isempty(GOMF_annotations(col))
    Master_Gaussian.GOMF_terms(protein_counter)=GOMF_annotations(col);
  else
    Master_Gaussian.GOMF_terms(protein_counter)={''};
  end
  
  %Write out GOBP annotation
  if ~isempty(GOBP_annotations(col))
    Master_Gaussian.GOBP_terms(protein_counter)=GOBP_annotations(col);
  else
    Master_Gaussian.GOBP_terms(protein_counter)={''};
  end
  
  %Write out CORUM annotation
  if ~isempty(CORUM_annotations(col))
    Master_Gaussian.Corum_terms(protein_counter,1)=CORUM_annotations(col);
  else
    Master_Gaussian.Corum_terms(protein_counter,1)={''};
  end
  
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Unpack GO/CORUM terms
Master_Gaussian.Unpacked_GOMF_terms_unpacked=cell(0,0);
Master_Gaussian.Unpacked_GOBP_terms_unpacked=cell(0,0);
Master_Gaussian.Unpacked_CORUM_terms_unpacked=cell(0,0);

for unique_complex_counter=1:number_unique_guassian
  
  %Unpack GOMF
  Unpacked_terms1=regexp(Master_Gaussian.GOMF_terms(unique_complex_counter),';','split');
  Unpacked_terms=Unpacked_terms1{:};
  
  for temp_count=1:length(Unpacked_terms)
    Master_Gaussian.Unpacked_GOMF_terms_unpacked(unique_complex_counter,temp_count)=Unpacked_terms(temp_count);
  end
  
  %Unpack GOBP
  Unpacked_terms1=regexp(Master_Gaussian.GOBP_terms(unique_complex_counter),';','split');
  Unpacked_terms=Unpacked_terms1{:};
  
  for temp_count=1:length(Unpacked_terms)
    Master_Gaussian.Unpacked_GOBP_terms_unpacked(unique_complex_counter,temp_count)=Unpacked_terms(temp_count);
  end
  
  %Unpack CORUM
  Unpacked_terms1=regexp(Master_Gaussian.Corum_terms(unique_complex_counter),';','split');
  Unpacked_terms=Unpacked_terms1{:};
  
  for temp_count=1:length(Unpacked_terms)
    Master_Gaussian.Unpacked_CORUM_terms_unpacked(unique_complex_counter,temp_count)=Unpacked_terms(temp_count);
  end
end

%Replace empty cell with NaN
Master_Gaussian.Unpacked_GOMF_terms_unpacked(cellfun(@isempty,  Master_Gaussian.Unpacked_GOMF_terms_unpacked))={'NaN'};
Master_Gaussian.Unpacked_GOBP_terms_unpacked(cellfun(@isempty, Master_Gaussian.Unpacked_GOBP_terms_unpacked))={'NaN'};
Master_Gaussian.Unpacked_CORUM_terms_unpacked(cellfun(@isempty, Master_Gaussian.Unpacked_CORUM_terms_unpacked))={'NaN'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for Global_counter=1:length(Interaction_file_names)
  
  %% Import interaction list
  fileID_A = fopen(Interaction_file_names{Global_counter},'r');
  Import_Interaction_70pc_1 = textscan(fileID_A, '%s%s%s%s%s%s%s%s%s%s%s%s%[^\n\r]', 'Delimiter', ',', 'HeaderLines' ,1, 'ReturnOnError', false);
  fclose(fileID_A);
  
  %Count number of interactions
  detected_interaction1=length(Import_Interaction_70pc_1{1,1});
  
  %Convert for cyto interaction to a list that can be used
  interaction_A.unique_interactions=Import_Interaction_70pc_1{1,1};
  interaction_A.proteinA =Import_Interaction_70pc_1{1,2};
  interaction_A.proteinB =Import_Interaction_70pc_1{1,3};
  interaction_A.centerA_formated_cyto=Import_Interaction_70pc_1{1,4};
  interaction_A.centerB_formated_cyto=Import_Interaction_70pc_1{1,5};
  interaction_A.Replicates_formated_cyto=Import_Interaction_70pc_1{1,6};
  interaction_A.DeltaHeight_formated_cyto=Import_Interaction_70pc_1{1,7};
  interaction_A.DeltaCenter_formated_cyto=Import_Interaction_70pc_1{1,8};
  interaction_A.Deltawidth_formated_cyto=Import_Interaction_70pc_1{1,9};
  interaction_A.DeltaEuc_formated_cyto=Import_Interaction_70pc_1{1,11};
  interaction_A.proteinInCorum=Import_Interaction_70pc_1{1,11};
  interaction_A.interactionInCorum=Import_Interaction_70pc_1{1,12};
  
  
  %Convert Centers into a numerical value to be used
  for num_convert=1:detected_interaction1
    Temp1=regexp(interaction_A.centerA_formated_cyto{num_convert},';','split');
    Temp2=regexp(interaction_A.centerB_formated_cyto{num_convert},';','split');
    
    for temp_counter=1:length(Temp1)
      Temp1(1,temp_counter)=Temp1(temp_counter);
    end
    
    for temp_counter=1:length(Temp2)
      Temp2(1,temp_counter)=Temp2(temp_counter);
    end
    
    %Copy first entry to use downstream
    interaction_A.centerA_first_value(num_convert,1)=str2num(Temp1{1,1});
    interaction_A.centerB_first_value(num_convert,1)=str2num(Temp2{1,1});
  end
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  %Create matrix to write data to
  Interaction_matrix=zeros(number_unique_guassian,number_unique_guassian);
  
  Can_not_find_counter=1;
  
  for unique_Gaus_counter=1:number_unique_guassian
    
    %Test is center of protein is between fraction 2 and 50
    
    if (Master_Gaussian.Centerofbestgaussian(unique_Gaus_counter)<50 && Master_Gaussian.Centerofbestgaussian(unique_Gaus_counter)>2)
      %look up protein in interation list
      if ~isempty(strmatch(Master_Gaussian.Proteinname(unique_Gaus_counter),interaction_A.proteinA, 'exact'))==1 |...
          ~isempty(strmatch(Master_Gaussian.Proteinname(unique_Gaus_counter),interaction_A.proteinB, 'exact'))==1
        
        %reset values
        Protein_pair_to_test=cell(0,0);
        test_values=zeros(1,1);
        test_values2=zeros(1,1);
        
        %look up times an interaction observed
        [position_in_proteinA]=ind2sub(size(interaction_A.proteinA), strmatch(Master_Gaussian.Proteinname(unique_Gaus_counter), interaction_A.proteinA, 'exact'));
        
        [position_in_proteinB]=ind2sub(size(interaction_A.proteinB), strmatch(Master_Gaussian.Proteinname(unique_Gaus_counter), interaction_A.proteinB, 'exact'));
        
        %add position together
        total_position=[position_in_proteinA' position_in_proteinB'];
        
        Center_to_test=[interaction_A.centerA_first_value(position_in_proteinA)' interaction_A.centerB_first_value(position_in_proteinB)'];
        Protein_pair_to_test=[interaction_A.proteinB(position_in_proteinA)' interaction_A.proteinA(position_in_proteinB)'];
        
        %test if any of the values are within two fractions
        test_values=(Center_to_test-Master_Gaussian.Centerofbestgaussian(unique_Gaus_counter));
        %remove values greater then two away
        total_position(abs(test_values)>2)=[];
        
        if ~isempty(total_position)==1
          %Determine the connectivity of the protein
          number_unique_interactions=length(total_position);
          
          for density_counter=1:number_unique_interactions
            %find position of protein pairs
            [position_of_interactor]=ind2sub(size(Master_Gaussian.Proteinname), strmatch(Protein_pair_to_test{density_counter}, Master_Gaussian.Proteinname, 'exact'));
            
            %If can not find proteins
            if isempty(position_of_interactor)==1
              
              Can_not_find_proteins(Can_not_find_counter)=Protein_pair_to_test(density_counter);
              Can_not_find_counter=Can_not_find_counter+1;
              
            else
              %Determine if interaction is within two fract
              test_values2=(Master_Gaussian.Centerofbestgaussian(unique_Gaus_counter)-Master_Gaussian.Centerofbestgaussian(position_of_interactor));
              position_of_interactor(abs(test_values2)>2)=[];
              
              if ~isempty(position_of_interactor)==1
                %Record data in matrix
                Interaction_matrix(unique_Gaus_counter,position_of_interactor)=1;
              end
              
            end
          end
        end
      end
    end
  end
  
  %Normalize interaction matrix
  total_unique_proteins_interacting_with=sum(Interaction_matrix');
  
  Normalized_Interaction_matrix=zeros(number_unique_guassian,number_unique_guassian);
  %divide
  for unique_Gaus_counter=1:number_unique_guassian
    if total_unique_proteins_interacting_with(unique_Gaus_counter)>0
      Normalized_Interaction_matrix(:,unique_Gaus_counter)=Interaction_matrix(unique_Gaus_counter,:)./total_unique_proteins_interacting_with(unique_Gaus_counter);
    end
  end
  
  %Save Normalized Matrix for MCL analysis
  Output_name=strcat(datadir,'Normalized_interaction_matrix_pc',num2str(Interaction_pc_values(Global_counter)),'.csv');
  csvwrite(Output_name,Normalized_Interaction_matrix);
  
  %attempt mcl
  mTemp=Normalized_Interaction_matrix;
  
  %Set varibales
  energy = 1;
  emax = 0.001;
  
  p = 2;
  minval = 0.000001;
  
  %set
  iteration=1;
  %Repeat till the variation is less then 10% OR 10 interations have been
  %evaluated, Note the limit on interations is due convergence always being
  %outcome from MCL
  while (energy > emax*1.10 | energy < emax*0.90) | iteration== 10
    
    %Set emax
    emax =energy;
    
    %write out values
    iteration
    energy
    
    % expand by multiplying m * m
    % this preserves column (or row) normalisation
    m2 = mTemp * mTemp;
    
    % inflation
    m2 = mTemp .^ p;
    
    % pruning
    m2(find(m2 < minval)) = 0;
    
    % normalisation
    dinv = diag(1./sum(m2));
    m2 = m2 * dinv;
    
    %Remove NaN with zeros
    m2(isnan(m2))=0;
    
    %Test results, uncomment if need to print out matrix
    %Output_name=strcat('MCL_Normalized_interaction_matrix_iteration',num2str(iteration),'.csv');
    %csvwrite(Output_name,m2);
    
    % calculate residual energy
    maxs = max(m2);
    sqsums = sum(m2 .^ 2);
    energy = max(maxs - sqsums);
    
    %set m2 to mTemp
    mTemp=m2;
    
    iteration=iteration+1;
    
  end % while e
  
  %Test results, uncomment if need to print out matrix
  Output_name=strcat(datadir, 'MCL_Normalized_interaction_matrix_iteration_pc',num2str(Interaction_pc_values(Global_counter)),'_',num2str(iteration-1),'.csv');
  csvwrite(Output_name,mTemp);
  
  %Determine how many unique combinations
  
  String_combinations=cell(1,1);
  Names_proteins_in_complexes=cell(1,1);
  Center_proteins_in_complexes=zeros(1,1);
  
  Counter_1=1;
  for unique_Gaus_counter=1:number_unique_guassian
    %check if atleast one measuerment is within row
    if ~sum(mTemp(unique_Gaus_counter,:))==0
      %Turn row into a string and find unique
      String_combinations{Counter_1}=num2str(mTemp(unique_Gaus_counter,:));
      %Save names
      Names_proteins_in_complexes(Counter_1)=Master_Gaussian.Proteinname(unique_Gaus_counter);
      %Save Center
      Center_proteins_in_complexes(Counter_1)=Master_Gaussian.Centerofbestgaussian(unique_Gaus_counter);
      %Add one to counter
      Counter_1=Counter_1+1;
    end
  end
  
  %Find unique combinations
  [Unique_combinations Unique_positions]=unique(String_combinations);
  Number_unique_complex=length(Unique_combinations);
  
  %Convert back to number
  for unique_complex_counter=1:Number_unique_complex
    Number_temp(unique_complex_counter,:)= str2num(Unique_combinations{unique_complex_counter});
  end
  
  %Write out names of proteins in complexes
  NR_Names_proteins_in_complexes=Names_proteins_in_complexes(Unique_positions);
  
  %Write out names of proteins in complexes
  NR_Center_proteins_in_complexes=Center_proteins_in_complexes(Unique_positions);
  
  %Create array to save values to
  Number_interaction_parters=zeros(Number_unique_complex,1);
  
  Output_name2=strcat(datadir,'Unique_Cluster_identified_pc',num2str(Interaction_pc_values(Global_counter)),'.txt');
  
  fid = fopen(Output_name2,'w');
  fprintf (fid,'%s\t%s\t%s\t%s\t%s\t%s\t\n',...
    'Protein Complex number','Number of interactors in complex','Center of central node','Protein identity of Central node','Centers of protein interactors','Protein interactors'); %Write Header
  %Build interaction network
  for unique_complex_counter=1:Number_unique_complex
    
    %Write out cluster number
    Cluster_name=strcat(num2str(unique_complex_counter));
    fprintf(fid,'%s	',Cluster_name);
    
    %Find non zero position
    position_of_interactor=find(Number_temp(unique_complex_counter,:));
    
    %Determine number of proteins to write out
    number_to_write_out=length(position_of_interactor);
    
    %Add one to interaction
    number_to_write_out=number_to_write_out+1;
    
    %Save number of interactions
    Number_interaction_parters(unique_complex_counter,1)=number_to_write_out;
    
    %Write out number of interactions within
    fprintf(fid,'%f	',(number_to_write_out));
    
    %Write out Center
    fprintf(fid,'%s	',num2str(NR_Center_proteins_in_complexes(unique_complex_counter)));
    
    %Write out Centers of interactors
    fprintf(fid,'%s	',strjoin(Master_Gaussian.Centerofbestgaussian(position_of_interactor),';'));
    
    %Write out Protein A
    fprintf(fid,'%s	',NR_Names_proteins_in_complexes{unique_complex_counter});
    
    
    for write_out=1:(number_to_write_out-1)
      fprintf(fid,'%s	',Master_Gaussian.Proteinname{position_of_interactor(write_out)});
    end
    
    fprintf(fid,'\n');
  end
  fclose(fid);
  
  %Write out stats file for clutstering
  mean_number_interactions=mean(Number_interaction_parters);
  median_number_interactions=median(Number_interaction_parters);
  Max_number_interactions=max(Number_interaction_parters);
  Min_number_interactions=min(Number_interaction_parters);
  STdev_interactions=std(Number_interaction_parters);
  
  Output_name3=strcat(datadir,'Summary_clustering_information_pc',num2str(Interaction_pc_values(Global_counter)),'.txt');
  
  fid = fopen(Output_name3,'w');
  fprintf (fid,'%s\t%s\t%s\t%s\t%s\t%s\t\n',...
    'Total number of cluster','Mean number of interactions','Median Number of interactions','Maximum number of interactions','Minimum Number of interactions','Standard derviation of interactions'); %Write Header
  fprintf(fid,'%f	',Number_unique_complex);
  fprintf(fid,'%f	',mean_number_interactions);
  fprintf(fid,'%f	',median_number_interactions);
  fprintf(fid,'%f	',Max_number_interactions);
  fprintf(fid,'%f	',Min_number_interactions);
  fprintf(fid,'%f	',STdev_interactions);
  fclose(fid);
  
  %Create figures
  
  %number of complexes over SEC
  Hist_array=zeros(fraction_number(2),7);
  
  for fraction_counter1= 1:Number_unique_complex
    %Define Center of Complex
    Center_of_gaussain= floor(NR_Center_proteins_in_complexes(fraction_counter1));
    %Ensure center is greater then zero
    if  Center_of_gaussain>0
      
      if Number_interaction_parters(fraction_counter1) >=1 && Number_interaction_parters(fraction_counter1) <=5 && Center_of_gaussain <=fraction_number(2)
        
        Hist_array(Center_of_gaussain,1)=Hist_array(Center_of_gaussain,1)+1;
        
      elseif Number_interaction_parters(fraction_counter1) >5 && Number_interaction_parters(fraction_counter1) <=10 && Center_of_gaussain <=fraction_number(2)
        
        Hist_array(Center_of_gaussain,2)=Hist_array(Center_of_gaussain,2)+1;
        
      elseif Number_interaction_parters(fraction_counter1) >10 && Number_interaction_parters(fraction_counter1) <=20 && Center_of_gaussain <=fraction_number(2)
        
        Hist_array(Center_of_gaussain,3)=Hist_array(Center_of_gaussain,3)+1;
      elseif Number_interaction_parters(fraction_counter1) >20 && Number_interaction_parters(fraction_counter1) <=30 && Center_of_gaussain <=fraction_number(2)
        
        Hist_array(Center_of_gaussain,4)=Hist_array(Center_of_gaussain,4)+1;
      elseif Number_interaction_parters(fraction_counter1) >30 && Number_interaction_parters(fraction_counter1) <=40 && Center_of_gaussain <=fraction_number(2)
        
        Hist_array(Center_of_gaussain,5)=Hist_array(Center_of_gaussain,5)+1;
      elseif Number_interaction_parters(fraction_counter1) >40 && Number_interaction_parters(fraction_counter1) <=50 && Center_of_gaussain <=fraction_number(2)
        
        Hist_array(Center_of_gaussain,6)=Hist_array(Center_of_gaussain,6)+1;
      elseif Number_interaction_parters(fraction_counter1) >50 && Center_of_gaussain <=fraction_number(2)
        
        Hist_array(Center_of_gaussain,7)=Hist_array(Center_of_gaussain,7)+1;
      end
    end
  end
  
  
  %Graph data as log2 scatter
  f1=figure;
  f6_figure=bar(1:fraction_number(2), [Hist_array(:,1) Hist_array(:,2) Hist_array(:,3) Hist_array(:,4) Hist_array(:,5) Hist_array(:,6) Hist_array(:,7)],0.6, 'stack');
  for k=1:7
    set(f6_figure(k),'facecolor', colour_to_use(k,:), 'EdgeColor', 'k' )
  end
  legend('1-5 protein members', '6-10 protein members', '11-20 protein members', '21-30 protein members','31-40 protein members', '41-50 protein members','>50 members','FontSize',8, 'Location', 'Best');
  xlim([0,fraction_number(2)+2]);
  title('Observed Complexes Across the Cytoplasm interactome');
  xlabel('Fractions');
  ylabel('# Complex clusters');
  Output_name4=strcat(figdir1,'Detected_complex_across_SEC_pc',num2str(Interaction_pc_values(Global_counter)),'.pdf');
  print(f1,'-dpdf','-r600', Output_name4)
  
  
  %Pie chart of complexes
  f2=figure;
  h=pie(sum(Hist_array))
  hp = findobj(h, 'Type', 'patch');
  set(hp(1), 'FaceColor', myC(2,:));
  set(hp(2), 'FaceColor', colour_to_use(2,:));
  set(hp(3), 'FaceColor', colour_to_use(4,:));
  set(hp(4), 'FaceColor', colour_to_use(7,:));
  set(hp(3), 'FaceColor', colour_to_use(6,:));
  set(hp(3), 'FaceColor', colour_to_use(5,:));
  legend('1-5 protein members', '6-10 protein members', '11-20 protein members', '21-30 protein members','31-40 protein members', '41-50 protein members','>50 members','FontSize',6, 'Location', 'southoutside')
  title('Number of protein members in Complexes cluster');
  Output_name5=strcat(figdir1,'Pie_chart_of_complexes_pc',num2str(Interaction_pc_values(Global_counter)),'.pdf');
  print(f2,'-dpdf','-r600', Output_name5);
  
  
  %Create arrya to write GO terms to
  Combined_GOMF_vaules=cell(Number_unique_complex,1);
  Combined_GOMF_terms=cell(Number_unique_complex,1);
  
  Max_terms=cell(Number_unique_complex,1);
  Max_values_array=zeros(Number_unique_complex,1);
  Max_values_pc_array=zeros(Number_unique_complex,1);
  
  Output_name2=strcat(datadir,'Unique_Cluster_identified_GOMF_terms_pc',num2str(Interaction_pc_values(Global_counter)),'.txt');
  
  fid = fopen(Output_name2,'w');
  fprintf (fid,'%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t\n',...
    'Protein Complex number','Number of interactors in complex','most common observed GOMF term in complex','Number of time common observed GOMF term in complex',...
    '% common observed GOMF term in complex','Observed GOMF terms','Number of GOMF terms observed',...
    'Center of central node','Protein identity of Central node','Centers of protein interactors','Protein interactors'); %Write Header
  %Build interaction network
  
  
  %Determine the most common observed GO Terms from complexes
  for unique_complex_counter=1:Number_unique_complex
    
    %Write out cluster number
    Cluster_name=strcat(num2str(unique_complex_counter));
    fprintf(fid,'%s	',Cluster_name);
    
    %Find non zero position
    position_of_interactor=find(Number_temp(unique_complex_counter,:));
    
    %Determine number of proteins to write out
    number_to_write_out=length(position_of_interactor);
    
    %Find Protein A in master list
    [Position_of_proteinA]=ind2sub(size(Master_Gaussian.Proteinname), strmatch(NR_Names_proteins_in_complexes{unique_complex_counter}, Master_Gaussian.Proteinname, 'exact'));
    
    
    %Copy Protein A terms to array
    Observed_terms=Master_Gaussian.Unpacked_GOMF_terms_unpacked(Position_of_proteinA(1),:);
    
    %Copy Protein B terms to array
    for write_out=1:number_to_write_out
      Observed_terms=[Observed_terms Master_Gaussian.Unpacked_GOMF_terms_unpacked(position_of_interactor(write_out),:)];
    end
    
    %Replace empty cells with NaN
    Observed_terms(cellfun(@isempty, Observed_terms))={'NaN'};
    
    %Determine unique terms
    Unique_terms_observed=unique(Observed_terms);
    
    %Remove NaN
    Unique_terms_observed(strmatch('NaN',Unique_terms_observed))=[];
    
    %Determine the lenght of array
    number_of_unique_observation=length(Unique_terms_observed);
    
    %Create array to count to
    Counter_array=zeros(1,number_of_unique_observation);
    
    %Count observation of Unique terms in Protein A
    for Term_counter=1:number_of_unique_observation
      for Term_counter2=1:length(Master_Gaussian.Unpacked_GOMF_terms_unpacked(Position_of_proteinA(1),:))
        if strcmp(Master_Gaussian.Unpacked_GOMF_terms_unpacked{Position_of_proteinA(1),Term_counter2},Unique_terms_observed{Term_counter})==1
          Counter_array(Term_counter)=  Counter_array(Term_counter)+1;
        end
      end
    end
    
    %Count observation of Unique terms in Protein B
    for write_out=1:number_to_write_out
      for Term_counter=1:number_of_unique_observation
        for Term_counter2=1:length(Master_Gaussian.Unpacked_GOMF_terms_unpacked(Position_of_proteinA(1),:))
          if strcmp(Master_Gaussian.Unpacked_GOMF_terms_unpacked{position_of_interactor(write_out),Term_counter2},Unique_terms_observed{Term_counter})==1
            Counter_array(Term_counter)=  Counter_array(Term_counter)+1;
          end
        end
      end
    end
    
    %Copy Count array to a string
    Combined_GOMF_vaules{unique_complex_counter}=strjoin(Counter_array,';');
    Combined_GOMF_terms{unique_complex_counter}=strjoin(Unique_terms_observed,';');
    
    if ~isempty(Counter_array)
      %Determine the GO term with the highest number of observation
      [value_of_max location_of_max]=max(Counter_array);
      
      %Save max values
      Max_terms{unique_complex_counter}=Unique_terms_observed(location_of_max);
      Max_values_array(unique_complex_counter)=value_of_max;
      Max_values_pc_array(unique_complex_counter)=value_of_max/(number_to_write_out+1);
    else
      Max_terms{unique_complex_counter}=NaN;
      Max_values_array(unique_complex_counter)=NaN;
      Max_values_pc_array(unique_complex_counter)=NaN;
    end
    
    
    %Write out number of interactions within
    fprintf(fid,'%f	',(number_to_write_out+1));
    
    %Write Max GOMF terms
    if ~isempty(Counter_array)==1
      fprintf(fid,'%s	', Unique_terms_observed{location_of_max});
      
      %Write out max number of times GOMF term observed
      fprintf(fid,'%f	',Max_values_array(unique_complex_counter));
      
      %Write out number of interactions within
      fprintf(fid,'%f	',Max_values_pc_array(unique_complex_counter));
    else
      fprintf(fid,'%s	', NaN);
      
      %Write out max number of times GOMF term observed
      fprintf(fid,'%f	',NaN);
      
      %Write out number of interactions within
      fprintf(fid,'%f	',NaN);
    end
    
    %Write List of GOMF terms
    fprintf(fid,'%s	',Combined_GOMF_terms{unique_complex_counter});
    
    %Write List of GOMF number
    fprintf(fid,'%s	',Combined_GOMF_vaules{unique_complex_counter});
    
    %Write out Center
    fprintf(fid,'%s	',num2str(NR_Center_proteins_in_complexes(unique_complex_counter)));
    
    %Write out Centers of interactors
    fprintf(fid,'%s	',strjoin(Master_Gaussian.Centerofbestgaussian(position_of_interactor),';'));
    
    %Write out Protein A
    fprintf(fid,'%s	',NR_Names_proteins_in_complexes{unique_complex_counter});
    
    %Write out Protein B
    for write_out=1:(number_to_write_out)
      fprintf(fid,'%s	',Master_Gaussian.Proteinname{position_of_interactor(write_out)});
    end
    
    fprintf(fid,'\n');
    
    
    
  end
  fclose(fid);
  
  
  
  
  Output_name2=strcat(datadir,'Unique_Cluster_identified_GOBP_terms_pc',num2str(Interaction_pc_values(Global_counter)),'.txt');
  
  fid = fopen(Output_name2,'w');
  fprintf (fid,'%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t\n',...
    'Protein Complex number','Number of interactors in complex','most common observed GOBP term in complex','Number of time common observed GOBP term in complex',...
    '% common observed GOBP term in complex','Observed GOBP terms','Number of GOBP terms observed',...
    'Center of central node','Protein identity of Central node','Centers of protein interactors','Protein interactors'); %Write Header
  %Build interaction network
  
  
  %Determine the most common observed GO Terms from complexes
  for unique_complex_counter=1:Number_unique_complex
    
    %Write out cluster number
    Cluster_name=strcat(num2str(unique_complex_counter));
    fprintf(fid,'%s	',Cluster_name);
    
    %Find non zero position
    position_of_interactor=find(Number_temp(unique_complex_counter,:));
    
    %Determine number of proteins to write out
    number_to_write_out=length(position_of_interactor);
    
    %Find Protein A in master list
    [Position_of_proteinA]=ind2sub(size(Master_Gaussian.Proteinname), strmatch(NR_Names_proteins_in_complexes{unique_complex_counter}, Master_Gaussian.Proteinname, 'exact'));
    
    
    %Copy Protein A terms to array
    Observed_terms=Master_Gaussian.Unpacked_GOBP_terms_unpacked(Position_of_proteinA(1),:);
    
    %Copy Protein B terms to array
    for write_out=1:number_to_write_out
      Observed_terms=[Observed_terms Master_Gaussian.Unpacked_GOBP_terms_unpacked(position_of_interactor(write_out),:)];
    end
    
    %Replace empty cells with NaN
    Observed_terms(cellfun(@isempty, Observed_terms))={'NaN'};
    
    %Determine unique terms
    Unique_terms_observed=unique(Observed_terms);
    
    %Remove NaN
    Unique_terms_observed(strmatch('NaN',Unique_terms_observed))=[];
    
    %Determine the lenght of array
    number_of_unique_observation=length(Unique_terms_observed);
    
    %Create array to count to
    Counter_array=zeros(1,number_of_unique_observation);
    
    %Count observation of Unique terms in Protein A
    for Term_counter=1:number_of_unique_observation
      for Term_counter2=1:length(Master_Gaussian.Unpacked_GOBP_terms_unpacked(Position_of_proteinA(1),:))
        if strcmp(Master_Gaussian.Unpacked_GOBP_terms_unpacked{Position_of_proteinA(1),Term_counter2},Unique_terms_observed{Term_counter})==1
          Counter_array(Term_counter)=  Counter_array(Term_counter)+1;
        end
      end
    end
    
    %Count observation of Unique terms in Protein B
    for write_out=1:number_to_write_out
      for Term_counter=1:number_of_unique_observation
        for Term_counter2=1:length(Master_Gaussian.Unpacked_GOBP_terms_unpacked(Position_of_proteinA(1),:))
          if strcmp(Master_Gaussian.Unpacked_GOBP_terms_unpacked{position_of_interactor(write_out),Term_counter2},Unique_terms_observed{Term_counter})==1
            Counter_array(Term_counter)=  Counter_array(Term_counter)+1;
          end
        end
      end
    end
    
    %Copy Count array to a string
    Combined_GOBP_vaules{unique_complex_counter}=strjoin(Counter_array,';');
    Combined_GOBP_terms{unique_complex_counter}=strjoin(Unique_terms_observed,';');
    
    if ~isempty(Counter_array)
      %Determine the GO term with the highest number of observation
      [value_of_max location_of_max]=max(Counter_array);
      
      %Save max values
      Max_terms_BP{unique_complex_counter}=Unique_terms_observed(location_of_max);
      Max_values_array_BP(unique_complex_counter)=value_of_max;
      Max_values_pc_array_BP(unique_complex_counter)=value_of_max/(number_to_write_out+1);
    else
      Max_terms_BP{unique_complex_counter}=NaN;
      Max_values_array_BP(unique_complex_counter)=NaN;
      Max_values_pc_array_BP(unique_complex_counter)=NaN;
    end
    
    
    %Write out number of interactions within
    fprintf(fid,'%f	',(number_to_write_out+1));
    
    %Write Max GOMF terms
    if ~isempty(Counter_array)==1
      fprintf(fid,'%s	', Unique_terms_observed{location_of_max});
      
      %Write out max number of times GOMF term observed
      fprintf(fid,'%f	',Max_values_array_BP(unique_complex_counter));
      
      %Write out number of interactions within
      fprintf(fid,'%f	',Max_values_pc_array_BP(unique_complex_counter));
    else
      fprintf(fid,'%s	', NaN);
      
      %Write out max number of times GOMF term observed
      fprintf(fid,'%f	',NaN);
      
      %Write out number of interactions within
      fprintf(fid,'%f	',NaN);
    end
    
    %Write List of GOMF terms
    fprintf(fid,'%s	',Combined_GOBP_terms{unique_complex_counter});
    
    %Write List of GOMF number
    fprintf(fid,'%s	',Combined_GOBP_vaules{unique_complex_counter});
    
    %Write out Center
    fprintf(fid,'%s	',num2str(NR_Center_proteins_in_complexes(unique_complex_counter)));
    
    %Write out Centers of interactors
    fprintf(fid,'%s	',strjoin(Master_Gaussian.Centerofbestgaussian(position_of_interactor),';'));
    
    %Write out Protein A
    fprintf(fid,'%s	',NR_Names_proteins_in_complexes{unique_complex_counter});
    
    %Write out Protein B
    for write_out=1:(number_to_write_out)
      fprintf(fid,'%s	',Master_Gaussian.Proteinname{position_of_interactor(write_out)});
    end
    
    fprintf(fid,'\n');
    
    
    
  end
  fclose(fid);
  
  
  Output_name2=strcat(datadir,'Unique_Cluster_identified_CORUM_terms_pc',num2str(Interaction_pc_values(Global_counter)),'.txt');
  
  fid = fopen(Output_name2,'w');
  fprintf (fid,'%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t\n',...
    'Protein Complex number','Number of interactors in complex','most common observed CORUM term in complex','Number of time common observed CORUM term in complex',...
    '% common observed CORUM term in complex','Observed CORUM terms','Number of CORUM terms observed',...
    'Center of central node','Protein identity of Central node','Centers of protein interactors','Protein interactors'); %Write Header
  %Build interaction network
  
  
  %Determine the most common observed GO Terms from complexes
  for unique_complex_counter=1:Number_unique_complex
    
    %Write out cluster number
    Cluster_name=strcat(num2str(unique_complex_counter));
    fprintf(fid,'%s	',Cluster_name);
    
    %Find non zero position
    position_of_interactor=find(Number_temp(unique_complex_counter,:));
    
    %Determine number of proteins to write out
    number_to_write_out=length(position_of_interactor);
    
    %Find Protein A in master list
    [Position_of_proteinA]=ind2sub(size(Master_Gaussian.Proteinname), strmatch(NR_Names_proteins_in_complexes{unique_complex_counter}, Master_Gaussian.Proteinname, 'exact'));
    
    
    %Copy Protein A terms to array
    Observed_terms=Master_Gaussian.Unpacked_CORUM_terms_unpacked(Position_of_proteinA(1),:);
    
    %Copy Protein B terms to array
    for write_out=1:number_to_write_out
      Observed_terms=[Observed_terms Master_Gaussian.Unpacked_CORUM_terms_unpacked(position_of_interactor(write_out),:)];
    end
    
    %Replace empty cells with NaN
    Observed_terms(cellfun(@isempty, Observed_terms))={'NaN'};
    
    %Determine unique terms
    Unique_terms_observed=unique(Observed_terms);
    
    %Remove NaN
    Unique_terms_observed(strmatch('NaN',Unique_terms_observed))=[];
    
    %Determine the lenght of array
    number_of_unique_observation=length(Unique_terms_observed);
    
    %Create array to count to
    Counter_array=zeros(1,number_of_unique_observation);
    
    %Count observation of Unique terms in Protein A
    for Term_counter=1:number_of_unique_observation
      for Term_counter2=1:length(Master_Gaussian.Unpacked_CORUM_terms_unpacked(Position_of_proteinA(1),:))
        if strcmp(Master_Gaussian.Unpacked_CORUM_terms_unpacked{Position_of_proteinA(1),Term_counter2},Unique_terms_observed{Term_counter})==1
          Counter_array(Term_counter)=  Counter_array(Term_counter)+1;
        end
      end
    end
    
    %Count observation of Unique terms in Protein B
    for write_out=1:number_to_write_out
      for Term_counter=1:number_of_unique_observation
        for Term_counter2=1:length(Master_Gaussian.Unpacked_CORUM_terms_unpacked(Position_of_proteinA(1),:))
          if strcmp(Master_Gaussian.Unpacked_CORUM_terms_unpacked{position_of_interactor(write_out),Term_counter2},Unique_terms_observed{Term_counter})==1
            Counter_array(Term_counter)=  Counter_array(Term_counter)+1;
          end
        end
      end
    end
    
    %Copy Count array to a string
    Combined_CORUM_vaules{unique_complex_counter}=strjoin(Counter_array,';');
    Combined_CORUM_terms{unique_complex_counter}=strjoin(Unique_terms_observed,';');
    
    if ~isempty(Counter_array)
      %Determine the GO term with the highest number of observation
      [value_of_max location_of_max]=max(Counter_array);
      
      %Save max values
      Max_terms_CORUM{unique_complex_counter}=Unique_terms_observed(location_of_max);
      Max_values_array_CORUM(unique_complex_counter)=value_of_max;
      Max_values_pc_array_CORUM(unique_complex_counter)=value_of_max/(number_to_write_out+1);
    else
      Max_terms_CORUM{unique_complex_counter}=NaN;
      Max_values_array_CORUM(unique_complex_counter)=NaN;
      Max_values_pc_array_CORUM(unique_complex_counter)=NaN;
    end
    
    
    %Write out number of interactions within
    fprintf(fid,'%f	',(number_to_write_out+1));
    
    %Write Max GOMF terms
    if ~isempty(Counter_array)==1
      fprintf(fid,'%s	', Unique_terms_observed{location_of_max});
      
      %Write out max number of times GOMF term observed
      fprintf(fid,'%f	',Max_values_array_CORUM(unique_complex_counter));
      
      %Write out number of interactions within
      fprintf(fid,'%f	',Max_values_pc_array_CORUM(unique_complex_counter));
    else
      fprintf(fid,'%s	', NaN);
      
      %Write out max number of times GOMF term observed
      fprintf(fid,'%f	',NaN);
      
      %Write out number of interactions within
      fprintf(fid,'%f	',NaN);
    end
    
    %Write List of GOMF terms
    fprintf(fid,'%s	',Combined_CORUM_terms{unique_complex_counter});
    
    %Write List of GOMF number
    fprintf(fid,'%s	',Combined_CORUM_vaules{unique_complex_counter});
    
    %Write out Center
    fprintf(fid,'%s	',num2str(NR_Center_proteins_in_complexes(unique_complex_counter)));
    
    %Write out Centers of interactors
    fprintf(fid,'%s	',strjoin(Master_Gaussian.Centerofbestgaussian(position_of_interactor),';'));
    
    %Write out Protein A
    fprintf(fid,'%s	',NR_Names_proteins_in_complexes{unique_complex_counter});
    
    %Write out Protein B
    for write_out=1:(number_to_write_out)
      fprintf(fid,'%s	',Master_Gaussian.Proteinname{position_of_interactor(write_out)});
    end
    
    fprintf(fid,'\n');
    
    
    
  end
  fclose(fid);
  
  
  close all
end

diary('off')



