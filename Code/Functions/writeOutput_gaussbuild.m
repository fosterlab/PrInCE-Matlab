%function writeOutput_gaussbuild(datadir,datadir1,datadir2,datadir3,datadir4,...
%  Coef,SSE,adjrsquare,Try_Fit,...
%  txt_MvsL,txt_HvsL,replicate,SEC_size_alignment,experimental_channels,...
%  cleandata,rawdata,tmp1,tmp2,...
%  protgausI)

% Handles all the output writing of Gauss_Build.m
%
% Makes these files:
%   1_1_HvsL_OutputGaus.csv
%   1_1_HvsL_Output_Chromatograms.csv                   % One
%   1_1_HvsL_OutputGaus_filtered_out.csv                % of
%   1_1_HvsL_Output_Chromatograms_filtered_out.csv      % these
%   1_1_MvsL_OutputGaus.csv                             % for
%   1_1_MvsL_Output_Chromatograms.csv                   % each
%   1_1_MvsL_OutputGaus_filtered_out.csv                % protein
%   1_1_MvsL_Output_Chromatograms_filtered_out.csv      % and Gaussian.
%   HvsL_Summary_Gausians_identifed.csv
%   MvsL_Summary_Gausians_identifed.csv
%   HvsL_Combined_OutputGaus.csv
%   HvsL_Summary_Gausians_for_individual_proteins.csv
%   MvsL_Combined_OutputGaus.csv
%   MvsL_Summary_Gausians_for_individual_proteins.csv
%   HvsL_Combined_OutputGaus_rep*.csv
%   MvsL_Combined_OutputGaus_rep*.csv
%
%   Combined_Chromatograms.csv
%   Combined_Chromatograms_filtered_out.csv
%   Combined_OutputGaus_filtered_out.csv
%   Summary_Proteins_with_Gausians.csv
%   Proteins_not_fitted_to_gaussian_HvsL.csv
%   Proteins_not_fitted_to_gaussian_MvsL.csv
%
% Adapted from Nichollas Scott's Gaus_build_24_1.m.
% Made by Greg Stacey on Nov 25 2015.


%% 0. Initialize

SEC_fit=polyfit(SEC_size_alignment(1,:),SEC_size_alignment(2,:),1);

Nchannels = size(Coef,1);
Nproteins = size(Coef,2);

% Make No_Gaus, i.e. the number of Gaussians selected for each protein
No_Gaus = zeros(1,Nproteins);
for ri = 1:Nproteins
  No_Gaus(ri) = length(Coef{1,ri})/3;
end


%% Make files for individual fits
% 1_1_HvsL_OutputGaus.csv
% 1_1_HvsL_Output_Chromatograms.csv
% 1_1_HvsL_OutputGaus_filtered_out.csv
% 1_1_HvsL_Output_Chromatograms_filtered_out.csv
% 1_1_MvsL_OutputGaus.csv
% 1_1_MvsL_Output_Chromatograms.csv
% 1_1_MvsL_OutputGaus_filtered_out.csv
% 1_1_MvsL_Output_Chromatograms_filtered_out.csv
disp('        Writing individual fit files...')

Gaussians_used_in_analysis_counter = nan(Nchannels, Nproteins);
Gaussians_excluded_from_analysis_counter = nan(Nchannels, Nproteins);

for ci = 1:Nchannels
  Experimental_channel = experimental_channels{ci};
  
  for ri = 1:Nproteins
    if Coef{ci,ri}==9999;
    else
      Num_Gaus=(numel(Coef{ci,ri}))/3;
      Gaus = reshape(Coef{ci,ri},3,Num_Gaus);
      Values_Considered_for_analysis=0;
      Values_Not_Considered_for_analysis=0;
      for i = 1:Num_Gaus
        Height = Gaus(1,i);
        Center = Gaus(2,i)-5;  %5 is because that are added 5 points before chrom
        Width = Gaus(3,i);
        Size_of_complex=SEC_fit(1)*Center+SEC_fit(2);
        %Test guassian properties
        if Height>=0.5 && Width>=1 && Center>1 && Center<Nfractions && ~isnan(cleandata{ci}(ri,round(Center+5)))
          %write out Gaussians for further analysis
          Values_Considered_for_analysis = Values_Considered_for_analysis+1; %Records the number of gaussians which will outputed
          
          %For Gausians
          fileName1A = [datadir3 num2str(ri),'_',mat2str(i),'_',Experimental_channel,'_OutputGaus.csv'];
          fid1A = fopen(fileName1A,'at'); % create the output file
          fprintf (fid1A,'%6.4f,%6.4f,%s,%6.4f,%6.4f,%6.4f,%6.4f,%6.4f,%6.4f,\n',...
            ri, replicate(ri) , txt_MvsL{ri+1}, Height, Center, Width, SSE(ci,ri), adjrsquare(ci,ri), round(Size_of_complex/10)*10);
          fclose(fid1A);
          
          %For Chromatograms
          fileName2A = [datadir1 num2str(ri),'_',mat2str(i),'_',Experimental_channel,'_Output_Chromatograms.csv'];
          fid2A = fopen(fileName2A,'at');
          fprintf(fid2A,'%6.4f,%6.4f,%s,',ri, replicate(ri), txt_MvsL{ri+1});
          fprintf(fid2A,'%6.4g,', cleandata{ci}(ri,:));
          fprintf(fid2A,'\n');
          fclose(fid2A);
        elseif Height<0.5 || Width<1 || Center<=1 || Center>=Nfractions || isnan(cleandata{ci}(ri,round(Center+5)))
          %write out Gaussians excluded from further analysis
          Values_Not_Considered_for_analysis=Values_Not_Considered_for_analysis+1; %Records the number of guassian which will outputed
          
          %For Gausians
          fileName1B = [datadir4 num2str(ri),'_',mat2str(i),'_',Experimental_channel,'_OutputGaus_filtered_out.csv'];
          fid2A = fopen(fileName1B,'at'); % create the output file
          fprintf (fid2A,'%6.4f,%6.4f,%s,%6.4f,%6.4f,%6.4f,%6.4f,%6.4f,%6.4f,\n',...
            ri, nan ,txt_MvsL{ri+1}, Height, Center, Width, SSE, adjrsquare, (round(Size_of_complex/10)*10));
          fclose(fid2A);
          
          %For Chromatograms
          fileName2B = [datadir2 'GaussBuild/' num2str(ri),'_',mat2str(i),'_',Experimental_channel,'_Output_Chromatograms_filtered_out.csv'];
          fid2B = fopen(fileName2B,'at');
          fprintf(fid2B,'%6.4f,%6.4f,%s,',ri, replicate(ri), txt_MvsL{ri+1});
          fprintf(fid2B,'%6.4g,', cleandata{ci}(ri,:));
          fprintf(fid2B,'\n');
          fclose(fid2B);
        end
      end
      Gaussians_used_in_analysis_counter(ci,ri) = Values_Considered_for_analysis;
      Gaussians_excluded_from_analysis_counter(ci,ri) = Values_Not_Considered_for_analysis;
    end
  end
end


%% MvsL_Summary_Gausians_identifed.csv
disp('        Writing MvsL_Summary_Gausians_identifed.csv...')
fn = strcat([datadir 'GaussBuild/MvsL_Summary_Gausians_identifed.csv']);

ci=1;
Mean_No_Quant_String = sum(~isnan(rawdata{ci}(:))) / size(rawdata{ci},1);
Mean_No_After_adding_one_missing_data = sum(tmp1{ci}(:)>0) / size(tmp1{ci},1);
Mean_No_After_filtering = sum(tmp2{ci}(:)>.05) / size(tmp2{ci},1);

fid3 = fopen(fn,'w');
fprintf (fid3,'%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n',...
  'Number_of_proteins','Mean_No_Quant_String', 'Mean_No_After_adding_one_missing_data',...
  'Mean_No_After_filtering', '0 Gausian', '1 Gausian', '2 Gausians',...
  '3 Gausians', '4 Gausians', '5 Gausians', 'No_try_fit');                      %Write Header

fprintf (fid3,'%6.4f,%6.4f,%6.4f,%6.4f,%6.4f,%6.4f,%6.4f,%6.4f,%6.4f,%6.4f,%6.4f\n',...
  Nproteins, Mean_No_Quant_String, Mean_No_After_adding_one_missing_data,Mean_No_After_filtering,...
  sum(No_Gaus==0 & Try_Fit(ci,:)==1), sum(No_Gaus==1 & Try_Fit(ci,:)==1), sum(No_Gaus==2 & Try_Fit(ci,:)==1),...
  sum(No_Gaus==3 & Try_Fit(ci,:)==1), sum(No_Gaus==4 & Try_Fit(ci,:)==1), sum(No_Gaus==5& Try_Fit(ci,:)==1), sum(Try_Fit(1,:)));
fclose(fid3);



%% HvsL_Summary_Gausians_identifed.csv
disp('        Writing HvsL_Summary_Gausians_identifed.csv...')
fn = strcat([datadir 'GaussBuild/HvsL_Summary_Gausians_identifed.csv']);

ci=2;
Mean_No_Quant_String = sum(~isnan(rawdata{ci}(:))) / size(rawdata{ci},1);
Mean_No_After_adding_one_missing_data = sum(tmp1{ci}(:)>0) / size(tmp1{ci},1);
Mean_No_After_filtering = sum(tmp2{ci}(:)>.05) / size(tmp2{ci},1);

fid3 = fopen(fn,'w');
fprintf (fid3,'%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n',...
  'Number_of_proteins','Mean_No_Quant_String', 'Mean_No_After_adding_one_missing_data',...
  'Mean_No_After_filtering', '0 Gausian', '1 Gausian', '2 Gausians',...
  '3 Gausians', '4 Gausians', '5 Gausians', 'No_try_fit');                      %Write Header

fprintf (fid3,'%6.4f,%6.4f,%6.4f,%6.4f,%6.4f,%6.4f,%6.4f,%6.4f,%6.4f,%6.4f,%6.4f\n',...
  Nproteins, Mean_No_Quant_String, Mean_No_After_adding_one_missing_data,Mean_No_After_filtering,...
  sum(No_Gaus==0 & Try_Fit(ci,:)==1), sum(No_Gaus==1 & Try_Fit(ci,:)==1), sum(No_Gaus==2 & Try_Fit(ci,:)==1),...
  sum(No_Gaus==3 & Try_Fit(ci,:)==1), sum(No_Gaus==4 & Try_Fit(ci,:)==1), sum(No_Gaus==5& Try_Fit(ci,:)==1), sum(Try_Fit(1,:)));
fclose(fid3);



%% Proteins_not_fitted_to_gaussian_MvsL.csv
disp('        Writing Proteins_not_fitted_to_gaussian_MvsL.csv...')
fn = strcat([datadir 'GaussBuild/Proteins_not_fitted_to_gaussian_MvsL.csv']);

ci=1;
fid5 = fopen(fn,'w');
for ri = 1:Nproteins
  if length(Coef{ci,ri})<3 && Try_Fit(ci,ri)==1
    fprintf(fid5,'%6.4f,%s,%6.4f,', replicate(ri),txt_MvsL{(ri+1)},replicate(ri));
    fprintf(fid5,'%6.4g,', rawdata{ci}(ri,:));
    fprintf(fid5,'\n');
  end
end
fclose(fid5);



%% Proteins_not_fitted_to_gaussian_HvsL.csv
disp('        Writing Proteins_not_fitted_to_gaussian_HvsL.csv...')
fn = strcat([datadir 'GaussBuild/Proteins_not_fitted_to_gaussian_HvsL.csv']);

ci=2;
fid5 = fopen(fn,'w');
for ri = 1:Nproteins
  if length(Coef{ci,ri})<3 && Try_Fit(ci,ri)==1
    fprintf(fid5,'%6.4f,%s,%6.4f,', replicate(ri),txt_HvsL{(ri+1)},replicate(ri));
    fprintf(fid5,'%6.4g,', rawdata{ci}(ri,:));
    fprintf(fid5,'\n');
  end
end
fclose(fid5);



%% MvsL_Summary_Gausians_for_individual_proteins.csv
disp('        Writing MvsL_Summary_Gausians_for_individual_proteins.csv...')
fn = strcat([datadir 'GaussBuild/MvsL_Summary_Gausians_for_individual_proteins.csv']);

ci=1;
fid4 = fopen(fn,'w');
fprintf (fid4,'%s,%s,%s,%s,%s\n',...
  'Protein_number', 'Gene_name', 'Number_of_Gausians_detected','Number_of_Gausians_within_defined_boundaries','Number_of_Gausians_filtered');               %Write Header
for ri=1:Nproteins
  fprintf(fid4, '%s,%s,%s,%s,%s\n', num2str(ri),...
    txt_MvsL{ri+1},...  % Protein_names
    num2str(No_Gaus(ri)),...
    num2str(Gaussians_used_in_analysis_counter(ci,ri)),...
    num2str(Gaussians_excluded_from_analysis_counter(ci,ri)));
end
fclose(fid4);



%% HvsL_Summary_Gausians_for_individual_proteins.csv
disp('        Writing HvsL_Summary_Gausians_for_individual_proteins.csv...')
fn = strcat([datadir 'GaussBuild/HvsL_Summary_Gausians_for_individual_proteins.csv']);

ci=1;
fid4 = fopen(fn,'w');
fprintf (fid4,'%s,%s,%s,%s,%s\n',...
  'Protein_number', 'Gene_name', 'Number_of_Gausians_detected','Number_of_Gausians_within_defined_boundaries','Number_of_Gausians_filtered');               %Write Header
for ri=1:Nproteins
  fprintf(fid4, '%s,%s,%s,%s,%s\n', num2str(ri),...
    txt_HvsL{ri+1},...  % Protein_names
    num2str(No_Gaus(ri)),...
    num2str(Gaussians_used_in_analysis_counter(ci,ri)),...
    num2str(Gaussians_excluded_from_analysis_counter(ci,ri)));
end
fclose(fid4);



%% MvsL_Combined_OutputGaus.csv
disp('        Writing MvsL_Combined_OutputGaus.csv...')
fn = strcat([datadir 'GaussBuild/MvsL_Combined_OutputGaus.csv']);

ci=1;
Combined_OutputGaus_length = size(protgausI{ci},1);

%Write out Combined process data to csv file
fid_combined_Gaus = fopen(fn,'w');
fprintf (fid_combined_Gaus,'%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n',... %header for OutputGaus output
  'Guassian index number','Protein_number','Replicate','Protein name', 'Height', 'Center', 'Width', 'SSE', 'adjrsquare','Complex size'); %Write Header
for kk = 1:Combined_OutputGaus_length
  ri = protgausI{ci}(kk,1);
  gn = protgausI{ci}(kk,2);
  rep = protgausI{ci}(kk,3);
  Height = Coef{ci,ri}((gn-1)*3 + 1);
  Center = Coef{ci,ri}((gn-1)*3 + 2);
  Width = Coef{ci,ri}((gn-1)*3 + 3);
  Size_of_complex=SEC_fit(1)*Center+SEC_fit(2);
  fprintf(fid_combined_Gaus,'%6.4f,%6.4f,%6.4f,',kk,ri,rep); % index information
  fprintf(fid_combined_Gaus,'%s,',txt_MvsL{ri+1}); % protein name
  fprintf(fid_combined_Gaus,'%6.4f,%6.4f,%6.4f,%6.4f,%6.4f,%6.4f,\n',...
    Height, Center, Width, SSE(ci,ri), adjrsquare(ci,ri), round(Size_of_complex/10)*10); % Gaussian fitting information
end
fclose(fid_combined_Gaus);



%% HvsL_Combined_OutputGaus.csv
disp('        Writing HvsL_Combined_OutputGaus.csv...')
fn = strcat([datadir 'GaussBuild/HvsL_Combined_OutputGaus.csv']);

ci=2;
Combined_OutputGaus_length = size(protgausI{ci},1);

%Write out Combined process data to csv file
fid_combined_Gaus = fopen(fn,'w');
fprintf (fid_combined_Gaus,'%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n',... %header for OutputGaus output
  'Guassian index number','Protein_number','Replicate','Protein name', 'Height', 'Center', 'Width', 'SSE', 'adjrsquare','Complex size'); %Write Header
for kk = 1:Combined_OutputGaus_length
  ri = protgausI{ci}(kk,1);
  gn = protgausI{ci}(kk,2);
  Height = Coef{ci,ri}((gn-1)*3 + 1);
  Center = Coef{ci,ri}((gn-1)*3 + 2);
  Width = Coef{ci,ri}((gn-1)*3 + 3);
  Size_of_complex=SEC_fit(1)*Center+SEC_fit(2);
  fprintf(fid_combined_Gaus,'%6.4f,%6.4f,%6.4f,',kk,ri,replicate(ri)); % index information
  fprintf(fid_combined_Gaus,'%s,',txt_HvsL{ri+1}); % protein name
  fprintf(fid_combined_Gaus,'%6.4f,%6.4f,%6.4f,%6.4f,%6.4f,%6.4f,\n',...
    Height, Center, Width, SSE(ci,ri), adjrsquare(ci,ri), round(Size_of_complex/10)*10); % Gaussian fitting information
end
fclose(fid_combined_Gaus);



%% MvsL_Combined_Chromatograms.csv
disp('        Writing MvsL_Combined_Chromatograms.csv...')
fn = strcat([datadir 'GaussBuild/MvsL_Combined_Chromatograms.csv']);

ci=1;
Combined_OutputGaus_length = size(protgausI{ci},1);

fid_combined_Chromatogram = fopen(fn,'w');
for kk = 1:Combined_OutputGaus_length
  ri = protgausI{ci}(kk,1);
  fprintf(fid_combined_Chromatogram,'%6.4f,%6.4f,%6.4f,',kk,ri,replicate(ri)); %Write out the index information
  fprintf(fid_combined_Chromatogram,'%s,',txt_MvsL{ri+1}); %Write out protein name
  fprintf(fid_combined_Chromatogram,'%6.4g,',cleandata{ci}(ri,5:end)); %Chromatogram information
  fprintf(fid_combined_Chromatogram,'\n');
end
fclose(fid_combined_Chromatogram);



%% HvsL_Combined_Chromatograms.csv
disp('        Writing HvsL_Combined_Chromatograms.csv...')
fn = strcat([datadir 'GaussBuild/HvsL_Combined_Chromatograms.csv']);

ci=2;
Combined_OutputGaus_length = size(protgausI{ci},1);

fid_combined_Chromatogram = fopen(fn,'w');
for kk = 1:Combined_OutputGaus_length
  ri = protgausI{ci}(kk,1);
  fprintf(fid_combined_Chromatogram,'%6.4f,%6.4f,%6.4f,',kk,ri,replicate(ri)); %Write out the index information
  fprintf(fid_combined_Chromatogram,'%s,',txt_MvsL{ri+1}); %Write out protein name
  fprintf(fid_combined_Chromatogram,'%6.4g,',cleandata{ci}(ri,5:end)); %Chromatogram information
  fprintf(fid_combined_Chromatogram,'\n');
end
fclose(fid_combined_Chromatogram);



%% MvsL_Combined_OutputGaus_rep*.csv
% Divide up MvsL_Combined_OutputGaus.csv
disp('        Writing MvsL_Combined_OutputGaus_rep*.csv and HvsL_Combined_OutputGaus_rep*.csv...')

Experimental_channels = {'MvsL' 'HvsL' 'HvsM'};
Unique_replicate = unique(replicate);

for ci = 1:2
  Experimental_channel = Experimental_channels{ci};
  for divider_counter1 = 1:length(Unique_replicate)
    %Create name of gaussian file to output divided gaus data to
    Process_Gaus_import_Name= [datadir 'GaussBuild/' Experimental_channel '_Combined_OutputGaus_rep' mat2str(divider_counter1) '.csv'];
    
    fid_processing= fopen(Process_Gaus_import_Name,'wt'); % create the output file with the header infromation
    fprintf (fid_processing,'%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n',... %header for OutputGaus output
      'Guassian index number', 'Protein_number','Replicate',...
      'Protein name', 'Height', 'Center', 'Width', 'SSE',...
      'adjrsquare', 'Complex Size'); %Write Header
    
    for kk = 1:Ngauss(ci)
      ri = protgausI{ci}(kk,1);
      if replicate(ri) == divider_counter1
        gn = protgausI{ci}(kk,2);
        Height = Coef{ci,ri}((gn-1)*3 + 1);
        Center = Coef{ci,ri}((gn-1)*3 + 2) - 5;
        Width = Coef{ci,ri}((gn-1)*3 + 3);
        Size_of_complex=SEC_fit(1)*Center+SEC_fit(2);
        fprintf(fid_processing,'%6.4f,%6.4f,%6.4f,',kk,ri,replicate(ri)); % index information
        fprintf(fid_processing,'%s,',txt_MvsL{ri+1}); % protein name
        fprintf(fid_processing,'%6.4f,%6.4f,%6.4f,%6.4f,%6.4f,%6.4f,\n',...
          Height, Center, Width, SSE(ci,ri), adjrsquare(ci,ri), round(Size_of_complex/10)*10); % Gaussian fitting information
      end
    end
    fclose(fid_processing);
  end
end


%% MvsL_Summary_Gausians_for_individual_proteins_rep*.csv
% Divide up MvsL_Summary_Gausians_for_individual_proteins.csv
disp('        Writing MvsL_Summary_Gausians_for_individual_proteins_rep*.csv and MvsL_Summary_Gausians_for_individual_proteins_rep*.csv...')

Experimental_channels = {'MvsL' 'HvsL' 'HvsM'};
Unique_replicate = unique(replicate);

for ci = 1:2
  Experimental_channel = Experimental_channels{ci};
  for divider_counter1=1:Unique_replicate
    %Create name of gaussian file to output divided gaus data to
    Process_Summary_gausian_info_Name= [datadir 'GaussBuild/' Experimental_channel,'_Summary_Gausians_for_individual_proteins_rep',mat2str(divider_counter1),'.csv'];
    
    fid_processing= fopen(Process_Summary_gausian_info_Name,'wt'); % create the output file with the header infromation
    fprintf (fid_processing,'%s,%s,%s,%s,%s\n',... %header for OutputGaus output
      'Protein_number', 'Protein_name','Number_of_Gausians_detected',...
      'Number_of_Gausians_within_defined_boundaries',...
      'Number_of_Gausians_filtered'); %Write Header
    
    for kk = 1:Ngauss(ci)
      ri = protgausI{ci}(kk,1);
      if replicate(ri) == divider_counter1
        fprintf(fid4, '%s,%s,%s,%s,%s\n', num2str(ri),...
          txt_MvsL{ri+1},...  % Protein_names
          num2str(No_Gaus(ri)),...
          num2str(Gaussians_used_in_analysis_counter(ci,ri)),...
          num2str(Gaussians_excluded_from_analysis_counter(ci,ri)));
      end
    end
    fclose(fid_processing);
  end
end



%% MvsL_Raw_data_maxquant_rep*.csv
% Divide up MvsL_Raw_data_maxquant.csv
disp('        Writing MvsL_Raw_data_maxquant_rep*.csv...')

ci = 1;

%Divid up Summary_Gausians_for_individual_proteins_rep
for divider_counter1=1:Unique_replicate
  
  %Create name of gaussian file to output divided gaus data for MvsL
  Process_Replicate_raw_data1= [datadir 'GaussBuild/' Experimental_channels{ci},'_Raw_data_maxquant_rep',mat2str(divider_counter1),'.csv'];
  
  fid_processing3= fopen(Process_Replicate_raw_data1,'wt'); % create the output file with the header infromation
  fprintf(fid_processing3,'%s, ', txt_MvsL{1,1:Nfractions+2}); %header for OutputGaus output
  fprintf(fid_processing3,'\n');
  
  for kk = 1:Ngauss(ci)
    ri = protgausI{ci}(kk,1);
    if replicate(ri) == divider_counter1
      fprintf(fid_processing3,'%s, %6.4f,', txt_MvsL{ri+1}, num_val_MvsL(ri,1));
      fprintf(fid_processing3,'%6.4g,',rawdata{ci}(ri,2:end)); %Chromatogram information
      fprintf(fid_processing3,'\n');
    end
  end
  fclose(fid_processing3);
end



%% HvsL_Raw_data_maxquant_rep*.csv
% Divide up HvsL_Raw_data_maxquant.csv
disp('        Writing HvsL_Raw_data_maxquant_rep*.csv...')

ci = 2;

%Divid up Summary_Gausians_for_individual_proteins_rep
for divider_counter1=1:Unique_replicate
  
  %Create name of gaussian file to output divided gaus data for MvsL
  Process_Replicate_raw_data1= [datadir 'GaussBuild/' Experimental_channels{ci},'_Raw_data_maxquant_rep',mat2str(divider_counter1),'.csv'];
  
  fid_processing3= fopen(Process_Replicate_raw_data1,'wt'); % create the output file with the header infromation
  fprintf(fid_processing3,'%s, ', txt_HvsL{1,1:Nfractions+2}); %header for OutputGaus output
  fprintf(fid_processing3,'\n');
  
  for kk = 1:Ngauss(ci)
    ri = protgausI{ci}(kk,1);
    if replicate(ri) == divider_counter1
      fprintf(fid_processing3,'%s, %6.4f,', txt_HvsL{ri+1}, num_val_MvsL(ri,1));
      fprintf(fid_processing3,'%6.4g,',rawdata{ci}(ri,2:end)); %Chromatogram information
      fprintf(fid_processing3,'\n');
    end
  end
  fclose(fid_processing3);
end




%% HvsM_Raw_data_maxquant_rep*.csv
% Divide up HvsM_Raw_data_maxquant.csv
if 0
disp('        Writing HvsM_Raw_data_maxquant_rep*.csv...')

ci = 3;
HvsM_filename=[datadir0 'Combined_replicates_2014_04_22_contaminates_removed_for_HvsM_scripts.xlsx'];
[num_val_HvsM,txt_HvsM] = xlsread(HvsM_filename); %Import file raw Maxqaunt output %#

%Divid up Summary_Gausians_for_individual_proteins_rep
for divider_counter1=1:Unique_replicate
  
  %Create name of gaussian file to output divided gaus data for MvsL
  Process_Replicate_raw_data1= strcat(Experimental_channels{ci},'_Raw_data_maxquant_rep',mat2str(divider_counter1),'.csv');
  
  fid_processing3= fopen(Process_Replicate_raw_data1,'wt'); % create the output file with the header infromation
  fprintf(fid_processing3,'%s, ', txt_MvsL{1,1:Nfractions+2}); %header for OutputGaus output
  fprintf(fid_processing3,'\n');
  
  for kk = 1:Ngauss(ci)
    ri = protgausI{ci}(kk,1);
    if replicate(ri) == divider_counter1
      fprintf(fid_processing3,'%s, %6.4f,', txt_HvsM{ri+1,1}, num_val_HvsM(ri,1));
      fprintf(fid_processing3,'%6.4g,',num_val_HvsM(ri,2:end)); %Chromatogram information
      fprintf(fid_processing3,'\n');
    end
  end
  fclose(fid_processing3);
end
end
