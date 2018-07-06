%WRITEOUTPUT_GAUSSBUILD Writes output tables for the PRINCE GaussBuild
%   module.




%% 0. Initialize

try
  SEC_fit=polyfit(SEC_size_alignment(1,:),SEC_size_alignment(2,:),1);
catch
  disp('Gauss_Build: writeOutput: SEC fitting failed. Size of Complex will be zero.')
  SEC_fit = [0 0];
end

Nchannels = size(Coef,1);
Nproteins = size(Coef,2);

% Make No_Gaus, i.e. the number of Gaussians selected for each protein
No_Gaus = zeros(1,Nproteins);
for ri = 1:Nproteins
  No_Gaus(ri) = length(Coef{1,ri})/3;
end


%% Make files for individual fits
% 1_1_*_OutputGaus.csv
% 1_1_*_Output_Chromatograms.csv
% 1_1_*_OutputGaus_filtered_out.csv
% 1_1_*_Output_Chromatograms_filtered_out.csv
%disp('        Writing individual fit files...')

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
          
        elseif Height<0.5 || Width<1 || Center<=1 || Center>=Nfractions || isnan(cleandata{ci}(ri,round(Center+5)))
          %write out Gaussians excluded from further analysis
          Values_Not_Considered_for_analysis=Values_Not_Considered_for_analysis+1; %Records the number of guassian which will outputed
          
        end
      end
      Gaussians_used_in_analysis_counter(ci,ri) = Values_Considered_for_analysis;
      Gaussians_excluded_from_analysis_counter(ci,ri) = Values_Not_Considered_for_analysis;
    end
  end
end




%% *_Summary_Gausians_for_individual_proteins.csv
%disp('        Writing *_Summary_Gausians_for_individual_proteins.csv...')
if 0
  for ci = 1:Nchannels
    fn = strcat([datadir user.silacratios{ci} '_Summary_Gausians_for_individual_proteins.csv']);
    fid4 = fopen(fn,'w');
    fprintf (fid4,'%s,%s,%s,\n',...
      'Protein_number', 'Gene_name', 'Number_of_Gausians_detected');               %Write Header
    for ri=1:Nproteins
      fprintf(fid4, '%s,%s,%s,\n', num2str(ri),...
        txt_val{ci}{ri+1},...  % Protein_names
        num2str(No_Gaus(ri)));
      %num2str(length(Coef{ci,ri})/3),...
      %num2str(Gaussians_excluded_from_analysis_counter(ci,ri)));
    end
    fclose(fid4);
  end
end



%% *_Combined_OutputGaus.csv
%disp('        Writing *_Combined_OutputGaus.csv...')

for ci = 1:Nchannels
  fn = strcat([datadir user.silacratios{ci} '_Combined_OutputGaus.csv']);
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
    fprintf(fid_combined_Gaus,'%s,',txt_val{ci}{ri+1}); % protein name
    fprintf(fid_combined_Gaus,'%6.4f,%6.4f,%6.4f,%6.4f,%6.4f,%6.4f,\n',...
      Height, Center, Width, SSE(ci,ri), adjrsquare(ci,ri), round(Size_of_complex/10)*10); % Gaussian fitting information
  end
  fclose(fid_combined_Gaus);
end





%% *_Combined_OutputGaus_rep*.csv
% Divide up MvsL_Combined_OutputGaus.csv
%disp('        Writing *_Combined_OutputGaus_rep*.csv...')

for ci = 1:Nchannels
  Experimental_channel = user.silacratios{ci};
  Unique_replicate = unique(replicate{ci});
  for divider_counter1 = 1:length(Unique_replicate)
    %Create name of gaussian file to output divided gaus data to
    Process_Gaus_import_Name= [maindir '/Output/tmp/' Experimental_channel '_Combined_OutputGaus_rep' mat2str(divider_counter1) '.csv'];
    
    fid_processing= fopen(Process_Gaus_import_Name,'wt'); % create the output file with the header infromation
    fprintf (fid_processing,'%s,%s,%s,%s,%s,%s,%s\n',... %header for OutputGaus output
      'Protein name', 'Height', 'Center', 'Width', 'SSE',...
      'adjrsquare', 'Complex Size'); %Write Header
    
    for kk = 1:Ngauss(ci)
      ri = protgausI{ci}(kk,1);
      if replicate{ci}(ri) == divider_counter1
        gn = protgausI{ci}(kk,2);
        Height = Coef{ci,ri}((gn-1)*3 + 1);
        Center = Coef{ci,ri}((gn-1)*3 + 2) - 5;
        Width = Coef{ci,ri}((gn-1)*3 + 3);
        Size_of_complex=SEC_fit(1)*Center+SEC_fit(2);
        fprintf(fid_processing,'%s,',txt_val{ci}{ri+1}); % protein name
        fprintf(fid_processing,'%6.4f,%6.4f,%6.4f,%6.4f,%6.4f,%6.4f,\n',...
          Height, Center, Width, SSE(ci,ri), adjrsquare(ci,ri), round(Size_of_complex/10)*10); % Gaussian fitting information
      end
    end
    fclose(fid_processing);
  end
end


%% *_Raw_data_maxquant_rep*.csv
% Divide up *_Raw_data_maxquant.csv
%disp('        Writing *_Raw_data_maxquant_rep*.csv...')

for ci = 1:Nchannels
  Unique_replicate = unique(replicate{ci});

  %Divid up Summary_Gausians_for_individual_proteins_rep
  for divider_counter1 = 1:length(Unique_replicate)
    
    %Create name of gaussian file to output divided gaus data for MvsL
    Process_Replicate_raw_data1= [maindir '/Output/tmp/' user.silacratios{ci},'_Raw_data_maxquant_rep',mat2str(divider_counter1),'.csv'];
    
    fid_processing3= fopen(Process_Replicate_raw_data1,'wt'); % create the output file with the header infromation
    fprintf(fid_processing3,'%s, ', txt_val{ci}{1,1:end}); %header for OutputGaus output
    fprintf(fid_processing3,'\n');
    
    for kk = 1:Ngauss(ci)
      ri = protgausI{ci}(kk,1);
      if replicate{ci}(ri) == divider_counter1
        fprintf(fid_processing3,'%s, %6.4f,', txt_val{ci}{ri+1}, replicate{ci}(ri));
        fprintf(fid_processing3,'%6.4g,',rawdata{ci}(ri,:)); %Chromatogram information
        fprintf(fid_processing3,'\n');
      end
    end
    fclose(fid_processing3);
  end
end

