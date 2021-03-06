%WRITEOUTPUT_ALIGNMENT Writes output tables for the PRINCE Alignment module.


New_Fraction_numbering = -5:fraction_number+5;
lenght_new_fraction = length(New_Fraction_numbering);
for ci = 1:length(Experimental_channels)
  Experimental_channel = Experimental_channels{ci};
  for rr = 1:Nreplicates
    Column_name = cell(lenght_new_fraction,1);
    for column_counter=1:lenght_new_fraction
      Column_name{column_counter}=strcat('Ratio_',Experimental_channel,'_',mat2str(New_Fraction_numbering(column_counter)));
    end
  end
end



%%        'Adjusted_',Experimental_channel,'_Raw_data_maxquant_rep',mat2str(alignment_counter),'.csv'
% disp('        Adjusted_*vsL_Raw_data_maxquant_rep*.csv')
% 
% for ci = 1:length(Experimental_channels)
%   Experimental_channel = Experimental_channels{ci};
%   for rr = 1:Nreplicates
%     Column_name = cell(lenght_new_fraction,1);
%     for column_counter=1:lenght_new_fraction
%       Column_name{column_counter}=strcat('Ratio_',Experimental_channel,'_',mat2str(New_Fraction_numbering(column_counter)));
%     end
%     
%     I = find(cleandata{1}(:,1) == rr);
%     
%     fid9_Name = [maindir '/Output/tmp/' 'Adjusted_' Experimental_channel '_Raw_data_maxquant_rep' mat2str(rr) '.csv'];
%     Column_header=repmat('%s,', 1, (lenght_new_fraction)+2);
%     
%     fid9 = fopen(fid9_Name,'at');
%     fprintf (fid9, [Column_header, '\n'],'Protein name', 'Replicate', Column_name{:});  %Write Header
%     for ri=1:length(I)
%       ii = I(ri);
%       
%       fprintf(fid9, '%s,', txt_val{ci}{ii+1});
%       fprintf(fid9,'%6.4f,', rr);
%       fprintf(fid9,'%6.4f,', adjusted_raw_data{ci}(ii,:));
%       fprintf(fid9,'\n');
%     end
%     fclose(fid9);
%     
%   end
% end


%%    'Adjusted_',Experimental_channel,'_Raw_for_ROC_analysis_rep',mat2str(alignment_counter),'.csv'

for ci = 1:length(Experimental_channels)
  Experimental_channel = Experimental_channels{ci};
  for rr = 1:Nreplicates
    
    Irawdata=zeros(1,1);
    %find location of protein to write out
    for ii = 1:size(Gaus_import{ci,rr}.textdata,1)
      %find location of protein to write out
      %location_Protein_in_Raw = ind2sub(length(Summary_gausian_infomration{ci,rr}.textdata(:,2)),...
      %  strmatch(Gaus_import{ci,rr}.textdata(Roc_counter+1,1),Summary_gausian_infomration{ci,rr}.textdata(:,2),'exact'));
      %locations_of_protein_data(Roc_counter) = location_Protein_in_Raw(1);
      protName = Gaus_import{ci,rr}.textdata{ii};
      I = find(ismember(txt_val{ci}(2:end),protName) & replicates==rr);
      Irawdata(ii) = I(1);
    end
    
    fid9B_Name = [maindir '/Output/tmp/' 'Adjusted_' Experimental_channel '_Raw_for_ROC_analysis_rep' mat2str(rr) '.csv'];
    Column_header=repmat('%s,', 1, ((lenght_new_fraction)+2));
    
    tmp1 = zeros(length(Irawdata),size(adjusted_raw_data{ci},2));
    
    fid9B = fopen(fid9B_Name,'at');
    fprintf (fid9B, [Column_header, '\n'],...
      'Protein name', 'Replicate', Column_name{:});  %Write Header
    for Roc_counter2 = 1:length(Irawdata)
      fprintf(fid9B, '%s,', txt_val{ci}{Irawdata(Roc_counter2)+1,1});
      fprintf(fid9B,'%6.4f,',rr);
      fprintf(fid9B,'%6.4f,', adjusted_raw_data{ci}(Irawdata(Roc_counter2),:));
      fprintf(fid9B,'\n');
      
      tmp1(Roc_counter2,:) = adjusted_raw_data{ci}(Irawdata(Roc_counter2),:);
    end
    fclose(fid9B);
    
  end
end



%%    fid9B_Name = strcat('Adjusted_HvsM_Raw_data_maxquant_rep',mat2str(alignment_counter),'.csv');



%%    fid6_Name = strcat('Adjusted_Chromatograms_vobose_rep',mat2str(alignment_counter),'_','.csv');



%%    fid7_Name = strcat('Adjusted_Combined_OutputGaus_rep',mat2str(alignment_counter),'.csv');

for ci = 1:length(Experimental_channels)
  Experimental_channel = Experimental_channels{ci};
  for rr = 1:Nreplicates
    
    fid7_Name = [maindir '/Output/tmp/'  'Adjusted_Combined_OutputGaus_' Experimental_channel '_rep' mat2str(rr) '.csv'];
    fid7 = fopen(fid7_Name,'at');
    fprintf (fid7,'%s,%s,%s,%s,%s,%s,%s\n',...
      'Protein name', 'Height', 'Center','Width','SSE','adjrsquare', 'Complex Size');  %Write Header
    for ri = 1:size(Adjusted_Gaus_import{ci,rr}.data,1)
      fprintf(fid7, '%s,', Adjusted_Gaus_import{ci,rr}.textdata{ri,1});
      fprintf(fid7,'%6.4f,', Adjusted_Gaus_import{ci,rr}.data(ri,:));
      fprintf(fid7,'\n');
    end
    fclose(fid7);
  end
end



%%  fid10_Name = strcat('Adjusted_',Experimental_channel,'_Combined_OutputGaus.csv');

for ci = 1:length(Experimental_channels)
  Experimental_channel = Experimental_channels{ci};
  
  fid7_Name = [maindir '/Output/tmp/' 'Adjusted_' Experimental_channel '_Combined_OutputGaus.csv'];
  fid = fopen(fid7_Name,'wt');
  fprintf (fid,'%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n',...
    'Guassian index number', 'Protein_number','Replicate',...
    'Protein name', 'Height', 'Center','Width','SSE','adjrsquare', 'Complex Size');  %Write Header
  
  cc = 0;
  for rr = 1:Nreplicates
    for combined_gaus_writeout = 1:size(Adjusted_Gaus_import{ci,rr}.data,1)
      fprintf(fid,'%6.4f,%6.4f,%6.4f,',cc,0,rr); % index information
      fprintf(fid, '%s,', Adjusted_Gaus_import{ci,rr}.textdata{combined_gaus_writeout,1});
      fprintf(fid,'%6.4f,', Adjusted_Gaus_import{ci,rr}.data(combined_gaus_writeout,:));
      fprintf(fid,'\n');
    end
  end
  fclose(fid);
end



%%  fid11_Name = strcat('Adjusted_',Experimental_channel,'_Raw_data_maxquant.csv');

for ci = 1:length(Experimental_channels)
  Experimental_channel = Experimental_channels{ci};
  
  fid11_Name = [maindir '/Output/tmp/' 'Adjusted_' Experimental_channel '_Raw_data_maxquant.csv'];
  fid11 = fopen(fid11_Name,'wt');
  %fprintf (fid11,[Column_header, '\n'], Title_import2{:});  %Write Header
  % protein name, replicate, data(:)
  fprintf(fid11,'%s,%s,%s\n','Protein name','Replicate','Ratios');
  
  for ii = 1:size(adjusted_raw_data{ci},1)
    fprintf(fid11, '%s,',  txt_val{ci}{ii+1,1});
    fprintf(fid11, '%6f,', replicates(ii));
    fprintf(fid11, '%f,',  adjusted_raw_data{ci}(ii,:));
    fprintf(fid11,'\n');
  end
  fclose(fid11);
end

