
%Figure showing changes in response to treatment

%Graph data as log2 scatter
figure
%log2_sorted=sort(normaliased_log_2_fold_comparsion_of_gaussian(:));
log2_sorted=sort(Combined_Gaussians.log2_normalised_gaussians(:));
subplot(2,1,1);hold on
P4 =scatter(1:length(log2_sorted(:)),log2_sorted(:),8,'fill');
title('Log2 changes between replicate: BN-PAGE Mitochondira Fas treated','FontSize', 10);
ylabel('Normalised Log2 ratio','FontSize', 10);
xlabel('Gaussians identified','FontSize', 10);
%Set limits based on observed values
%ylim([log2_sorted(1)*1.05,log2_sorted(Total_number_of_unique_gaussians)*1.05]);
%Set limits based on user defined cut off
plot([-100 length(log2_sorted(:))+100],[1 1],':');
plot([-100 length(log2_sorted(:))+100],[-1 -1],':');
xlim([-10,sum(~isnan(log2_sorted(:)))+10]);
ylim([-6,6]);
grid on

%plot distribution of protein chnages across the fractions
%define colours
subplot(2,1,2);
f6_figure=bar(1:Nfraction, [Hist_array(:,1) Hist_array(:,2) Hist_array(:,3)],0.6, 'stack');
for k=1:3
  set(f6_figure(k),'facecolor', myC(k,:), 'EdgeColor', 'k' )
end
legend('No change', 'Increase', 'Decrease', 'FontSize',8, 'Location', 'Best');
xlim([-1,frac2+1]);
title('Changes in Gaussians observed across fractions');
xlabel('Fractions');
ylabel('# of Guassian centers');
saveas(gcf, [figdir 'Comparison/All_gaussians_changes_observed.png']);



% for figure_counter=1:replicate_num
%
%   %Determine Threshold to use for determine signficant change, use 1.96*stdev
%   for writeout_counter2= 1:length(Unique_protein_names)
%     percentage_derivation(writeout_counter2)= log2(Finalised_Master_Gaussian_list.Summed_Gaussian_areaMvsL(writeout_counter2,figure_counter)/...
%       Finalised_Master_Gaussian_list.Summed_Gaussian_areaHvsL(writeout_counter2,figure_counter));
%     if isinf(percentage_derivation(writeout_counter2))==1
%       percentage_derivation(writeout_counter2)=NaN;
%     end
%   end
%
%   %Remove nan values
%   percentage_derivation(isnan(percentage_derivation))=[];
%
%   %calculate stdev
%   std_of_area_measurment=std(percentage_derivation);
%   if isnan(std_of_area_measurment)
%     threshold_of_change=max(abs(percentage_derivation));
%   else
%     threshold_of_change=1.96*std_of_area_measurment;
%   end
%
%   %Figure showing coverage within PCP-SEC-SILAC experiments
%   %%Count how many proteins change across the observed SEC fraction
%   bin_number_hist_MvsL= 0.05:0.05:ceil(max(max(Finalised_Master_Gaussian_list.Summed_Gaussian_areaMvsL)));
%   bin_number_hist_MvsL_plus1= 0.05:0.05:(ceil(max(max(Finalised_Master_Gaussian_list.Summed_Gaussian_areaMvsL)))+0.05);
%   bin_number_hist2_MvsL=length(bin_number_hist_MvsL);
%
%   Hist_array3A=zeros(Dimension_of_master_gaussian_list(1),(bin_number_hist2_MvsL+1));
%   Hist_array3B=zeros(Dimension_of_master_gaussian_list(1),(bin_number_hist2_MvsL+1));
%   Hist_array3C=zeros(Dimension_of_master_gaussian_list(1),(bin_number_hist2_MvsL+1));
%   Hist_array3D=zeros(Dimension_of_master_gaussian_list(1),(bin_number_hist2_MvsL+1));
%
%   %Count the number of gaussian detected vs the area
%   for hist2_counter1= 1:(bin_number_hist2_MvsL)+1
%     for writeout_counter2= 1:Dimension_of_master_gaussian_list(1)
%       %determine bin values
%       Coverage_Bin_number=ceil(Finalised_Master_Gaussian_list.Summed_Gaussian_areaMvsL(writeout_counter2,figure_counter)/0.05);
%       %determine if the fold change is within the bin value
%       if  hist2_counter1-1==Coverage_Bin_number
%         HvsL_values=Finalised_Master_Gaussian_list.Summed_Gaussian_areaHvsL(writeout_counter2,figure_counter);
%         MvsL_values=Finalised_Master_Gaussian_list.Summed_Gaussian_areaMvsL(writeout_counter2,figure_counter);
%         %determine bin of the total observed isotopologue 'volumn' in MvsL
%         Hist_array3A(writeout_counter2,hist2_counter1)=...
%           Hist_array3A(writeout_counter2,hist2_counter1) + 1;
%         %determine if the 'volumn' is different from HvsL
%         if log2(MvsL_values/HvsL_values)<=threshold_of_change &...
%             log2(MvsL_values/HvsL_values)>=-threshold_of_change
%           Hist_array3B(writeout_counter2,hist2_counter1)=...
%             Hist_array3B(writeout_counter2,hist2_counter1) + 1;
%         elseif log2(MvsL_values/HvsL_values)>threshold_of_change
%           Hist_array3C(writeout_counter2,hist2_counter1)=...
%             Hist_array3C(writeout_counter2,hist2_counter1) + 1;
%         elseif log2(MvsL_values/HvsL_values)<-threshold_of_change
%           Hist_array3D(writeout_counter2,hist2_counter1)=...
%             Hist_array3D(writeout_counter2,hist2_counter1) + 1;
%         end
%       end
%     end
%   end
%
%   %Figure showing coverage within PCP-SEC-SILAC experiments
%   %%Count how many proteins change across the observed SEC fraction
%   bin_number_hist_HvsL= 0.05:0.05:ceil(max(max(Finalised_Master_Gaussian_list.Summed_Gaussian_areaHvsL)));
%   bin_number_hist_HvsL_plus1= 0.05:0.05:(ceil(max(max(Finalised_Master_Gaussian_list.Summed_Gaussian_areaHvsL)))+0.05);
%   bin_number_hist2_HvsL=length(bin_number_hist_HvsL);
%
%   Hist_array4A=zeros(Dimension_of_master_gaussian_list(1),(bin_number_hist2_HvsL+1));
%   Hist_array4B=zeros(Dimension_of_master_gaussian_list(1),(bin_number_hist2_HvsL+1));
%   Hist_array4C=zeros(Dimension_of_master_gaussian_list(1),(bin_number_hist2_HvsL+1));
%   Hist_array4D=zeros(Dimension_of_master_gaussian_list(1),(bin_number_hist2_HvsL+1));
%
%   %Count the number of gaussian detected vs the area
%   for hist2_counter1= 1:(bin_number_hist2_HvsL)+1
%     for writeout_counter2= 1:Dimension_of_master_gaussian_list(1)
%       %determine bin values
%       Coverage_Bin_number=ceil(Finalised_Master_Gaussian_list.Summed_Gaussian_areaHvsL(writeout_counter2,figure_counter)/0.05);
%       %determine if the fold change is within the bin value
%       if  hist2_counter1-1==Coverage_Bin_number
%         HvsL_values=Finalised_Master_Gaussian_list.Summed_Gaussian_areaHvsL(writeout_counter2,figure_counter);
%         MvsL_values=Finalised_Master_Gaussian_list.Summed_Gaussian_areaMvsL(writeout_counter2,figure_counter);
%         %determine bin of the total observed isotopologue 'volumn' in MvsL
%         Hist_array4A(writeout_counter2,hist2_counter1)=...
%           Hist_array4A(writeout_counter2,hist2_counter1) + 1;
%         %determine if the 'volumn' is different from HvsL
%         if log2(HvsL_values/MvsL_values)<=threshold_of_change &...
%             log2(HvsL_values/MvsL_values)>=-threshold_of_change
%           Hist_array4B(writeout_counter2,hist2_counter1)=...
%             Hist_array4B(writeout_counter2,hist2_counter1) + 1;
%         elseif log2(HvsL_values/MvsL_values)>threshold_of_change
%           Hist_array4C(writeout_counter2,hist2_counter1)=...
%             Hist_array4C(writeout_counter2,hist2_counter1) + 1;
%         elseif log2(HvsL_values/MvsL_values)<-threshold_of_change
%           Hist_array4D(writeout_counter2,hist2_counter1)=...
%             Hist_array4D(writeout_counter2,hist2_counter1) + 1;
%         end
%       end
%     end
%   end
%
%   %add value to a single array to grpah
%   Isotopologue_areaA(:,1)=sum(Hist_array3A)';
%   Isotopologue_areaA(:,2)=sum(Hist_array3B)';
%   Isotopologue_areaA(:,3)=sum(Hist_array3C)';
%   Isotopologue_areaA(:,4)=sum(Hist_array3D)';
%
%   Isotopologue_areaB(:,1)=sum(Hist_array4A)';
%   Isotopologue_areaB(:,2)=sum(Hist_array4B)';
%   Isotopologue_areaB(:,3)=sum(Hist_array4C)';
%   Isotopologue_areaB(:,4)=sum(Hist_array4D)';
%
%   try
%     %Find bin which contains the max
%     [fitted_valuesMvsL fitted_statsMvsL] = fit(bin_number_hist_MvsL_plus1.',Isotopologue_areaA(:,1),'gauss1');
%   catch
%     fitted_valuesMvsL.a1=NaN;
%     fitted_valuesMvsL.b1=NaN;
%     fitted_valuesMvsL.c1=NaN;
%     fitted_statsMvsL.rsquare=NaN;
%   end
%
%   %Generate Gaussian for MvsL and HvsL
%   for test_value=1:(bin_number_hist2_MvsL+1)
%     %MvsL
%     gaussian_fit_MvsL_bin(test_value)=(fitted_valuesMvsL.a1*exp(-((0.05*(test_value)- fitted_valuesMvsL.b1)...
%       /fitted_valuesMvsL.c1).^2));
%   end
%
%   try
%     %Find bin which contains the max
%     [fitted_valuesHvsL fitted_statsHvsL] = fit(bin_number_hist_HvsL_plus1.',Isotopologue_areaB(:,1),'gauss3');
%   catch
%     fitted_valuesHvsL.a1=NaN;
%     fitted_valuesHvsL.a2=NaN;
%     fitted_valuesHvsL.a3=NaN;
%     fitted_valuesHvsL.b1=NaN;
%     fitted_valuesHvsL.b2=NaN;
%     fitted_valuesHvsL.b3=NaN;
%     fitted_valuesHvsL.c1=NaN;
%     fitted_valuesHvsL.c2=NaN;
%     fitted_valuesHvsL.c3=NaN;
%     fitted_statsHvsL.rsquare=NaN;
%   end
%
%   for test_value=1:(bin_number_hist2_HvsL+1)
%     %HvsL
%     gaussian_fit_HvsL_bin1(test_value)=(fitted_valuesHvsL.a1*exp(-((0.05*(test_value)- fitted_valuesHvsL.b1)...
%       /fitted_valuesHvsL.c1).^2));
%     gaussian_fit_HvsL_bin2(test_value)=(fitted_valuesHvsL.a2*exp(-((0.05*(test_value)- fitted_valuesHvsL.b2)...
%       /fitted_valuesHvsL.c2).^2));
%     gaussian_fit_HvsL_bin3(test_value)=(fitted_valuesHvsL.a3*exp(-((0.05*(test_value)- fitted_valuesHvsL.b3)...
%       /fitted_valuesHvsL.c3).^2));
%   end
%
%
%   %Figure showing coverage within PCP-SEC-SILAC experiments
%   %%Count how many proteins change across the observed SEC fraction
%   bin_number_hist_SEC=0.05:0.05:ceil(max(max(Finalised_Master_Gaussian_list.SEC_coverageMvsL)));
%   bin_number_hist_SEC_plus1= 0.05:0.05:(ceil(max(max(Finalised_Master_Gaussian_list.SEC_coverageMvsL)))+0.05);
%   bin_number_hist2_SEC=length(bin_number_hist_SEC);
%
%   Hist_array5A_SEC=zeros(Dimension_of_master_gaussian_list(1),(bin_number_hist2_SEC+1));
%   Hist_array5B_SEC=zeros(Dimension_of_master_gaussian_list(1),(bin_number_hist2_SEC+1));
%   Hist_array5C_SEC=zeros(Dimension_of_master_gaussian_list(1),(bin_number_hist2_SEC+1));
%   Hist_array5D_SEC=zeros(Dimension_of_master_gaussian_list(1),(bin_number_hist2_SEC+1));
%
%   %Count the number of gaussian detected vs the area
%   for hist2_counter1= 1:(bin_number_hist2_SEC)+1
%     for writeout_counter2= 1:Dimension_of_master_gaussian_list(1)
%       %determine bin values
%       Coverage_Bin_number=ceil(Finalised_Master_Gaussian_list.SEC_coverageMvsL(writeout_counter2,figure_counter)/0.05);
%       %determine if the fold change is within the bin value
%       if  hist2_counter1-1==Coverage_Bin_number
%         HvsL_values=Finalised_Master_Gaussian_list.Summed_Gaussian_areaHvsL(writeout_counter2,figure_counter);
%         MvsL_values=Finalised_Master_Gaussian_list.Summed_Gaussian_areaMvsL(writeout_counter2,figure_counter);
%         %determine bin of the total observed isotopologue 'volumn' in MvsL
%         Hist_array5A_SEC(writeout_counter2,hist2_counter1)=...
%           Hist_array5A_SEC(writeout_counter2,hist2_counter1)+ 1;
%         %determine if the 'volumn' is different from HvsL
%         if log2(MvsL_values/HvsL_values)<=threshold_of_change &...
%             log2(MvsL_values/HvsL_values)>=-threshold_of_change
%           Hist_array5B_SEC(writeout_counter2,hist2_counter1)=...
%             Hist_array5B_SEC(writeout_counter2,hist2_counter1) + 1;
%         elseif log2(MvsL_values/HvsL_values)>threshold_of_change
%           Hist_array5C_SEC(writeout_counter2,hist2_counter1)=...
%             Hist_array5C_SEC(writeout_counter2,hist2_counter1) + 1;
%         elseif log2(MvsL_values/HvsL_values)<-threshold_of_change
%           Hist_array5D_SEC(writeout_counter2,hist2_counter1)=...
%             Hist_array5D_SEC(writeout_counter2,hist2_counter1) + 1;
%         end
%       end
%     end
%   end
%
%   Hist_array6A_SEC=zeros(Dimension_of_master_gaussian_list(1),(bin_number_hist2_SEC+1));
%   Hist_array6B_SEC=zeros(Dimension_of_master_gaussian_list(1),(bin_number_hist2_SEC+1));
%   Hist_array6C_SEC=zeros(Dimension_of_master_gaussian_list(1),(bin_number_hist2_SEC+1));
%   Hist_array6D_SEC=zeros(Dimension_of_master_gaussian_list(1),(bin_number_hist2_SEC+1));
%
%   %Count the number of gaussian detected vs the area
%   for hist2_counter1= 1:(bin_number_hist2_SEC)+1
%     for writeout_counter2= 1:Dimension_of_master_gaussian_list(1)
%       %determine bin values
%       Coverage_Bin_number=ceil(Finalised_Master_Gaussian_list.SEC_coverageHvsL(writeout_counter2,figure_counter)/0.05);
%       %determine if the fold change is within the bin value
%       if  hist2_counter1-1==Coverage_Bin_number
%         HvsL_values=Finalised_Master_Gaussian_list.Summed_Gaussian_areaHvsL(writeout_counter2,figure_counter);
%         MvsL_values=Finalised_Master_Gaussian_list.Summed_Gaussian_areaMvsL(writeout_counter2,figure_counter);
%         %determine bin of the total observed isotopologue 'volumn' in MvsL
%         Hist_array6A_SEC(writeout_counter2,hist2_counter1)=...
%           Hist_array6A_SEC(writeout_counter2,hist2_counter1)+ 1;
%         %determine if the 'volumn' is different from HvsL
%         if log2(HvsL_values/MvsL_values)<=threshold_of_change &...
%             log2(HvsL_values/MvsL_values)>=-threshold_of_change
%           Hist_array6B_SEC(writeout_counter2,hist2_counter1)=...
%             Hist_array6B_SEC(writeout_counter2,hist2_counter1) + 1;
%         elseif  log2(HvsL_values/MvsL_values)>threshold_of_change
%           Hist_array6C_SEC(writeout_counter2,hist2_counter1)=...
%             Hist_array6C_SEC(writeout_counter2,hist2_counter1) + 1;
%         elseif log2(HvsL_values/MvsL_values)<-threshold_of_change
%           Hist_array6D_SEC(writeout_counter2,hist2_counter1)=...
%             Hist_array6D_SEC(writeout_counter2,hist2_counter1) + 1;
%         end
%       end
%     end
%   end
%
%   %add value to a single array to grpah
%   Coverage_areaA(:,1)=sum(Hist_array5A_SEC)';
%   Coverage_areaA(:,2)=sum(Hist_array5B_SEC)';
%   Coverage_areaA(:,3)=sum(Hist_array5C_SEC)';
%   Coverage_areaA(:,4)=sum(Hist_array5D_SEC)';
%
%   Coverage_areaB(:,1)=sum(Hist_array6A_SEC)';
%   Coverage_areaB(:,2)=sum(Hist_array6B_SEC)';
%   Coverage_areaB(:,3)=sum(Hist_array6C_SEC)';
%   Coverage_areaB(:,4)=sum(Hist_array6D_SEC)';
%
%
%   %Compare coverage to observed guassian within replicates
%   %MvsL
%   Number_Gaussian_MvsL=zeros(Dimension_of_master_gaussian_list(1),1);
%
%   for writeout_counter2= 1:Dimension_of_master_gaussian_list(1)
%     %Find position of protein in MvsL array
%     [position_inMvsL]=ind2sub(length(MvsL_Gaussians.Protein_name),...
%       strmatch(Finalised_Master_Gaussian_list.Protein_name(writeout_counter2), MvsL_Gaussians.Protein_name, 'exact'));
%     %determine the length of position_inMvsL
%     number_observation=length(position_inMvsL);
%
%     %Find how many observations were recorded within each isotopologue channel
%     for counter=1:number_observation
%       if MvsL_Gaussians.Replicate(position_inMvsL(counter))==figure_counter
%         temp_value=MvsL_Gaussians.Unique_identifier(position_inMvsL(counter),:);
%         Number_Gaussian_MvsL(writeout_counter2)=sum(~cellfun('isempty',temp_value));
%       end
%     end
%   end
%
%   %HvsL
%   Number_Gaussian_HvsL=zeros(Dimension_of_master_gaussian_list(1),1);
%
%   for writeout_counter2= 1:Dimension_of_master_gaussian_list(1)
%     %Find position of protein in MvsL array
%     [position_inHvsL]=ind2sub(length(HvsL_Gaussians.Protein_name),...
%       strmatch(Finalised_Master_Gaussian_list.Protein_name(writeout_counter2), HvsL_Gaussians.Protein_name, 'exact'));
%     %determine the length of position_inMvsL
%     number_observation=length(position_inHvsL);
%
%     %Find how many observations were recorded within each isotopologue channel
%     for counter=1:number_observation
%       if HvsL_Gaussians.Replicate(position_inHvsL(counter))==figure_counter
%         temp_value=HvsL_Gaussians.Unique_identifier(position_inHvsL(counter),:);
%         Number_Gaussian_HvsL(writeout_counter2)=sum(~cellfun('isempty',temp_value));
%       end
%     end
%   end
%
%   %Count how many of each Guassian were observed for each bin
%   Hist_array7A=zeros(bin_number_hist2_MvsL,6);
%   Hist_array7B=zeros(bin_number_hist2_SEC,6);
%   Hist_array7C=zeros(bin_number_hist2_HvsL,6);
%   Hist_array7D=zeros(bin_number_hist2_SEC,6);
%
%   %Count the number of gaussian detected vs the area based on isotopologue area
%   for hist2_counter1= 1:(bin_number_hist2_MvsL)+1
%     for writeout_counter2= 1:Dimension_of_master_gaussian_list(1)
%       %determine bin values
%       Coverage_Bin_number=ceil(Finalised_Master_Gaussian_list.Summed_Gaussian_areaMvsL(writeout_counter2,figure_counter)/0.05);
%       %determine if the fold change is within the bin value
%       if  hist2_counter1-1==Coverage_Bin_number
%         %Create matrix couting how Gaussian were detected
%         Hist_array7A(hist2_counter1, (Number_Gaussian_MvsL(writeout_counter2)+1))=...
%           Hist_array7A(hist2_counter1, (Number_Gaussian_MvsL(writeout_counter2)+1))+1;
%       end
%     end
%   end
%
%   %Count the number of gaussian detected vs the area based on isotopologue area
%   for hist2_counter1= 1:(bin_number_hist2_HvsL)+1
%     for writeout_counter2= 1:Dimension_of_master_gaussian_list(1)
%       %determine bin values
%       Coverage_Bin_number=ceil(Finalised_Master_Gaussian_list.Summed_Gaussian_areaHvsL(writeout_counter2,figure_counter)/0.05);
%       %determine if the fold change is within the bin value
%       if  hist2_counter1-1==Coverage_Bin_number
%         %Create matrix couting how Gaussian were detected
%         Hist_array7C(hist2_counter1, (Number_Gaussian_HvsL(writeout_counter2)+1))=...
%           Hist_array7C(hist2_counter1, (Number_Gaussian_HvsL(writeout_counter2)+1))+1;
%       end
%     end
%   end
%
%   %Count the number of gaussian detected vs the area based on SEC coverage
%   for hist2_counter1= 1:(bin_number_hist2_SEC)+1
%     for writeout_counter2= 1:Dimension_of_master_gaussian_list(1)
%       %determine bin values
%       Coverage_Bin_number=ceil(Finalised_Master_Gaussian_list.SEC_coverageMvsL(writeout_counter2,figure_counter)/0.05);
%       %determine if the fold change is within the bin value
%       if  hist2_counter1-1==Coverage_Bin_number
%         %Create matrix couting how Gaussian were detected
%         Hist_array7B(hist2_counter1, (Number_Gaussian_MvsL(writeout_counter2)+1))=...
%           Hist_array7B(hist2_counter1, (Number_Gaussian_MvsL(writeout_counter2)+1))+1;
%       end
%     end
%   end
%
%   %Count the number of gaussian detected vs the area based on SEC coverage
%   for hist2_counter1= 1:(bin_number_hist2_SEC)+1
%     for writeout_counter2= 1:Dimension_of_master_gaussian_list(1)
%       %determine bin values
%       Coverage_Bin_number=ceil(Finalised_Master_Gaussian_list.SEC_coverageHvsL(writeout_counter2,figure_counter)/0.05);
%       %determine if the fold change is within the bin value
%       if  hist2_counter1-1==Coverage_Bin_number
%         %Create matrix couting how Gaussian were detected
%         Hist_array7D(hist2_counter1, (Number_Gaussian_HvsL(writeout_counter2)+1))=...
%           Hist_array7D(hist2_counter1, (Number_Gaussian_HvsL(writeout_counter2)+1))+1;
%       end
%     end
%   end
%
%   %write out figures of coverage
%   Text_for_figureMvsL1=strcat('Apex: ', mat2str(round(fitted_valuesMvsL.b1*100)/100));
%   Text_for_figureMvsL2=strcat('r-square: ', mat2str(round(fitted_statsMvsL.rsquare*100)/100));
%   Text_for_figureHvsL1=strcat('Apex1: ', mat2str(round(fitted_valuesHvsL.b1*100)/100));
%   Text_for_figureHvsL2=strcat('Apex2: ', mat2str(round(fitted_valuesHvsL.b2*100)/100));
%   Text_for_figureHvsL3=strcat('Apex3: ', mat2str(round(fitted_valuesHvsL.b3*100)/100));
%   Text_for_figureHvsL4=strcat('r-square: ', mat2str(round(fitted_statsHvsL.rsquare*100)/100));
%
%   %Write out histogram of complex coverage within sample
%   figure
%   f7=subplot(2,1,1);
%   f7_1_figure=bar(bin_number_hist_MvsL_plus1,Isotopologue_areaA(:,1), 0.6, 'stack');
%   set(f7_1_figure(1),'facecolor', myC(1,:), 'EdgeColor', 'k' );
%   hold on
%   plot(bin_number_hist_MvsL_plus1.',gaussian_fit_MvsL_bin,'linewidth', 2, 'Color', colour_to_use(3,:));
%   try
%     xlim([0,max(Finalised_Master_Gaussian_list.Summed_Gaussian_areaHvsL(:,figure_counter))*1.1])
%   catch
%   end
%   try
%     ylim([0,max(Isotopologue_areaA(:,1))*1.1]);
%   catch
%   end
%   title('Total observed protein PCP-SILAC coverage within the Lys4Arg6 sample','FontSize', 12);
%   ylabel('Number of proteins','FontSize', 8);
%   xlabel('Coverage (expected area under the curve vs observed)','FontSize', 8);
%   try
%     text(max(Finalised_Master_Gaussian_list.Summed_Gaussian_areaHvsL(:,figure_counter))*0.9,max(Isotopologue_areaA(:,1)), Text_for_figureMvsL1, 'FontSize', 8);
%     text(max(Finalised_Master_Gaussian_list.Summed_Gaussian_areaHvsL(:,figure_counter))*0.9,max(Isotopologue_areaA(:,1))*0.90, Text_for_figureMvsL2, 'FontSize', 8);
%   catch
%   end
%
%   f7=subplot(2,1,2);
%   f7_2_figure=bar(bin_number_hist_HvsL_plus1,Isotopologue_areaB(:,1) , 0.6, 'stack');
%   set(f7_2_figure(1),'facecolor', myC(1,:), 'EdgeColor', 'k' );
%   hold on
%   fit_fig1=plot(bin_number_hist_HvsL_plus1.',gaussian_fit_HvsL_bin1,'linewidth', 2, 'Color', colour_to_use(6,:));
%   hold on
%   fit_fig2=plot(bin_number_hist_HvsL_plus1.',gaussian_fit_HvsL_bin2,'linewidth', 2,'Color', colour_to_use(3,:));
%   hold on
%   fit_fig3=plot(bin_number_hist_HvsL_plus1.',gaussian_fit_HvsL_bin3,'linewidth', 2,'Color', colour_to_use(5,:));
%   try
%     xlim([0,max(Finalised_Master_Gaussian_list.Summed_Gaussian_areaHvsL(:,figure_counter))*1.1])
%   catch
%   end
%   try
%     ylim([0,max(Isotopologue_areaB(:,1))*1.1]);
%   catch
%   end
%   title('Total observed protein PCP-SILAC coverage within the Lys8Arg10 samples','FontSize', 12);
%   ylabel('Number of proteins','FontSize', 8);
%   xlabel('Coverage (expected area under the curve vs observed)','FontSize', 8);
%   try
%     text(max(Finalised_Master_Gaussian_list.Summed_Gaussian_areaHvsL(:,figure_counter))*0.9,max(Isotopologue_areaB(:,1)), Text_for_figureHvsL1, 'FontSize', 8);
%     text(max(Finalised_Master_Gaussian_list.Summed_Gaussian_areaHvsL(:,figure_counter))*0.9,max(Isotopologue_areaB(:,1))*0.90, Text_for_figureHvsL2, 'FontSize', 8);
%     text(max(Finalised_Master_Gaussian_list.Summed_Gaussian_areaHvsL(:,figure_counter))*0.9,max(Isotopologue_areaB(:,1))*0.80, Text_for_figureHvsL3, 'FontSize', 8);
%     text(max(Finalised_Master_Gaussian_list.Summed_Gaussian_areaHvsL(:,figure_counter))*0.9,max(Isotopologue_areaB(:,1))*0.70, Text_for_figureHvsL4, 'FontSize', 8);
%   catch
%   end
%
%   % save Histogram
%   Save_name_plot=[figdir 'Comparison/Isotopologue_coverage_Histogram_replicate',mat2str(figure_counter),'.png'];
%   saveas(f7, Save_name_plot);
%
%   figure
%   f8=subplot(2,1,1);
%   f8_1_figure=bar(bin_number_hist_SEC_plus1,Coverage_areaA(:,1), 0.6, 'stack');
%   set(f8_1_figure(1),'facecolor', myC(1,:), 'EdgeColor', 'k' );
%   xlim([0,1]);
%   ylim([0,max(Coverage_areaA(:,1))*1.1]);
%   title('Total observed protein coverage across SEC fractions within the Lys4Arg6 sample','FontSize', 12);
%   ylabel('Number of proteins','FontSize', 8);
%   xlabel('Coverage (% of fracations analyzed)','FontSize', 8);
%   % Convert y-axis values to percentage values by multiplication
%   a=[cellstr(num2str(get(gca,'xtick')'*100))];
%   % Create a vector of '%' signs
%   pct = char(ones(size(a,1),1)*'%');
%   % Append the '%' signs after the percentage values
%   new_xticks = [char(a),pct];
%   % 'Reflect the changes on the plot
%   set(gca,'xticklabel',new_xticks);
%
%   f8=subplot(2,1,2);
%   f8_2_figure=bar(bin_number_hist_SEC_plus1,Coverage_areaB(:,1), 0.6, 'stack');
%   set(f8_2_figure(1),'facecolor', myC(1,:), 'EdgeColor', 'k' );
%   xlim([0,1]);
%   ylim([0,max(Coverage_areaA(:,1))*1.1]);
%   title('Total observed protein coverage across SEC fractions within the Lys8Arg10 sample','FontSize', 12);
%   ylabel('Number of proteins','FontSize', 8);
%   xlabel('Coverage (% of fracations analyzed)','FontSize', 8);
%   % Convert y-axis values to percentage values by multiplication
%   a=[cellstr(num2str(get(gca,'xtick')'*100))];
%   % Create a vector of '%' signs
%   pct = char(ones(size(a,1),1)*'%');
%   % Append the '%' signs after the percentage values
%   new_xticks = [char(a),pct];
%   % 'Reflect the changes on the plot
%   set(gca,'xticklabel',new_xticks);
%
%   % save Histogram
%   Save_name_plot=[figdir 'Comparison/SEC_fraction_coverage_Histogram_replicate',mat2str(figure_counter),'.png'];
%   saveas(f8, Save_name_plot);
%
%   %Create figure of observed changes across Isotopologue_coverage
%   figure
%   f9=subplot(2,1,1);
%   f9_1_figure=bar(bin_number_hist_MvsL_plus1,[Isotopologue_areaA(:,4),Isotopologue_areaA(:,2),Isotopologue_areaA(:,3)] , 0.6, 'stack');
%   order_of_columns=[1 3 6];
%   for k=1:3
%     set(f9_1_figure(k),'facecolor', colour_to_use(order_of_columns(k),:), 'EdgeColor', 'k' );
%   end
%   legend('Decrease','No change','Increase');
%   xlim([0,max(Finalised_Master_Gaussian_list.Summed_Gaussian_areaHvsL(:,figure_counter))*1.1])
%   title('Changes (95%), PCP-SILAC coverage within experiments (Lys4Arg6 sample)','FontSize', 12);
%   ylabel('Number of proteins','FontSize', 8);
%   xlabel('Coverage (expected area under the curve vs observed)','FontSize', 8);
%
%   f9=subplot(2,1,2)
%   f9_2_figure=bar(bin_number_hist_SEC_plus1,[Coverage_areaA(:,4),Coverage_areaA(:,2),Coverage_areaA(:,3)], 0.6, 'stack');
%   for k=1:3
%     set(f9_2_figure(k),'facecolor', colour_to_use(order_of_columns(k),:), 'EdgeColor', 'k' );
%   end
%   xlim([0,1]);
%   ylim([0,max(Coverage_areaA(:,1))*1.1]);
%   title('Changes (95%), coverage across SEC fractions within experiments (Lys4Arg6 sample)','FontSize', 12);
%   ylabel('Number of proteins','FontSize', 8);
%   xlabel('Coverage (% of fracations analyzed)','FontSize', 8);
%   % Convert y-axis values to percentage values by multiplication
%   a=[cellstr(num2str(get(gca,'xtick')'*100))];
%   % Create a vector of '%' signs
%   pct = char(ones(size(a,1),1)*'%');
%   % Append the '%' signs after the percentage values
%   new_xticks = [char(a),pct];
%   % 'Reflect the changes on the plot
%   set(gca,'xticklabel',new_xticks);
%
%   % save Histogram
%   Save_name_plot=[figdir 'Comparison/Changes_Histogram_Lys4Arg6_replicate',mat2str(figure_counter),'.png'];
%   saveas(f9, Save_name_plot);
%
%   %Create figure of observed changes across Isotopologue_coverage
%   figure
%   f9B=subplot(2,1,1);
%   f9B_1_figure=bar(bin_number_hist_HvsL_plus1,[Isotopologue_areaB(:,4),Isotopologue_areaB(:,2),Isotopologue_areaB(:,3)] , 0.6, 'stack');
%   order_of_columns=[1 3 6];
%   for k=1:3
%     set(f9B_1_figure(k),'facecolor', colour_to_use(order_of_columns(k),:), 'EdgeColor', 'k' );
%   end
%   legend('Decrease','No change','Increase');
%   xlim([0,max(Finalised_Master_Gaussian_list.Summed_Gaussian_areaHvsL(:,figure_counter))*1.1])
%   title('Changes (95%), PCP-SILAC coverage within experiments (Lys8Arg10 sample)','FontSize', 12);
%   ylabel('Number of proteins','FontSize', 8);
%   xlabel('Coverage (expected area under the curve vs observed)','FontSize', 8);
%
%   f9B=subplot(2,1,2);
%   f9B_2_figure=bar(bin_number_hist_SEC_plus1,[Coverage_areaB(:,4),Coverage_areaB(:,2),Coverage_areaB(:,3)], 0.6, 'stack');
%   for k=1:3
%     set(f9B_2_figure(k),'facecolor', colour_to_use(order_of_columns(k),:), 'EdgeColor', 'k' );
%   end
%   xlim([0,1]);
%   ylim([0,max(Coverage_areaB(:,1))*1.1]);
%   title('Changes (95%), coverage across SEC fractions within experiments (Lys8Arg10 sample)','FontSize', 12);
%   ylabel('Number of proteins','FontSize', 8);
%   xlabel('Coverage (% of fracations analyzed)','FontSize', 8);
%   % Convert y-axis values to percentage values by multiplication
%   a=[cellstr(num2str(get(gca,'xtick')'*100))];
%   % Create a vector of '%' signs
%   pct = char(ones(size(a,1),1)*'%');
%   % Append the '%' signs after the percentage values
%   new_xticks = [char(a),pct];
%   % 'Reflect the changes on the plot
%   set(gca,'xticklabel',new_xticks);
%
%   % save Histogram
%   Save_name_plot=[figdir 'Comparison/Changes_Histogram_Lys8Arg10_replicate ',mat2str(figure_counter),'.png'];
%   saveas(f9B, Save_name_plot);
%
%   figure
%   f10=subplot(2,1,1)
%   f10_1_figure=bar(bin_number_hist_MvsL,[Hist_array7A(:,1),Hist_array7A(:,2),Hist_array7A(:,3),...
%     Hist_array7A(:,4),Hist_array7A(:,5),Hist_array7A(:,6)], 0.6, 'stack');
%   order_of_columns=[1 3 6];
%   for k=1:6
%     set(f10_1_figure(k),'facecolor', colour_to_use((k),:), 'EdgeColor', 'k' );
%   end
%   legend('zero','one','two','three','four','five','FontSize',8, 'Location', 'EastOutside');
%   xlim([0,max(Finalised_Master_Gaussian_list.Summed_Gaussian_areaHvsL(:,figure_counter))*1.1])
%   title('Number of fitted Gaussian compared to PCP-SILAC coverage','FontSize', 12);
%   ylabel('Number of proteins','FontSize', 8);
%   xlabel('Coverage (expected area under the curve vs observed)','FontSize', 8);
%
%   f10=subplot(2,1,2);
%   f10_2_figure=bar(bin_number_hist_SEC,[Hist_array7B(:,1),Hist_array7B(:,2),Hist_array7B(:,3),...
%     Hist_array7B(:,4),Hist_array7B(:,5),Hist_array7B(:,6)], 0.6, 'stack');
%   for k=1:6
%     set(f10_2_figure(k),'facecolor', colour_to_use((k),:), 'EdgeColor', 'k' );
%   end
%   xlim([0,1]);
%   ylim([0,max(Coverage_areaA(:,1))*1.1]);
%   legend('zero','one','two','three','four','five','FontSize',8, 'Location', 'EastOutside');
%   title('Number of fitted Gaussian compared to SEC fractions coverage','FontSize', 12);
%   ylabel('Number of proteins','FontSize', 8);
%   xlabel('Coverage (% of fracations analyzed)','FontSize', 8);
%   % Convert y-axis values to percentage values by multiplication
%   a=[cellstr(num2str(get(gca,'xtick')'*100))];
%   % Create a vector of '%' signs
%   pct = char(ones(size(a,1),1)*'%');
%   % Append the '%' signs after the percentage values
%   new_xticks = [char(a),pct];
%   % 'Reflect the changes on the plot
%   set(gca,'xticklabel',new_xticks);
%
%   % save Histogram
%   Save_name_plot=[figdir 'Comparison/Gaussians_Histogram_Lys4Arg6_replicate ',mat2str(figure_counter),'.png'];
%   saveas(f10, Save_name_plot);
%
%   figure
%   f10B=subplot(2,1,1);
%   f10B_1_figure=bar(bin_number_hist_HvsL,[Hist_array7C(:,1),Hist_array7C(:,2),Hist_array7C(:,3),...
%     Hist_array7C(:,4),Hist_array7C(:,5),Hist_array7C(:,6)], 0.6, 'stack');
%   order_of_columns=[1 3 6];
%   for k=1:6
%     set(f10B_1_figure(k),'facecolor', colour_to_use((k),:), 'EdgeColor', 'k' );
%   end
%   legend('zero','one','two','three','four','five','FontSize',8, 'Location', 'EastOutside');
%   xlim([0,max(Finalised_Master_Gaussian_list.Summed_Gaussian_areaHvsL(:,figure_counter))*1.1])
%   title('Number of fitted Gaussian compared to PCP-SILAC coverage','FontSize', 12);
%   ylabel('Number of proteins','FontSize', 8);
%   xlabel('Coverage (expected area under the curve vs observed)','FontSize', 8);
%
%   f10B=subplot(2,1,2);
%   f10B_2_figure=bar(bin_number_hist_SEC,[Hist_array7D(:,1),Hist_array7D(:,2),Hist_array7D(:,3),...
%     Hist_array7D(:,4),Hist_array7D(:,5),Hist_array7D(:,6)], 0.6, 'stack');
%   for k=1:6
%     set(f10B_2_figure(k),'facecolor', colour_to_use((k),:), 'EdgeColor', 'k' );
%   end
%   xlim([0,1]);
%   ylim([0,max(Coverage_areaB(:,1))*1.1]);
%   legend('zero','one','two','three','four','five','FontSize',8, 'Location', 'EastOutside');
%   title('Number of fitted Gaussian compared to SEC fractions coverage','FontSize', 12);
%   ylabel('Number of proteins','FontSize', 8);
%   xlabel('Coverage (% of fracations analyzed)','FontSize', 8);
%   % Convert y-axis values to percentage values by multiplication
%   a=[cellstr(num2str(get(gca,'xtick')'*100))];
%   % Create a vector of '%' signs
%   pct = char(ones(size(a,1),1)*'%');
%   % Append the '%' signs after the percentage values
%   new_xticks = [char(a),pct];
%   % 'Reflect the changes on the plot
%   set(gca,'xticklabel',new_xticks);
%
%   % save Histogram
%   Save_name_plot=[figdir 'Comparison/Gaussians_Histogram_Lys8Arg10_replicate ',mat2str(figure_counter),'.png'];
%   saveas(f10B, Save_name_plot);
%
%   figure
%   f11=subplot(2,1,1);
%   f11_1_figure=bar(bin_number_hist_HvsL_plus1,Isotopologue_areaB(:,1), 0.6, 'stack');
%   set(f11_1_figure(1),'facecolor', myC(1,:), 'EdgeColor', 'k' );
%   hold on
%   fit_fig1=plot(bin_number_hist_HvsL_plus1.',gaussian_fit_HvsL_bin1,'linewidth', 2, 'Color', colour_to_use(6,:));
%   hold on
%   fit_fig2=plot(bin_number_hist_HvsL_plus1.',gaussian_fit_HvsL_bin2,'linewidth', 2,'Color', colour_to_use(3,:));
%   hold on
%   fit_fig3=plot(bin_number_hist_HvsL_plus1.',gaussian_fit_HvsL_bin3,'linewidth', 2,'Color', colour_to_use(5,:));
%   try
%     xlim([0,max(Finalised_Master_Gaussian_list.Summed_Gaussian_areaHvsL(:,figure_counter))*1.1])
%   catch
%   end
%   try
%     ylim([0,max(Isotopologue_areaB(:,1))*1.1]);
%   catch
%   end
%   title('Total observed protein PCP-SILAC coverage within the Lys8Arg10 samples','FontSize', 12);
%   ylabel('Number of proteins','FontSize', 8);
%   xlabel('Coverage (expected area under the curve vs observed)','FontSize', 8);
%   try
%     text(max(Finalised_Master_Gaussian_list.Summed_Gaussian_areaHvsL(:,figure_counter))*0.9,max(Isotopologue_areaB(:,1)), Text_for_figureHvsL1, 'FontSize', 8);
%     text(max(Finalised_Master_Gaussian_list.Summed_Gaussian_areaHvsL(:,figure_counter))*0.9,max(Isotopologue_areaB(:,1))*0.90, Text_for_figureHvsL2, 'FontSize', 8);
%     text(max(Finalised_Master_Gaussian_list.Summed_Gaussian_areaHvsL(:,figure_counter))*0.9,max(Isotopologue_areaB(:,1))*0.80, Text_for_figureHvsL3, 'FontSize', 8);
%     text(max(Finalised_Master_Gaussian_list.Summed_Gaussian_areaHvsL(:,figure_counter))*0.9,max(Isotopologue_areaB(:,1))*0.70, Text_for_figureHvsL4, 'FontSize', 8);
%   catch
%   end
%
%   f11=subplot(2,1,2)
%   f11_1_figure=bar(bin_number_hist_HvsL_plus1,[Isotopologue_areaB(:,4),Isotopologue_areaB(:,2),Isotopologue_areaB(:,3)] , 0.6, 'stack');
%   order_of_columns=[1 3 6];
%   for k=1:3
%     set(f11_1_figure(k),'facecolor', colour_to_use(order_of_columns(k),:), 'EdgeColor', 'k' );
%   end
%   legend('Decrease','No change','Increase');
%   xlim([0,max(Finalised_Master_Gaussian_list.Summed_Gaussian_areaHvsL(:,figure_counter))*1.1])
%   title('Changes (95%), PCP-SILAC coverage within experiments (Lys8Arg10 sample)','FontSize', 12);
%   ylabel('Number of proteins','FontSize', 8);
%   xlabel('Coverage (expected area under the curve vs observed)','FontSize', 8);
%
%   %save Histogram
%   Save_name_plot=[figdir 'Comparison/Comparsion_Histogram_Lys8Arg10_replicate ',mat2str(figure_counter),'.png'];
%   saveas(f11, Save_name_plot);
% end
%
% %%Count how many guassian change in each fraction
% Hist_array4_MvsL=zeros(fraction_to_plot,3);
% Hist_array4_HvsL=zeros(fraction_to_plot,3);
%
% for fraction_counter1= 1:Dimension_of_master_gaussian_list(1)
%   for fraction_counter2= 1:nnz(Finalised_Master_Gaussian_list.Center(fraction_counter1,:))
%     %Find Center
%     Center_of_gaussain= floor(Finalised_Master_Gaussian_list.Center(fraction_counter1,fraction_counter2));
%     %Determine if Center is greater then zero
%     if  Center_of_gaussain>0
%
%       if Finalised_Master_Gaussian_list.Average_fold_change(fraction_counter1,fraction_counter2) >0 &&...
%           Center_of_gaussain <=fraction_to_plot && isequal(Finalised_Master_Gaussian_list.Satifies_Adjusted_pvalue(fraction_counter1,fraction_counter2),{'+'})
%
%         Hist_array4_MvsL(Center_of_gaussain,1)=Hist_array4_MvsL(Center_of_gaussain,1)+1;
%         Hist_array4_MvsL(Center_of_gaussain,2)=Hist_array4_MvsL(Center_of_gaussain,2)+1;
%       elseif Finalised_Master_Gaussian_list.Average_fold_change(fraction_counter1,fraction_counter2) <0 &&...
%           Center_of_gaussain <=fraction_to_plot && isequal(Finalised_Master_Gaussian_list.Satifies_Adjusted_pvalue(fraction_counter1,fraction_counter2),{'+'})
%
%         Hist_array4_MvsL(Center_of_gaussain,1)=Hist_array4_MvsL(Center_of_gaussain,1)+1;
%         Hist_array4_MvsL(Center_of_gaussain,3)=Hist_array4_MvsL(Center_of_gaussain,3)+1;
%       elseif Center_of_gaussain <=fraction_to_plot
%
%         Hist_array4_MvsL(Center_of_gaussain,1)=Hist_array4_MvsL(Center_of_gaussain,1)+1;
%       end
%     end
%   end
% end

%Create array to plot
plotting_counter=1;

for writeout_counter1= 1:length(Unique_protein_names)
  for writeout_counter2= 1:nnz(Finalised_Master_Gaussian_list.Center(writeout_counter1,:))
    data_scatter_plot(plotting_counter,1)=Finalised_Master_Gaussian_list.foldChange(writeout_counter1,writeout_counter2);
    data_scatter_plot(plotting_counter,2)= 1;%Finalised_Master_Gaussian_list.Stdev_fold_change(writeout_counter1,writeout_counter2);
    if isequal(Finalised_Master_Gaussian_list.Satifies_Adjusted_pvalue(writeout_counter1,writeout_counter2),{'+'})
      data_scatter_plot(plotting_counter,3)=1;
    else
      data_scatter_plot(plotting_counter,3)=0;
    end
    plotting_counter=plotting_counter+1;
  end
end

%Sort arrya for plot and divided
[unique_log2_sorted1,unique_log2_sorted2]=sort(data_scatter_plot(:,1));
data_scatter_plot=[data_scatter_plot(unique_log2_sorted2,1) data_scatter_plot(unique_log2_sorted2,2) data_scatter_plot(unique_log2_sorted2,3)];

%Create array
length_values_plot=length(data_scatter_plot);
Multiple_measurement_Adjusted_pvalue_corrected=zeros(length_values_plot,2);

%Work out mean stdev
work_out_mean=data_scatter_plot(:,2);
work_out_mean(work_out_mean==0)=[];
mean_stdev=mean(work_out_mean);

%Divide dataset to colour code individual changes
Gaussian_ploted=0;
Gaussian_not_ploted=1;
for counter= 1:length(data_scatter_plot)
  if data_scatter_plot(counter,3) == 1 && (data_scatter_plot(counter,1) >0.2 | data_scatter_plot(counter,1) <-0.2 )
    Multiple_measurement_Adjusted_pvalue_corrected(counter,1)=data_scatter_plot(counter,1);
    Multiple_measurement_Adjusted_pvalue_corrected(counter,2)=data_scatter_plot(counter,2);
    Gaussian_ploted=1+Gaussian_ploted;
  elseif ~(data_scatter_plot(counter,3) == 1 && (data_scatter_plot(counter,1) >0.2 | data_scatter_plot(counter,1) <-0.2 ))
    Multiple_measurement_Adjusted_pvalue_corrected(counter,1)=NaN;
    Multiple_measurement_Adjusted_pvalue_corrected(counter,2)=NaN;
    Gaussian_not_ploted=1+Gaussian_not_ploted;
  end
end

%Text for figure
Text_for_figure1=strcat('Mean std dev: ', mat2str((round(mean_stdev*1000))/1000));
Text_for_figure2=strcat('Number of Guassian below bonferroni adjusted p-value: ', mat2str(Gaussian_ploted));

%Graph data as log2 scatter
figure;
f_all=subplot(2,1,1); hold on
P4A =scatter([1:length_values_plot],data_scatter_plot(:,1),5,'fill', 'markerfacecolor',colour_to_use(1,:));
P4B =scatter([1:length_values_plot],Multiple_measurement_Adjusted_pvalue_corrected(:,1), 5,'fill', 'markerfacecolor',colour_to_use(6,:));
title('Log2 changes in Gaussians: Cyotplasmic- Fas treated','FontSize', 12);
ylabel('Log2 ratio','FontSize', 10);
xlabel('Gaussians identified','FontSize', 10);
%ylim([data_scatter_plot(1,1)*1.05,data_scatter_plot(end,1)*1.05]);
xlim([-100,length_values_plot+100]);
plot([-100,length_values_plot+100],[1 1],':');
plot([-100,length_values_plot+100],[-1 -1],':');
text(100, 3.5, Text_for_figure1, 'FontSize', 8);
text(100, 3.0, Text_for_figure2, 'FontSize', 8);

%plot distribution of protein chnages across the fractions
%define colours
subplot(2,1,2);
%f_all_figure=bar(1:fraction_to_plot, [Hist_array4_MvsL(:,1) Hist_array4_MvsL(:,2) Hist_array4_MvsL(:,3)],0.6, 'stack');
f_all_figure=bar(1:Nfraction, [Hist_array(:,1) Hist_array(:,2) Hist_array(:,3)],0.6, 'stack');
for k=1:3
  set(f_all_figure(k),'facecolor', myC(k,:), 'EdgeColor', 'k' )
end
legend('No change', 'Increase', 'Decrease', 'FontSize',8, 'Location', 'Best');
xlim([-2,frac2+2]);
title('Changes in Gaussians observed across fractions');
xlabel('Fractions');
ylabel('# of Guassian centers');
saveas(f_all, [figdir 'Comparison/Unique_gaussian_changes_observed.png']);

%Create figures for all proteins showing overlay of
%Create a file to write pdf name out to
List_of_pdf=cell(length(Unique_protein_names)*2,2);
List_of_pdf_counter=1;

%%

if user.fastcomparison==0
  
  %Determine which Gaussians change by comparing fitted data with raw data
  for Gaussian_counter1 = 1:length(Unique_protein_names)
    
    %Protein being plotted
    Protein_to_plot=Finalised_Master_Gaussian_list.Protein_name{Gaussian_counter1};
    
    % Find this protein in the raw data
    Iraw = find(strcmp(Protein_to_plot,txt_val{1}(:,2)));
    
    
    if isempty(Iraw)
      disp([Protein_to_plot ' not found (not sure why!)'])
      continue
    end
    
    % Get the data to plot
    data2plot = zeros(Nchannels,replicate_num,length(frac1:frac2));
    for jj = 1:replicate_num
      for ii = 1:Nchannels
        tmp = num_val{ii}(Iraw(jj),frac1:frac2);
        tmp(isnan(tmp)) = 0;
        data2plot(ii,jj,:) = tmp;
      end
    end
    
    figure
    
    % Plot the Gaussians
    subplot(replicate_num+2,1,1), hold on
    %Count the number of gaussians detected
    number_of_gaussian_to_plot= nnz(Finalised_Master_Gaussian_list.Center(Gaussian_counter1,:));
    for hold_on_counter=1:number_of_gaussian_to_plot
      Center_output=Finalised_Master_Gaussian_list.Center(Gaussian_counter1, hold_on_counter);
      Height_output=Finalised_Master_Gaussian_list.Height(Gaussian_counter1, hold_on_counter);
      Width_output=Finalised_Master_Gaussian_list.Width(Gaussian_counter1, hold_on_counter);
      Fitted_gaus=1:0.1:frac2;
      number_of_data_points=length(Fitted_gaus);
      gaus_colour_set=colour_to_use(hold_on_counter,:);
      %     for fill_in_counter=1:number_of_data_points
      %       hold on  %Graph Gaus fill area
      %       y1 =  Height_output*exp(-((Fitted_gaus(fill_in_counter)-Center_output)/Width_output).^2);
      %       patch([Fitted_gaus(fill_in_counter) Fitted_gaus(fill_in_counter)], [0 y1], 'w','EdgeColor',gaus_colour_set(:),'EdgeAlpha',0.2,'LineWidth',2);
      %     end
      y1 =  Height_output*exp(-((Fitted_gaus-Center_output)/Width_output).^2);
      patch(Fitted_gaus([1 1:end end]), [0 y1 0], gaus_colour_set(:)','EdgeColor',gaus_colour_set(:),'FaceAlpha',0.2,'LineWidth',2);
      plot(Fitted_gaus,y1,'Color','black','LineWidth',1);
    end
    y = ylim;
    s = ['R^2 = ' num2str(Finalised_Master_Gaussian_list.adjrsquare(Gaussian_counter1,1))];
    text(1,y(1) + diff(y)*1.02,s)
    ylim([y(1) y(2)*1.1])
    xlim([0,frac2+1])
    title_name_plot=strcat('Gaussian curves identified across replicates of :',Protein_to_plot);
    title(title_name_plot,'FontSize', 12);
    ylabel('Isotopologue ratio','FontSize', 10);
    xlabel('Fractions','FontSize', 10);
    
    % Plot the raw data
    for jj = 1:replicate_num
      subplot(replicate_num+2,1,jj+1),hold on
      for ii = 1:Nchannels
        plot((frac1:frac2)-frac1+1,squeeze(data2plot(ii,jj,:)),'s',...
          'MarkerFaceColor',colour_to_use(2,:),'MarkerSize',4);
      end
    end
    legend(user.silacratios,'location','best')
    for jj = 1:replicate_num
      subplot(replicate_num+2,1,jj+1),hold on
      for ii = 1:Nchannels
        plot((frac1:frac2)-frac1+1,squeeze(data2plot(ii,jj,:)),...
          'Color', [0.5 0.5 0.5],'LineWidth',0.5);
      end
    end
    xlim([0,frac2+1])
    
    % Plot the relative changes
    subplot(replicate_num+2,1,replicate_num+2),hold on
    Replicate_normalised_raw_data = zeros(replicate_num,number_of_gaussian_to_plot);
    for ii = 1:replicate_num
      Replicate_normalised_raw_data(ii,:) = Finalised_Master_Gaussian_list.foldChange_byreplicate(Gaussian_counter1,1:number_of_gaussian_to_plot,ii);
    end
    Names_for_bar_graph = cell(1*number_of_gaussian_to_plot,1);
    for jj = 1:number_of_gaussian_to_plot
      %Values_for_bar_graph(ii,1) = Replicate_normalised_raw_data(ii,index_minimum);
      Names_for_bar_graph{jj} = ['G_' mat2str(jj)];% '_R_' num2str(ii)];
    end
    xbar = linspace(-1/4,1/4,replicate_num);
    if replicate_num==1
      xbar = 0;
    end
    for ii = 1:number_of_gaussian_to_plot
      b1 = bar(ii+xbar,Replicate_normalised_raw_data(:,ii));
      set(b1,'facecolor',colour_to_use(ii,:))
    end
    y = ylim;
    y2(1) = min([-1 y(1)]);
    y2(2) = max([1 y(2)]);
    for ii =1:number_of_gaussian_to_plot-1
      plot([1 1]*(ii+0.5),y2,'--r')
    end
    title_name_bar=strcat('Log2 ratio of treated to untreated at gaussian apex of :',Protein_to_plot);
    plot([0 replicate_num+1],[1 1],':k','LineWidth',1);
    plot([0 replicate_num+1],[-1 -1],':k','LineWidth',1);
    title(title_name_bar,'FontSize', 12);
    ylabel('Log2 ratio','FontSize', 10);
    %xlabel('Fractions','FontSize', 10);
    set(gca, 'XTickLabel',Names_for_bar_graph, 'XTick',1:numel(Names_for_bar_graph),'FontSize', 12);
    xlim([0.5 number_of_gaussian_to_plot+0.5]);
    
    Save_name_replicates=[figdir 'Comparison/ProteinGaussianMaps/' mat2str(Gaussian_counter1),'_1_PCP_SEC_Profiles of_',Protein_to_plot,'.png'];
    %List_of_pdf{List_of_pdf_counter,1}=Save_name_replicates;
    %List_of_pdf{List_of_pdf_counter,2}=List_of_pdf_counter;
    %List_of_pdf_counter=List_of_pdf_counter+1;
    saveas(gcf, Save_name_replicates);
    %print('-dpdf', '-r600', Save_name_replicates);
    close 'all';
    
  end
  
end


%%
%if user.fastcomparison == 0
if 0
  
  for Gaussian_counter1= 1:length(Unique_protein_names)
    % Determine if Gaussians should be plotted
    
    %Count the number of gaussians detected
    number_of_gaussian_to_plot= nnz(Finalised_Master_Gaussian_list.Center(Gaussian_counter1,:));
    
    %create counter
    hold_on_counter1=1;
    
    %rest varibles
    Center_test=[];
    Height_test=[];
    Width_test=[];
    
    %determine if Center is less then fraction_to_plot
    for hold_on_counter=1:number_of_gaussian_to_plot
      if ~(Finalised_Master_Gaussian_list.Center(Gaussian_counter1, hold_on_counter)>frac2-2) % minus two add for comsetics
        Center_test(hold_on_counter1)=Finalised_Master_Gaussian_list.Center(Gaussian_counter1, hold_on_counter);
        Height_test(hold_on_counter1)=Finalised_Master_Gaussian_list.Height(Gaussian_counter1, hold_on_counter);
        Width_test(hold_on_counter1)=Finalised_Master_Gaussian_list.Width(Gaussian_counter1, hold_on_counter);
        hold_on_counter1=hold_on_counter1+1;
      end
    end
    
    number_of_gaussian_to_plot=length(Center_test);
    
    if ~isempty(Center_test)
      %Plot Quantiation_of proteins
      figure
      f9 = subplot(2,1,1);
      for hold_on_counter=1:number_of_gaussian_to_plot
        hold on  %Graph Gaus
        Center_output=Center_test(hold_on_counter);
        Height_output=Height_test(hold_on_counter);
        Width_output=Width_test(hold_on_counter);
        Fitted_gaus=1:0.1:frac2;
        number_of_data_points=length(Fitted_gaus);
        gaus_colour_set=colour_to_use(hold_on_counter,:);
        for fill_in_counter=1:number_of_data_points
          hold on  %Graph Gaus fill area
          y1 =  Height_output*exp(-((Fitted_gaus(fill_in_counter)-Center_output)/Width_output).^2);
          patch([Fitted_gaus(fill_in_counter) Fitted_gaus(fill_in_counter)], [0 y1], 'w','EdgeColor',gaus_colour_set(:),'EdgeAlpha',0.2,'LineWidth',2);
        end
        y1 =  Height_output*exp(-((Fitted_gaus-Center_output)/Width_output).^2);
        P1 = plot(Fitted_gaus,y1);
        set(P1,'Color','black','LineWidth',1);
        xlim([0,frac2]);
      end
      title_name_plot=strcat('Gaussian curves identified across replicates of :',Protein_to_plot);
      title(title_name_plot,'FontSize', 12);
      ylabel('Isotopologue ratio','FontSize', 10);
      xlabel('Fractions','FontSize', 10);
      
      
      f9 = subplot(2,1,2);hold on
      %Format log2 values into ascending order to plot
      %Values_for_bar_graph=zeros(number_of_gaussian_to_plot*replicate_num,1);
      %Names_for_bar_graph=cell(number_of_gaussian_to_plot*replicate_num,1);
      %bar_position=[1:replicate_num:(number_of_gaussian_to_plot*replicate_num)];
      
      Center_to_plot =Center_test(:);
      Center_counter1=1;
      Gaussian_number_counter1=1;
      
      %Copy value to matrix to manipulate
      %     Replicate1_normalised_raw_data=Finalised_Master_Gaussian_list.normaliased_log2_comparsion_replicate_1(Gaussian_counter1,1:number_of_gaussian_to_plot);
      %     Replicate2_normalised_raw_data=Finalised_Master_Gaussian_list.normaliased_log2_comparsion_replicate_2(Gaussian_counter1,1:number_of_gaussian_to_plot);
      %     Replicate3_normalised_raw_data=Finalised_Master_Gaussian_list.normaliased_log2_comparsion_replicate_3(Gaussian_counter1,1:number_of_gaussian_to_plot);
      %     Replicate4_normalised_raw_data=Finalised_Master_Gaussian_list.normaliased_log2_comparsion_replicate_4(Gaussian_counter1,1:number_of_gaussian_to_plot);
      
      
      %  while ~isempty(Center_to_plot)
      %Find minimum center to plot
      %[~, index_minimum]=min(Center_to_plot);
      Replicate_normalised_raw_data = zeros(replicate_num,number_of_gaussian_to_plot);
      for ii = 1:replicate_num
        Replicate_normalised_raw_data(ii,:) = Finalised_Master_Gaussian_list.foldChange_byreplicate(Gaussian_counter1,1:number_of_gaussian_to_plot,ii);
      end
      Names_for_bar_graph = cell(1*number_of_gaussian_to_plot,1);
      for jj = 1:number_of_gaussian_to_plot
        %Values_for_bar_graph(ii,1) = Replicate_normalised_raw_data(ii,index_minimum);
        Names_for_bar_graph{jj} = ['G_' mat2str(jj)];% '_R_' num2str(ii)];
      end
      
      %       if replicate_num >=2
      %         %replicate 2
      %         Values_for_bar_graph(bar_position(Center_counter1)+1,1)=Replicate2_normalised_raw_data(index_minimum);
      %         Bar_bin_name_replicate_2=strcat('G_',mat2str(Gaussian_number_counter1),'_R_2');
      %         Names_for_bar_graph{bar_position(Center_counter1)+1,1}=Bar_bin_name_replicate_2;
      %       end
      %
      %       if replicate_num >=3
      %         %replicate 3
      %         Values_for_bar_graph(bar_position(Center_counter1)+2,1)=Replicate3_normalised_raw_data(index_minimum);
      %         Bar_bin_name_replicate_3=strcat('G_',mat2str(Gaussian_number_counter1),'_R_3');
      %         Names_for_bar_graph{bar_position(Center_counter1)+2,1}=Bar_bin_name_replicate_3;
      %       end
      %
      %       if replicate_num ==4
      %         %replicate 4
      %         Values_for_bar_graph(bar_position(Center_counter1)+3,1)=Replicate4_normalised_raw_data(index_minimum);
      %         Bar_bin_name_replicate_4=strcat('G_',mat2str(Gaussian_number_counter1),'_R_4');
      %         Names_for_bar_graph{bar_position(Center_counter1)+3,1}=Bar_bin_name_replicate_4;
      %       end
      %
      % %Remove values of data already plotted
      % Replicate1_normalised_raw_data(index_minimum)=[];
      % Replicate2_normalised_raw_data(index_minimum)=[];
      % Replicate3_normalised_raw_data(index_minimum)=[];
      % Replicate4_normalised_raw_data(index_minimum)=[];
      % Center_to_plot(index_minimum)=[];
      %
      %       Center_counter1=Center_counter1+1;
      %       Gaussian_number_counter1=Gaussian_number_counter1+1;
      %  end
      
      
      %P8= bar(Values_for_bar_graph);
      xbar = linspace(-1/4,1/4,replicate_num);
      if replicate_num==1
        xbar = 0;
      end
      for ii = 1:number_of_gaussian_to_plot
        b1 = bar(ii+xbar,Replicate_normalised_raw_data(:,ii));
        set(b1,'facecolor',colour_to_use(ii,:))
        %       for jj = 1:replicate_num
        %         y = Replicate_normalised_raw_data(jj,ii);
        %         text(ii+xbar(jj)-0.025,y + 0.2*sign(y),['Rep' num2str(jj)]);
        %       end
      end
      y = ylim;
      y(1) = min([-1 y(1)]);
      y(2) = max([1 y(2)]);
      for ii =1:number_of_gaussian_to_plot-1
        plot([1 1]*(ii+0.5),y,'--r')
      end
      title_name_bar=strcat('Log2 ratio of treated to untreated at gaussian apex of :',Protein_to_plot);
      plot([0 replicate_num+1],[1 1],':k','LineWidth',1);
      plot([0 replicate_num+1],[-1 -1],':k','LineWidth',1);
      title(title_name_bar,'FontSize', 12);
      ylabel('Log2 ratio','FontSize', 10);
      %xlabel('Fractions','FontSize', 10);
      set(gca, 'XTickLabel',Names_for_bar_graph, 'XTick',1:numel(Names_for_bar_graph),'FontSize', 12);
      xlim([0.5 number_of_gaussian_to_plot+0.5]);
      
      %Save image
      Save_name_plot=[figdir 'Comparison/ProteinGaussianMaps/' mat2str(Gaussian_counter1),'_2_PCP_SEC_Profiles of of_',Protein_to_plot,'.png'];
      List_of_pdf{List_of_pdf_counter,1}=Save_name_plot;
      List_of_pdf{List_of_pdf_counter,2}=List_of_pdf_counter;
      List_of_pdf_counter=List_of_pdf_counter+1;
      saveas(f9, Save_name_plot);
      %print('-dpdf', '-r600', Save_name_plot);
      
      close 'all';
    end
    
  end
end


%%
if 0 % this is unfinished
  
  % 4. Make pie chart figure
  
  %Create array to store which replicate each protein
  Proteins_observed_in_replicates=zeros(length(Unique_protein_names), 8);
  
  for ii = 1:Nchannels
    for count_shared_guassians3 = 1:length(Unique_protein_names)
      internal_location_of_protein_of_interest_MvsL = find(strmatch(Unique_protein_names(count_shared_guassians3), GaussSummary(ii).Protein_name, 'exact'));
      %Record which repliates the protein of interest was seen in MvsL replicates
      Number_of_times_seen_MvsL=length(internal_location_of_protein_of_interest_MvsL);
      for Number_of_times_seen_counter1= 1:Number_of_times_seen_MvsL
        if GaussSummary(ii).Replicate(internal_location_of_protein_of_interest_MvsL(Number_of_times_seen_counter1)) == 1;
          Proteins_observed_in_replicates(count_shared_guassians3,1)=1;
        elseif GaussSummary(ii).Replicate(internal_location_of_protein_of_interest_MvsL(Number_of_times_seen_counter1)) == 2;
          Proteins_observed_in_replicates(count_shared_guassians3,2)=1;
        elseif GaussSummary(ii).Replicate(internal_location_of_protein_of_interest_MvsL(Number_of_times_seen_counter1)) == 3;
          Proteins_observed_in_replicates(count_shared_guassians3,3)=1;
        elseif GaussSummary(ii).Replicate(internal_location_of_protein_of_interest_MvsL(Number_of_times_seen_counter1)) == 4;
          Proteins_observed_in_replicates(count_shared_guassians3,4)=1;
        elseif GaussSummary(ii).Replicate(internal_location_of_protein_of_interest_MvsL(Number_of_times_seen_counter1)) == 0;
          Proteins_observed_in_replicates(count_shared_guassians3,4)=0;
          Proteins_observed_in_replicates(count_shared_guassians3,3)=0;
          Proteins_observed_in_replicates(count_shared_guassians3,2)=0;
          Proteins_observed_in_replicates(count_shared_guassians3,1)=0;
        end
      end
    end
    %Write Protein information to Protein structure array
    Protein_information.Protein_name=Unique_protein_names;
    Protein_information.Observed_MvsL_replicate1= Proteins_observed_in_replicates(:,1);
    Protein_information.Observed_MvsL_replicate2= Proteins_observed_in_replicates(:,2);
    Protein_information.Observed_MvsL_replicate3= Proteins_observed_in_replicates(:,3);
    Protein_information.Observed_MvsL_replicate4= Proteins_observed_in_replicates(:,4);
    Protein_information.Observed_HvsL_replicate1= Proteins_observed_in_replicates(:,5);
    Protein_information.Observed_HvsL_replicate2= Proteins_observed_in_replicates(:,6);
    Protein_information.Observed_HvsL_replicate3= Proteins_observed_in_replicates(:,7);
    Protein_information.Observed_HvsL_replicate4= Proteins_observed_in_replicates(:,8);
    
  end
  
  
  %Number of times observed in specific channels
  Protein_information.Total_number_channels_observed= (sum(Proteins_observed_in_replicates(:,:)'))';
  Protein_information.Observed_in_MvsL= (sum(Proteins_observed_in_replicates(:,1:4)')>0)';
  Protein_information.Observed_in_HvsL= (sum(Proteins_observed_in_replicates(:,5:8)')>0)';
  Protein_information.Number_Observed_in_MvsL= (sum(Proteins_observed_in_replicates(:,1:4)'))';
  Protein_information.Number_Observed_in_HvsL= (sum(Proteins_observed_in_replicates(:,5:8)'))';
  Protein_information.Observed_only_HvsL= ((Protein_information.Observed_in_HvsL>0) & (Protein_information.Observed_in_MvsL==0));
  Protein_information.Observed_only_MvsL= ((Protein_information.Observed_in_MvsL>0) & (Protein_information.Observed_in_HvsL==0));
  Protein_information.Observed_both= ((Protein_information.Observed_in_MvsL>0) & (Protein_information.Observed_in_HvsL>0));
  
  %Number of time observed in replicates
  Protein_information.Replicate1=Protein_information.Observed_MvsL_replicate1 & Protein_information.Observed_HvsL_replicate1;
  Protein_information.Replicate2=Protein_information.Observed_MvsL_replicate2 & Protein_information.Observed_HvsL_replicate2;
  Protein_information.Replicate3=Protein_information.Observed_MvsL_replicate3 & Protein_information.Observed_HvsL_replicate3;
  Protein_information.Replicate4=Protein_information.Observed_MvsL_replicate4 & Protein_information.Observed_HvsL_replicate4;
  
  %Number of times observed between replicates
  Protein_information.Observed_in_zero_HvsL= (Protein_information.Number_Observed_in_HvsL==0);
  Protein_information.Observed_in_zero_MvsL= (Protein_information.Number_Observed_in_MvsL==0);
  Protein_information.Observed_in_one_HvsL= (Protein_information.Number_Observed_in_HvsL==1);
  Protein_information.Observed_in_one_MvsL= (Protein_information.Number_Observed_in_MvsL==1);
  Protein_information.Observed_in_two_HvsL= (Protein_information.Number_Observed_in_HvsL==2);
  Protein_information.Observed_in_two_MvsL= (Protein_information.Number_Observed_in_MvsL==2);
  Protein_information.Observed_in_three_HvsL= (Protein_information.Number_Observed_in_HvsL==3);
  Protein_information.Observed_in_three_MvsL= (Protein_information.Number_Observed_in_MvsL==3);
  Protein_information.Observed_in_four_HvsL= (Protein_information.Number_Observed_in_HvsL==4);
  Protein_information.Observed_in_four_MvsL= (Protein_information.Number_Observed_in_MvsL==4);
  
  %Number of times observed between replicates and channels
  %HvsL
  Protein_information.Observed_in_one_HvsL_both_channels= Protein_information.Observed_in_one_HvsL & Protein_information.Observed_both;
  Protein_information.Observed_in_two_HvsL_both_channels= Protein_information.Observed_in_two_HvsL & Protein_information.Observed_both;
  Protein_information.Observed_in_three_HvsL_both_channels= Protein_information.Observed_in_three_HvsL & Protein_information.Observed_both;
  Protein_information.Observed_in_four_HvsL_both_channels= Protein_information.Observed_in_four_HvsL & Protein_information.Observed_both;
  %MvsL
  Protein_information.Observed_in_one_MvsL_both_channels= Protein_information.Observed_in_one_MvsL & Protein_information.Observed_both;
  Protein_information.Observed_in_two_MvsL_both_channels= Protein_information.Observed_in_two_MvsL & Protein_information.Observed_both;
  Protein_information.Observed_in_three_MvsL_both_channels= Protein_information.Observed_in_three_MvsL & Protein_information.Observed_both;
  Protein_information.Observed_in_four_MvsL_both_channels= Protein_information.Observed_in_four_MvsL & Protein_information.Observed_both;
  %HvsL
  Protein_information.Observed_in_zero_MvsL_HvsL_channels= Protein_information.Observed_in_zero_MvsL & Protein_information.Observed_only_HvsL;
  Protein_information.Observed_in_one_HvsL_HvsL_channels= Protein_information.Observed_in_one_HvsL & Protein_information.Observed_only_HvsL;
  Protein_information.Observed_in_two_HvsL_HvsL_channels= Protein_information.Observed_in_two_HvsL & Protein_information.Observed_only_HvsL;
  Protein_information.Observed_in_three_HvsL_HvsL_channels= Protein_information.Observed_in_three_HvsL & Protein_information.Observed_only_HvsL;
  Protein_information.Observed_in_four_HvsL_HvsL_channels= Protein_information.Observed_in_four_HvsL & Protein_information.Observed_only_HvsL;
  %MvsL
  Protein_information.Observed_in_zero_HvsL_MvsL_channels= Protein_information.Observed_in_zero_HvsL & Protein_information.Observed_only_MvsL;
  Protein_information.Observed_in_one_MvsL_MvsL_channels= Protein_information.Observed_in_one_MvsL & Protein_information.Observed_only_MvsL;
  Protein_information.Observed_in_two_MvsL_MvsL_channels= Protein_information.Observed_in_two_MvsL & Protein_information.Observed_only_MvsL;
  Protein_information.Observed_in_three_MvsL_MvsL_channels= Protein_information.Observed_in_three_MvsL & Protein_information.Observed_only_MvsL;
  Protein_information.Observed_in_four_MvsL_MvsL_channels= Protein_information.Observed_in_four_MvsL & Protein_information.Observed_only_MvsL;
  
  %Determine overlap of proteins identifications between replicates
  
  %Compare how many proteins were observed in each channel (use for Figure 1)
  pie_figure1(1,1)=nnz(Protein_information.Total_number_channels_observed(Protein_information.Total_number_channels_observed==1));
  pie_figure1(1,2)=nnz(Protein_information.Total_number_channels_observed(Protein_information.Total_number_channels_observed==2));
  pie_figure1(1,3)=nnz(Protein_information.Total_number_channels_observed(Protein_information.Total_number_channels_observed==3));
  pie_figure1(1,4)=nnz(Protein_information.Total_number_channels_observed(Protein_information.Total_number_channels_observed==4));
  pie_figure1(1,5)=nnz(Protein_information.Total_number_channels_observed(Protein_information.Total_number_channels_observed==5));
  pie_figure1(1,6)=nnz(Protein_information.Total_number_channels_observed(Protein_information.Total_number_channels_observed==6));
  pie_figure1(1,7)=nnz(Protein_information.Total_number_channels_observed(Protein_information.Total_number_channels_observed==7));
  pie_figure1(1,8)=nnz(Protein_information.Total_number_channels_observed(Protein_information.Total_number_channels_observed==8));
  
  %Comparing what channels the proteins were observed in across all
  %replicates
  pie_figure2(1,1)=nnz(Protein_information.Observed_both);
  pie_figure2(1,2)=nnz(Protein_information.Observed_only_HvsL);
  pie_figure2(1,3)=nnz(Protein_information.Observed_only_MvsL);
  
  %Compare what the number of times observed between replicates (Figure 3)
  pie_figure3(1,1)=nnz(Protein_information.Observed_in_zero_MvsL);
  pie_figure3(1,3)=nnz(Protein_information.Observed_in_one_MvsL);
  pie_figure3(1,5)=nnz(Protein_information.Observed_in_two_MvsL);
  pie_figure3(1,7)=nnz(Protein_information.Observed_in_three_MvsL);
  pie_figure3(1,9)=nnz(Protein_information.Observed_in_four_MvsL);
  pie_figure3(1,2)=nnz(Protein_information.Observed_in_zero_HvsL);
  pie_figure3(1,4)=nnz(Protein_information.Observed_in_one_HvsL);
  pie_figure3(1,6)=nnz(Protein_information.Observed_in_two_HvsL);
  pie_figure3(1,8)=nnz(Protein_information.Observed_in_three_HvsL);
  pie_figure3(1,10)=nnz(Protein_information.Observed_in_four_HvsL);
  
  %Compare the Number of times observed between replicates and channels
  %(Figure 4)
  pie_figure4(1,1)=nnz(Protein_information.Observed_in_zero_HvsL_MvsL_channels);
  pie_figure4(1,3)=nnz(Protein_information.Observed_in_one_HvsL_HvsL_channels);
  pie_figure4(1,5)=nnz(Protein_information.Observed_in_two_HvsL_HvsL_channels);
  pie_figure4(1,2)=nnz(Protein_information.Observed_in_zero_MvsL_HvsL_channels);
  pie_figure4(1,4)=nnz(Protein_information.Observed_in_one_MvsL_MvsL_channels);
  pie_figure4(1,6)=nnz(Protein_information.Observed_in_two_MvsL_MvsL_channels);
  pie_figure4(1,7)=nnz(Protein_information.Observed_in_three_HvsL_HvsL_channels);
  pie_figure4(1,8)=nnz(Protein_information.Observed_in_three_MvsL_MvsL_channels);
  pie_figure4(1,9)=nnz(Protein_information.Observed_in_four_HvsL_HvsL_channels);
  pie_figure4(1,10)=nnz(Protein_information.Observed_in_four_MvsL_MvsL_channels);
  
  %(Figure 5)
  pie_figure5(1,1)=nnz(Protein_information.Observed_in_one_HvsL_both_channels);
  pie_figure5(1,3)=nnz(Protein_information.Observed_in_two_HvsL_both_channels);
  pie_figure5(1,5)=nnz(Protein_information.Observed_in_three_HvsL_both_channels);
  pie_figure5(1,7)=nnz(Protein_information.Observed_in_four_HvsL_both_channels);
  pie_figure5(1,2)=nnz(Protein_information.Observed_in_one_MvsL_both_channels);
  pie_figure5(1,4)=nnz(Protein_information.Observed_in_two_MvsL_both_channels);
  pie_figure5(1,6)=nnz(Protein_information.Observed_in_three_MvsL_both_channels);
  pie_figure5(1,8)=nnz(Protein_information.Observed_in_four_MvsL_both_channels);
  
  
  %Compare what the number of times proteins observed between replicates (Figure 3)
  pie_figure6(1,1)=nnz(Protein_information.Replicate1);
  pie_figure6(1,2)=nnz(Protein_information.Replicate2);
  pie_figure6(1,3)=nnz(Protein_information.Replicate3);
  pie_figure6(1,4)=nnz(Protein_information.Replicate4);
  
  % #output
  %write out list of which proteins were observed (as gaussians in each replicate)
  fid_Proteins_in_each_rep = fopen([datadir 'Comparison/Protein_gaussian_observed_in_each_replicate.csv'],'w');
  fprintf (fid_Proteins_in_each_rep,'%s,%s,%s,%s,%s,%s,%s,%s,%s,\n',...
    'Protein_name','Replicate 1 MvsL channel','Replicate 2 MvsL channel','Replicate 3 MvsL channel','Replicate 4 MvsL channel','Replicate 1 HvsL channel','Replicate 2 HvsL channel','Replicate 3 HvsL channel','Replicate 4 HvsL channel'); %Write Header
  for write_out_Proteins_in_rep = 1:length(Unique_protein_names)
    fprintf(fid_Proteins_in_each_rep,'%s,', Unique_protein_names{write_out_Proteins_in_rep});
    fprintf(fid_Proteins_in_each_rep,'%6.4g,', Proteins_observed_in_replicates(write_out_Proteins_in_rep,:));
    fprintf(fid_Proteins_in_each_rep,'\n');
  end
  fclose(fid_Proteins_in_each_rep);
  
  %write out analysis of observation
  fid_Proteins_in_each_rep2 = fopen([datadir 'Comparison/Protein_observed_in_each_replicate.csv'],'w');
  fprintf (fid_Proteins_in_each_rep2,'%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,\n',...
    'Protein_Observed_in_both_channels','Protein_Observed_only_in_MvsL_channel','Protein_Observed_only_in_HvsL_channel',...
    'Protein_Observed_in_MvsL_once_channels','Protein_Observed_in_MvsL_twice_channel','Protein_Observedin_MvsL_three_times_channel','Protein_Observedin_MvsL_four_times_channel',...
    'Protein_Observed_in_HvsL_once_channels','Protein_Observed_in_HvsL_twice_channel','Protein_Observedin_HvsL_three_times_channel','Protein_Observedin_HvsL_four_times_channel'); %Write Header
  fprintf(fid_Proteins_in_each_rep2,'%6.4f,', pie_figure2(1,1));
  fprintf(fid_Proteins_in_each_rep2,'%6.4f,', pie_figure2(1,3));
  fprintf(fid_Proteins_in_each_rep2,'%6.4f,', pie_figure2(1,2));
  fprintf(fid_Proteins_in_each_rep2,'%6.4f,', pie_figure3(1,3));
  fprintf(fid_Proteins_in_each_rep2,'%6.4f,', pie_figure3(1,5));
  fprintf(fid_Proteins_in_each_rep2,'%6.4f,', pie_figure3(1,7));
  fprintf(fid_Proteins_in_each_rep2,'%6.4f,', pie_figure3(1,9));
  fprintf(fid_Proteins_in_each_rep2,'%6.4f,', pie_figure3(1,4));
  fprintf(fid_Proteins_in_each_rep2,'%6.4f,', pie_figure3(1,6));
  fprintf(fid_Proteins_in_each_rep2,'%6.4f,', pie_figure3(1,8));
  fprintf(fid_Proteins_in_each_rep2,'%6.4f,', pie_figure3(1,10));
  fclose(fid_Proteins_in_each_rep2);
  
  %write out analysis of observation
  fid_Proteins_in_each_rep3 = fopen([datadir 'Comparison/Protein_observed_in_each_replicate_and_channels.csv'],'w');
  fprintf (fid_Proteins_in_each_rep3,'%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,\n',...
    'Protein_Observed_only_in_MvsL_channel','Protein_Observed_only_in_HvsL_channel','Protein_Observed_only_in_HvsL_channel_once',...
    'Protein_Observed_only_in_HvsL_channel_twice','Protein_Observed_only_in_HvsL_channel_three_times','Protein_Observed_only_in_HvsL_channel_four_times',...
    'Protein_Observed_only_in_MvsL_channel_once','Protein_Observed_only_in_MvsL_channesl_twice','Protein_Observed_only_in_MvsL_channel_three_times',...
    'Protein_Observed_only_in_MvsL_channel_four_times','Protein_Observed_Both_channesl_once_HvsL','Protein_Observed_Both_channesl_twice_HvsL',...
    'Protein_Observed_Both_channesl_three_times_HvsL','Protein_Observed_Both_channesl_four_times_HvsL','Protein_Observed_Both_channesl_once_MvsL',...
    'Protein_Observed_Both_channesl_twice_MvsL','Protein_Observed_Both_channesl_three_times_MvsL','Protein_Observed_Both_channesl_four_times_MvsL'); %Write Header
  fprintf(fid_Proteins_in_each_rep3,'%6.4f,', pie_figure4(1,1));
  fprintf(fid_Proteins_in_each_rep3,'%6.4f,', pie_figure4(1,2));
  fprintf(fid_Proteins_in_each_rep3,'%6.4f,', pie_figure4(1,3));
  fprintf(fid_Proteins_in_each_rep3,'%6.4f,', pie_figure4(1,5));
  fprintf(fid_Proteins_in_each_rep3,'%6.4f,', pie_figure4(1,7));
  fprintf(fid_Proteins_in_each_rep3,'%6.4f,', pie_figure4(1,9));
  fprintf(fid_Proteins_in_each_rep3,'%6.4f,', pie_figure4(1,4));
  fprintf(fid_Proteins_in_each_rep3,'%6.4f,', pie_figure4(1,6));
  fprintf(fid_Proteins_in_each_rep3,'%6.4f,', pie_figure4(1,8));
  fprintf(fid_Proteins_in_each_rep3,'%6.4f,', pie_figure4(1,10));
  fprintf(fid_Proteins_in_each_rep3,'%6.4f,', pie_figure5(1,1));
  fprintf(fid_Proteins_in_each_rep3,'%6.4f,', pie_figure5(1,3));
  fprintf(fid_Proteins_in_each_rep3,'%6.4f,', pie_figure5(1,5));
  fprintf(fid_Proteins_in_each_rep3,'%6.4f,', pie_figure5(1,7));
  fprintf(fid_Proteins_in_each_rep3,'%6.4f,', pie_figure5(1,2));
  fprintf(fid_Proteins_in_each_rep3,'%6.4f,', pie_figure5(1,4));
  fprintf(fid_Proteins_in_each_rep3,'%6.4f,', pie_figure5(1,6));
  fprintf(fid_Proteins_in_each_rep3,'%6.4f,', pie_figure5(1,8));
  fclose(fid_Proteins_in_each_rep3);
  
  
  %Figure parameters, Figure 1
  %Note: Option are present to also output a visulation on the overlap
  %between replicates and isotopologue channels but this has found to be
  %largely uninformative
  f1 = figure;
  P(1)=subplot(1,1,1);
  Figure1= pie(pie_figure1);
  hp = findobj(P(1), 'Type', 'patch');
  figure_pos = get(P(1),'position');
  
  %create varible for legend
  combinedstrings=cell(length(hp),1);
  
  for colour_count=1:length(hp)
    set(hp(colour_count), 'FaceColor', colour_to_use(colour_count,:)); % set colour according to defined palate
    combinedstrings{colour_count,1}=mat2str(colour_count);
  end
  
  set(P(1), 'position',[figure_pos(1)*1.6 figure_pos(2)*1.0 figure_pos(3:4)]);
  
  %Create label
  hText = findobj(Figure1,'Type','text'); % text handles
  percentValues = get(hText,'String'); % percent values
  
  %Position text
  textPositions_cell = get(hText,'Position'); % cell array
  textPositions = cell2mat(textPositions_cell); % numeric array
  textPositions(:,1) = textPositions(:,1)*1.2; % add offset
  %set(hText,'Position',num2cell(textPositions,[3,2]));
  for jj = 1:size(hText,1)
    set(hText(jj),'Position',textPositions(jj,:));
  end
  
  %# add legend
  Legend_1_position = legend(P(1), combinedstrings);
  Legend_pos = get(Legend_1_position,'position');
  set(Legend_1_position, 'position',[0.05 0.25 0.15 0.5]);
  
  saveas(f1, [figdir 'Comparison/Protein_analysis.png']);
  
  tt = toc;
  tt = round(tt/10)*10;
  fprintf('  ...  %.2f seconds\n',tt)
  
end