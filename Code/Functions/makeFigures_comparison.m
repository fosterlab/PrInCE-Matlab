
colour_to_use=[0.254 0.411 0.882 %Colour: Royal Blue
  205/255 92/255 92/255 %Colour: Indian Red
  106/255 90/255 205/255 %Colour: Slate Blue
  135/255 206/255 250/255 %Colour: Sky Blue
  0.28 0.82 0.8 %Colour: Medium Turquoise
  178/255 34/255 34/255 %Colour: Firebrick
  255/255 69/255 0 %Colour: Orange Red
  244/255 238/255 224/255 %Colour: Honeydew 2
  34/255 139/255 34/255 %Colour: Forest Green
  222/255 184/255 135/255 %Colour: Burlywood
  186/255 85/255 211/255 %Colour: Medium Orchid
  219/255 112/255 147/255]; %Colour: Pale Violet Red


%% Log2 changes, normalized
figure
log2_sorted = sort(Combined_Gaussians.log2_normalised_gaussians(:));
subplot(2,1,1);hold on
scatter(1:length(log2_sorted(:)),log2_sorted(:),8,'fill');
title('Log2 changes, normalized (mean = 0)','FontSize', 10);
ylabel('Normalised Log2 ratio','FontSize', 10);
xlabel('Gaussians identified','FontSize', 10);
%Set limits based on observed values
%ylim([log2_sorted(1)*1.05,log2_sorted(Total_number_of_unique_gaussians)*1.05]);
%Set limits based on user defined cut off
plot([-100 length(log2_sorted(:))+100],[1 1],'--','color', [.7 .7 0]);
plot([-100 length(log2_sorted(:))+100],[-1 -1],'--r');
xlim([-10,sum(~isnan(log2_sorted(:)))+10]);
y1 = min(log2_sorted)*1.1;
y2 = max(log2_sorted)*1.1;
ylim([min([-1.1 y1]) max([1.1 y2])]);
grid on

%plot distribution of protein chnages across the fractions
% Count how many guassian change in each fraction
Hist_array=zeros(Nfraction,3);
for ii= 1:number_of_proteins
  %Define Center of Gaussian
  Center_of_gaussain= floor(Combined_Gaussians.Center(ii));
  %Ensure center is greater then zero
  
  if Combined_Gaussians.log2_normalised_gaussians(ii) >=1 && Center_of_gaussain <=fraction_to_plot
    
    Hist_array(Center_of_gaussain,1)=Hist_array(Center_of_gaussain,1)+1;
    Hist_array(Center_of_gaussain,2)=Hist_array(Center_of_gaussain,2)+1;
  elseif Combined_Gaussians.log2_normalised_gaussians(ii) <=-1 && Center_of_gaussain <=fraction_to_plot
    
    Hist_array(Center_of_gaussain,1)=Hist_array(Center_of_gaussain,1)+1;
    Hist_array(Center_of_gaussain,3)=Hist_array(Center_of_gaussain,3)+1;
  elseif Combined_Gaussians.log2_normalised_gaussians(ii) >=-1....
      && Combined_Gaussians.log2_normalised_gaussians(ii) <=1....
      && Center_of_gaussain <=fraction_to_plot
    
    Hist_array(Center_of_gaussain,1)=Hist_array(Center_of_gaussain,1)+1;
  end
end
subplot(2,1,2);
f6_figure=bar(1:Nfraction, [Hist_array(:,1) Hist_array(:,2) Hist_array(:,3)],0.6, 'stack');
for k=1:3
  set(f6_figure(k),'facecolor', myC(k,:), 'EdgeColor', 'k' )
end
legend('No change', 'Increase', 'Decrease', 'Location', 'Best');
xlim([-1,frac2+1]);
title('Changes in Gaussians observed across fractions');
xlabel('Fractions');
ylabel('# of Guassian centers');
saveas(gcf, [figdir 'Comparison/Fold_changes_Gaussians_normalized.png']);



%% Log2 changes, non normalized
figure
log2_sorted = sort(Combined_Gaussians.log2_of_gaussians(:));
subplot(2,1,1);hold on
scatter(1:length(log2_sorted(:)),log2_sorted(:),8,'fill');
title('Log2 changes','FontSize', 10);
ylabel('Log2 ratio','FontSize', 10);
xlabel('Gaussians identified','FontSize', 10);
%Set limits based on observed values
%ylim([log2_sorted(1)*1.05,log2_sorted(Total_number_of_unique_gaussians)*1.05]);
%Set limits based on user defined cut off
plot([-100 length(log2_sorted(:))+100],[1 1],'--','color', [.7 .7 0]);
plot([-100 length(log2_sorted(:))+100],[-1 -1],'--r');
xlim([-10,sum(~isnan(log2_sorted(:)))+10]);
y1 = min(log2_sorted)*1.1;
y2 = max(log2_sorted)*1.1;
ylim([min([-1.1 y1]) max([1.1 y2])]);
grid on

%plot distribution of protein chnages across the fractions
%plot distribution of protein chnages across the fractions
% Count how many guassian change in each fraction
Hist_array=zeros(Nfraction,3);
for ii= 1:number_of_proteins
  %Define Center of Gaussian
  Center_of_gaussain= floor(Combined_Gaussians.Center(ii));
  %Ensure center is greater then zero
  
  if Combined_Gaussians.log2_of_gaussians(ii) >=1 && Center_of_gaussain <=fraction_to_plot
    
    Hist_array(Center_of_gaussain,1)=Hist_array(Center_of_gaussain,1)+1;
    Hist_array(Center_of_gaussain,2)=Hist_array(Center_of_gaussain,2)+1;
  elseif Combined_Gaussians.log2_of_gaussians(ii) <=-1 && Center_of_gaussain <=fraction_to_plot
    
    Hist_array(Center_of_gaussain,1)=Hist_array(Center_of_gaussain,1)+1;
    Hist_array(Center_of_gaussain,3)=Hist_array(Center_of_gaussain,3)+1;
  elseif Combined_Gaussians.log2_of_gaussians(ii) >=-1....
      && Combined_Gaussians.log2_of_gaussians(ii) <=1....
      && Center_of_gaussain <=fraction_to_plot
    
    Hist_array(Center_of_gaussain,1)=Hist_array(Center_of_gaussain,1)+1;
  end
end
subplot(2,1,2);
f6_figure=bar(1:Nfraction, [Hist_array(:,1) Hist_array(:,2) Hist_array(:,3)],0.6, 'stack');
for k=1:3
  set(f6_figure(k),'facecolor', myC(k,:), 'EdgeColor', 'k' )
end
legend('No change', 'Increase', 'Decrease', 'Location', 'Best');
xlim([-1,frac2+1]);
title('Changes in Gaussians observed across fractions');
xlabel('Fractions');
ylabel('# of Guassian centers');
saveas(gcf, [figdir 'Comparison/Fold_changes_Gaussians.png']);



%% Individual protein plot 1
% One plot per protein in each replicate, i.e. can make multiple plots of the same protein.

cc = 0;
if user.fastcomparison==0
  for ii = 1:length(Finalised_Master_Gaussian_list.Protein_name)
    for rep = 1:user.Nreplicate
      cc = cc+1
      
      % get replicate + protein name
      protName = Finalised_Master_Gaussian_list.Protein_name{ii,1};
      %rep_protName = Finalised_Master_Gaussian_list.Replicate_Protein_identifier{ii};
      rep_protName = [num2str(rep) '_' protName]
      
      % Find this replicate + protein in Combined_Gaussians
      Icg = find(Combined_Gaussians.Replicate == rep & strcmp(protName,Combined_Gaussians.Protein_name));
      [C_comp,Iccomp] = sort(Combined_Gaussians.Center(Icg));
      
      % Get channel names and indices
      chan_den = Combined_Gaussians.denominatorChannnel{Icg};
      chan_num = Combined_Gaussians.numeratorChannnel{Icg};
      Iden = find(ismember(user.silacratios,chan_den));
      Inum = find(ismember(user.silacratios,chan_num));
      
      % get raw data
      Iraw = find(ismember(txt_val{Iden}(:,1), rep_protName));
      raw_num = num_val{Inum}(Iraw,2:end);
      raw_den = num_val{Iden}(Iraw,2:end);
      
      % get gaussian params
      Igd_den = find(ismember(GaussData{Iden}(:,2),rep_protName));
      hd = (GaussData{Iden}(Igd_den,7));
      cd = (GaussData{Iden}(Igd_den,8));
      wd = (GaussData{Iden}(Igd_den,9));
      gp_den = zeros(length(hd),3);
      for jj = 1:length(hd)
        gp_den(jj,1) = str2num(hd{jj});
        gp_den(jj,2) = str2num(cd{jj});
        gp_den(jj,3) = str2num(wd{jj});
      end
      Igd_num = find(ismember(GaussData{Inum}(:,2),rep_protName));
      hn = (GaussData{Inum}(Igd_num,7));
      cn = (GaussData{Inum}(Igd_num,8));
      wn = (GaussData{Inum}(Igd_num,9));
      gp_num = zeros(length(hn),3);
      for jj = 1:length(hn)
        gp_num(jj,1) = str2num(hn{jj});
        gp_num(jj,2) = str2num(cn{jj});
        gp_num(jj,3) = str2num(wn{jj});
      end
      
      % make gaussian curves
      xfit = linspace(0,length(raw_num)+1,101);
      yfit_den = zeros(size(xfit));
      for jj = 1:length(cd)
        yfit_den = yfit_den + gp_den(jj,1)*exp(-((xfit-gp_den(jj,2))/gp_den(jj,3)).^2);
      end
      xfit = linspace(0,length(raw_num)+1,101);
      yfit_num = zeros(size(xfit));
      for jj = 1:length(cn)
        yfit_num = yfit_num + gp_num(jj,1)*exp(-((xfit-gp_num(jj,2))/gp_num(jj,3)).^2);
      end
      
      % get fold changes for each gaussian
      fold_raw = Combined_Gaussians.log2_of_gaussians(Icg);
      fold_norm = Combined_Gaussians.log2_normalised_gaussians(Icg);
      
      figure
      subplot(3,1,1), hold on % raw chromatograms
      xraw = 1:length(raw_num);
      scatter(xraw,raw_num,20,colour_to_use(2,:),'filled')
      scatter(xraw,raw_den,20,colour_to_use(1,:),'filled')
      plot(xraw,raw_num,'color',colour_to_use(2,:));
      plot(xraw,raw_den,'color',colour_to_use(1,:));
      xlim([xraw(1) xraw(end)])
      y = ylim;
      if y(1)>0
        y(1) = 0.01;
      end
      ylim(y)
      for jj = 1:length(C_comp)
        plot([1 1]*C_comp(jj), y, '--k')
        text(C_comp(jj)+0.1, y(1) + diff(y)*.8, num2str(jj))
      end
      xlabel('Fraction','fontsize',10)
      ylabel('Isotopologue ratio','FontSize', 10);
      legend(user.silacratios([Inum Iden]),'location','best')
      title('Raw chromatogram','fontsize',10)
      ax = axis;
      
      subplot(3,1,2), hold on % fit gaussians, numerator channel
      plot(xfit,yfit_num,'color',colour_to_use(2,:));
      plot(xfit,yfit_den,'color',colour_to_use(1,:));
      xlim([xraw(1) xraw(end)])
      for jj = 1:length(C_comp)
        plot([1 1]*C_comp(jj), y, '--k')
        text(C_comp(jj)+0.1, y(1) + diff(y)*.8, num2str(jj))
      end
      xlabel('Fraction','fontsize',10)
      ylabel('Isotopologue ratio','FontSize', 10);
      legend(user.silacratios([Inum Iden]),'location','best')
      title('Fitted Gaussians','fontsize',10)
      axis(ax)
      
      subplot(3,2,5), hold on
      bar(fold_raw(Iccomp))
      y = ylim;
      y(1) = min([-1.1 y(1)]);
      y(2) = max([1.1 y(2)]);
      ylim(y)
      I = find(isnan(fold_raw(Iccomp)));
      for jj = 1:length(I)
        patch([-.4 -.4 .4 .4]+I(jj),y([1 2 2 1]),[-1 -1 -1 -1],[.9 .9 .9],'edgecolor',[.9 .9 .9])
      end
      ylabel('Fold change','fontsize',10)
      xlabel('Gaussian number','fontsize',10)
      
      subplot(3,2,6), hold on
      bar(fold_norm(Iccomp))
      y = ylim;
      y(1) = min([-1.1 y(1)]);
      y(2) = max([1.1 y(2)]);
      ylim(y)
      I = find(isnan(fold_norm(Iccomp)));
      for jj = 1:length(I)
        patch([-.4 -.4 .4 .4]+I(jj),y([1 2 2 1]),[-1 -1 -1 -1],[.9 .9 .9],'edgecolor',[.9 .9 .9])
      end
      ylabel('Norml. fold change','fontsize',10)
      xlabel('Gaussian number','fontsize',10)
      
      sf = [figdir 'Comparison/IndividualProteins/' rep_protName '_foldchange.png'];
      saveas(gcf, sf);
      close all
      
    end
  end
end





%%

if 0
  
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
    title_name_bar=strcat('Normalized log2 ratio of treated to untreated at gaussian apex of :',Protein_to_plot);
    plot([0 replicate_num+1],[1 1],':k','LineWidth',1);
    plot([0 replicate_num+1],[-1 -1],':k','LineWidth',1);
    title(title_name_bar,'FontSize', 12);
    ylabel('Log2 ratio','FontSize', 10);
    set(gca, 'XTickLabel',Names_for_bar_graph, 'XTick',1:numel(Names_for_bar_graph),'FontSize', 12);
    xlim([0.5 number_of_gaussian_to_plot+0.5]);
    
    Save_name_replicates=[figdir 'Comparison/ProteinGaussianMaps/' mat2str(Gaussian_counter1),'_1_PCP_SEC_Profiles of_',Protein_to_plot,'.png'];
    saveas(gcf, Save_name_replicates);
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
