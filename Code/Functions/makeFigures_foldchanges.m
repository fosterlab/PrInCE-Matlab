%MAKEFIGURES_FOLDCHANGES Makes figures for the PRINCE FoldChanges module.

colour_to_use=[0.254 0.411 0.882 %Colour: Royal Blue
  205/255 92/255 92/255 %Colour: Indian Red
  0.28 0.82 0.8 %Colour: Medium Turquoise
  135/255 206/255 250/255 %Colour: Sky Blue
  222/255 184/255 135/255 %Colour: Burlywood
  178/255 34/255 34/255 %Colour: Firebrick
  255/255 69/255 0 %Colour: Orange Red
  106/255 90/255 205/255 %Colour: Slate Blue
  244/255 238/255 224/255 %Colour: Honeydew 2
  34/255 139/255 34/255 %Colour: Forest Green
  186/255 85/255 211/255 %Colour: Medium Orchid
  219/255 112/255 147/255]; %Colour: Pale Violet Red




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
      if isempty(Icg)
        continue
      end
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
      
      sf = [figdir '/IndividualProteins/' rep_protName '_foldchange.png'];
      saveas(gcf, sf);
      close all
      
    end
  end
end




%% Log2 changes, normalized + non-normalized
% Use Finalised_Master_Gaussian_list

figure, hold on
log2_sorted = sort(Finalised_Master_Gaussian_list.foldChange_normalized(:));
title('Log2 changes','FontSize', 10);
ylabel('Normalised Log2 ratio','FontSize', 10);
xlabel('Gaussians identified','FontSize', 10);
scatter(1:length(log2_sorted(:)),log2_sorted(:),8,[27,158,119]/255,'fill');

log2_sorted = sort(Finalised_Master_Gaussian_list.foldChange(:));
title('Log2 changes','FontSize', 10);
ylabel('Log2 ratio','FontSize', 10);
xlabel('Gaussians identified','FontSize', 10);
scatter(1:length(log2_sorted(:)),log2_sorted(:),8,[217,95,2]/255,'fill');

plot([-100 length(log2_sorted(:))+100],[1 1],'--','color', [.7 .7 0]);
plot([-100 length(log2_sorted(:))+100],[-1 -1],'--r');
xlim([-10,sum(~isnan(log2_sorted(:)))+10]);
y1 = min(log2_sorted)*1.1;
y2 = max(log2_sorted)*1.1;
ylim([min([-1.1 y1]) max([1.1 y2])]);
grid on

legend('Normalized (median = 0)','Raw','location','northwest')

saveas(gcf, [figdir '/Fold_changes_Gaussians.png']);



%% Histograms of changes per fraction

x = 1:fraction_to_plot;
ttest_sig = nan(length(x),3);
ranksum_sig = nan(length(x),3);
fold_sig = nan(length(x),3);
for xi = 1:length(x)
  I = Finalised_Master_Gaussian_list.Center>x(xi)-.5 & Finalised_Master_Gaussian_list.Center<x(xi)+.5;
  D1 = Finalised_Master_Gaussian_list.Ttest_p_values(I);
  D2 = Finalised_Master_Gaussian_list.MWW_p_values(I);
  D3 = Finalised_Master_Gaussian_list.foldChange_normalized(I);
  
  % significant t-test
  ttest_sig(xi,2) = sum(D1(:)<.05 & D3(:)<0);
  ttest_sig(xi,3) = sum(D1(:)<.05 & D3(:)>0);
  ttest_sig(xi,1) = sum(I(:)) - sum(ttest_sig(xi,2:3));
  
  % significant t-test
  ranksum_sig(xi,2) = sum(D2(:)<.05 & D3(:)<0);
  ranksum_sig(xi,3) = sum(D2(:)<.05 & D3(:)>0);
  ranksum_sig(xi,1) = sum(I(:)) - sum(ranksum_sig(xi,2:3));
  
  % "significant" fold change
  fold_sig(xi,2) = sum(D3(:)<-1);
  fold_sig(xi,3) = sum(D3(:)>1);
  fold_sig(xi,1) = sum(I(:)) - sum(fold_sig(xi,2:3));
end

figure
subplot(2,1,1)
b1 = bar(fold_sig,.6,'stacked');
for kk = 1:length(b1)
  set(b1(kk),'facecolor', colour_to_use(kk,:), 'EdgeColor', colour_to_use(kk,:) )
end
legend('No change', 'Decrease', 'Increase', 'Location', 'northeast');
xlim([-1,frac2+3]);
title('Change detected (normalized fold > 2)');
xlabel('Fractions');
ylabel('# of Guassian centers');

subplot(2,1,2)
b2 = bar(ttest_sig,.6,'stacked');
for kk = 1:length(b1)
  set(b2(kk),'facecolor', colour_to_use(kk,:), 'EdgeColor', colour_to_use(kk,:))
end
legend('No change', 'Decrease', 'Increase', 'Location', 'northeast');
xlim([-1,frac2+3]);
title('Change detected, t-test (p<0.05)');
xlabel('Fractions');
ylabel('# of Guassian centers');

saveas(gcf, [figdir '/Histogram_changes_per_fraction.png']);


