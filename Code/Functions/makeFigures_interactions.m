%MAKEFIGURES_INTERACTIONS Makes figures for the PRINCE Interactions module.

precPlot = num2str(round(desiredPrecision(di)*100));


%% How many interactions were found in at least N replicates?

try
  sf = [figdir 'Number_interactions_per_channel.png'];
  
  jj = 1;
  PA = Precision_array;
  
  bar_x = 1:(number_of_channels);
  if size(Precision_array,1)==1
    bar_y = Precision_array;
    bar_x(2) = 0;
    bar_y(2,:) = 0;
    x0 = [.5 1.5];
  else
    x0 = [bar_x(1)-1 bar_x(end)+1];
    bar_x = 1:(number_of_channels);
    bar_y = Precision_array;
  end
  prec_fig = bar_y(:,2) ./ (bar_y(:,1) + bar_y(:,2));
  
  myC= [30/255 144/255 255/255
    255/255 215/255 0/255
    178/255 34/255 34/255
    193/255 205/255 193/255];
  
  figure
  subplot(2,1,1)
  bar(bar_x, bar_y, 0.6, 'stack');
  y = ylim;
  for ii = 1:length(prec_fig)
    if isnan(prec_fig(ii))
      s = '--';
    else
      s = [num2str(round(prec_fig(ii)*1000)/10) '%'];
    end
    texty = sum(bar_y(ii,:));
    text(ii-0.25,texty+diff(y)*.04,s)
  end
  legend('Inter-complex (FP)', 'Intra-complex (TP)',  'Novel');
  xlim(x0)
  grid on
  ylabel('Number of interactions','FontSize', 8);
  xlabel('N conditions','FontSize', 8);
  title('Number of interactions in AT LEAST N conditions','FontSize', 12);
  
  % How many interactions were found in each channel?
  kk = 2;
  if kk == 1 % global precision
    I = interaction_final.global;
  else % replicate precision
    I = ones(size(interaction_final.global));
  end
  [Ia,Ib] = find(interaction_final.replicate_numbers>0 & repmat(I,1,size(interaction_final.replicate_numbers,2))>0 );
  N_intperrep = zeros(length(Ia), 1);
  itype = zeros(length(Ia), 1);
  for ii = 1:length(Ia)
    N_intperrep(ii) = interaction_final.replicate_numbers(Ia(ii),Ib(ii));
    if interaction_final.proteinInCorum(Ia(ii)) && interaction_final.interactionInCorum(Ia(ii))
      itype(ii) = 1;
    elseif interaction_final.proteinInCorum(Ia(ii)) && ~interaction_final.interactionInCorum(Ia(ii))
      itype(ii) = 2;
    else
      itype(ii) = 3;
    end
  end
  
  x = 1:number_of_replicates*number_of_channels;
  h_tp = hist(N_intperrep(itype==1),x);
  h_fp = hist(N_intperrep(itype==2),x);
  h_nov = hist(N_intperrep(itype==3),x);
  hh = [h_fp; h_tp; h_nov];
  xlim_flag = 0;
  if length(x)==1
    x(2) = x(1)+1;
    hh(:,2) = zeros(length(hh),1);
    xlim_flag = 1;
  end
  
  xlab = cell(number_of_replicates*number_of_channels,1);
  for ii = 1:number_of_replicates*number_of_channels
    sr = nan(length(user.silacratios),1);
    for jj = 1:length(user.silacratios)
      tmp = findstr(user.silacratios{jj},GaussIn{ii});
      if ~isempty(tmp)
        sr(jj) = tmp;
      end
    end
    Irep = findstr('rep',GaussIn{ii});
    rep = GaussIn{ii}(Irep+3);
    xlab{ii} = ['Rep' rep '-' user.silacratios{sr>0}];
  end
  
  subplot(2,1,2),hold on
  f2 = bar(x,hh', 0.6,'stacked');
  xlim(x0)
  ylabel('Number of interactions','FontSize', 10);
  xlabel('Replicate - Isotopologue channel','FontSize', 10);
  title('Number of interactions in EACH condition','FontSize', 12);
  set(gca,'xtick',x,'xticklabel',xlab,'fontsize',9)
  grid on
  
  saveas(gcf, sf);
catch
  warning('Failed to plot Number_interactions_per_channel.')
end


%% Final precision-recall, ROC curves

try
  [~,I] = sort(interaction_final.precisionDropoutavg(:,1),'descend');
  tmp = interaction_final.precisionDropoutavg(:,1);
  tmp(tmp==0) = nan;
  xi = 1:length(tmp);
  
  figure,hold on
  plot(xi,tmp(I),'k')
  grid on
  ylabel('Precision, TP/(TP+FP)','fontsize',12)
  xlabel('Number of interactions, TP+FP','fontsize',12)
  % Save figure
  set(gcf,'paperunits','inches','paperposition',[.25 2.5 9 9])
  sf=[figdir 'Precision_vs_NumberOfInteractions'];
  saveas(gcf, sf, 'png');
  
  
  ntp = cumsum(interaction_final.proteinInCorum(I)==1 & interaction_final.interactionInCorum(I)==1);
  figure,hold on
  plot(ntp / all_positives,tmp(I),'k')
  grid on
  ylabel('Precision, TP/(TP+FP)','fontsize',12)
  xlabel('Recall, TP/(TP+FN)','fontsize',12)
  % Save figure
  set(gcf,'paperunits','inches','paperposition',[.25 2.5 9 9])
  sf=[figdir 'Precision_vs_Recall'];
  saveas(gcf, sf, 'png');
catch
  warning('Failed to plot Precision_vs_Recall.')
end



%% Plot best 200 interactions

f1 = [user.maindir 'Output/Figures/Interactions/BestInteractions/'];
if not(exist(f1,'dir'))
  mkdir(f1)
end

% get chromatograms and proteins
chroms = cell(number_of_channels,length(rep2channel));
prots = cell(number_of_channels,length(rep2channel));
nfrac = 0;
for channel_counter = 1:number_of_channels
  replicatesThisChannel = find(rep2channel == channel_counter);
  for replicate_counter = replicatesThisChannel
    sf = [maindir '/Output/tmp/' 'data_rep' num2str(replicate_counter) '_chan' num2str(channel_counter) '.mat'];
    load(sf)
    chroms{channel_counter,replicate_counter} = Chromatograms_raw;
    prots{channel_counter,replicate_counter} = Protein;
    nfrac = max([nfrac size(Chromatograms_raw,2)]);
    clear TP_Matrix possibleInts Protein inverse_self Chromatograms Chromatograms_raw Dist
  end
end

[tmp,I] = sort(interaction_final.precisionDropout,'descend');
nplots = min([200  length(I)]);
h = figure;
for ii = 1:nplots
  try
    intI = I(ii);
    chans = interaction_final.channel{intI};
    protA = interaction_final.proteinA{intI};
    protB = interaction_final.proteinB{intI};
    
    subplots = [length(user.silacratios) number_of_replicates];
    
    % Select the figure and clear it
    set(0, 'CurrentFigure', h);
    clf reset;
    cc = 0;
    for jj = 1:length(user.silacratios)
      for kk = 1:length(rep2channel)
        if isempty(prots{jj,kk})
          continue;
        end
        cc = cc+1;
        Ia = find(ismember(prots{jj,kk}.Isoform,protA));
        Ib = find(ismember(prots{jj,kk}.Isoform,protB));
        
        subplot(subplots(1),subplots(2),cc), hold on
        plot(1:length(chroms{jj,kk}(Ia,:)),chroms{jj,kk}(Ia,:),'k','linewidth',2)
        plot(1:length(chroms{jj,kk}(Ib,:)),chroms{jj,kk}(Ib,:),'color',[.75 .75 .75],'linewidth',2)
        legend(protA,protB,'location','best')
        title([user.silacratios{jj} ' - replicate ' num2str(mod(kk-1,number_of_replicates)+1)])
        ylabel('Protein amount')
        xlabel('Fraction number')
        xlim([0 nfrac+1])
      end
    end
    ax = axis;
    s = ['Prec=' num2str(round(interaction_final.precisionDropout(intI)*100)) '%'];
    text(ax(1)+(ax(2)-ax(1))*.2,ax(3)+(ax(4)-ax(3))*.5,s)
    % Save figure
    sf=[f1 'Interactions_' num2str(ii) '_' protA '_' protB];
    set(gcf,'paperunits','inches','paperposition',[.25 2.5 9 7],'units','inches','position',[.25 2.5 12 7])
    saveas(gcf, sf, 'png');
  catch ME
    disp(['Interactions.m: Interaction plot ', num2str(ii), ' :', ME.message]);
  end
end

