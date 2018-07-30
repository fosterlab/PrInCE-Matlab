%MAKEFIGURES_ALIGNMENT Makes figures for the PRINCE Alignment module.

%colsr = {[1 0 0 ] [0 1 0] [0 0 1]};
colsr = rand(Nreplicates,3);


%% Bar alignment figures
% use pfit from Alignment.m
try
  Image=figure;
  
  for kk = 1:2
    subplot(2,1,kk)
    hold on
    
    ci = 1;
    x = 1:55;
    
    x_rep = nan(Nreplicates,length(x));
    for rr = 1:Nreplicates
      if kk == 1
        x_rep(rr,:) = x;
      else
        b = pfit(ci,rr,1); % intercept
        m = pfit(ci,rr,2); % slope
        x_rep(rr,:) = x*m + b;
      end
    end
    
    % make initial patches for legend
    for rr = 1:Nreplicates
      patch([-1 -1 -2 -2]*.01 - 100,[-2 -1 -1 -2]*.01 -100,colsr(rr,:))
    end
    
    for rr = 1:Nreplicates
      for pp = 1:length(x)-1
        if mod(pp,5) == 0;
          cc = [.5 .5 .5];
        else
          cc = colsr(rr,:);
        end
        patch([x_rep(rr,pp) x_rep(rr,pp) x_rep(rr,pp+1) x_rep(rr,pp+1)],[rr rr+1 rr+1 rr],cc)
      end
      patch([x_rep(rr,end) x_rep(rr,end) x_rep(rr,end)+1 x_rep(rr,end)+1],[rr rr+1 rr+1 rr],cc)
    end
    
    %Figure details
    axis([-3 60 0.5 Nreplicates+1.5]);
    if kk == 1
      title('Before alignment','FontSize', 12);
    else
      title('After alignment','FontSize', 12);
    end
    ylabel('Replicates','FontSize', 12);
    xlabel('Fractions','FontSize', 12);
    sleg = cell(Nreplicates,1);
    for jj = 1:Nreplicates
      sleg{jj} = ['Replicate ' num2str(jj)];
    end
    legend(sleg, 'Orientation', 'horizontal','Location', 'South')
    
    set(gca,'YTickLabel',[],'ytick',[]);
    
    %generate image name
    Image_name = [figdir 'Bar_alignment_image.png'];
    
    %Save image
    saveas(Image, Image_name)
  end
catch
  warning('Failed to plot Bar_alignment_image.')
end



%% Delta_C scatter, rep vs rep_align

try
  for ci = 1:Nchannels
    figure
    RR = replicate_to_align_against(ci);
    
    for ri = 1:Nreplicates
      subplot(Nreplicates,1,ri),hold on
      plot([0 55],[0 55],'--r')
      overlap = intersect(summerised_names_G1{ci,ri},summerised_names_G1{ci,RR});
      try overlap([1 2]) = [];end
      Ia = find(ismember(Gaus_import{ci,ri}.textdata(:,1),overlap));
      Ib = find(ismember(Gaus_import{ci,RR}.textdata(:,1),overlap));
      x = Gaus_import{ci,ri}.data(Ia,2);
      y = Gaus_import{ci,RR}.data(Ib,2);
      scatter(x,y,20,abs(x-y),'filled')
      xlabel(['replicate ' num2str(ri)])
      ylabel('alignment replicate')
    end
    
    Image_name=[figdir 'Scatter_allreps_' Experimental_channels{ci} '.png'];
    set(gcf,'paperunits','inches','paperposition',[0.1 0.1 6 2+2*Nreplicates])
    saveas(gcf, Image_name, 'png');
  end
catch
  warning('Failed to plot Scatter_allreps_.')
end

