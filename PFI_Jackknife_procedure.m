% PFI_Jackknife_procedure
% confirm latency differences for F1 and F2 (during PFI).


% clear all; 
close all;
%%
%need to load the PFI data:
% copy-paste of BP_SSVEP script (PFI and Catch).

%background dirs
% addpath('/Users/MattDavidson/Desktop/SSVEP-Wills')
% cd('/Users/matthewDavidson/Desktop/SSVEP-feedbackproject/AA_ANALYZED DATA')
cd('/Users/mDavidson/Desktop/SSVEP-feedbackproject/AA_ANALYZED DATA')

basefol=pwd;
clearvars -except basefol allppants LOO*
dbstop if error
getelocs
cd(basefol)
cd('EEG')
cd('GFX_EEG-RESS')
%
%     cd('newplots-MD')
colsare={'b' , 'k','b','k','b','g','m'}; % blue for tg, black for BG.
hzare={'TG (f1)' , 'BG (f2)','TG (2f1)','BG (2f2)','TG (3f1)','TGBG','IM (f2-f1)'}; % blue for tg, black for BG.
linetypes={'--', '-', [],[],'--'};
legendis=[];
lgc=1;
%     clf
figure(1); hold on;
clearvars yl;

useMiller=0;
useinterp=1;

%% first as sanity check, plot the group results
[Overall_rawsamp, OverallSigdiff]=deal([]);

% all positive direction, for (disap-reapp) SNR increases., compares significant points to this reference,

hzdir = [1 1 1]; 
% to avoid retaining spurious clusters.
for hzis=1:3 %for each frequency
    
    switch hzis
        case 1
            load('GFX_PFIperformance_withSNR_15_min0_RESS')
            %rename for next part of script
            storeacrossPpant_onsetSNR15=storeacrossPpant_onsetSNR;
            storeacrossPpant_offsetSNR15=storeacrossPpant_offsetSNR;
            col='b';
        case 2
            
            load('GFX_PFIperformance_withSNR_20_min0_RESS')
            %rename for next part of script
            storeacrossPpant_onsetSNR20=storeacrossPpant_onsetSNR;
            storeacrossPpant_offsetSNR20=storeacrossPpant_offsetSNR;
            col='k';
            
             case 3
            
            load('GFX_PFIperformance_withSNR_5_min0_RESS')
            %rename for next part of script
            storeacrossPpant_onsetSNR5=storeacrossPpant_onsetSNR;
            storeacrossPpant_offsetSNR5=storeacrossPpant_offsetSNR;
            col='m';
    end
    
    %%        Jackknife leave one out
    
    useONSETs = squeeze(nanmean(storeacrossPpant_onsetSNR,2)); %participant averages.
    useOFFSETs = squeeze(nanmean(storeacrossPpant_offsetSNR,2));
    
    ttestdata=[];
    if useinterp==1
        % resample to estimate the latencies at ms level:
        
        %find nearest SNR/ closest time point:
        %interpolate to increase resolution and find closest point (ms)
        %
        samps = timeidDYN(1):1/1000:timeidDYN(end);
        %
        data1= pchip(timeidDYN, useONSETs,samps);
        data2= pchip(timeidDYN, useOFFSETs,samps);
        
    else
        data1=useONSETs;
        data2=useOFFSETs;
    end
    
    %which method?
    if useMiller==1 % calculate using 0.5*peak.
        tmp_Monset = squeeze(mean(data1,1));
        tmp_Moffset = squeeze(mean(data2,1));
        %zero this trace:
%         tmp_use = tmp_Monset - tmp_Monset(1);
        %or use difference\
        tmp_use = tmp_Monset-tmp_Moffset;
        tmp_use = tmp_use - tmp_use(1);
        peakat = 0.9* max(tmp_use); % 90% criterion for response locked ERPs (see Miller et al, 1998 psychophys).
        
        
        % now find closest point
        diffV= abs(tmp_use-peakat);
        if useinterp==1
        [~,Xmark] = find(diffV<.005);
         %we want minimum crossover point:
        Xmarkmin=min(Xmark);
        else
            Xmarkmin = dsearchn(diffV(1:7)', min(diffV));
        end
        
       
        if useinterp==1
            STC = samps(Xmarkmin);
        else
            STC= timeidDYN(Xmarkmin);
        end
        OverallSigdiff(hzis) = STC;
        
        Overall_rawsamp(hzis) = Xmarkmin;
    else
        % calculate first sig time point:
        
        %use orig data.
        ttestdata(1,:,:)= data1;
        ttestdata(2,:,:)= data2;
        
        %              %check for sig
        pvals=zeros(1,size(ttestdata,3));
        tvals=zeros(1,size(ttestdata,3));
        for itime = 1:size(ttestdata,3)
            
            [h,pvals(itime),~,stat]=ttest(ttestdata(1,:,itime), ttestdata(2,:,itime));
                        
            tvals(itime)= stat.tstat;
            if pvals(itime) < .05 && tvals(itime) <0 && hzdir(hzis)>0 %incorrect 
                pvals(itime) = nan;
            elseif pvals(itime) < .05 &&  tvals(itime) >0 && hzdir(hzis)<0 %incorrect
                pvals(itime)=nan;
            end
        end
        sigs=find(pvals<.05);
        %find max cluster.
        % check for cluster of sigs:
        vect1 = diff(sigs);
        v1 = (vect1(:)==1);
        d = diff(v1);
        if ~isempty(d)
            clusterSTandEND= [find([v1(1);d]==1) find([d;-v1(end)]==-1)];
            [~,maxClust] = max(clusterSTandEND(:,2)-clusterSTandEND(:,1));
            
            % first point in largest cluster
            minSIG=sigs(clusterSTandEND(maxClust,1));
            if useinterp==1
                STC=samps(minSIG);
            else
                STC=timeidDYN(minSIG);
            end
        else
            STC =nan; % no sig cluster this iteration
        end
                 OverallSigdiff(hzis) = STC;
        
    end
    
    
    
    
    %%%%%% Sanity check plotting
    % new mean SNR:
    tmp_Monset = squeeze(nanmean(data1,1));
    tmp_Moffset = squeeze(nanmean(data2,1));
    
    if useMiller~=1
        figure(1);
        %             subplot(211)
        if useinterp~=1
        plot(timeidDYN, tmp_Monset, col); hold on;
        plot(timeidDYN, tmp_Moffset, col); hold on;
        else
            plot(samps, tmp_Monset, col); hold on;
        plot(samps, tmp_Moffset, col); hold on;
        end
        plot([STC, STC], ylim, ['-'], 'color', 'r');
        title(['min sig at ' num2str(STC)])
    else
        figure(1);
        if useinterp==1
            plot(samps, tmp_use, 'k')
            
        else
            plot(timeidDYN, tmp_use, 'k')
        end
        hold on;
        plot([STC,STC], ylim, ['k:']);
        
    end
    
end
%% Now repeat the analysis in a Leave one out subsample.


    %preallocate times....
    [LOO_firstSIG,LOO_firstrawsamp]=deal(zeros(size(useONSETs,1),2));
    
    
    nppantsA=1:size(useONSETs,1); % use all participants
    
    hzdir = [.001, .001, .001]; % all positive for (disap-reap) SNR.
    for iJK = 1:length(nppantsA)
        
        %remove this participant for averaging.
        useppants= nppantsA;
        useppants(iJK)=[];
        
        
        for hzis=1:3 %for each frequency
            
            %             if iJK==1 % just load the first time.
            switch hzis
                case 1
                    storeacrossPpant_onsetSNR=storeacrossPpant_onsetSNR15;
                    storeacrossPpant_offsetSNR=storeacrossPpant_offsetSNR15;
                    col='b';
                case 2
                    
                    storeacrossPpant_onsetSNR=storeacrossPpant_onsetSNR20;
                    storeacrossPpant_offsetSNR=storeacrossPpant_offsetSNR20;
                    col='k';
                     case 3
                       storeacrossPpant_onsetSNR=storeacrossPpant_onsetSNR5;
                    storeacrossPpant_offsetSNR=storeacrossPpant_offsetSNR5;
                    col='m';
                    
                    
            end
            
            %%        Jackknife leave one out
            %
            % need to compute cross-over points of average timecourses, using leave-one-out procedure.
            
            
            useONSETs = squeeze(nanmean(storeacrossPpant_onsetSNR,2)); %participant averages.
            useOFFSETs = squeeze(nanmean(storeacrossPpant_offsetSNR,2));
            
            %Compute first significat point for this subset of ppants:
            
            subuseONSETs= squeeze(useONSETs(useppants,:));
            subuseOFFSETs= squeeze(useOFFSETs(useppants,:));
            %
            if useinterp==1
                % resample to estimate the latencies at ms level:
                
                %find nearest SNR/ closest time point:
                %interpolate to increase resolution and find closest point (ms)
                %
                samps = timeidDYN(1):1/1000:timeidDYN(end);
                %
                data1= pchip(timeidDYN, subuseONSETs,samps);
                data2= pchip(timeidDYN, subuseOFFSETs,samps);
                xlabis = samps;
            else
                data1=subuseONSETs;
                data2=subuseOFFSETs;
                xlabis = timeidDYN;
            end
            
            %which method?
            if useMiller==1 % calculate using 0.5*peak.
                
                tmp_Monset = squeeze(mean(data1,1));
                tmp_Moffset = squeeze(mean(data2,1));
                %zero this trace:
%                 tmp_use = tmp_Monset - tmp_Monset(1);
                tmp_use = tmp_Monset - tmp_Moffset;
                tmp_use=tmp_use - tmp_use(1);
                
                peakat = 0.9* max(tmp_use);
                
                % now find closest point
                diffV= abs(tmp_use-peakat);
                if useinterp==1
                    [~,Xmark] = find(diffV<.005);
                    %we want minimum crossover point:
                    Xmarkmin=min(Xmark);
                else
                    Xmarkmin = dsearchn(diffV(1:7)', min(diffV));
                end
                
                
                if useinterp==1
                    STC = samps(Xmarkmin);
                    xlabis = samps;
                else
                    STC= timeidDYN(Xmarkmin);
                    xlabis=timeidDYN;
                end
                
                LOO_firstSIG(iJK,hzis)= STC;
                LOO_firstrawsamp(iJK, hzis)=Xmarkmin;
            else
                % calculate first sig time point:
                ttestdata=[];
                %use orig data.
                ttestdata(1,:,:)= data1;
                ttestdata(2,:,:)= data2;
                
                %              %check for sig
                pvals=zeros(1,size(ttestdata,3));
                tvals=zeros(1,size(ttestdata,3));
                for itime = 1:size(ttestdata,3)
                    
                    [h,pvals(itime),~,stat]=ttest(ttestdata(1,:,itime), ttestdata(2,:,itime));
                    
                    
                    tvals(itime)= stat.tstat;
                          
                    % make sure we only keep significant changes in the
                    % right direction, to avoid spurious early clusters.
                    if pvals(itime)<.05 &&  tvals(itime)< hzdir(hzis); % wrong direction.
                        pvals(itime) = nan;
                    end
                end
                sigs=find(pvals<.05);
                %find max cluster.
                % check for cluster of sigs:
                vect1 = diff(sigs);
                v1 = (vect1(:)==1);
                d = diff(v1);
               if ~isempty(d)
                clusterSTandEND= [find([v1(1);d]==1) find([d;-v1(end)]==-1)];
                [~,maxClust] = max(clusterSTandEND(:,2)-clusterSTandEND(:,1));
                
                % first point in largest cluster
                minSIG=sigs(clusterSTandEND(maxClust,1));
                if useinterp==1
                    STC=samps(minSIG);
                else
                    STC=timeidDYN(minSIG);
                end
                else
                    STC =nan; % no sig cluster this iteration
                end
                LOO_firstSIG(iJK,hzis)= STC;
            end
            
            %
            
            %% plot also the mean SNR, and first SIG for sanity check.
            figure(2); subplot(4,4, iJK)
            if useMiller==1
                plot(xlabis, tmp_use)
                hold on;
                plot([STC,STC], ylim, 'color', col);
            else
                plot(xlabis, squeeze(mean(data1,1))); hold on
                plot(xlabis, squeeze(mean(data2,1)));
                plot([STC, STC], ylim, ['-'], 'color', 'r');
            end
            %%
            
        end
                title(['diff is ' num2str(LOO_firstSIG(iJK,2)-LOO_firstSIG(iJK,1))]);
        
        
        
    end
    %% SUMMARY STATS, if using Miller Method
    if useMiller==1
    %From Miller et al., Psychophysiology, 1998. The l-o-o variability
    %can be used to estimate standard error, using eqn (2) of their paper:
    
    %se= sqrt((n-1/n) * sum((J - K)^2)); where:
    % J = estimated value in subsample
    % K = average of values in subsample
    % N=subsample size
    
    N=size(LOO_firstSIG,1);
    % we are interested in the difference in our sample beween
    % f1 and f2, when slopes reached 0.5 max amp. This difference:
    
       
%     J= LOO_firstSIG(:,2)-LOO_firstSIG(:,1); %  F2 and F1.
% D= diff(OverallSigdiff);  %calculated  above    

J= LOO_firstrawsamp(:,1)-LOO_firstrawsamp(:,2); %  F2 and F1.
D= Overall_rawsamp(1)-Overall_rawsamp(2)    
    K= mean(J);
    
    StErrHzfromJK =  sqrt(((N-1)/N)  *   sum((J-K).^2));
    
    %note that using this se estimate, we can compute inferential
    %statistics.
    
    %eqn 3 from Miller et al., 1998;
    %jack-knifed t-stat = D/se;
    
    %D=difference in overall sample (true difference in
    %sig), se is computed above.
    
    
    tstatJK= D/StErrHzfromJK;

    pval = tcdf(tstatJK, N-1)
    end
    
      %% plot ditribution of results
    
    
    figure(3); clf
    set(gcf, 'units', 'normalized', 'position', [0 .35 .35 .6])
    subplot(211)
    usecols = [0,0,1; 0,0,0; 1,0,1]; % b'k'm'.
    legh=[];
    for ihz = 1:3
    h1=raincloud_plot(LOO_firstSIG(:,ihz),'box_on', 1, 'color', usecols(ihz,:), 'alpha', 0.5,...
        'box_dodge', 1, 'box_dodge_amount', .15, 'dot_dodge_amount', .15,...
        'box_col_match', 1);
    %increase scatter size:
    h1{2}.SizeData=20;
    h1{2}.LineWidth=5;
    legh(ihz) = h1{1};
    end
    
    axis tight
    set(gca, 'fontsize', 20, 'ytick', []);
    xlim([-1.5 0])
    xlabel('time from button press [s]')    
    title({['First significant SNR points during PFI,'];['Jackknife (N-1) subsamples']})
    legend([legh(1) legh(2) legh(3)], 'Target (f1)', 'Surround (f2)', 'IM (f2-f1)', 'autoupdate', 'off')
    set(gcf, 'color', 'w')
    xlim([-1.5 1.5])
    hold on; plot([0 0 ], ylim, ':', 'linew', 3, 'col', 'k')
    %
    subplot(212);
    h1=raincloud_plot((LOO_firstSIG(:,2)-LOO_firstSIG(:,1)),'box_on', 1, 'color', [1,0,0], 'alpha', 0.5,...
        'box_dodge', 1, 'box_dodge_amount', .15, 'dot_dodge_amount', .15,...
        'box_col_match', 0);
    %increase scatter size:
    h1{2}.SizeData=100;
    
    h1{2}.LineWidth=2;
    
    title('Latency difference (f2-f1), each subsample')
    axis tight
    set(gca, 'fontsize', 20, 'ytick', []);
    xlim([-1.5 0])
    hold on; plot([0 0], ylim)
    set(gca, 'fontsize', 20);
    xlabel('seconds')
    xlim([-1.5 1.5])
 hold on; plot([0 0 ], ylim, ':', 'linew', 3, 'col', 'k')

 

%% fig 2. 
% prevBoxplotfigure; % for back compatability.
figure(6); clf % arrange columns as f1, IM , f2 (1,3,2)
set(gcf, 'units', 'normalized', 'position', [0 .55 .4 .4]);
tmp = LOO_firstSIG;
tmp(:,2) = LOO_firstSIG(:,3); %LOO 3 was ihz. LOO 2 was f2.
tmp(:,3) = LOO_firstSIG(:,2);

SCcol = {'b', 'm', 'k'};
nppants = size(LOO_firstSIG,1);
keepscatter_x= []; % for linking the individual data points afterwards.
keepscatter_y= []; % for linking the individual data points afterwards.
%

for ibox = 1:3%:3%1:3 %f1,im,f2. 

    
    %define median, and interquartile range.   
    % overlay individual data points. 
     useD =squeeze(tmp(:,ibox));
     %plot Median using a solid line.
     medval = nanmedian(useD);     
     

%% plot an arrow instead.
startp = [ibox+.35, medval]; % x,y
endp = [ibox-.35, medval]; % x,y
pa=plot_arrow(ibox-.35, medval, ibox+.35, medval, ...
    'linewidth', 3,'Color', SCcol{ibox},...
    'facecolor', SCcol{ibox}, 'edgecolor', SCcol{ibox}, 'headwidth', .15);


 hold on
 %%
 % add STD patch.
 pstd =nanstd(useD);
 xv = [ibox-.15;ibox-.15; ibox+.15; ibox+.15];
 yv = [medval-pstd; medval+pstd;medval+pstd; medval-pstd];
%  cv = [0; 0; 1; 1]; %shows increasing gradient (disap > reap);
 ptch= patch(xv, yv, SCcol{ibox});
 ptch.EdgeColor = 'w';%SCcol{ibox};
 ptch.FaceAlpha= .1;
 %add boxplot lines.
 %width
 plot([ibox, ibox], [medval-pstd,medval+pstd], 'linew', 1, 'color', SCcol{ibox});
 %whiskers
 plot([ibox-.15, ibox+.15], [medval-pstd,medval-pstd], 'linew', 2, 'color', SCcol{ibox});
 plot([ibox-.15, ibox+.15], [medval+pstd,medval+pstd], 'linew', 2, 'color', SCcol{ibox});
     %% now adjust scatter points, to offset slightly.     
     useppants = sum(~isnan(useD)); % adjust so that median bar is centred.
     ph = repmat(ibox, [useppants,1]); % y position
     offsets = linspace(-.4, .4, useppants);
    if ibox==2
        offsets = linspace(-.2, .2, useppants);
    end
         ph = ph+ offsets';
     
     % before plotting, reorder chronologically (by order at f2)
%      tmp2=LOO_firstSIG(:,2); 
%      [~, sID] = sort(tmp2,'ascend');
     [~, sID] = sort(useD,'ascend');
     chron = useD(sID);
     
     
     
%    accomodate for the IM (ibox =2), having missing values:     
     if ibox==2 % this is the IM, so continue only with the relevant points
     % fill the new order with correct placements.
     pht=ph;
     ind = ~isnan(chron); % find timepoints with an IM
     nant = nan(1,length(chron)); % pre array
     nant(ind) = ph(1:sum(ind)); % fill new array with timings, where SNR exists
     ph = nant; % continue
     end
     
     %% plot with and without markers, to show no IM cases.
        
     %using the correct order (sID), determine which entries have IM and
     %which do not, to change scatter appearance.
     
     IM = LOO_firstSIG(sID,3);
     noIM = find(isnan(IM));
     wIM = find(~isnan(IM));
     
     
 sc=scatter(ph(wIM), chron(wIM)); % with IM case
 sc.LineWidth = 1;
 sc.SizeData = 25;
% %  sc.MarkerFaceColor = SCcol{ibox};
  sc.MarkerFaceColor = [1,1,1];
 sc.MarkerEdgeColor = SCcol{ibox};
 
 hold on
%%     
 sc=scatter(ph(noIM), chron(noIM)); % no IM case
 sc.LineWidth = 1;
 sc.SizeData = 25;
%  sc.MarkerFaceColor = [1,1,1];
 sc.MarkerFaceColor = SCcol{ibox};
 sc.MarkerEdgeColor = SCcol{ibox};
% 
%  hold on
%  sc=scatter(ph(wIM), chron(wIM)); % with IM case
%  sc.LineWidth = 1;
%  sc.SizeData = 25;
%  sc.MarkerFaceColor = [1 1 1];
%  sc.MarkerEdgeColor = SCcol{ibox};
%  
%
 hold on % we need to reorganise to correctly link with individual lines.
 %at the moment, the yaxis position is reset per freq. we want to reset to
 %the order shown, starting with order of f2.
      tmp2=LOO_firstSIG(:,2); % f2 order
     [~, fID] = sort(tmp2,'ascend');
 keepscatter_y(:,ibox) = useD(fID); % we now have matched the SNR together.
 % to retrieve the correct ycoords per SNR point, we need to revert from
 % chronological order.
 % revert:
 [~,rev]=sort(sID, 'ascend');
 %fix ph
 phn = ph(rev);
 %now can store, as the ph matched to SNR value
  keepscatter_x(:,ibox) = phn(fID);
end
 %%
 set(gca, 'Xtick', 1:3, 'Xticklabels', {'Target (f1)' 'IM (f2-f1)', 'Surround (f2)'})
 hold on; plot(xlim, [0, 0], ['k',':'], 'linew', 2);
 %%
 % add connecting lines between f2 and f1.
for isubs = 1:size(keepscatter_x,1)    
    hold on
    %these xvals
    xvals = keepscatter_x(isubs,:); % row.
    yvals = keepscatter_y(isubs,:); % row.
    
    %sort by timing?
      if any(isnan(yvals))
          xvals=xvals([1,3]);
          yvals=yvals([1,3]);
          linest = '-';
          linew=1;
      else
          linest = '-';
          linew=1;
      end
      
   
  
    line(xvals,yvals, 'Color', [.8 .8 .8],'linestyle', linest);%, 'linew', linew);
    

    
end
%
view([90 -90]);
    axis ij    
    set(gca, 'ydir', 'normal', 'fontsize', 25);
 shg
 ylabel('time from button press (s)')
 ylim([-1.5 1.5]);
 set(gcf, 'color', 'w')
 

