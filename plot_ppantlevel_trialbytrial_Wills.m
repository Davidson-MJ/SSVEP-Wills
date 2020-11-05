


clear all
dbstop if error
cd('/Users/MDavidson/Desktop/SSVEP-feedbackproject/AA_ANALYZED DATA/EEG')
%
dirs =dir([pwd filesep '*_*' 'EEG']);
basefol=pwd;


%

basefol=pwd;

counter=1;
close all
for ifol=[9 16];
        for ihz=5%:6; %TGs1 TGs 2, 20hz 40 hz.
            switch ihz                                
                case 5
                    usehz=20;
                    iloc=1;
                case 6
                    usehz=40;
                    iloc=1; %doesnt matter, always the same RESS filter.
            end
            
            
            
            icounter=1;
            
            cd(basefol)
            cd(dirs(ifol).name)
            
            load(['ppant_PFI_Epoched_RESS'])
            
            for itimezero = 1%:2
                if itimezero==1
                    
                    %append them all. TARG-> Disappearing (more buttons
                    %pressed).

                    datatouse = cat(1, ress_PFI_0_1_Hz(ihz).ressTS,ress_PFI_1_2_Hz(ihz).ressTS,ress_PFI_2_3_Hz(ihz).ressTS, ress_PFI_3_4_Hz(ihz).ressTS);
%                     RTstouse = [durs0_1'; durs1_2'; durs2_3'];
                    BPstouse = cat(1, ress_PFI_0_1_Hz(ihz).BPs,ress_PFI_1_2_Hz(ihz).BPs,ress_PFI_2_3_Hz(ihz).BPs,ress_PFI_3_4_Hz(ihz).BPs);
                    ctype = 'PFI increase';
                    bsrem = [-3 -1]; %seconds
                    
                else
                    
                    datatouse = cat(1, ress_PFI_1_0_Hz(ihz).ressTS,ress_PFI_2_1_Hz(ihz).ressTS,ress_PFI_3_2_Hz(ihz).ressTS,ress_PFI_4_3_Hz(ihz).ressTS);
%                     RTstouse = [durs1_0'; durs2_1'; durs3_2'];
                    BPstouse= cat(1, ress_PFI_1_0_Hz(ihz).BPs,ress_PFI_2_1_Hz(ihz).BPs,ress_PFI_3_2_Hz(ihz).BPs,ress_PFI_4_3_Hz(ihz).BPs);
                    ctype = 'PFI decrease';
                    
                    
                    bsrem = [1 3]; %seconds
                    
                    
                end
                
                
                
                %plot the topo for sanity check: (pre disap)
                windcheck= [-3 -0.1];
%                 tidx=dsearchn(timeidDYN', [windcheck]');
                
                
                
                %restrict to certain PFI duration?
                allt= 1:size(datatouse,1);
                %                 baddurs = find(RTstouse<30);
                %remove those trials.
                %                 allt(baddurs)=[];
                
                
                
               
                
                %%
                %%
                c=colormap('viridis');
                %space out colormap accordingly (5 colours for this data)
                indices = ceil(linspace(1, size(c,1),5));
                % round up
                
                %first row is lowest color
                colorsblock=[c(indices(1),:); ...
                    c(indices(2),:); ...
                    c(indices(3),:); ...
                    c(indices(4),:);
                    c(indices(5),:)];
                
                
                
%                  colorsblock=[0.9932    0.9062    0.1439];
                colormap(colorsblock)
%                 c=colormap('viridis')
                shg
                %%
                    
                %plot BP for comparison
                figure(1);
                if counter==1
                subplot(2,2, [1 3])
                else
                    subplot(2,2, [2])
                end
                %sort trial index by longest RT first.
                %                 [sortedRTs, cid] = sort(RTstouse(allt), 'descend');
                
                %sort by accum PFI in window after onset.
                if itimezero==1
                    accumPeriod = sum(BPstouse(:,180:end),2);
                else
                    accumPeriod = sum(BPstouse(:,1:180),2);
                    
                end
                
                [checkBP, cid] = sort(accumPeriod, 'descend');
                
                
                % MOVE nan to bottom of order (not top)
                nanind= find(isnan(checkBP));
                %new start point
                cid=[cid(length(nanind)+1:end) ; cid(nanind)];
                %                     checkBP=[checkBP(length(nanind)+1:end); checkBP(nanind)]
                
                %rearrange.
                BPstouse=BPstouse(cid,:);
                datais=squeeze(datatouse(cid,:));
%                 RTstouse=RTstouse(cid);
shg
                %

                tBP=-3:1/60:3;
                imagesc(-3:1/60:3,1:length(cid),BPstouse)
                if ifol==16
                                c=colorbar;
                ylabel(c, 'buttons pressed')
                set(c, 'Ytick', [0:4], 'Location', 'SouthOutside')
                %change location.
                end
                %                 set(gca, 'ytick', 1:length(cid), 'yticklabel', round(sortedRTs./60,2))
                ylabel('PFI events')
                xlabel('Time [sec]')
                hold on
                plot([0 0] , ylim, ['w:'], 'linewidth', 3)
                title({['Participant ' num2str(ifol)]})
                
                set(gca, 'fontsize', 25, 'Ytick', [0:20:size(BPstouse,1)])
            end
        end
        shg
        set(gcf, 'color', 'w');
        figure(2)
        subplot(1,2,counter)
        colormap('viridis')
       
        %plot smoothed ver
         % smooth across trials
                smoothedBP=zeros(size(BPstouse));
                for itime=1:size(BPstouse,2)
                    smoothedBP(:,itime)= smooth(BPstouse(:,itime),15);
                end
                
         imagesc(-3:1/60:3,1:length(cid),smoothedBP)
                %                 set(gca, 'ytick', 1:length(cid), 'yticklabel', round(sortedRTs./60,2))
                ylabel({['Normalized'];['trial count']})
                xlabel('Time [sec]')
                hold on
                plot([0 0] , ylim, ['w:'], 'linewidth', 3)
%                 title({['Subject ' num2str(ifol)]})
                
                set(gca, 'fontsize', 25, 'Ytick', [])
                set(gcf, 'color', 'w')
                 if ifol==16
            
       
                c=colorbar;
                ylabel(c, 'buttons pressed')
                set(c, 'Ytick', [0:4])
        end
                counter=counter+1;
end
%%
set(gcf, 'color', 'w')