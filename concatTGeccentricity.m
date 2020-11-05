cd('/Users/MDavidson/Desktop/SSVEP-feedbackproject/AA_ANALYZED DATA/Behaviour')
load('Raw_Behavioral_AllPP.mat', 'AllTrialsAllPPBottomLeft'); 
%%
targX=[];
ntrials=48;
nppants = size(AllTrialsAllPPBottomLeft,1)/ntrials;
for ippant = 1:nppants
    
    targX(ippant) = AllTrialsAllPPBottomLeft(ippant*48,12);
end
%%
%note that the target eccentricity is actually th hypotenuse, of value
%extracted.
targXn = round(sqrt(2*(targX.*targX)));
disp([[1:nppants]',targX', targXn'])
%%
figure(1);
clf

%convert to DVA:
%steps (in pixels, along diagonal, were roughly:
pixels =[141, 156, 170, 184, 198, 212, 226, 240,255,269];
%dva =

% FOR MBI EEG DOWNSTAIRS:
display.dist=50;
display.width = 51;
display.resolution = 1920;
TGdva= zeros(size(targXn));
tmpd=[];
for ippant = 1:nppants
pix = targXn(ippant); %pixel distance, diagonally from centre to centre of TGs.

TGdva(ippant)=pix2angle(display,pix);
% keep the total dva range, for fixing the histogram below:
if ippant<length(pixels)
tmpd(ippant) = pix2angle(display, pixels(ippant));
end
end
%%
tmpd=round(tmpd,1);
fixrange = 4.3:.4:9;
clf
TGdva = round(TGdva, 1);
%clean up for plot
fixval = find(~ismember(TGdva, fixrange));
%%
for ifix=fixval
    badv = TGdva(ifix);
    %find nearest for replacement.
    repl= nearest(fixrange, badv);
    TGdva(ifix) = fixrange(repl);
end

%%
clf
 xvec=[3.9:0.4:8.7];
hg=histogram(TGdva, length(xvec));
hg.FaceColor = 'b';
hold on
hg.BinWidth = .4;
set(gca, 'xtick',xvec);
%%
tmp=TGdva;
tmp(:) = nan;
tmp([3,5,8]) = TGdva(1,[3,5,8]);
%add removed ppants in different colour.
hr = histogram(tmp, ceil(nppants));
hr.EdgeColor='r';
hr.FaceColor = hg.FaceColor;
hr.LineWidth = 4;
hr.BinWidth = hg.BinWidth;
hr.BinLimits= hg.BinLimits;
ylim([.01 10])


shg

%%
xlabel({['Calibrated target eccentricity'];['from center (\circ visual angle)']});
ylabel('Participant count')
set(gcf, 'color', 'w');
set(gca, 'fontsize', 25);
xlim([3.5 9]);
set(gca, 'xtick',4.3:.8:8.7);
legend(hr, 'excluded from sample');
% %%
% hg=histogram(round(TGdva,2), 25);
% hg.FaceColor = 'b';
% xlabel({['Calibrated target eccentricity (d.v.a.)']});
% ylabel('# participants')
% set(gcf, 'color', 'w');
% %
% set(gca, 'fontsize', 25, 'Xminortick', 'on','xtick', 3:.5:7);
% shg
%%
cd ../Figures
cd('Behavioural Bar Plots')
%%
print('-dpng', 'Target eccentricity per participant')