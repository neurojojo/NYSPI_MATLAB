%APO (0 Glutamate):
dir_apo = 'F:\\Disk E\\#BackupMicData\\2017 May-Jun\\170511_WA-MH_SNAPf-mGlu2-Cy3AC-Cy5AC_w&wo Glu_Live FRET\\333 nM Cy3AC_666 nM Cy5AC_PassX_100 ngmL tet_50 PCA PCD\\Track_Opt28-TW5_52418\\LT20FreeDiff 82719\\'

%15 uM Glutamate:
glu_15_dir = 'F:\\Disk E\\#BackupMicData\\2017 Jul-Dec\\170901_WA_SNAPfast-mGlu2_Cy3AC-Cy5AC_Live FRET glutamate\\333 nM Cy3AC_666 nM Cy5AC_Pass2_100 ngmL tet_50 PCA PCD_15 uM L-Glutamate\\Track_Opt28-TW5_080919\\LT20\\'

%100 uM Glutamate:
glu_100_dir = 'F:\\Disk E\\#BackupMicData\\2017 Jul-Dec\\170816_WA_SNAPf-mGlu2_Cy3AC-Cy5AC_Live FRET\\333 nM Cy3AC_555 nM Cy5AC_Pass14_100 ngmL tet_50 PCA PCD_100 uM L-Glu\\Track_Opt28-TW5_111518\\LT20\\'

%Antagonist (10 uM LY341495):
antagonist_dir = 'F:\\Disk E\\#BackupMicData\\2017 Jul-Dec\\170816_WA_SNAPf-mGlu2_Cy3AC-Cy5AC_Live FRET\\333 nM Cy3AC_555 nM Cy5AC_Pass14_100 ngmL tet_50 PCA PCD_10 uM LY341495\\Track_Opt28-TW5_082519\\LT20\\'

%% My analysis will compute background and remove interferences

for N = 1:6;
    myfret1 = fretAnalysisObject(sprintf('%s\\#0%iCh2\\',dir_apo,N),'Apo');
    myfret1.calculateBackground();
    myfret1.getInterferences(); 
    myfret1.calculateFret('free',1);
    myfret1.calculateFret('immobile',1);
    myfret1.saveFret();
    myfret1={};
end
close all


for N = 1:6;
    myfret1 = fretAnalysisObject(sprintf('%s\\#0%iCh2\\',glu_15_dir,N),'Glu_15');
    myfret1.calculateBackground();
    myfret1.getInterferences(); 
    myfret1.calculateFret('free',1);
    myfret1.calculateFret('immobile',1);
    myfret1.saveFret();
    myfret1={};
end
close all


for N = 1:6;
    myfret1 = fretAnalysisObject(sprintf('%s\\#0%iCh2\\',glu_100_dir,N),'Glu_100');
    myfret1.calculateBackground();
    myfret1.getInterferences(); 
    myfret1.calculateFret('free',1);
    myfret1.calculateFret('immobile',1);
    myfret1.saveFret();
    myfret1={};
end
close all


for N = 1:6;
    myfret1 = fretAnalysisObject(sprintf('%s\\#0%iCh2\\',antagonist_dir,N),'Antagonist');
    myfret1.calculateBackground();
    myfret1.getInterferences(); 
    myfret1.calculateFret('free',1);
    myfret1.calculateFret('immobile',1);
    myfret1.saveFret();
    myfret1={};
end
close all
%% Figure

alldata = [cellfun( @(x) mean( x.fretTraces.Fret.Calculated_free ), myfret1, 'ErrorHandler', @(x,y) nan )';...
cellfun( @(x) mean( x.fretTraces.Fret.Calculated_free ), myfret2, 'ErrorHandler', @(x,y) nan )';...
cellfun( @(x) mean( x.fretTraces.Fret.Calculated_free ), myfret3, 'ErrorHandler', @(x,y) nan )';...
cellfun( @(x) mean( x.fretTraces.Fret.Calculated_free ), myfret4, 'ErrorHandler', @(x,y) nan )'];
lbls = [ 1*ones(6,1); 3*ones(6,1); 2*ones(6,1); 4*ones(6,1) ];

figure('color','w');
h=boxplot( alldata, lbls ); hold on;

scatter( lbls, alldata , 64, 'Jitter', 'on', 'JitterAmount', 0.1, 'MarkerEdgeColor','k', 'MarkerFaceColor', 'm', 'MarkerFaceAlpha', 0.5 )


all_lines = findobj(ax,'Type','Line');arrayfun( @(x) set(x,'LineStyle','-','Color','k','LineWidth',1), all_lines )
all_lines = findobj(gca,'Type','Line');arrayfun( @(x) set(x,'LineStyle','-','Color','k','LineWidth',1), all_lines )
myboxes = findobj(gca,'Tag','Box')arrayfun( @(box) patch( box.XData, box.YData, 'm', 'FaceAlpha', 0.5), myboxes(1:5) )
myboxes = findobj(gca,'Tag','Box'); arrayfun( @(box) patch( box.XData, box.YData, 'm', 'FaceAlpha', 0.5), myboxes(1:5) )
%% FREE

thisfretset = myfret4;

collect_Fret_values = arrayfun( @(N) thisfretset{N}.fretTraces.Fret.Calculated_free', [1:6] , 'UniformOutput', false, 'ErrorHandler', @(x,y) nan);

bins = [0:0.025:1];

results = cell2mat( cellfun( @(cell_) (histc( cell_, bins )./numel(cell_))', collect_Fret_values, 'UniformOutput', false ) )
figure('color','w'); 


counts = histc(cell2mat(collect_Fret_values),bins);
b=bar( bins, mean(results,2) ); hold on;

title( sprintf('%s - All cells (free)',regexprep(thisfretset{N}.experimentName,'_',' ')), 'color', 'k', 'fontweight', 'normal', 'fontsize', 24); box off;

hold on;
mygm = fitgmdist( cell2mat(collect_Fret_values)',1 );
counts = pdf( mygm, bins' ) * mean(diff(bins));
line( bins', counts, 'color', 'k', 'linewidth', 3 )
errorbar( bins, mean(results,2), std(results,[],2)/sqrt(6), 'linewidth', 2 )
set(gca,'color',[1,1,1],'xcolor',[0,0,0],'ycolor',[0,0,0],'linewidth',4,'fontsize',18,'TickDir','out');
saveas(gcf, sprintf('%s_summary_free.png',thisfretset{N}.experimentName) )

%IMMOBILE

thisfretset = myfret4;

collect_Fret_values = arrayfun( @(N) thisfretset{N}.fretTraces.Fret.Calculated_immobile', [1:6] , 'UniformOutput', false, 'ErrorHandler', @(x,y) nan);

bins = [0:0.025:1];

results = cell2mat( cellfun( @(cell_) (histc( cell_, bins )./numel(cell_))', collect_Fret_values, 'UniformOutput', false ) )
figure('color','w'); 


counts = histc(cell2mat(collect_Fret_values),bins);
b=bar( bins, mean(results,2) ); hold on;

title( sprintf('%s - All cells (immobile)',regexprep(thisfretset{N}.experimentName,'_',' ')), 'color', 'k', 'fontweight', 'normal', 'fontsize', 24); box off;

hold on;
mygm = fitgmdist( cell2mat(collect_Fret_values)',1 );
counts = pdf( mygm, bins' ) * mean(diff(bins));
line( bins', counts, 'color', 'k', 'linewidth', 3 )
errorbar( bins, mean(results,2), std(results,[],2)/sqrt(6), 'linewidth', 2 )
set(gca,'color',[1,1,1],'xcolor',[0,0,0],'ycolor',[0,0,0],'linewidth',4,'fontsize',18,'TickDir','out');
saveas(gcf, sprintf('%s_summary_immobile.png',thisfretset{N}.experimentName) )



%%
experiment=1

ch1_bleach_first = arrayfun( @(x) (myfret{experiment}.fretTraces.Ch2.traceMetadata(x).endOfTrace-myfret{experiment}.fretTraces.Ch1.traceMetadata(x).endOfTrace)>10, [1:myfret{experiment}.Ntracks] );
ch1_post_bleach = arrayfun( @(N) nanmedian( myfret{experiment}.fretTraces.Ch1.int_clean(N,[myfret{experiment}.fretTraces.Ch1.traceMetadata(N).endOfTrace:myfret{experiment}.fretTraces.Ch1.traceMetadata(N).endOfTrace+100]) ), [1:myfret{experiment}.Ntracks], 'ErrorHandler', @(x,y) nan ) 
ch2_post_bleach = arrayfun( @(N) nanmedian( myfret{experiment}.fretTraces.Ch2.int_clean(N,[myfret{experiment}.fretTraces.Ch1.traceMetadata(N).endOfTrace:myfret{experiment}.fretTraces.Ch1.traceMetadata(N).endOfTrace+100]) ), [1:myfret{experiment}.Ntracks], 'ErrorHandler', @(x,y) nan ) 

mytbl = table( [1:numel(ch1_post_bleach)]', ch1_bleach_first', ch1_post_bleach', ch2_post_bleach', sign((ch2_post_bleach'-ch1_post_bleach')), ((ch2_post_bleach'-ch1_post_bleach')), ...
    'VariableNames', {'Idx','Ch1_first','Ch1_pb','Ch2_pb','Ch2_over_Ch1','Ch2_sub_Ch1'})
mytbl = sortrows(mytbl,'Ch2_over_Ch1'); 
mytbl_sub = mytbl( mytbl.Ch1_first==true, : );
mytbl_sub = mytbl_sub( and( mytbl_sub.Ch2_pb > 0, mytbl_sub.Ch2_over_Ch1 > 0) , : )
mytbl_sub = sortrows(mytbl_sub,'Ch2_sub_Ch1')


%%
N=233
figure; plot( myfret{experiment}.fretTraces.Ch1.int_clean(N,:) ); hold on; plot( myfret{experiment}.fretTraces.Ch2.int_clean(N,:) );

xlim([0,myfret{experiment}.fretTraces.Ch2.traceMetadata(N).endOfTrace+200])
title( sprintf('%i %i', myfret{experiment}.fretTraces.Ch1.traceMetadata(N).endOfTrace, myfret{experiment}.fretTraces.Ch2.traceMetadata(N).endOfTrace) )
