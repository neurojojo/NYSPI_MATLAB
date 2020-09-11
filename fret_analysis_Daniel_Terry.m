mytable = table()
myfretobjects = dir('*.mat');

fret_calc = @(x,y,mask) mask.*(x./(x+y));

for experiment = 1:numel( myfretobjects )
    load( myfretobjects(experiment).name )
    
    imm_cnf = or( obj.fretTraces.Diff.Ch1.idlImb, obj.fretTraces.Diff.Ch1.idlCnf );
    imm_cnf = and( imm_cnf, obj.traceInterference>0.85 );
    
    fre = and( obj.fretTraces.Diff.Ch1.idlFre, obj.traceInterference>0.85 );
    fre_fret = fret_calc( obj.Ch1.adjusted_intensity, obj.Ch2.adjusted_intensity, fre );
    [ch1_mean,ch1_sd] = deal( nanmean(log(obj.Ch1.adjusted_intensity(obj.Ch1.adjusted_intensity>0))),...
        nanstd(log(obj.Ch1.adjusted_intensity(obj.Ch1.adjusted_intensity>0))) );
    
    %lower_bound = exp( ch1_mean - ch1_sd );
        
    for i = 1:obj.Ntracks
        locs = regexp( sprintf('%i',imm_cnf(i,:)), '1{1,}' );
        lens = cellfun( @numel, regexp( sprintf('%i',imm_cnf(i,:)), '1{1,}', 'match' ) );
        
        tmp = arrayfun( @(x) locs(x):locs(x)+lens(x)-1, [1:numel(locs)], 'UniformOutput', false )'
        
        ch1 = arrayfun( @(x) obj.Ch1.adjusted_intensity(i,x{1}) , tmp  , 'UniformOutput' , false);
        ch2 = arrayfun( @(x) obj.Ch2.adjusted_intensity(i,x{1}) , tmp  , 'UniformOutput' , false);

        mytable = [mytable; table( repmat( {obj.experimentName}, numel(locs), 1), ...
            repmat(str2double(obj.cellnum{1}),numel(locs),1), ...
            i*ones(numel(locs),1), locs', lens',ch1,ch2, ...
            'VariableNames', {'Experiment','Cell','Molecule','Frame','Length','Ch1_Intensity','Ch2_Intensity'}) ]
    end    
end

%%


figure; hist( cell2mat( arrayfun(@(x) mytable(x,:).Ch1_Intensity{1}./(mytable(x,:).Ch1_Intensity{1}+mytable(x,:).Ch2_Intensity{1}),...
    [1:size(mytable,1)],'UniformOutput',false) ), 50 )