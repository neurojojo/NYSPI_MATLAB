function output = moveImageData( originalDir )

    channel_dirs = dir( sprintf('%s\\#*Ch2',originalDir) );
    
    for i = 1:numel(channel_dirs)
        
        lookhere = fullfile(channel_dirs(i).folder,channel_dirs(i).name);
        % files and directories to move %
        all_files = dir( sprintf('%s\\**\\**', lookhere ) );
        
        % Create the directory
        mkdir( regexprep( lookhere, 'F:\\', 'C:\\Users\\meszaro\\Documents\\' ) );
        % Add imageData directory as well
        mkdir( sprintf('%s\\ImageData\\',regexprep( lookhere, 'F:\\', 'C:\\Users\\meszaro\\Documents\\' )) );
        
        % Avoid traces files at all cost and all directories that might
        % contain them
        not_traces_files = find(arrayfun( @(x) and( ~x.isdir, isempty(strfind( x.name, 'traces' )) ), all_files )==1);
        
        arrayfun( @(x) copyfile( fullfile( all_files(x).folder, all_files(x).name ), ...
                    regexprep( fullfile( all_files(x).folder, all_files(x).name ), 'F:\\', 'C:\\Users\\meszaro\\Documents\\' ) ), ...
                    not_traces_files )
        
        output = arrayfun( @(x) {sprintf( '%s', fullfile( all_files(x).folder, all_files(x).name )),...
            regexprep( fullfile( all_files(x).folder, all_files(x).name ), 'F:\\', 'C:\\Users\\meszaro\\Documents\\' )}, ...
                    not_traces_files, 'UniformOutput', false );
    end
    
end