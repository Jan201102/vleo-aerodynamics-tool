function lut = load_lut(lut_file)
    % LOAD_LUT Load a lookup table from a file in .csv format.
    % 
    % Parameters:
    %   lut_file (string): loouptable file path relative from the folder this function is called from.
    %
    % Returns:
    %   lut (struct): The loaded lookup table.

    [parent_folder,~,~] = fileparts(matlab.desktop.editor.getActiveFilename);
    display(parent_folder)
    lut_path = fullfile(parent_folder,lut_file);
    if ~isfile(lut_path)
        error("Look-up table file not found. Please check the path: %s", lut_path);
    end
    lut = readmatrix(lut_path);
end

    