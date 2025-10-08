function data_tab = generate_input_tab(number_Module, Module_name, filename, storage_location)
% number_Module,Module_name,filename
% number_Module = 7;
% Module_name = ["SH2" "RH2" "SH1" "RH1" "EVAP" "HP_ECO" "LP_ECO"];
% filename = 'C:\Users\tbrunaue\Documents\GitHub\HEX_recalc\input\1_Test_table.xlsx';

% addpath("C:\Users\tbrunaue\Documents\GitHub\HEX_recalc\input")
filename =string(storage_location)+"\"+string(filename);

dataLines = [2, Inf];

%% Set up the Import Options and import the data
opts = spreadsheetImportOptions("NumVariables", number_Module+2);

% Specify sheet and range
% opts.Sheet = sheetName;
opts.DataRange = dataLines(1, :);

% Specify column names and types
opts.VariableNames = ["variable", "unit", Module_name];
opts.SelectedVariableNames = ["variable", "unit", Module_name];
opts.VariableTypes = ["string", "string", repmat("double",1,number_Module)];

% Specify variable properties
% opts = setvaropts(opts, ["Var1", "Var2"], "WhitespaceRule", "preserve");
% opts = setvaropts(opts, ["Var1", "Var2"], "EmptyFieldRule", "auto");

% Import the data
data_tab = readtable(filename, opts, "UseExcel", false,ReadRowNames=true);

for idx = 2:size(dataLines, 1)
    opts.DataRange = dataLines(idx, :);
    tb = readtable(workbookFile, opts, "UseExcel", false);
    Testtable = [Testtable; tb]; %#ok<AGROW>
end

end