function generate_input_sheet(number_Module,Module_name,filename,storage_location)
%Input
% number_Module = 7;
% Module_name = ["SH2" "RH2" "SH1" "RH1" "EVAP" "HP_ECO" "LP_ECO"];
% filename = 'C:\Users\tbrunaue\Documents\GitHub\HEX_recalc\input\1_Test_table.xlsx';

%function
filename = string(storage_location)+"\"+string(filename);
T = table('Size',[14,number_Module+2],'VariableTypes',['string' "string" repmat("double",1,number_Module)]);
T.Properties.VariableNames = ["variable" "unit" Module_name];
col = [{"da"; "s"; "K"; "nCells"; "mDot1"; "mDot3"; "nTube"; "l_Tube"; "l_cavity"; "l_bed"; "w_channel"; "T1_start"; "T3_start"; "p_start"; "Qdot_start"}...
{"m"; "m"; "m"; "-";  "kg/s"; "kg/s"; "-"; "m"; "m"; "m"; "m"; "K"; "K"; "Pa"; "W"}];
T(:,1:2) = col;

writetable(T,filename)

hExcel = actxserver('Excel.Application');
hExcel.Visible = 1;
hWorkbook = hExcel.Workbooks.Open(filename);
hWorksheet = hWorkbook.Sheets.Item(1);
% hWorksheet.Columns.Item(1).columnWidth = 10; %first column
hWorksheet.Columns.columnWidth = 12; %first column
hWorkbook.Save
release(hExcel.Workbooks)
delete(hExcel); 
end
