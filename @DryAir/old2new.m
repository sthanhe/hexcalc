vars={'rho','h','s','c_p','c_v','beta','w','lambda','eta','ny','a','Pr'};

tables=cell(1,length(vars));
for i=1:length(vars)
    tab=readtable('dryAirOld.xls','Sheet',vars{i},'Range','A:Y');

    isstr=arrayfun(@(x) isa(tab{1,x},'cell'),1:size(tab,2));
    tab{:,isstr}=strrep(tab{:,isstr},'','-');
    tab{:,isstr}=strrep(tab{:,isstr},',','.');

    tabnew=table('Size',size(tab),'VariableTypes',repmat({'double'},1,size(tab,2)));
    tabnew(:,~isstr)=tab(:,~isstr);
    tabnew{:,isstr}=str2double(tab{:,isstr});

    T=tabnew{1,2:end}+273;
    T=arrayfun(@(x) num2str(T(x)),1:length(T),'UniformOutput',false);
    T=cellfun(@(x,y) [x,y],repmat({'T'},1,length(T)),T,'UniformOutput',false);
    tabnew.Properties.VariableNames=[{'p'},T];

    tabnew(1,:)=[];

    writetable(tabnew,'dryAir.xls','WriteMode','overwritesheet','Sheet',vars{i});
end