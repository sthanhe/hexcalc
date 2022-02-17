function val=tppAir(prop,T)

persistent dryAir
if isempty(dryAir)
    dryAirStruct=load('dryAirTable.mat');
    dryAir=dryAirStruct.tab;
end


try
    propVal=dryAir{:,prop};
catch ME
    if strcmp(ME.identifier,'MATLAB:table:UnrecognizedVarName')
        varName=strsplit(ME.message);
        varName=varName{end};
        error(['Unknown or unsupported physical property: ',varName]);
    else
        rethrow(ME);
    end
end


val=interp1(dryAir.T,propVal,T);

end