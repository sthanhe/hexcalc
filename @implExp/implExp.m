classdef implExp
    %CC-By Stefan Thanheiser
    methods(Static)
        function sz=size(varargin)
            dims=max(cellfun(@(x) ndims(x),varargin));
            szCell=cellfun(@(x) implExp.padsz(dims,x),varargin,'UniformOutput',false);
            szs=cell2mat(szCell');

            %all dimensions not equal to 1 must be the same
            not1=szs~=1;
            check=arrayfun(@(x) all(szs(not1(:,x),x)==szs(find(not1(:,x),1),x)),1:size(szs,2));

            if all(check)
                sz=max(szs,[],1);
                sz(any(szs==0,1))=0;
            else
                throw(MException('implExp:incompatibleArrays','Arrays have incompatible sizes'));
            end
        end


        function varargout=normalize(sz,varargin)
            n=prod(sz);
            if n~=0
                pad=@(x) implExp.padsz(numel(sz),x);
                varargout=cellfun(@(x) reshape(repmat(x,sz./pad(x)),1,n),varargin,'UniformOutput',false);
            else
                varargout=repmat({NaN(sz)},1,length(varargin));
            end
        end
    end
    
    
    methods(Static, Access=private)
        function padsz=padsz(dims,x)
            padsz=[size(x),ones(1,dims-ndims(x))];
        end
    end
end




