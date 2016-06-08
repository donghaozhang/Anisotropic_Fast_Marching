%  varargout=ifftshiftedcoormatrix(dimension)
% The dimension is a vector specifying the size of the returned coordinate
% matrices. The number of output argument is equals to the dimensionality
% of the vector "dimension". All the dimension is starting from "1"
function varargout=ifftshiftedcoormatrix(dimension)
dim=length(dimension);
p = floor(dimension/2);

    for i=1:dim
        a=([p(i)+1:dimension(i) 1:p(i)])-p(i)-1;
        reshapepara=ones(1,dim);
        reshapepara(i)=dimension(i);
        A=reshape(a, reshapepara);
        repmatpara=dimension;
        repmatpara(i)=1;
        varargout{i}=repmat(A, repmatpara);
    end
end