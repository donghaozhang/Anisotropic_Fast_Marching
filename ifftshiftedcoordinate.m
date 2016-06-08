function output=ifftshiftedcoordinate(dimension, dimindex, pixelspacing)
    dim=length(dimension);
    p = floor(dimension/2);
    a=single([p(dimindex)+1:dimension(dimindex) 1:p(dimindex)])-p(dimindex)-1;
    a=a/pixelspacing(dimindex)/dimension(dimindex);
    reshapepara=ones(1,dim, 'single');
    reshapepara(dimindex)=dimension(dimindex);
    
    a=reshape(a, reshapepara);
    repmatpara=dimension;
    repmatpara(dimindex)=1;
    output=repmat(a, repmatpara);
end