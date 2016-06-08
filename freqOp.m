function result=freqOp(freq, marginwidth)
    result=freq(marginwidth(1)+1:size(freq,1)-marginwidth(1), marginwidth(2)+1:size(freq,2)-marginwidth(2), marginwidth(3)+1:size(freq,3)-marginwidth(3));
end