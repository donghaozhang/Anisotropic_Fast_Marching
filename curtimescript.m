curtime = clock;
timelist = fix(curtime);
timestring = [];
for i = 1 : numel(timelist) 
    curstring = num2str(timelist(i));
    timestring = [timestring curstring];
end