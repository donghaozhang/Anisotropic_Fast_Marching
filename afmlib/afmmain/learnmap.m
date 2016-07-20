keySet =   [3];
valueSet = [327.2];
mapObj = containers.Map(keySet,valueSet)
% mapObj = containers.Map
mapObj(1) = 100.0;
keys(mapObj)
mapObj(7) = 200.0;
keys(mapObj)
mapObj(5) = 50.0;
keys(mapObj)
mapObj(2) = 540.0;
keys(mapObj)
values(mapObj)
mapObj(4) = 540.0;
keys(mapObj)
values(mapObj)
remove(mapObj, 5);
keysmapObj = keys(mapObj)
firstkey = keysmapObj{1}
remove(mapObj, firstkey); 
keys(mapObj)
values(mapObj)





