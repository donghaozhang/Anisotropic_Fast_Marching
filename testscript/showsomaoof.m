load('mat\eigvoof.mat');
eigvoofone = eigvoof(:,:,:,1);
eigvooftwo = eigvoof(:,:,:,2);
eigvoofthree = eigvoof(:,:,:,3);
eigsum = eigvoofone + eigvooftwo + eigvoofthree;
I = eigsum < -1;
showbibox(I);
