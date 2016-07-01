I = [1 2 3 9 10; 
    4 2 1 9 2;
    7 3 1 5 1;
    3 3 11 1 2;
    5 2 4 1 9;]
Fx=zeros(size(I),class(I));
Fy=zeros(size(I),class(I));
Fz=zeros(size(I),class(I));

J=zeros(size(I)+2,class(I));
J(:,:)=max(I(:));
J(2:end-1,2:end-1)=I;
Ne=[-1 -1; -1  0; -1  1; 0 -1; 0  1; 1 -1;  1  0; 1  1];
for i=1:length(Ne);
   In=J(2+Ne(i,1):end-1+Ne(i,1),2+Ne(i,2):end-1+Ne(i,2));
   size(In)
   check = In<I
   I(check)= In(check);
   D=Ne(i,:) 
   D=D./sqrt(sum(D.^2));
   Fx(check)= D(1) 
   Fy(check)= D(2);
end