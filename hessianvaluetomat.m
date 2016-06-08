function hessianmat = hessianvaluetomat(vecI)
    % convert [d11 d12 d13 d22 d23 d33] into matrix
    % result matix is [d11 d12 d13
    %                  d12 d22 d23
    %                  d31 d32 d33] 
    %
    hessianmat = zeros(3,3);
    hessianmat = [vecI(1),vecI(2),vecI(3);vecI(2),vecI(4),vecI(5);vecI(3),vecI(5),vecI(6);];
end