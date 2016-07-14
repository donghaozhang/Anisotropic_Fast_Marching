	clc; clear all; close all;
    trXF = zeros([8, 3]);
    trXF = [ [0,0,1];[0, 0, -1];[0,0,1];[0,0,-1];[0,0,1];[0, 0, -1];[0,0,1];[0,0,-1] ];
    trYF = zeros([8, 3]);
	trYF = [[1,0,0];[1,0,0];[1,0,0];[1,0,0];[-1,0,0];[-1,0,0];[-1,0,0];[-1,0,0] ];
	trZF = zeros([8, 3]);
    trZF = [[0,1,0];[0,1,0];[0,-1,0];[0,-1,0];[0,1,0];[0,1,0];[0,-1,0];[0,-1,0]];
	% float a[3], b[3], c[3], dir1, temp;
    trX = zeros([8, 3]); trY = zeros([8, 3]); trZ = zeros([8, 3]);
for k = 1 : 8	
		a(1) = -1*trXF(k,1);
		a(2) = -1*trYF(k,1);
		a(3) = -1*trZF(k,1);
      
		b(1) = -1*trXF(k,2);
		b(2) = -1*trYF(k,2);
		b(3) = -1*trZF(k,2);
      
		c(1) = -1*trXF(k,3);
		c(2) = -1*trYF(k,3);
		c(3) = -1*trZF(k,3);

		dir1 = c(1)*(a(2)*b(3) - a(3)*b(2)) - c(2)*(a(1)*b(3)-a(3)*b(1)) + c(3)*(a(1)*b(2)-a(2)*b(1));
        if (dir1 > 0)
            temp = a(1);
			a(1) = b(1);
			b(1) = temp;
			temp = a(2);
			a(2) = b(2);
			b(2) = temp;
			temp = a(3);
			a(3) = b(3);
			b(3) = temp;
			trX(k,1) = -a(1);
			trY(k,1) = -a(2);
			trZ(k,1) = -a(3);
			trX(k,2) = -b(1);
			trY(k,2) = -b(2);
			trZ(k,2) = -b(3);
			trX(k,3) = -c(1);
			trY(k,3) = -c(2);
			trZ(k,3) = -c(3);
        else
            trX(k,1) = -a(1);
			trY(k,1) = -a(2);
			trZ(k,1) = -a(3);
			trX(k,2) = -b(1);
			trY(k,2) = -b(2);
			trZ(k,2) = -b(3);
			trX(k,3) = -c(1);
			trY(k,3) = -c(2);
			trZ(k,3) = -c(3);
        end
end
X1=[];Y1=[];Z1=[];
ColList='ymcrgbk'; %no w letter because you won't see the white line 
pcol = [255,0,0];
hold on
for k = 1 : 8
    col=[rand, rand, rand];
    X1=[];Y1=[];Z1=[];
    trX_tmp = trX(k, :);
    trY_tmp = trY(k, :);
    trZ_tmp = trZ(k, :);
    X1 = [X1; trX_tmp(1); trY_tmp(1); trZ_tmp(1); trX_tmp(1)];
    Y1 = [Y1; trX_tmp(2); trY_tmp(2); trZ_tmp(2); trX_tmp(2)];
    Z1 = [Z1; trX_tmp(3); trY_tmp(3); trZ_tmp(3); trX_tmp(3)];
    plot3(X1, Y1, Z1, 'color', col)
    pause(2)
end 

hold off



% Points = rand(50,3);
% points = [1 0 0; 0 1 0; 0 0 1; -1 0 0];
% dt = delaunayTriangulation(trX);
% faceColor  = [0.6875 0.8750 0.8984];
% figure
% tetramesh(dt,'FaceColor',faceColor,'FaceAlpha',0.3);