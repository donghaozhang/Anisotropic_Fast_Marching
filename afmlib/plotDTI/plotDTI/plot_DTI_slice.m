function plot_DTI_slice(D, zslice, bI)
    delta = 1;
    nx = size(D, 1);
    ny = size(D, 2);
    nz = size(D, 3);
    figure
    hold on
    constrain_z = constrain(zslice, 1, nz);
    for i = 1 : nx
        for j = 1 : ny
            for k = constrain_z
                if bI(i,j,k) == 1
                    T_vec = squeeze(D(i,j,k,:));
                    hessianmat = hessianvaluetomat(T_vec);
                    [v,l] = eig(hessianmat);
                    [X,Y,Z] = ellipsoid(0,0,0,l(1,1),l(2,2),l(3,3),10);
                    sz = size(X);
                    for x=1:sz(1)
                        for y = 1 : sz(2)
                            A = [X(x,y) Y(x,y) Z(x,y)]';
                            A = v * A;
                            X(x,y) = A(1); Y(x,y) = A(2); Z(x,y) = A(3);
                        end
                    end
                    X = X + (i-1) * delta * 2;
                    Y = Y + (j-1) * delta * 2;
                    Z = Z + (k-1) * delta * 2;
                    h(i) = surf(X,Y,Z);
                    % fprintf('x: %d, y: %d, z: %d\n', i, j, k);
                end
            end
        end
    end   
    axis equal
    view([0 90]);
    set(gca,'GridLineStyle','none')
    set(gca,'XTick',[])
    set(gca,'YTick',[])
    set(gca,'ZTick',[])
    shading interp
    colormap([0.8 0.8 0.8])
    lighting phong
    light('Position',[0 0 1],'Style','infinite','Color',[ 1.000 0.584 0.000]);
    hold off
end