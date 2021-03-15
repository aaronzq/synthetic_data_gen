%% generate 3D tubulins (don't know if it exists :))
%% ==============user-defined parameters starts=============================

imgNum = 2;                   % the number of generated images

gardenSize = [500, 500, 80];   % a garden to grow tubulins (x, y and z, um)
cropSize = [300,300, 60];      % crop a sub-region of the garden
brightness = 60;               % (uint8)

tubulinNum = 60;               % number of tubulins
walkLength = 1800;             % how many step of the growing process
dk_lat = 0.02/2.5;             % stiffness: higher, more flexible (lateral)
dk_axi = 0.01/2.5;             % stiffness: higher, more flexible (axial)
v0 = 1;                        % grow speed (um/step): higher, faster

%% ==============user-defined parameters ends=============================

nTry = 1;
while nTry  <= imgNum
    try
        trajectory = zeros(walkLength+1, 3, tubulinNum);
        % grow some tubulins from x axis
        for i = 1 : round(tubulinNum/3)
            x0 = [gardenSize(1)*rand(1) 0 0];
            k0_lat = 0.25 + 0.05*(rand(1)-0.5)*2;
            k0_axi = 0.25 + 0.05*(rand(1)-0.5)*2;
            kmin_lat = 0.2;
            kmax_lat = 0.3;
            kmin_axi = 0.2;
            kmax_axi = 0.3;
            [trajectory(:,1,i), trajectory(:,2,i), trajectory(:,3,i)] = ...
                gen_tubulin3d(walkLength, dk_lat, dk_axi, k0_lat, k0_axi, kmin_lat, kmax_lat,...
                kmin_axi, kmax_axi, v0, x0);
        end
        % grow some tubulins from y axis
        for i = round(tubulinNum/3)+1 : 2*round(tubulinNum/3)
            x0 = [0 gardenSize(2)*rand(1) 0];
            k0_lat = 0.05*(rand(1)-0.5)*2;
            k0_axi = 0.25 + 0.05*(rand(1)-0.5)*2;
            kmin_lat = -0.1;
            kmax_lat = 0.1;
            kmin_axi = 0.2;
            kmax_axi = 0.3;
            [trajectory(:,1,i), trajectory(:,2,i), trajectory(:,3,i)] = ...
                gen_tubulin3d(walkLength, dk_lat, dk_axi, k0_lat, k0_axi, kmin_lat, kmax_lat,...
                kmin_axi, kmax_axi, v0, x0);
        end
        % grow some tubulins from z axis
        for i = 2*round(tubulinNum/3)+1 : tubulinNum
            x0 = [0 0 gardenSize(3)*rand(1)];
            k0_lat = 0.125 + 0.05*(rand(1)-0.5)*2;
            k0_axi = 0.25 + 0.05*(rand(1)-0.5)*2;
            kmin_lat = 0;
            kmax_lat = 0.2;
            kmin_axi = 0.2;
            kmax_axi = 0.3;
            [trajectory(:,1,i), trajectory(:,2,i), trajectory(:,3,i)] = ...
                gen_tubulin3d(walkLength, dk_lat, dk_axi, k0_lat, k0_axi, kmin_lat, kmax_lat,...
                kmin_axi, kmax_axi, v0, x0);    
        end
    %     figure;
    %     for i = 1 : tubulinNum    
    %         scatter3(trajectory(:,1,i), trajectory(:,2,i), trajectory(:,3,i), 10, 'k','fill'); hold on;
    %     end
    %     hold off;
    %     axis equal;
        img = zeros(cropSize);
        for i = 1 : tubulinNum
            for p = 1 : walkLength
                tX = trajectory(p, 1, i);
                tY = trajectory(p, 2, i);
                tZ = trajectory(p, 3, i);

                if ( (tX > (gardenSize(1)-cropSize(1))/2 + 1) && (tX < (gardenSize(1)+cropSize(1))/2 - 1) && ...
                        (tY > (gardenSize(2)-cropSize(2))/2 + 1) && (tY < (gardenSize(2)+cropSize(2))/2 - 1) && ...
                        (tZ > (gardenSize(3)-cropSize(3))/2 + 1) && (tZ < (gardenSize(3)+cropSize(3))/2 - 1) )
                    iX = tX - (gardenSize(1)-cropSize(1))/2 ;
                    iY = tY - (gardenSize(2)-cropSize(2))/2 ;
                    iZ = tZ - (gardenSize(3)-cropSize(3))/2 ;            
                    img(floor(iY):ceil(iY), floor(iX):ceil(iX), floor(iZ):ceil(iZ)) = ...
                        img(floor(iY):ceil(iY), floor(iX):ceil(iX), floor(iZ):ceil(iZ)) + brightness;
                end
            end
        end
        out = flipud(img);
        out = imgaussfilt3(out,1.2);
        out = uint8(out);
        out_noise = imnoise(out,'gaussian',1/255,0.000005);

        saveName = sprintf('./tubulins3d%d.tif', nTry);
        saveName2 = sprintf('./tubulins3d_noise%d.tif', nTry);
        write3d(out, saveName);
        write3d(out_noise, saveName2);
        
        disp(['Synthetic Tubulins Image ' num2str(nTry) ' ...done.']);  
        nTry = nTry + 1;
    catch
%         disp(['Synthetic Tubulins Image ' num2str(nTry) ' ...skipped.']);  
    end
end
disp('Done.');  

