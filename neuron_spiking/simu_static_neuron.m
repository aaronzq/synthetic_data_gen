%% generate 3D time sequence of synthetic spiking neurons
%% ==============user-defined parameters starts============================= 
%  dimension of the volume cube
height = 150;                         % px
width  = 150;                         % px
depth  = 76;                          % px
voxel_size   = [1.81, 4];             % um, [xy, z]

%  shape of the neuron
bead_size    = 10;                     % um
aniso_factor = 2;                     % z elongation factor

% spike signals, tune the fire rate and sampling rate for your prefered
% spike train
fr = 0.1;
samplingRate = 5;
tDecay = 1.2;
signalNoise = 0;                   % 0~1, normally distributed noise added to signal

%  statistics of the sequence
n_beads = 75;                        % the number of neurons
n_time = 6000;                         % the length of the sequence
Max_intensity = 100;                 % (in 16 bit)
addGaussianNoise = true;                  % add image gaussian noise at final step
gaussian_mean = 0.003;
gaussian_var = 0;

savePath = './Results';
saveLabel = 'Neuron75';
dirLabel = 'Label';
dirDynamic = 'Dynamic';
%% ==============user-defined parameters ends============================= 

if ~exist(savePath, 'dir')
    mkdir(savePath);
end
if ~exist(fullfile(savePath,saveLabel), 'dir')
    mkdir(fullfile(savePath,saveLabel));
end
if ~exist(fullfile(savePath,saveLabel,dirLabel), 'dir')
    mkdir(fullfile(savePath,saveLabel,dirLabel));
end
if ~exist(fullfile(savePath,saveLabel,dirDynamic), 'dir')
    mkdir(fullfile(savePath,saveLabel,dirDynamic));
end

bead_size = [bead_size, bead_size * aniso_factor];
bead_size = bead_size ./ voxel_size;
px_xy = floor(bead_size(1) + 0.5);
px_z = floor(bead_size(2) + 0.5);
bead_unit_inten = generateBead3D([px_xy, px_xy, px_z] * 2 + 1, [px_xy, px_z] / 4, 1);         % an ellipsoid neuron

bead_label_inten = bead_unit_inten;                                                           
bead_label_inten(bead_label_inten>1/2.718) = 1; bead_label_inten(bead_label_inten<=1/2.718) = 0;      % a mask for the neuron

beads_posi = rand([n_beads, 3]);                                                              % random positions of the neurons
beads_posi = beads_posi * 0.8 + 0.10;                                                         % constraint the neuron distribution in the inner volume
beads_posi = [beads_posi(:, 1) * height, beads_posi(:, 2) * width, beads_posi(:, 3) * depth]; % (height width depth)


figure;
scatter3(beads_posi(:,1)*voxel_size(1), beads_posi(:,2)*voxel_size(1), beads_posi(:,3)*voxel_size(2), 'filled'); 
axis equal
set(gca, 'Projection', 'orthographic');
set(gca, 'GridColor', 'k');
set(gca, 'lineWidth', 1);
set(gca, 'GridAlpha', 0.3);
set(gcf,'Color','w');
grid on
title('The distribution of the synthetic neurons');
    
[spikeMat, tVec] = poissonSpikeGen(fr, n_time, n_beads, samplingRate);     % generate signals from poisson process
figure;
plotRaster(spikeMat, tVec);
xlabel('time(s)')
ylabel('neuron id')
title('The spike signals of the synthetic neurons');
decayKernel = expDecay(tDecay, samplingRate);                              % decay kernel to emulate the calicum signal

bead_dynamic_inten = zeros(n_time, n_beads);
for n=1:size(spikeMat,1)
    bead_dynamic_inten(:,n) = conv(double(spikeMat(n,:)), decayKernel, 'same');
    noise = signalNoise*randn(n_time,1);
    bead_dynamic_inten(:,n) = abs(bead_dynamic_inten(:,n) + noise);
    bead_dynamic_inten(:,n) = bead_dynamic_inten(:,n) * Max_intensity;
end


% generate label for neuron, each neuron is labeled with a different pixel
% intensity
volume_label = zeros(height, width, depth);
for n = 1 : n_beads
    x = floor(beads_posi(n, 2));
    y = floor(beads_posi(n, 1));
    z = floor(beads_posi(n, 3));
    brightness = n;
    bead_label = bead_label_inten * brightness;
    radius = [px_xy, px_z];
    %% place the center of the bead to point(y, x, z)
    for i = y - radius(1) : y + radius(1)
        for j = x - radius(1) : x + radius(1)
            for k = z - radius(2) : z + radius(2)
                if i > 0 && j > 0 && k > 0 && i <= height && j <= width && k <= depth
                    volume_label(i, j ,k) = max(volume_label(i, j ,k), bead_label(i - (y - radius(1)) + 1,  j - (x - radius(1)) + 1, k - (z - radius(2)) + 1));
                end
            end
        end
    end

end
volume_label = uint16(volume_label);
write3d(volume_label, fullfile(savePath,saveLabel,dirLabel, 'neuron_label.tif'), 16);

% generate the volume sequence
for t = 1:n_time
    t
    volume_dynamic = zeros(height, width, depth);
    for n = 1 : n_beads
        x = floor(beads_posi(n, 2));
        y = floor(beads_posi(n, 1));
        z = floor(beads_posi(n, 3));

        brightness =   bead_dynamic_inten (t,n);
        bead_dynamic = bead_unit_inten * brightness;

        radius = [px_xy, px_z];
        %% place the center of the bead to point(y, x, z)
        for i = y - radius(1) : y + radius(1)
            for j = x - radius(1) : x + radius(1)
                for k = z - radius(2) : z + radius(2)
                    if i > 0 && j > 0 && k > 0 && i <= height && j <= width && k <= depth
                        volume_dynamic(i, j ,k) = volume_dynamic(i, j ,k) + bead_dynamic(i - (y - radius(1)) + 1,  j - (x - radius(1)) + 1, k - (z - radius(2)) + 1);
                    end
                end
            end
        end

    end
    
    volume_dynamic = uint16(volume_dynamic);
    if addGaussianNoise
        volume_dynamic_noise = imnoise(volume_dynamic, 'gaussian', gaussian_mean, gaussian_var);   % add more noise to the volume
    else
        volume_dynamic_noise = volume_dynamic;
    end
    
    write3d(volume_dynamic_noise, fullfile(savePath,saveLabel, dirDynamic, sprintf('neuron_dynamic_%03d.tif', t)), 16);
    
    
end

save(fullfile(savePath,saveLabel,'neuron_signal_gt.mat'), 'bead_dynamic_inten', 'beads_posi', '-v7.3');
disp('Done.');
