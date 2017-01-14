% (c) 2013-2015 Miika Aittala, Jaakko Lehtinen, Tim Weyrich, Aalto 
% University, University College London. This code is released under the 
% Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International 
% license (http://creativecommons.org/licenses/by-nc-sa/4.0/).

function res = tex_compute_transport_cuda(fdata, tpar)
    rng(123456);

    if nargin < 2
        tpar = struct;
        tpar.iters = 180;
        tpar.desc.w = 16;
    end

    res = struct;    
    
	tpar.desc.w = 2*tpar.desc.w;    
    img_guide = fdata.img_double.guide;
    img_flash = fdata.img_double.flash;

    H = size(img_guide,1);
    W = size(img_guide,2);
    S = [H W];
    M = S(1) ./ fdata.par.B(1);
    
    %% Form the gradient image and operators
    
    alpha = 96/M*1;
    
    Dx = conv2d_mat([-1 1], [M M], 'replicate'); 
    Dy = conv2d_mat([1;-1], [M M], 'replicate');

    
    %% CUDA setup
    blocksize = 512;
    dsor = 1;
    blocknum = M^2/blocksize;
    kernel_feature = ...
            parallel.gpu.CUDAKernel('cuda_feature.ptx','cuda_feature.cu');
	kernel_feature.ThreadBlockSize = [blocksize 1 1];
	kernel_feature.GridSize = [blocknum/dsor dsor 1];
    
    kernel_feature.ThreadBlockSize
    kernel_feature.GridSize
    
    %% Pre-processing...

    WW = tpar.desc.w;
    N = 8;
    
    img_guide_pad = padarray(img_guide, [WW WW], 'symmetric');
    GG = exp(-linspace(-2,2, ceil(min(W,H)/4)).^2);
    GG = GG / sum(GG);

    img_guide_pad_low = imfilter(img_guide_pad,GG,'same','symmetric');    
    img_guide_pad_low = imfilter(img_guide_pad_low,GG','same','symmetric');    
    img_guide_pad = img_guide_pad - img_guide_pad_low;

    G2 = fspecial('gaussian',[WW*16 WW*16], WW*4);
    img_sq = imfilter(abs(img_guide_pad).^2,GG,'same','symmetric');
    img_sq = imfilter(img_sq,GG','same','symmetric');
    std(img_guide_pad(:))
    img_guide_pad = img_guide_pad ./ sqrt(img_sq);
    clf
    std(img_guide_pad(:))
    
    img_guide_std = img_guide_pad(WW+1:end-WW, WW+1:end-WW, :);

    imagec((0.2*img_guide_pad+0.5));%-min(img_guide_pad(:))))
    
    res.img_guide_std = img_guide_std;
    
    %% Compute the descriptors
     

    desc = [];
    for c = 1:3
        img_guide_pad_ = img_guide_pad(:,:,c);
        desc = cat(3,  ...
            desc, ...
            brief_blur(img_guide_pad_, WW/2, 4, WW/8), ...
            brief_blur(img_guide_pad_, WW,   3, WW/4), ...  
            brief_blur(img_guide_pad_, 2,    1, 0));

        
    end
    desc = desc(WW+1:end-WW, WW+1:end-WW,:);
    
    fdata.feat.desc = desc;

    
    %% Cut out the master tile 
    % (notice we used to call them "blocks" prior writing the paper...)
    % (in fact, "main block" == "master tile" here)

    dataname = fdata.name;
	out_fn = strcat(fdata.path, 'out/transport.mat')

	moff = round(fdata.par.main_offs .* [H W])
    
    nblocks = [H/M W/M]
    
    desc_gpu = gpuArray(desc);
    mainblock_desc_gpu = desc_gpu(moff(1)+1:moff(1)+M, moff(2)+1:moff(2)+M, :);

    mainblock_desc = desc(moff(1)+1:moff(1)+M, moff(2)+1:moff(2)+M, :);

	mainblock_guide = img_guide(moff(1)+1:moff(1)+M, moff(2)+1:moff(2)+M, :);
	mainblock_guide_std = img_guide_std(moff(1)+1:moff(1)+M, moff(2)+1:moff(2)+M, :);
    mainblock_guide = (mainblock_guide);
 	mainblock_flash = img_flash(moff(1)+1:moff(1)+M, moff(2)+1:moff(2)+M, :);
    mainblock_flash = (mainblock_flash);
    
    %% Set up buffers etc
    
    res.transport = cell(nblocks(1), nblocks(2));
  
    res.recons1 = zeros(M,M,3,nblocks(1),nblocks(2));
    res.recons2 = zeros(M,M,3,nblocks(1),nblocks(2));
    res.mask = zeros(M,M,nblocks(1),nblocks(2));
    res.recons3 = zeros(M,M,3,nblocks(1),nblocks(2));
    res.recons4 = zeros(M,M,3,nblocks(1),nblocks(2));
    res.P_b2m = zeros(M,M,nblocks(1),nblocks(2));
    res.P_m2b = zeros(M,M,nblocks(1),nblocks(2));

    % GPU buffers
    g_mins_idx = zeros([blocknum+1 1], 'uint32', 'gpuArray');
    g_mins = zeros([blocknum 1], 'single', 'gpuArray');
    g_fdists = zeros([M^2 1], 'single', 'gpuArray');
    mins_idx = zeros([M^2 1], 'uint32');
    mins = zeros([M^2 1], 'single');

    
    % Gradient reconstruction setup
    A = [alpha*speye(M^2); Dx; Dy];
    ATA = A'*A;

    MIC = [];
    % michol sometimes fails (not sure if ever here), but just to
    % be sure...
    try 
        MIC = ichol(ATA,struct('michol','on'));
    catch err
        err
        MIC = speye(size(ATA,1));
    end
    
    
    for bj = 1:1:nblocks(2)
        for bi =1:1:nblocks(1)
            %figure
            tic
            COMPUTING_BLOCK = [bi bj nblocks]
            block_desc = desc((bi-1)*M+1:bi*M, (bj-1)*M+1:bj*M, :);
            block_desc_gpu = desc_gpu((bi-1)*M+1:bi*M, (bj-1)*M+1:bj*M, :);

            block_guide = img_guide((bi-1)*M+1:bi*M, (bj-1)*M+1:bj*M, :);
            block_guide_std = img_guide_std((bi-1)*M+1:bi*M, (bj-1)*M+1:bj*M, :);
         	block_flash = img_flash((bi-1)*M+1:bi*M, (bj-1)*M+1:bj*M, :);

            % Compute the transport plans from both master to target, and
            % vice versa
            
            tic
            % block->main
            g_P = zeros([M^2 1], 'uint32', 'gpuArray');
            for i = 1:M^2
                [g_P, g_mins_idx, g_mins, g_fdists] = feval(kernel_feature, g_P, g_mins_idx, g_mins, g_fdists, mainblock_desc_gpu, block_desc_gpu, mainblock_guide_std, block_guide_std, M^2, N, 0.01, i-1);  % Note: -1 in index for CUDA
                if mod(i, round(M^2/10)) == 0
%                    printf('.');
                end
            end                        
            T = speye(M^2);
            P_b2m = gather(g_P);
            T = T(gather(g_P),:);%
            toc
            tic
            % main->block
            g_P = zeros([M^2 1], 'uint32', 'gpuArray');
            for i = 1:M^2
                [g_P, g_mins_idx, g_mins, g_fdists] = feval(kernel_feature, g_P, g_mins_idx, g_mins, g_fdists, block_desc_gpu, mainblock_desc_gpu, block_guide_std, mainblock_guide_std, M^2, N, 0, i-1);  % Note: -1 in index for CUDA
                if mod(i, round(M^2/10)) == 0
       %             printf('.');
                end
            end                        
            T2 = speye(M^2);
            T2 = T2(gather(g_P),:)';%
            P_m2b = gather(g_P);
            toc
            
            % Done.
            % Now rearrange the samples according to the plan, and solve
            % the gradient reconstruction problem. Most of the code below
            % is old junk that is no longer really needed for anything.
            %
            % This all happens in a very messy piece of code where we we show some
            % stuff on screen but also compute some useful stuff for later
            % use...
            
            block_guide = (block_guide);
            block_flash = (block_flash);
            clf;

            disp('Poisson');

            colormap gray;
            
            subplot(4,3,1);
            %imagesc(mainblock_guide*64);
            imagec(mainblock_guide);
            
            subplot(4,3,2);

            imagec(block_guide);
            
            subplot(4,3,4);

            recon = reshape(unit_rows((T)) * reshape(block_guide, [], 3), size(block_guide));
            imagec(recon);
            
            subplot(4,3,7);
            recon = reshape(unit_rows(T) * reshape(block_flash, [], 3), size(block_flash));
            res.recons3(:,:,:,bi,bj) = recon;
                        
            recon = reshape(unit_rows(T) * reshape(remove_linear_log(block_flash), [], 3), size(block_flash));
            res.recons1_old(:,:,:,bi,bj) = recon;
            
            F = reshape(remove_linear_log(block_flash), [], 3);
            TFG = cat(1, alpha*T*F, T*Dx*F, T*Dy*F);

            b = TFG;
            ATb = A'*b;
            
            recon = zeros(M^2,3);
            for v = 1:3
                recon(:,v) = pcg(ATA, ATb(:,v), 1e-5,200,MIC,MIC');
            end
            recon = max(0,min(1,reshape(recon, [M M 3])));
            
            % This is the image that will be used in fitting later!
            res.recons1(:,:,:,bi,bj) = recon;

            imagec(recon);
            
            subplot(4,3,5);
            recon = reshape(unit_rows((T)') * reshape(mainblock_guide, [], 3), size(mainblock_guide));
            imagec(recon);

            subplot(4,3,8);
            recon = reshape(unit_rows((T)') * reshape(remove_linear_log(mainblock_flash), [], 3), size(mainblock_flash));
            res.recons2(:,:,:,bi,bj) = recon;
            imagec(recon);%(:,:,2));

            subplot(4,3,12);
            imagesc(reshape(sum(T,1), [M M]));
            res.mask(:,:,bi,bj) = reshape(sum(T,1), [M M]);
            

            subplot(4,3,10);
            
            imagec(remove_linear_log(mainblock_flash));
            subplot(4,3,11);
            
            imagec(remove_linear_log(block_flash));
            
            s = size(block_guide);
            map = zeros(size(block_guide,1), size(block_guide,2), 3);
            TT = sparse(size(T,1), size(T,2));

            recon_guide = zeros(size(block_flash));
            recon_flash = zeros(size(block_flash));
            
            
            for j = 1:size(T2,2)
                t = T2(:,j);

                sj = floor((j-1)/s(1));
                si = j - sj*s(1);
                sj = sj+1;
                %si = si+1;

                tidx = find(t == max(t),1);
                TT(tidx,j) = 1;
                tj = floor((tidx-1)/s(1));
                ti = tidx - tj*s(1);
                tj = tj+1;

                map(si,sj,1) = ti;
                map(si,sj,2) = tj;

                map(si,sj,1) = ti-si;
                map(si,sj,2) = tj-sj;

                recon_guide(si,sj,:) = mainblock_guide(ti,tj,:);
                recon_flash(si,sj,:) = mainblock_flash(ti,tj,:);

            end

            subplot(4,3,3);
            imagec(map/s(1)+0.5);
            
            MB = reshape(remove_linear_log(mainblock_flash), [], 3);
            BF = reshape(remove_linear_log(block_flash), [], 3);
            recon_guide = reshape(MB(P_m2b,:), size(mainblock_guide));
            recon_flash = reshape(BF(P_b2m,:), size(mainblock_guide));
            res.recons4(:,:,:,bi,bj) = recon_flash;
            
            % Save the transport maps for reverse transport (we do it here
            % already as all the necessary setup is done, even though it 
            % will only be used much later)
            res.P_b2m(:,:,bi,bj) = reshape(P_b2m, [M M]);
            res.P_m2b(:,:,bi,bj) = reshape(P_m2b, [M M]);
            
            
            subplot(4,3,6);
            imagec(recon_guide);    
            
            subplot(4,3,9);
            imagec(recon_flash);    
            
            drawnow;
        end
        
    end
	save(out_fn,'res');

end
