function write_video(A, filename, filetype, options)

% write_video(A, filename, filetype, options)
%
%   write_video(A, filename, filetype, options);
%
% Copyright (c) 2016 Gabriel Peyre

options.null = 0;
switch filetype
    case 'gif'
        q = 128; % color map size
        res = @(A)reshape(A, [size(A,1) size(A,2) 1 size(A,3)]);
        imwrite(rescale(res(A),1,q), gray(q), [filename '.gif'], 'DelayTime',1/20, 'loopcount', Inf);
    case {'mp4' 'avi'}
        meth = 'Indexed AVI';
        meth = 'MPEG-4';
        vidObj = VideoWriter(filename, meth);
        vidObj.FrameRate = 100; 
        vidObj.Quality = getoptions(options, 'quality', 30);
        % vidObj.VideoFormat = getoptions(options, 'format', 'Grayscale');
        % vidObj.CompressionRatio = getoptions(options, 'compression', 20);
        open(vidObj);
        for i=1:size(A,3)
            writeVideo(vidObj,A(:,:,i));
        end
        close(vidObj);
end