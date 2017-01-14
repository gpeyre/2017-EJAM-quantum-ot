% (c) 2013-2015 Miika Aittala, Jaakko Lehtinen, Tim Weyrich, Aalto 
% University, University College London. This code is released under the 
% Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International 
% license (http://creativecommons.org/licenses/by-nc-sa/4.0/).

function tex_render_rotate(l, g, mul, fname)
    writerObj = [];
    
    if nargin < 4
        writerObj = VideoWriter('rotate.avi');
    else
        writerObj = VideoWriter(fname);
    end
    writerObj.Quality = 95;
    open(writerObj);
    
    S = size(l);
    nv = S(3);
    S = S(1:2);
    
    asp = S(2)/S(1);
    
    C0 = cat(3, -repmat(linspace(-1,1,S(1))', [1 S(2)]), repmat(linspace(-asp,asp,S(2)), [S(1) 1]));

    l = reshape(l, prod(S(1:2)), nv)';

    for t = linspace(0,2*pi,100)
        R = tex_render(mul*(vec(C0(:,:,2))' + 0.8*sin(t)), mul*(vec(C0(:,:,1))' + 0.8*cos(t)), l, g);
        img = permute(reshape(R, [3 S]), [2 3 1]);
        imagec(img);
        drawnow;

        writeVideo(writerObj,sqrt(max(0,min(1,img))));
    end
    
    close(writerObj);
    
end