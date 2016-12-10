function [idx,p] = selectPoint(V)

% idx=selectPoint(figure, mesh)

datacursormode on
dcm_obj = datacursormode(gcf); 

fprintf(1,'Press any key when you have selected your vertex.\n');
pause;
a = getCursorInfo(dcm_obj);
p = a.Position;

[~,idx] = min( (V(1,:)-p(1)).^2 + (V(2,:)-p(2)).^2 + (V(3,:)-p(3)).^2 );


% searchResult = (mesh.vertices==repmat(p,size(mesh.vertices,1),1));
% searchResult = searchResult(:,1)&searchResult(:,2)&searchResult(:,3);

% idx = find(searchResult);

end