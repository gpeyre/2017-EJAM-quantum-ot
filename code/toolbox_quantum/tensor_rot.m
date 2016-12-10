function U = tensor_rot(Theta)

U = zeros(length(Theta), 2, 2);

U(:,1,1) = cos(Theta);   U(:,1,2) = sin(Theta); 
U(:,2,1) = -sin(Theta);  U(:,2,2) = cos(Theta); 

end