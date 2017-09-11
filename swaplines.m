function [ mat ] = swaplines( mat,i,j )
% swap corresponding lines in given matrix

temp = mat(j,:);
mat(j,:) = mat(i,:);
mat(i,:)= temp;

end

