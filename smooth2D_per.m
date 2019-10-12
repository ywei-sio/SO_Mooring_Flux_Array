function [As]=smooth2D_per(A,dx,dy,lat)
[m,n]=size(A);% periodic in 1st dimension.
As=A;
% missval=-1.e34;

for i=1:m
    for j=1:n
        % %      ix1=max(1,i-dx);ix2=min(m,i+dx);
        % %      iy1=max(1,j-dy);iy2=min(n,j+dy);
        % %      a=A(ix1:ix2,iy1:iy2);b=a(:);
        
        dx0=floor(dx/cos(lat(j)/90*pi/2));
        ix1=mod(i-dx0-1:i-1,m)+1;ix2=mod(i:i+dx0-1,m)+1;
        iy1=max(1,j-dy);iy2=min(n,j+dy);
        a=A([ix1,ix2],iy1:iy2);b=a(:);
        
        b(isnan(b))=[];
        b(abs(b)>10000)=[];
        
        if isnan(A(i,j)) ||abs(A(i,j))>10000
            As(i,j)=nan;
        else
            As(i,j)=mean(b);
            
        end
    end
end
