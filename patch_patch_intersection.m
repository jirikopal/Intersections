function [ uvt,pt,conv ] = patch_patch_intersection( uvt, ui,vi,u2i,v2i,u_knots, v_knots, X,us_knots, vs_knots, Xs, nit ,m,o,k,l)
%returns coordinetes of the intersection of the patches (if exist), one
%parameter must be fixed
pt =zeros(3,1);
conv =0;
tol = 1e-6; % in x,y,z
tol2 = 1e-4; % in u,v

if length(u2i) == 1
    [u2f, ~] = splinebasevec(us_knots,u2i,0,m);
end
if length(v2i) == 1
    [v2f, ~] = splinebasevec(vs_knots,v2i,0,o);
end

for i=1:nit
    [uf, ~] = splinebasevec(u_knots,uvt(1),0,k);
    [vf, ~] = splinebasevec(v_knots,uvt(2),0,l);
    [ufd, ~] = splinebasevec(u_knots,uvt(1),1,k);
    [vfd, ~] = splinebasevec(v_knots,uvt(2),1,l);
      
    if length(u2i) == 1
        [v2f, ~] = splinebasevec(vs_knots,uvt(3),0,o);
        [v2fd, ~] = splinebasevec(vs_knots,uvt(3),1,o);
        dXYZp2 = (kron(v2fd',u2f')*Xs)';
    end
    if length(v2i) == 1
        [u2f, ~] = splinebasevec(us_knots,uvt(3),0,m);
        [u2fd, ~] = splinebasevec(us_knots,uvt(3),1,m);
        dXYZp2 = (kron(v2f',u2fd')*Xs)';
    end

    dXYZu1 = (kron(vf',ufd')*X)'; %
    dXYZv1 = (kron(vfd',uf')*X)'; %
    
    J = [dXYZu1 dXYZv1 -dXYZp2];
    
    deltaXYZ = (kron(vf',uf') * X)' - (kron(v2f',u2f') * Xs)';
    %dist0 = norm(deltaXYZ);
    
%     alpha = 1;
%     for i=1:10 
%         uvt2 = uvt- alpha * J\deltaXYZ;
%         [~,uvt2] = rangetest(uvt2,ui,vi,u2i,v2i,0.0);
%         dist = get_delta(u_knots, v_knots,us_knots, vs_knots,uvt2,u2i,v2i,k,l,m,o,X,Xs);
%         if dist0 < dist
%             alpha = alpha/2;
%         else
%             uvt = uvt2;
%             break
%         end
%     end

    uvt = uvt- J\deltaXYZ;
    [test,uvt] = rangetest(uvt,ui,vi,u2i,v2i,0.0);    

end

[test,uvt] = rangetest(uvt,ui,vi,u2i,v2i,tol2);
if test == 1
   dist = get_delta(u_knots, v_knots,us_knots, vs_knots,uvt,u2i,v2i,k,l,m,o,X,Xs);
end

if test == 1
    if dist <= tol
        pt =kron(vf',uf')*X;
        conv =1;
    end 
else
    uvt = zeros(3,1);
end

end

