close all

% spline_surf_vec
 %intersections
 
  plot3(patch_bound_Xs,patch_bound_Ys,5*ones(us_n_basf-1,vs_n_basf-1),'LineWidth',4)
 %plot3(patch_bound_X,patch_bound_Y,5*ones(u_n_basf-1,v_n_basf-1),'LineWidth',4)
 hold on
 plot3(patch_bound_Ys,patch_bound_Xs,5*ones(us_n_basf-1,vs_n_basf-1),'LineWidth',4)
 %plot3(patch_bound_Y,patch_bound_X,5*ones(u_n_basf-1,v_n_basf-1),'LineWidth',4)
 
 plotresultsvec(us_knots, vs_knots,P0s,P1s,P2s,P3s,Xs,zs,Errs)
 plotresultsvec(u_knots, v_knots,P0,P1,P2,P3,X,z,Errs)

 for i=1:n_points
      plot3(point(i,5),point(i,6),point(i,7),'k.','MarkerSize',50);
 end
 
% [uvu2v2(1:2)',u_knot(h),uvu2v2(3),ptl(1),ptl(2),ptl(3),k,l,m,o];
%pointx = point;
pointb = [point];


%%% sort intersections by patches (1st patch)

 sp_i = (pointb(:,10)-1)*(us_n_intervs-1) + pointb(:,11);
 sp_i2 = (pointb(:,8)-1)*(u_n_intervs-1) + pointb(:,9);

% [ssp_i idx] = sort(sp_i);
% 
% m = length(sp_i); %m-out
% 
% 
% 
% pointb = pointb(idx,:);
% 
% %%% Compute intersectioned patches
% patch_int(1,2) = 1;
% patch_int(1,1) = ssp_i(1);
% %a = sp_i(1);
% different = 1;
% for i =2:n_points
%     if ssp_i(i) == patch_int(different,1)%a
%         patch_int(different,2) = patch_int(different,2)  +1;
%         continue
%     else
%         %a = sp_i(i);
%         different = different +1;
%         patch_int(different,1) = ssp_i(i);
%         patch_int(different,2) = 1;
%     end
% end
%  
% 
% different;



%return






%%% detect intersection point types 
% -1 - interion
%  0 - global boundary
%  1 - patch boundary

coinf = zeros(n_points,2);
%out = zeros(different,1);
point_type = zeros(n_points,1); % 1D

for j=1:n_points
        coinf(j,:) = boundary_point(pointb(j,3:4),pointb(j,10:11),us_knots,vs_knots);
         %type of interrsection (u,v) % 2D
        
%         if coinf(j,1) > -1 || coinf(j,2) > -1  % number of boundary points for patch
%             out(j) = out(j) +1;
%         end 
        
        % define intersection type
        if coinf(j,1) == -1 && coinf(j,2) == -1  % internal 
            point_type(j) = -1;
        elseif   coinf(j,1) == 0 || coinf(j,2) == 0 % global boundary
            point_type(j) = 0;
        elseif   coinf(j,1) == 1 && coinf(j,2) == 1 % patch internal boundary (patch edge)
            point_type(j) = 11;
        elseif   coinf(j,1) == 1 || coinf(j,2) == 1 % patch internal boundary
            point_type(j) = 1;
        end
        
        
end

 

% number of outputs



% % Sort points in patch
%   offset = 0;  
%   for j=1:different
%       
%       
%       patch_points= offset:offset+patch_int(j,2);
%       
%       boundg = find(point_type(patch_points) == 0)
%       boundi = find(point_type(patch_points) == 1)
%       bound = [boundg , boundi];
%       
%       l = lenght(bound)
%       if l >2
%           disp('more then 2 boundary points')
%       end
%       
%       dist = zeros(patch_int(j,2));
%       
%       for k = 1: offset+patch_int(j,2)
%           for i = 1: offset+patch_int(j,2)
%               
%               dist(i) = norm(pointsb(patch_points(k),1:2)- pointsb(patch_points(i),1:2))
%               i+offset
%               
%           end
%           
%       end
%       offset = offset + patch_int(j,2);
%   end
%

% Sort points in surface (create components)

  n_components = 0;

  % until 0 poits
  
  
  offset = 1;  
  pointb = [pointb, coinf, point_type, sp_i];
  matrixpoints = ones(n_points,1); 
  point_component = ones(n_points,1); 
  [a,b] = size(pointb);
  splinepoint = pointb;
  
  bound = find(pointb(offset:n_points,14) == 0); % find boundary point
  
  while(isempty(bound) == 0) % tests if there exist not yet connected point on the boundary of the surface 
      n_components = n_components +1;
      comp_bounds(n_components,1) = offset;
      
      
      for i=offset:n_points-1 % while "end of the component  - podminka"
          patch = pointb(i,15);
          id_point = find(pointb(i+1:n_points,15) == patch)+i;
          n_patch_points = length(id_point);
          
          if n_patch_points == 0 % "end of the component" - other
              %???
              break
          elseif n_patch_points == 1 %two points on patch
              pointb = swaplines(pointb,i+1,id_point );
              offset  = offset +1;
          else % more then two points on patch
              for l=1:n_patch_points
                  dist = 1/eps * ones(n_patch_points,1);
                  for k=l:n_patch_points
                      if id_point(k) ~= 0
                      dist(k) = norm(pointb(i+l-1,1:2) - pointb(id_point(k),1:2));
                      end
                  end
                  [~, b] = min(dist);
                  pointb = swaplines(pointb,i+l,id_point(b));
                  id_point(b) = 0;
              end
              offset = offset + n_patch_points;
              break
          end
               id = get_neigbour(pointb(offset-1,:))
          % connect new patch
          % test last point (zero flag)
      end

      
    comp_bounds(n_components,2) = offset; %% offset + 1
    bound = find(pointb(offset:n_points,14) == 0); % find next boundary point (if exist)
  end
  
      
      
      
      
%        %current patch ID
%        % +2
%       if isempty(id_out) == 1
%           break
%           %get_neigbour
%       else
%           splinepoint = swaplines(splinepoint,offset+2,offset+bound(1)); % indexy
%       end
%       
%       
%       
%       splinepoint = swaplines(splinepoint,offset+1,offset+bound(1));
%       
%       start = offset+1;
%       comp_bounds(n_components,1) = start;
%       
%       for j=offset+1:n_points-1
%           dist = 1/eps * ones(n_points,1);
%           for k=j+1:n_points
%               dist(k) = norm(splinepoint(j,1:2) - splinepoint(k,1:2));
%           end
%           [a b] = min(dist);
%           splinepoint = swaplines(splinepoint,b,j+1);
%           offset = offset +1;
%           if (splinepoint(j+1,14)==0)
%               break
%           end
%       end
%       comp_bounds(n_components,2) = offset;
%       
%        bound = find(pointb(offset+1:n_points,14) == 0)
%       pause
%   end      
      
      %%%% Remove duplicite points ()


%       for i=start:offset
%           %       if splinepoint(i,14) == 1
%           %           a = 1;
%           %       end
%           if splinepoint(i,14)*a == 1
%               matrixpoints(i) = 0;
%               a = 0;
%           elseif splinepoint(i,14) == 1
%               a =1;
%           else
%               a = 0;
%           end
%       end      
      
      
%       a = 0;
%       for i=start:offset
%           %       if splinepoint(i,14) == 1
%           %           a = 1;
%           %       end
%           if splinepoint(i,14)*a == 1
%               matrixpoints(i) = 0;
%               a = 0;
%           elseif splinepoint(i,14) == 1
%               a =1;
%           else
%               a = 0;
%           end
%       end

  
% Interpolate points
%  
% 
% pointb = pointb(1:m-out,:)
% 
% mm = m - out;
% 
% %return
% 
mr = sum(matrixpoints);
XYZ = zeros(mr,3);
uv = zeros(mr,2);
uvs = zeros(mr,2);
deltas = zeros(mr,1);
j = 0;

for i=1:m
    if matrixpoints(i)==1
        j = j+1;
        XYZ(j,:) =splinepoint(i,5:7);
        uv = splinepoint(i,1:2);
        uvs = splinepoint(i,3:4);
        if j~=1
            deltas(j) = norm(XYZ(j,:)-XYZ(j-1,:));
        end
    end
end

Deltas = zeros(mr,1);

for i=2:mr
    Deltas(i)=sum(deltas(1:i));
end

pointspersegment=4;
nbasf = ceil(mr/pointspersegment);
t_knots = get_knot_vector(nbasf+2);
t = Deltas/Deltas(mr);
X = zeros(mr,nbasf+2);
for i=1:mr
    
    [tf, ~] = splinebasevec(t_knots,t(i),0);
    X(i,:) = tf;
end

sol = X\XYZ;
sol1 = X\uv;
sol2 = X\uvs;
X*sol - XYZ;

%X'*X\X'*XYZ

ns = 100;
td = linspace(0,1,ns);
for i=1:ns
    PT(i,:)=splinebasevec(t_knots,td(i),0)'*sol;
    
end
plot3(PT(:,1),PT(:,2),PT(:,3),'r','LineWidth',3);
%
%
%
%
%
%