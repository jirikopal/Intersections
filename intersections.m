clc

% spline_surf_vec
 close all
% % %
 plotresultsvec(u_knots, v_knots,P0,P1,P2,P3,X,z,Err)
 hold on
 plotresultsvec(us_knots, vs_knots,P0s,P1s,P2s,P3s,Xs,zs,Errs)
hold on


%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Compute intersections
%%%%%%%%%%%%%%%%%%%%%%%%%

% u_n_basf = length(u_knots)-3;
% v_n_basf = length(v_knots)-3;
% us_n_basf = length(us_knots)-3;
% vs_n_basf = length(vs_knots)-3;

u_n_intervs = u_n_basf - 2;
v_n_intervs = v_n_basf - 2;
us_n_intervs = us_n_basf - 2;
vs_n_intervs = vs_n_basf - 2;

u_n_grid  = u_n_basf+1;
v_n_grid  = v_n_basf+1;
us_n_grid  = us_n_basf+1;
vs_n_grid  = vs_n_basf+1;


% Compute X & Y grid intersections

%%% Grid centers
[ X_coor,Y_coor] = compute_control_points( u_knots, v_knots, P0,P1,P2,P3);
[ Xs_coor,Ys_coor] = compute_control_points( us_knots, vs_knots, P0s,P1s,P2s,P3s);

%%% Compute Bounding Boxes
Z_coor  = vec2mat( z,u_n_basf,v_n_basf ); 
Zs_coor  = vec2mat( zs,us_n_basf,vs_n_basf );

[patch_bound_X,patch_bound_Y] = compute_patch_edges( u_knots, v_knots, P0,P1,P2,P3);
[patch_bound_Xs,patch_bound_Ys] = compute_patch_edges( us_knots, vs_knots, P0s,P1s,P2s,P3s);

% [ BB_X,BB_Y,BB_Z ] = compute_bounding_box( X_coor,Y_coor,Z_coor, u_n_intervs,v_n_intervs);
% [ BB_Xs,BB_Ys,BB_Zs ] = compute_bounding_box( Xs_coor,Ys_coor,Zs_coor, us_n_intervs,vs_n_intervs);

%%% Bonding boxes intersections
[isec,n_isec] = bounding_boxes_intersection( patch_bound_X,patch_bound_Y,Z_coor,patch_bound_Xs,patch_bound_Ys,Zs_coor);
[isec2,n_isec2] = bounding_boxes_intersection( patch_bound_Xs,patch_bound_Ys,Zs_coor,patch_bound_X,patch_bound_Y,Z_coor);


x = X_coor(:);
y = Y_coor(:);
xs = Xs_coor(:);
ys = Ys_coor(:);

Xs = [xs ys zs];

X = [x y z];
nt = 2; % number of main lines

nit = 5;
    
[point,  n_points,ninter] = get_intersection(X,Xs,u_n_intervs,v_n_intervs,u_knots,v_knots,...
     us_n_intervs,vs_n_intervs,us_knots,vs_knots,isec,n_isec, nit,nt);
    
[point2,  n_points2,ninter2] = get_intersection(Xs,X,us_n_intervs,vs_n_intervs,us_knots,vs_knots,...
    u_n_intervs,v_n_intervs,u_knots,v_knots,isec2,n_isec2, nit,nt);

  for i=1:n_points
      plot3(point(i,5),point(i,6),point(i,7),'r.','MarkerSize',50)
  end

 %for i=1:n_points2
 %    plot3(point2(i,5),point2(i,6),point2(i,7),'b.','MarkerSize',50)
 %end

%t = toc

%cas = sum(sum(ninter))/t


 plot3(patch_bound_Xs,patch_bound_Ys,5*ones(us_n_basf-1,vs_n_basf-1),'LineWidth',4)
 plot3(patch_bound_Ys,patch_bound_Xs,5*ones(us_n_basf-1,vs_n_basf-1),'LineWidth',4)
 %plot3(patch_bound_X,patch_bound_Y,5*ones(u_n_basf-1,v_n_basf-1),'LineWidth',4)
 %plot3(patch_bound_Y,patch_bound_X,5*ones(u_n_basf-1,v_n_basf-1),'LineWidth',4)
