% Write GMSH's geo-files

% This defines the domain.
  psi_min = 0.69; % suggested value = 0.69_dp
  psi_max = 1.20; % suggested value = 1.20_dp
    V_min = 0.50; % suggested value = 0.50_dp
% This defines the grid size.
       nx = 30;   % the number of nodes in psi-direction
       ny = 60;   % the number of nodes in V-direction

% Get boundary coordinates
[x,y] = mesh_ringlebFlow(psi_min, psi_max, V_min, nx, ny);

% Element Type
%EType = 'tri';
EType = 'qua';

% Element size;
l_boundary = 0.1;

% 1. Open geo-file
fileID = fopen(['ringlebFlow_',EType,'.geo'],'w');

% 2. Mesh parameters
fprintf(fileID,'// Set mesh Physical parameters\n');
fprintf(fileID,'l_boundary = %g;\n',l_boundary);
fprintf(fileID,'\n');

% 3. Set boundary splines
fprintf(fileID,'// Boundary points\n');
pID = 0;

    % Bottom
    for i = 1:nx
        pID = pID + 1;
        fprintf(fileID,'p%d = newp; Point(p%d) = {%g,%g,0,l_boundary};\n',pID,pID,x(i,1),y(i,1));
    end
    
    % Right
    for j = 1:ny
        pID = pID + 1;
        fprintf(fileID,'p%d = newp; Point(p%d) = {%g,%g,0,l_boundary};\n',pID,pID,x(nx+1,j),y(nx+1,j));
    end
    
    % Top
    for i = nx+1:-1:2
        pID = pID + 1;
        fprintf(fileID,'p%d = newp; Point(p%d) = {%g,%g,0,l_boundary};\n',pID,pID,x(i,ny+1),y(i,ny+1));
    end
    
    % Left
    for j = ny+1:-1:2
        pID = pID + 1;
        fprintf(fileID,'p%d = newp; Point(p%d) = {%g,%g,0,l_boundary};\n',pID,pID,x(1,j),y(1,j));
    end

fprintf(fileID,'\n');

% 4. Build boundary curves
fprintf(fileID,'// Build boundary curves\n');
fprintf(fileID,'bsp1 = newreg; BSpline(bsp1) = {');
for j = 1:nx
    fprintf(fileID,'p%d, ',j);
end
fprintf(fileID,'p%d};\n',nx+1);
fprintf(fileID,'bsp2 = newreg; BSpline(bsp2) = {');
for i = 1:ny
    fprintf(fileID,'p%d, ',nx+i);
end
fprintf(fileID,'p%d};\n',nx+ny+1);
fprintf(fileID,'bsp3 = newreg; BSpline(bsp3) = {');
for j = 1:nx
    fprintf(fileID,'p%d, ',nx+ny+j);
end
fprintf(fileID,'p%d};\n',nx+ny+nx+1);
fprintf(fileID,'bsp4 = newreg; BSpline(bsp4) = {');
for i = 1:ny
    fprintf(fileID,'p%d, ',nx+ny+nx+i);
end
fprintf(fileID,'p%d};\n\n',1);

% Build surface from curves
fprintf(fileID,'// Build surface from curves\n');
fprintf(fileID,'l1 = newl; Line Loop(l1) = {bsp1, bsp2, bsp3, bsp4};\n');
fprintf(fileID,'s1 = news; Plane Surface(s1) = {l1};\n');
fprintf(fileID,'\n');

% If quad-type elements requested
if strcmp(EType,'qua')
    fprintf(fileID,'Recombine Surface {:};\n\n');
end

% Define Physical elements
fprintf(fileID,'Physical Surface (1) = {s1};\n');
fprintf(fileID,'Physical Curve (1) = {bsp1};\n');
fprintf(fileID,'Physical Curve (2) = {bsp2};\n');
fprintf(fileID,'Physical Curve (3) = {bsp3};\n');
fprintf(fileID,'Physical Curve (4) = {bsp4};\n\n');

% Build mesh
fprintf(fileID,'Mesh 2;\n\n');

% Save format
fprintf(fileID,'// Save format\n');
fprintf(fileID,'Mesh.MshFileVersion = 2.2;\n');
fprintf(fileID,'Mesh.Format = 1;\n\n');

% Close file
fclose(fileID);
