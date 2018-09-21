function draw_normals(Mesh)

    if(~isfield(Mesh,'vertex_normals'))
        disp('Mesh does not have a "vertex normal" field. Calculating normals...');
        Mesh.vertex_normals = calc_normals(Mesh);
    end

    n = size(Mesh.vertex_normals,1);
    v = Mesh.vertices;
    vn = Mesh.vertex_normals;
    Mesh = rmfield(Mesh,'vertex_normals');
    disp_meshes(Mesh);
    hold on;
    quiver3(v(:,1),v(:,2),v(:,3),v(:,1)+4*vn(:,1),v(:,2)+4*vn(:,2),v(:,3)+4*vn(:,3),'LineWidth',1);
    
    figure();
    disp_meshes(Mesh);
    vn2 = vertexNormal(Mesh.vertices,Mesh.faces);
    quiver3(v(:,1),v(:,2),v(:,3),v(:,1)+4*vn2(:,1),v(:,2)+4*vn2(:,2),v(:,3)+4*vn2(:,3),'LineWidth',1);
%     for i = 1:n
%         quiver3(v(i,1),v(i,2),v(i,3),v(i,1)+vn(i,1),v(i,2)+vn(i,2),v(i,3)+vn(i,3));
%     end
    drawnow;
end