pointSource = [8827 2607 4875 4470 6708 6743 6359 3496 1606 6136 9979 5244 5673 5978 7295];

[U, SI_source, SI_target] = getSpinImageCorrespondences(Source,Target,pointSource);

figure();
v = Source.vertices;
scatter3(v(:,1),v(:,2),v(:,3),3); hold on;
scatter3(v(pointsSource,1),v(pointsSource,2),v(pointsSource,3),50,'r','filled');
axis equal;

figure();
v = Target.vertices;
scatter3(v(:,1),v(:,2),v(:,3),3); hold on;
scatter3(v(U,1),v(U,2),v(U,3),50,'r','filled');
axis equal;
