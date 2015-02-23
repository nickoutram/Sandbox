
D = dir('GIS Polygons/Sampled Reefs/*shp');

for i = 1:26
    SHPD = m_shaperead(['GIS Polygons/Sampled Reefs/' D(i).name(1:end-4)]);
    
    Ix = SHPD.ncst{1}(:,1);
    Iy = SHPD.ncst{1}(:,2);
    polyarea(Ix,Iy)
    disp(i);
end

% 
% 
% D = dir('Cvag_data_Philippines/projected shapefiles/*shp');
%     figure(1), clf, hold on; L = 0;
% for i = 1:29
%     plot(S(i).X,S(i).Y,'.')
%     L = L + length(S(i).Y);
%     disp(i);
% end
%     
%     
% for j = 1:26
%     SHPD = m_shaperead(['Cvag_data_Philippines/projected shapefiles/' D(j).name(1:end-4)]);
% 
%     
%     %     Ix = SHPD.ncst{j}(:,1)./1.1e5 - 1.485;
%     %     Iy = SHPD.ncst{j}(:,2)./1.1e5 - 0.15;
%         Ix = SHPD.ncst{1}(:,1);
%         Iy = SHPD.ncst{1}(:,2);
%         plot(Ix,Iy,'g','linewidth',2)
% 
%     disp('REEF = ');
%     disp(j);
% end
