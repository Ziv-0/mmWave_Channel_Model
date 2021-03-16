function pos_vector = generate_Kusers_pos(K)
% 产生K个用户的用户位置
% pos_vector K*3([x,y,z])
%半径分布10:100之间的均匀分布
%角度分布(-90到90)
pos_vector = zeros(K,3);
for i = 1:K
    theta = rand*pi-pi/2;
    r = 90*rand + 10;
    x = r*cos(theta);
    y = r*sin(theta);
    z = rand*0.5+0.7;
    pos_vector(i,:) = [x,y,z];
end
end

