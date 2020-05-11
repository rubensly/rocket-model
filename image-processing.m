% Script to determine distance between Rush Rhees and the Interfaith Chapel
clear all;close all;clc
 
% Read image and calculate conversion factor
img = imread('UofR_quad_edit.jpg');
imshow(img)
[x y] = getpts;
d1 = sqrt((x(2,1)-x(1,1))^2+(y(2,1)-y(1,1))^2);
% Calculate scaling factor (Walkpath width is 2.31 m)
f = 2.31/d1;
 
% Threshold to 96%
i96 = im2bw(img,0.96);
% Find centroid of all connected regions
c = regionprops(i96,'centroid');
cend = cat(1, c.Centroid);
 
% Find surface areas and sort in descending order
a = regionprops(i96,'area');
regionArea = cat(1, a.Area);
[val idx] = sort(regionArea,'descend');
 
% Sort centroids in the same order as areas
ranked_cent = zeros(length(cend),1);
xcentroid = cend(:,1);ycentroid = cend(:,2);
ranked_cent = [xcentroid(idx),ycentroid(idx)];
 
% Plot 2 largest centroids representing library and chapel
imshow(i96)
hold on
plot(ranked_cent(1:2,1),ranked_cent(1:2,2),'b*')
 
% Calculate pixel distance between desired centroids
d2 = sqrt((ranked_cent(2,1)-ranked_cent(1,1))^2+(ranked_cent(2,2)-ranked_cent(1,2))^2);
 
% Convert pixel distance to meters & feet
distance_meters = d2*f
distance_feet = d2*f*3.28084
