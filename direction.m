function [ dir ] = direction( x1,x2,y1,y2,o2 )
% This function calculates the direction (horizontal & zenith angles)

dir = atan2((y1-y2),(x1-x2))-o2;

end