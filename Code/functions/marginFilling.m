function img = marginFilling(org, r, val)

% comments

% org: original image
% r: filling radius
% val: values of filled regions

org(1:r,:) = val;
org(end-r:end,:) = val;
org(:,1:r) = val;
org(:,end-r:end) = val;

img = org;