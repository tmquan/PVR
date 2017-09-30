function u = bndcrop(v, sz)
    u = zeros(sz, class(v));
    u(end/2-size(u,1)/2+1:end/2+size(u,1)/2, end/2-size(u,2)/2+1:end/2+size(u,2)/2, end/2-size(u,3)/2+1:end/2+size(u,3)/2, :) = v(1:size(u,1), 1:size(u,2), 1:size(u,3), :);
return