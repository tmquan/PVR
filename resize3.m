function u = resize3(v, sz)
    u = zeros(sz, class(v));
    %u(1:size(u,1), 1:size(u,2), 1:size(u,3), :) = v(1:size(u,1), 1:size(u,2), 1:size(u,3), :);
    %[xu, yu, zu, ~] = sz;
    xu = sz(1);
    yu = sz(2);
    zu = sz(3);
    ku = sz(4);
    for k=1:ku
        size(u(:,:,:,k))
        size(v(:,:,:,k))
        u(:,:,:,k) = interp2(v(:,:,:,k), xu, yu);
    end
return