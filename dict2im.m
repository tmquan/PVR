function img = dict2im(D)
    [dimy, dimx, dimz, dimk] = size(D);
    padsize = 1;
    padval  = 1;
    
    data = D;
    data = data - min(data(:));
    data = data / max(data(:));
    
    data = permute(data, [4, 1, 2, 3]);
    %data = padarray(data, [0, padsize, padsize, 0], padval, 'both');
    data = data(:,:,:,1); % Take the first slice if 3D
    
    
    data = squeeze(data); %size(data)
    %figure; imagesc(data); axis equal off; colormap gray; drawnow;
    data = padarray(data, [0, padsize, padsize], padval, 'both');
    % size(data)
    
    % Force the number of filters to be square
    n = ceil(sqrt(dimk));
    % for i=1:n
    %     for j=1:n
    %         k = sub2ind([n, n], i, j);
    %         if k>dimk 
    %             continue
    %         end
    %         C{i, j} = squeeze(data(k,:,:));
    %     end
    % end
    % size(C)
    % img = cell2mat(C);
    img = zeros([data(2)*n, data(3)*n]);
    for i=1:n
        for j=1:n
            k = sub2ind([n, n], i, j);
            if k>dimk 
                continue
            end
            img((13*(i-1)+1):(13*i), ...
                (13*(j-1)+1):(13*j)) = squeeze(data(k,:,:));
        end
    end

end