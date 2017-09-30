function [] = imwritetif(mat, filename)
    dim = size(mat)
	dim = [dim 1 1];
	for k=1:dim(3)
		if (k==1)
			imwrite(uint8(mat(:, :, k)'), filename, 'tif', 'Compression','none');
		else
			imwrite(uint8(mat(:, :, k)'), filename, 'tif', 'WriteMode', 'append', 'Compression','none');
		end
	end
end