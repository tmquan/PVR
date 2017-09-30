function [s, d, y, gy] = saveMaps(S, G, Y, opt, prefix, folder)
	fprintf('Saving to maps\n');

    %Create folders
% 	if( isempty(folder))
% 		folder = 'maps/'
% 	end
% 	if(exist(folder, 'dir'))
% 		rmdir(folder, 's');
%     end
% 	mkdir(folder);

    elemSize = opt.elemSize;
	%dictSize = size(D); % TODO
	blobSize = opt.blobSize;
	dictSize = opt.dictSize;
	numMaps = blobSize(end);
    
    G = bndcrop(G, dictSize);
    G = reshape(G, dictSize);
    size(G)
    Y = reshape(Y, blobSize);
	%% Save raw map S
	S = 255*scale1(abs(S));
	
    imwritetif(S, [folder prefix '.omap']);

	%% Calculate response maps
    Gtmp = zeros(size(Y));
    Gtmp(1:dictSize(1), 1:dictSize(2), 1:dictSize(3),:) = G;
	Gf  = fft3(Gtmp); size(Gf);
	Yf  = fft3(Y); size(Yf);
	GYf = bsxfun(@times, Gf, Yf); 
	GY  = ifft3(GYf);
    
    tempS=zeros(blobSize(1),blobSize(2),blobSize(3));

   % GY(GY(:)<0)=0;
    
    tempS=sum(GY(:,:,:,:),4);
    
    minS=min(tempS(:))/(numMaps);
    maxS=max(tempS(:))-min(tempS(:));
    
    GY = GY-minS;
    GY = GY/maxS;
    GY = 255*GY;

	%% Calculate the feature maps
	Sf  = fft3(S);
	size(Sf)
	size(Gf)
	GSf = bsxfun(@times, Gf, Sf); 
	GS  = ifft3(GSf);
    %% Calculate the atom map
    D  = bndcrop(G, dictSize);
    
	Y   = 255*scale1(Y); %Scale the intensity
	%GY  = 255*scale1(GY); %Scale the intensity
	GS  = 255*scale1(abs(GS)); %Scale the intensity
    D   = 255*scale1(D); %Scale the intensity
    
	s = S;
    d = D;
	y = Y;
	gy = GY;
	%gs = GS;
    
    %% Write the maps
 	for n=1:numMaps
 		imwritetif(Y(:,:,:,n) , [folder prefix num2str(n, '%03d') '.smap']);	% Sparse maps
 		imwritetif(GY(:,:,:,n), [folder prefix num2str(n+numMaps, '%03d') '.rmap']);	% Response maps
 		imwritetif(GS(:,:,:,n), [folder prefix num2str(n+numMaps*2, '%03d') '.fmap']);	% Feature maps
        imwritetif(D(:,:,:,n),  [folder prefix num2str(n, '%03d') '.amap']);	% Atom maps
 	end


return