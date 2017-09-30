function res = ecsc_cpu(D0, S0, plan, isTrainingDictionary, folder, original)
    %% If we want to train the dictionary
    if isempty(isTrainingDictionary)
        isTrainingDictionary = 1;
    end
    
    %% Parameters extractions
    elemSize = plan.elemSize;
    dataSize = plan.dataSize;
    atomSize = plan.atomSize;
    dictSize = plan.dictSize;
    blobSize = plan.blobSize;

    numAtoms = dictSize(4);
    % plan.elemSize = [128, 128,  1,   1];
    % plan.dataSize = [128, 128,  1, 512]; % For example
    % plan.atomSize = [ 11,  11,  1,   1];
    % plan.dictSize = [ 11,  11,  1, 100];
    % plan.blobSize = [128, 128,  1, 100];
    % plan.iterSize = [128, 128,  1,  16];


    gNx = (prod(blobSize));
    gNd = (prod(blobSize));


    glambda = (plan.lambda.Value);
    grho    = (plan.rho.Value);
    gsigma  = (plan.sigma.Value);

    %% Operators here
    %% Mean removal and normalisation projections
    Pzmn    = @(x) bsxfun(@minus,   x, mean(mean(mean(x,1),2),3));
    % Pzmn    = @(x) bsxfun(@minus,   x, mean(mean(mean(mean(x,1),2),3),4));
    % Pnrm    = @(x) bsxfun(@rdivide, x, sqrt(sum(sum(sum(sum(x.^2,1),2),3),4)));
    Pnrm    = @(x) bsxfun(@rdivide, x, sqrt(sum(sum(sum(x.^2,1),2),3)));  
    %Pnrm    = @(x) bsxfun(@rdivide, x, (sum(sum(sum(abs(x).^1,1),2),3)));

    %% Projection of filter to full image size and its transpose
    % (zero-pad and crop respectively)
    Pzp     = @(x) zeropad(x, blobSize);
    PzpT    = @(x) bndcrop(x, dictSize);

    %% Projection of dictionary filters onto constraint set
    Pcn     = @(x) Pnrm(Pzp((PzpT(x))));  

    %% Memory reservation
    gS0     = S0; 
    %gS0     = reshape(gS0, dataSize);
    gD0     = D0;
    gD0     = Pnrm(gD0);

    grx = Inf;
    gsx = Inf;
    grd = Inf;
    gsd = Inf;
    geprix = 0;
    geduax = 0;
    geprid = 0;
    geduad = 0;

    gX      = zeros(blobSize ,'single');
    gY      = zeros(blobSize,'single');
    gYprv   = gY;
    gXf     = zeros(blobSize,'single');
    gYf     = zeros(blobSize,'single');

    % gS      = gS0; 
    %gSf     = zeros(dataSize);

    gD      = zeros(blobSize,'single');
    gG      = zeros(blobSize,'single');
    gGprv   = zeros(blobSize,'single');

    gD      = zeros(blobSize,'single');
    gG      = Pzp(gD); % Zero pad the dictionary
    gGprv   = gG;

    gDf     = zeros(blobSize,'single');
    gGf     = zeros(blobSize,'single');

    gU      = zeros(blobSize,'single');
    gH      = zeros(blobSize,'single');

    gGf     = zeros(blobSize,'single');
    gGf     = fft3(gG);
    
    % Temporary buffers
    gGSf    = zeros(blobSize,'single');
    gYSf    = zeros(blobSize,'single');


    %% Set up algorithm parameters and initialise variables
    res  = struct('itstat', [], 'plan', plan);
    %% Main loops
    k = 1;
    tstart = tic;
    while k <= plan.MaxIter && (grx > geprix | gsx > geduax | ...
                                grd > geprid | gsd > geduad),
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Permutation here
       
        for n=randperm(size(gS0, 4)) %1:size(gS0,4)
            gS = gS0(:,:,:,n);
            
            if isTrainingDictionary
                % gS =  gpuArray(imrotate(gS, 360*rand(1,1), 'crop', 'bilinear')) ;
                %gS = permute(gS, [randperm(2), 3, 4]);
                % size(gS0)
                % r = randi([0 8],1,1);
                % switch r
                    % case 0
                        % gS = (gS0);
                    % case 1
                        % gS = rot90(gS0, 1);
                    % case 2
                        % gS = rot90(gS0, 2);
                    % case 3
                        % gS = rot90(gS0, 3);
                    % case 4
                        % gS = rot90(gS0, 4);
                    % case 5
                        % gS = fliplr(gS0);
                    % case 6
                        % gS = flipud(gS0);
                    % case 7
                        % gS = gS0';
                    % otherwise
                        % gS = gS0;
                % end
                figure(3); imagesc(gS(:,:,ceil(end/2))); axis equal off; colormap gray; drawnow; 
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% Compute the signal in DFT domain
            gSf  = fft3(gS); 
            %% Extract the atom iteration
            % for iter = 1:numIters
                % chunk = 1:blobSize(4);
                % march = chunk+(iter-1)*iterSize(4); % marching through the dictionary
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                gD      = gD0; %(:,:,:,march);
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                gG      = Pzp(gD); % Zero pad the dictionary, PARTIALLY
                gGf     = fft3(gG);
                size(gGf)
                size(gSf)
                gGSf    = bsxfun(@times, conj(gGf), gSf); 

                %% Solve X subproblem
                gXf  = solvedbi_sm(gGf, grho, gGSf + grho*fft3(gY-gU)); 
                gX   = ifft3(gXf); 
                gXr  = gX; %relaxation
				
					
                %% Solve Y subproblem
                gY   = shrink(gXr + gU, (glambda/grho)*plan.weight); % Adjust threshold 
                if k<95
                    if original==1
                        if isTrainingDictionary
                            if mod(k,2)==0
                                idx = randperm(numAtoms, ceil(0.126*numAtoms));
                                gY(:,:,:,idx) = 0;
                            end
                            if mod(k,2)==1
                                 for kk=1:numAtoms
                                    t_response(kk)=sum(sum(sum(gY(:,:,:,kk)~=0)));
                                end
                                total_response=sum(t_response(:));
                                current_response=0;
                                for kk=1:floor(0.125*numAtoms)
                                    [t1 t2]=max(t_response);
                                    gY(:,:,:,t2)=0;%gY(:,:,:,t2).*(k/100);
                                    t_response(t2)=0;%t_response(t2).*(k/100);
                                end
                            end
                        end
                    end
                    if original==0
                        total_response=0;
                        t_response=zeros(numAtoms,1);
                        if isTrainingDictionary
                            for kk=1:numAtoms
                                t_response(kk)=sum(sum(sum(gY(:,:,:,kk)~=0)));
                            end
                            total_response=sum(t_response(:));
                            current_response=0;
                            for kk=1:floor(0.125*numAtoms)
                                [t1 t2]=max(t_response);
                                gY(:,:,:,t2)=0;%gY(:,:,:,t2).*(k/100);
                                t_response(t2)=0;%t_response(t2).*(k/100);
                                %current_response=sum(t_response(:));
                                %if(current_response~=0)
                                %    gY=gY(:,:,:,:).*(total_response/current_response);
                                %end

                            end
                        end
                    end
                end
                % gT   = mean(gY,4);
                % for k=1:numAtoms
                %     gY(:,:,:,k) = gT;
                % end
                % gY(gY<0) = 0;
                % idx = randperm(numAtoms);
                %gY(:,:,:,:) = gY(:,:,:,idx);
                % gY = reshape(gY, [blobSize(1)*blobSize(3)*sqrt(numAtoms), blobSize(2)*blobSize(3)*sqrt(numAtoms)]);
                % gY = histeq(real(gY));
                % gY = reshape(gY, blobSize);

                gYf  = fft3(gY);
                % gYf = reshape(gYf, [blobSize(1)*blobSize(3)*sqrt(numAtoms), blobSize(2)*blobSize(3)*sqrt(numAtoms)]);
                % gYf = histeq(abs(gYf));
                % gYf = reshape(gYf, blobSize);

                % size(gYf)
                % size(gSf)
                % size(bsxfun(@times, conj(gYf), gSf))
                % gYSf = sum(bsxfun(@times, co  nj(gYf), gSf), 4);
                gYSf = (bsxfun(@times, conj(gYf), gSf));
                %% Solve U subproblem
                gU = gU + gXr - gY;
                
                %% Update params 
                gnX = norm(gX(:)); gnY = norm(gY(:)); gnU = norm(gU(:));
                grx = norm(vec(gX - gY))/max(gnX,gnY);
                gsx = norm(vec(gYprv - gY))/gnU;
                
                geprix = sqrt(gNx)*plan.AbsStopTol/max(gnX,gnY)+plan.RelStopTol;
                geduax = sqrt(gNx)*plan.AbsStopTol/(grho*gnU)+plan.RelStopTol;

                if plan.rho.Auto,
                    if k ~= 1 && mod(k, plan.rho.AutoPeriod) == 0,
                        if plan.rho.AutoScaling,
                            grhomlt = sqrt(grx/gsx);
                            if grhomlt < 1, grhomlt = 1/grhomlt; end
                            if grhomlt > plan.rho.Scaling, grhomlt = plan.rho.Scaling; end
                        else
                            grhomlt = plan.rho.Scaling;
                        end
                        grsf = 1;
                        if grx > plan.rho.RsdlRatio*gsx, grsf = grhomlt; end
                        if gsx > plan.rho.RsdlRatio*grx, grsf = 1/grhomlt; end
                        grho = grsf*grho;plan.lambda.Value		= 0.01; %10; 1; 0.1; 0.01; 0.001; 
                        gU = gU/grsf;
                    end
                end

                %% Record information
                gYprv = gY;

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                if isTrainingDictionary
                    %% Solve D subproblem
                    % size(gYSf)
                    % size(gG)
                    gDf  = solvedbi_sm(gYf, gsigma, gYSf + gsigma*fft3(gG - gH));
                    %gXf  = solvedbi_sm(gGf, grho, gGSf + grho*fft3(gY-gU)); 
                    gD   = ifft3(gDf);
                    gDr  = gD;

                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %% Solve G subproblem
                    gG   = Pcn(gDr + gH);
					% idx = randperm(numAtoms, 2);
					% gG(:,:,:,idx) = 0;
                    % gG =  PzpT(gG);
                   
                    % gG = abs(gG);
                    % gG = Pcn(gG);
                    % gG(gG<0) = 0;
                    % gG(:,:,:,:) = gG(:,:,:,idx);
                    % G = gather(gG);
                    % for d=1:size(gG,4)
                    %     G(1:size(D0,1),1:size(D0,2),:,d) = imrotate(G(1:size(D0,1),1:size(D0,2),:,d),360*rand(1,1), 'crop', 'bilinear') ;
                    %     % gG(1:size(D0,1),1:size(D0,2),d) = imrotate(gG(1:size(D0,1),1:size(D0,2),d), 360*rand(1,1), 'crop', 'bilinear') ;
                    %     % gG(1:size(D0,1),1:size(D0,2),d) = gG(1:size(D0,1),1:size(D0,2),d)';
                    %     % gG(1:size(D0,1),1:size(D0,2),d) = imrotate(gG(1:size(D0,1),1:size(D0,2),d), 90*randi(3), 'crop') ;
                    % end
                    % gG = gpuArray(G);
                    % gG = Pcn(gG);


                    %% Solve H subproblem
                    gH = gH + gDr - gG;

                    %% Update params    
                    gnD = norm(gD(:)); gnG = norm(gG(:)); gnH = norm(gH(:));
                    grd = norm(vec(gD - gG))/max(gnD,gnG);
                    gsd = norm(vec(gGprv - gG))/gnH;


                    geprid = sqrt(gNd)*plan.AbsStopTol/max(gnD,gnG)+plan.RelStopTol;
                    geduad = sqrt(gNd)*plan.AbsStopTol/(gsigma*gnH)+plan.RelStopTol;
                    
                    if plan.sigma.Auto,
                        if k ~= 1 && mod(k, plan.sigma.AutoPeriod) == 0,
                            if plan.sigma.AutoScaling,
                                gsigmlt = sqrt(grd/gsd);
                                if gsigmlt < 1, gsigmlt = 1/gsigmlt; end
                                if gsigmlt > plan.sigma.Scaling, gsigmlt = plan.sigma.Scaling; end
                            else
                                gsigmlt = plan.sigma.Scaling;
                            end
                            gssf = 1;
                            if grd > plan.sigma.RsdlRatio*gsd, gssf = gsigmlt; end
                            if gsd > plan.sigma.RsdlRatio*grd, gssf = 1/gsigmlt; end
                            gsigma = gssf*gsigma;
                            gH = gH/gssf;
                        end
                    end
                    %% Record information
                    gGprv = gG;
                end
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %% Collect information
                % Compute l1 norm of Y
                gJl1 = sum(abs(vec( gY)));
                % Compute measure of D constraint violation

                if isTrainingDictionary
                   gJcn = norm(vec(Pcn(gD) - gD));
                   %gJcn = vec(Pcn(gD) - gD);
                end
                % Compute data fidelity term in Fourier domain (note normalisation)
                gJdf = sum(vec(abs(sum(bsxfun(@times,gGf,gYf),4)-gSf).^2))/(2*prod(blobSize));
                gJfn = gJdf + glambda*gJl1
				k
                % Record and display iteration details
                tk = toc(tstart);
                res.itstat = [res.itstat;...
                    [k gather(gJfn) gather(gJdf) gather(gJl1) gather(grx) gather(gsx)...
                    gather(grd) gather(gsd) gather(geprix) gather(geduax) gather(geprid)...
                    gather(geduad) gather(grho) gather(gsigma) tk]];
                f6=figure(6);
                plot(res.itstat(:,2));
                xlabel('Iterations');
                ylabel('Functional value');drawnow;
               % tt=sprintf('learn1/%d.png',k);
               % saveas(f6,tt);


                %% Debug
                %G = gather(PzpT(gG));
                %figure(5);
                %tmp = squeeze(G(:,:,1,:));
                %imdisp(dict2img(tmp)); drawnow;

                %% Update D partially
                gD0 = PzpT(gG);
                % [~,idx] = sort(mean(mean(mean(gD0,1),2),3), 'descend');
                % gD0 = gD0(:,:,:,idx);
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % end % End chunk
            %% Debug
            D0 = gather(gD0);
            % [~,idx] = sort(mean(mean(mean(D0,1),2),3), 'ascend');
            % D0 = D0(:,:,:,idx);
            % size(D0)
            f5=figure(5);
            imagesc(tiledict(squeeze(D0(:,:,floor(dictSize(3)/2),:)))); axis equal off; colormap gray; drawnow;
            
            f7=figure(7);
            imagesc(tiledict(squeeze(D0(:,:,2,:)))); axis equal off; colormap gray; drawnow;

            f8=figure(8);
            imagesc(tiledict(squeeze(D0(:,:,dictSize(3)-2,:)))); axis equal off; colormap gray; drawnow;

 %           tt=sprintf('%s%d.png',k);
 %           saveas(f5,tt);
            
        end % End for
        tempGY=ifft3(fft3(gG).*fft3(gY));
        tempS=zeros(blobSize(1),blobSize(2));

        
      %  tempS(1:blobSize(1),1:blobSize(2))=sum(tempGY(1:blobSize(1),1:blobSize(2),ceil(end/2),:));
        
        tempS=sum(tempGY(:,:,ceil(end/2),:),4);
        
        figure(8);
        imagesc(tempS);axis equal off; colormap gray; drawnow;
        % imdisp(tiledict(squeeze(gD0))); axis equal off; colormap gray; drawnow;
        %% Update iterations
        k = k+1;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end %% End main loop

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Collect the output
    gGY = ifft3(fft3(gG).*fft3(gY));
    gGS = ifft3(bsxfun(@times, fft3(gG), fft3(gS)));
    res.G  = gather(gG);  
    res.Y  = gather(gY);  
    res.GY = gather(gGY);  
    res.GS = gather(gGS);  
    
    if isTrainingDictionary
        if( isempty(folder))
            folder = 'maps/'
        end
        if(exist(folder, 'dir'))
            rmdir(folder, 's');
        end
        mkdir(folder);

        for kkk=1:dictSize(3)
            f5=figure(5);
            imagesc(tiledict(squeeze(D0(:,:,kkk,:)))); axis equal off; colormap gray; drawnow;
            saveas(f5,[folder 'dictionary' num2str(kkk, '%02d') '.png']);
        end
        saveas(f6,[folder 'energyfunction.png']);
    end
end