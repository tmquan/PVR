function [output] = SumR(x,y,z,num, prefix, folder,save_file_name)
    
    for k=1:num
        Mapss(:,:,:,k)=imreadtif([folder prefix num2str(k+num, '%03d') '.rmap']);
    end
    
                   % t_response=zeros(num,1);

                   % for kk=1:num
                   %     t_response(kk)=sum(sum(sum(Mapss(:,:,:,kk))));
                   % end
                   % for kk=1:3
                   %     [t1 t2]=max(t_response);
                   %     Mapss(:,:,:,t2)=0;
                   %     t_response(t2)=0;
                   % end


    output=zeros(x,y,z);
    for z_i=1:z
        for y_i=1:y
            for x_i=1:x
                output(x_i,y_i,z_i)=sum(Mapss(x_i,y_i,z_i,:));
            end
        end
    end
    

	if(exist(['data_' save_file_name], 'dir'))
		rmdir(['data_' save_file_name], 's');
    end
	mkdir(['data_' save_file_name]);

   imwritetif(output,['data_' save_file_name '/' save_file_name '_1.tif']);
return