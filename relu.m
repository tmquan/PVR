function y=relu(a11)
    y=gpuArray.zeros(size(a11),'single');
    [x11,y11,z11,k11]=size(a11);
    for x111=1:x11
        for y111=1:y11
            for z111=1:z11
                for k111=1:k11
                    if a11(x111,y111,z111,k111)<0
                        y(x111,y111,z111,k111)=0;
                    else
                        y(x111,y111,z111,k111)=a11(x111,y111,z111,k111);
                    end
                end
            end
        end
    end
return 