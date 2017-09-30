% make_map_gpu(100,90,160,8,8,8,24,0.01,'data_tooth/','maps_tooth_8_8_8_24/',1);
% SumR(100,90,160,24,'image','maps_tooth_8_8_8_24/','tooth8sum');
% make_map_gpu(100,90,160,16,16,16,24,0.01,'data_tooth8sum/','maps_tooth_16_16_16_24/',1);
% SumR(100,90,160,24,'image','maps_tooth_16_16_16_24/','tooth16sum');
% make_map_gpu(100,90,160,24,24,24,24,0.01,'data_tooth16sum/','maps_tooth_24_24_24_24/',1);
% SumR(100,90,160,24,'image','maps_tooth_24_24_24_24/','tooth24sum');



make_map(256,256,256,8,8,8,24,0.01,'data_kiwi/','maps_kiwi_8_8_8_24/',1);
SumR(256,256,256,24,'image','maps_kiwi_8_8_8_24/','kiwi8');
make_map(256,256,256,16,16,16,24,0.01,'data_kiwi8/','maps_kiwi_16_16_16_24/',1);
SumR(256,256,256,24,'image','maps_kiwi_16_16_16_24/','kiwi16');
make_map(256,256,256,24,24,24,24,0.01,'data_kiwi16/','maps_kiwi_24_24_24_24/',1);
SumR(256,256,256,24,'image','maps_kiwi_24_24_24_24/','kiwi24');



% make_map(256,256,256,8,8,8,24,0.01,'data_bonsai/','maps_bonsai_8_8_8_24/',1);
% SumR(256,256,256,24,'image','maps_bonsai_8_8_8_24/','bonsai8');
% make_map(256,256,256,16,16,16,24,0.01,'data_bonsai8/','maps_bonsai_16_16_16_24/',1);
% SumR(256,256,256,24,'image','maps_bonsai_16_16_16_24/','bonsai16');
% make_map(256,256,256,24,24,24,24,0.01,'data_bonsai16/','maps_bonsai_24_24_24_24/',1);




% make_map(256,256,128,5,5,5,24,0.01,'data_ldCT/','maps_ldCT_5_5_5_24_random_method/',1);
% SumR(256,256,128,24,'image','maps_ldCT_5_5_5_24_random_method/','ldCT5_random_method');
% make_map(256,256,128,11,11,11,24,0.01,'data_ldCT/','maps_ldCT_11_11_11_24_random_method/',1);
% SumR(256,256,128,24,'image','maps_ldCT_11_11_11_24_random_method/','ldCT11_random_method');
% make_map(256,256,128,21,21,21,24,0.01,'data_ldCT/','maps_ldCT_21_21_21_24_random_method/',1);
% SumR(256,256,128,24,'image','maps_ldCT_21_21_21_24_random_method/','ldCT21_random_method');

