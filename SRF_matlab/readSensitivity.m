close all
filename = '../data/srf/sen_sub';
nthread = 4;
SX = 128;
SY = 128;
SZ = 128;
sen = zeros(4,SX*SY*SZ);
for i = 1:nthread
    fid = fopen([filename,num2str(i-1)],'rb');
    a = fread(fid,'double');
    fclose(fid);
end
sen = reshape(sen,SX,SY,SZ);

senXY = sen(:,:,SZ/2);
figure,
imagesc(senXY);axis equal; colorbar;
senXZ = squeeze(sen(:,SY/2,:));
figure,axis equal; colorbar;
imagesc(senXZ);axis equal; colorbar;
senYZ = squeeze(sen(SZ/2,:,:));
figure,axis equal; colorbar;
imagesc(senYZ);axis equal; colorbar;

% filename = '../data/source';
% fid = fopen(filename,'rb');
% source = fread(fid, 'double');
% source = reshape(source,3,128*128*128);
% plot3(source(1,:),source(2,:),source(3,:),'.');
% fclose(fid);