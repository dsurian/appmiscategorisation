function [ix_cat,featuresAll,ixAllOutliers ] = createSynthetic(numbMemberPerGroup,dimMemberPerGroup)

ixAllOutliers = [];
ixAll = [];
ixStart = 1;
featuresAll = [];
ix_cat = [];
ixOutlier_cat = [];
ix1_all = [];
ix2_all = [];
ix3_all = [];
settingParams;

numbOutlier = ceil(prcogr * numbMemberPerGroup);
numbTopic = 3;
fprintf('True number of topics: %d\n',numbTopic);
numbCategory = 5;

for nc=1:numbCategory
    features1 = [];
    a = randi(rt,numbMemberPerGroup,dimMemberPerGroup);
    b = randi(rs,numbMemberPerGroup,dimMemberPerGroup);
    c = zeros(numbMemberPerGroup,dimMemberPerGroup);
    features1 = [ a b c ];
    featuresAll = [ featuresAll; features1 ];
    ix1 = [ixStart:ixStart+size(features1,1)-1]';
    ixAll = [ ixAll; ix1 ];
    ix1_all = [ ix1_all; ix1 ];
    ix_cat = [ ix_cat; ix1 (nc*ones(size(ix1,1),1)) ];
    ixStart = size(ixAll,1) + 1;
    
    features2 = [];
    a = randi(rt,numbMemberPerGroup,dimMemberPerGroup);
    b = randi(rs,numbMemberPerGroup,dimMemberPerGroup);
    c = zeros(numbMemberPerGroup,dimMemberPerGroup);
    features2 = [ b a c ];
    featuresAll = [ featuresAll; features2 ];
    ix2 = [ixStart:ixStart+size(features2,1)-1]';
    ixAll = [ ixAll; ix2 ];
    ix2_all = [ ix2_all; ix2 ];
    ix_cat = [ ix_cat; ix2 (nc*ones(size(ix2,1),1)) ];
    ixStart = size(ixAll,1) + 1;
    
    features3 = [];
    a = randi(rn,numbOutlier,dimMemberPerGroup);
    b = randi(rs,numbOutlier,dimMemberPerGroup);
    
    c = zeros(numbOutlier,dimMemberPerGroup);
    d = randi(ro,numbOutlier,dimMemberPerGroup);
    
    numbnoises = ceil(xn1 * size(c,1) * size(c,2));
    ixnoises = randperm(size(c,1) * size(c,2),numbnoises);
    noise = randi(rn,1,numbnoises);
    c(ixnoises) = noise;
    features3 = [ c c d ];    
        
    featuresAll = [ featuresAll; features3 ];
    ix3 = [ixStart:ixStart+size(features3,1)-1]';
    ixAll = [ ixAll; ix3 ];
    ix3_all = [ ix3_all; ix3 ];
    ix_cat = [ ix_cat; ix3 (nc*ones(size(ix3,1),1)) ];
    ixOutlier_cat = [ ixOutlier_cat; ix3 (nc*ones(size(ix3,1),1)) ];
    
    ixAllOutliers = [ ixAllOutliers; [ixStart:ixStart+size(features3,1)-1]' ];
    ixStart = size(ixAll,1) + 1;
end

numbOutliersTotal = length(ixAllOutliers);
featuresAll_copy = featuresAll;

%fprintf('Numb data true: %d, numb outliers: %d\n',length(ixAll)-numbOutliersTotal,numbOutliersTotal);

featuresAll = featuresAll_copy;

% Inject noises to all
% Pick randomly the vectors, inject to them
numbnoises = ceil(xn * size(featuresAll,1));
ixnoises = randperm(size(featuresAll,1),numbnoises);
noise = randi(rn,1,1);

toInject = featuresAll(ixnoises,:);
ixzeros = find(toInject==0);
x = randperm(length(ixzeros),ceil(xn * length(ixzeros)));
toInject(ixzeros(x)) = noise;
featuresAll(ixnoises,:) = toInject;

ix_cat_ori = ix_cat;         % Contains information per category
ix_cat = [ [1:size(featuresAll,1)]' ones(size(featuresAll,1),1) ];

figure(1);
imagesc(featuresAll_copy);
colormap(flipud(gray));

figure(2);
imagesc(featuresAll);
colormap(flipud(gray));

return