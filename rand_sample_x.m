function initx = rand_sample_x(N, D, xlimit)
% randomly sample N points in the space specified by xlimit
initx = zeros(N,D);

% parfor j=1:N
% initx(j,:) = rand(1,D).*(xlimit(2,:)-xlimit(1,:)) + xlimit(1,:);
% end

Ne = 100; % Number of points in each dimension delta =  (x_max-x_min)/100
for i = 1:D
    X_all(:,i) = linspace(xlimit(1,i),xlimit(2,i),Ne);
end
rxt = randi(Ne,N,D);
for i = 1: N
    for j = 1:D
        initx(i,j) = X_all(rxt(i,j),j); % Each column of this is a random vector for training
    end
end



end
% initx=sort(initx,1);