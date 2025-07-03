%%%%% Calculate area of nucleus from segmented area
% requires tom_toolbox package

% Make a mask for the nucleus in imod using sculpt and warp function
% 3 contours are enough, other slices can be intrapolated
% save as "nuc" modelfile, object should be set to "open"

% Pixelsize in um
pixelsize = 0.00138;
DirList  = dir('/PATH/WITH/TOMOGRAM-FOLDERS/');



% Get folders with tomos
% Extract only the folders, not regular files.
DirList  = DirList([DirList.isdir]);  % Folders only
% Get rid of first two folders: dot and dot dot.
DirList = DirList(3:end);

%exclude bad tomos /wrong folders
excl = [4,7,8];
DirList(excl) =[];

%%
% Loop through tomo-folders

for z =1:size(DirList,1)

% Convert to txt-file / nucleus
unix(['model2point ' DirList(z,1).folder '/' DirList(z,1).name '/nuc ' DirList(z,1).folder '/' DirList(z,1).name '/imod_seg/nuc.txt']);
data = tdfread([DirList(z,1).folder '/' DirList(z,1).name '/nuc.txt'],'\t');

% Search for tomogram & get filename
     try
         dinfo = dir([DirList(z,1).folder '/' DirList(z,1).name '/imod_seg/*','.rec']);
         [~,file,suffix] = fileparts(dinfo(1).name);
     catch
         dinfo = dir([DirList(z,1).folder '/' DirList(z,1).name '/imod_seg/*','_rec.mrc']);
         [~,file,suffix] = fileparts(dinfo(1).name);
     end
 
% Load tomo and show snapshot
tomo = tom_mrcread([DirList(z,1).folder '/' DirList(z,1).name '/imod_seg/' file suffix]);
tomo = tomo.Value;

% If tomo is rotated:
%tomo = permute(tomo,[1,3,2]);

figure(1)
imagesc(tomo(:,:,round(end/2,0))); colormap('gray');
pause(1)
x_dim = size(tomo,1);
y_dim = size(tomo,2);
z_dim = size(tomo,3);

% TESTCASE fixed paths
data = tdfread(['Y:TEMP\Matthias\nuc.txt'],'\t');
name = fieldnames(data);
pointlist = data.(name{1});
pointlist = pointlist(:,[1,3,2]);
x_dim = 1024;
y_dim = 1440;
z_dim = 336;
counter =1;
counter2 =1;
flag = 0;

% find out min z value and max z value
z_min = min(pointlist(:,3));
z_max = max(pointlist(:,3));

for k =1:(z_max-z_min + 1)
   interpolated_list{k}.z = z_min +k-1;
   interpolated_list{k}.xy = [];
end



%% Convert imod mask
for k=1:size(pointlist,1)
    % restrict z to volumesize
    if pointlist(k,3)>z_dim
        pointlist(k,3) = z_dim;
    elseif pointlist(k,3)< 1
        pointlist(k,3) = 1;
    end
    
    % save z height of previous coordinate
    if k~= size(pointlist,1)
        z_height_next = pointlist(k+1,3);
    elseif k == size(pointlist,1)
        z_height_next = -1;
    end
    
    z_height = pointlist(k,3);
    
    if z_height_next == z_height
        
        % restrict coordinates to volume dimensions
        if pointlist(k,1)>x_dim
            pointlist(k,1) = x_dim;
        elseif pointlist(k,1)<1
            pointlist(k,1) = 1;
        end

        if pointlist(k,2)>y_dim
            pointlist(k,2) = y_dim;
        elseif pointlist(k,2)<1
            pointlist(k,2) = 1;
        end
    
     
    % get all coordinates from one plane
    z_coord(counter,1) = pointlist(k,1);
    z_coord(counter,2) = pointlist(k,2);
    counter = counter +1;
    
    elseif z_height_next ~= z_height
        % find index of respective slice to insert into
        % interpolate_list
        for m = 1:size(interpolated_list,2)
            if z_height == interpolated_list{1,m}.z
                index = m;
            end
        end
            
        for l = 1: size(z_coord,1)
            
            % when last l reached compare to l=1
            if l == size(z_coord,1)
                dist1 = abs(z_coord(l,1) - z_coord(1,1));
                dist2 = abs(z_coord(l,2) - z_coord(1,2));
                flag = 1;
            else
                % find out distance between points
                dist1 = abs(z_coord(l,1) - z_coord(l+1,1));
                dist2 = abs(z_coord(l,2) - z_coord(l+1,2));
            end
            
            % get the bigger distance for interpolating
            if dist1 >= dist2
                dist_max = dist1;
            elseif dist1 < dist2
                dist_max = dist2;
            end
            
            if dist_max > 30
                disp('warning')
            end
            
            % if distance bigger than 10 interpolate between points
            if dist_max > 10
                if flag == 0
                    x = [round(linspace(z_coord(l,1),z_coord(l+1,1),dist_max));
                         round(linspace(z_coord(l,2),z_coord(l+1,2),dist_max))];
                    interpolated_list{index}.xy = [interpolated_list{index}.xy; x'];
                elseif flag == 1
                    x = [round(linspace(z_coord(l,1),z_coord(1,1),dist_max));
                         round(linspace(z_coord(l,2),z_coord(1,2),dist_max))];
                    interpolated_list{index}.xy = [interpolated_list{index}.xy; x'];
                    flag = 0;
                end
            else
                interpolated_list{index}.xy = [interpolated_list{index}.xy; z_coord(l,:)];
            end

        end
        counter =1;
        disp(['slice ' num2str(counter2) '  of ' num2str(z_max-z_min + 1) ' done']);
        counter2 = counter2 + 1;
        clear z_coord;
        
    end
    
   
end
%% Remove pixels at edge from structure

% Remove datapoints at edge of tomogram & round values
for k=1:size(interpolated_list,2)
    for l = size(interpolated_list{1,k}.xy,1):-1:1
        if interpolated_list{1,k}.xy(l,1)  <= 1 || interpolated_list{1,k}.xy(l,1) >= x_dim || interpolated_list{1,k}.xy(l,2)  <= 1 || interpolated_list{1,k}.xy(l,2) >= y_dim
            interpolated_list{1,k}.xy(l,:) = [];
        else
            interpolated_list{1,k}.xy(l,1) = round(interpolated_list{1,k}.xy(l,1));
            interpolated_list{1,k}.xy(l,2) = round(interpolated_list{1,k}.xy(l,2));
        end
    end
 end

%scatter(interpolated_list{1,1}.xy(:,1),interpolated_list{1,1}.xy(:,2));
figure(2)
scatter3(pointlist(:,1),pointlist(:,2),pointlist(:,3));
pause(3)

%% Calculate surface slice by slice

for k=1:size(interpolated_list,2)
    
    d = diff([interpolated_list{1,k}.xy(1:end,1),interpolated_list{1,k}.xy(1:end,2)]);
    
    
    
    d_abs = abs(d);
    d_mean = mean(d_abs);
    
    for l = 1:size(d,1)
        if d_abs(l,1) > 20 * d_mean(1) || d_abs(l,2) > 20 * d_mean(2)
            disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%');
            disp('WARNING, detected value higher than 20x mean distance between points');
            disp(['Value is    ' num2str(d_abs(l,:))])
            disp('Point will be remove but needs to be double checked');
            disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%');
            figure(4);plot(d)
            
            d(l,1) = 0; 
            d(l,2) = 0;
        end
    end
            
    figure(3)
    axis([0 1028 0 1028])
    scatter(interpolated_list{1,k}.xy(1:end,1),interpolated_list{1,k}.xy(1:end,2));
    axis([0 1028 0 1028])
    pause(0.01)
    length = sum(sqrt(sum(d.*d,2))); % length in pixels
    slice_length(k) = length * pixelsize; % length of contour times pixelsize(um)
    
end



area(z) = sum(slice_length(:)) * pixelsize; % length of contour muliplied by height

disp(['%%%%%%%%%%%%% running number ' num2str(z)])
disp(['Tomo-dimensions ' num2str(x_dim) ' ' num2str(y_dim) ' ' num2str(z_dim)])
disp(['Measured area is ' num2str(round(area(z),2)) ' um^2 from ' num2str(size(interpolated_list,2)) ' slices'])
disp(['%%%%%%%%%%%%%'])


clearvars -except z area DirList

end
