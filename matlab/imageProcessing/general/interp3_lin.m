% Function to replace interp2 to be applied for image processing. The
% interpolation is always applied using the lower resolution matrix to
% compute the weights of the interpolation. This means that a high
% resolution coordinates is interpolated into the four nearest neighbours
% of the low res matrix.

% TODO: Needs to be verified for the case where the new sampling implies a
% higher matrix in x but not in y. It's implemented separately, but there
% can be some issues.

function image_out = interp3_lin(X_coord_in, Y_coord_in, Z_coord_in, image_in, X_coord_out, Y_coord_out, Z_coord_out)
    
    outOfImageValues = 0;
    
    % Steps in each grid:
    x_step_in = X_coord_in(1, 2, 1)-X_coord_in(1, 1, 1);
    y_step_in = Y_coord_in(2, 1, 1)-Y_coord_in(1, 1, 1);
    z_step_in = Z_coord_in(1, 1, 2)-Z_coord_in(1, 1, 1);
    x_step_out = X_coord_out(1, 2, 1)-X_coord_out(1, 1, 1);
    y_step_out = Y_coord_out(2, 1, 1)-Y_coord_out(1, 1, 1);
    z_step_out = Z_coord_out(1, 1, 2)-Z_coord_out(1, 1, 1);
    
    % Use the higher step for the interpolation:
    if (x_step_in > x_step_out)
        x_step = x_step_in; % low resolution step
        X_coord_hr = X_coord_out; % high resolution matrix
        X_coord_lr = X_coord_in;
        
    else
        x_step = x_step_out;
        X_coord_hr = X_coord_in;
        X_coord_lr = X_coord_out;
    end
    if (y_step_in > y_step_out)
        y_step = y_step_in;
        Y_coord_hr = Y_coord_out;
        Y_coord_lr = Y_coord_in;
    else
        y_step = y_step_out;
        Y_coord_hr = Y_coord_in;
        Y_coord_lr = Y_coord_out;
    end
    if (z_step_in > z_step_out) % Needs to change in the future to z coordinate
        z_step = z_step_in;
        Z_coord_hr = Z_coord_out;
        Z_coord_lr = Z_coord_in;
    else
        z_step = z_step_out;
        Z_coord_hr = Z_coord_in;
        Z_coord_lr = Z_coord_out;
    end
    
    % Left upper point:
    x_min = X_coord_lr(1,1,1);
    y_min = Y_coord_lr(1,1,1);
    z_min = Z_coord_lr(1,1,1);
    
    % x, in pixel coordinates:
    x_new_coord = (X_coord_hr-x_min)./x_step+1; % the new coordinates are in the lr matrix.
    % y, in pixel coordiantes:
    y_new_coord = (Y_coord_hr-y_min)./y_step+1;
    % z, in pixel coordiantes:
    z_new_coord = (Z_coord_hr-z_min)./z_step+1;
    % Don't take into account those that are more tha one pixel out:
    indicesOut = (x_new_coord < 0) | (y_new_coord < 0) | (z_new_coord < 0) | (x_new_coord >= (size(X_coord_lr,2)+1)) | (y_new_coord >= (size(Y_coord_lr,1)+1)) | (z_new_coord >= (size(Z_coord_lr,3)+1));
    x_new_coord(indicesOut) = [];
    y_new_coord(indicesOut) = [];
    z_new_coord(indicesOut) = [];
    
    % Now interpolate in x:
    x_new_coord_0 = floor(x_new_coord);
    x_new_coord_1 = x_new_coord_0 + 1;
    x_weight_0 = 1 - (x_new_coord - x_new_coord_0);
    x_weight_1 = x_new_coord - x_new_coord_0;
    % Now interpoalte in y:
    y_new_coord_0 = floor(y_new_coord);
    y_new_coord_1 = y_new_coord_0 + 1;
    y_weight_0 = 1 - (y_new_coord - y_new_coord_0);
    y_weight_1 = y_new_coord - y_new_coord_0;
    % Now interpoalte in z:
    z_new_coord_0 = floor(z_new_coord);
    z_new_coord_1 = z_new_coord_0 + 1;
    z_weight_0 = 1 - (z_new_coord - z_new_coord_0);
    z_weight_1 = z_new_coord - z_new_coord_0;
    
    % Indices out of the image, set the same value as in the edge (we previously removed the ones that were completely out, but here we have the ones that are in the edge for interpoaltion:
    indices_out_edges = x_new_coord_0<1; % pixel on the ouside side of the edge.
    x_new_coord_0(indices_out_edges) = 1; % use the pixel inside the edge.
    indices_out_edges = y_new_coord_0<1; % pixel on the ouside side of the edge.
    y_new_coord_0(indices_out_edges) = 1; % use the pixel inside the edge.
    indices_out_edges = z_new_coord_0<1; % pixel on the ouside side of the edge.
    z_new_coord_0(indices_out_edges) = 1; % use the pixel inside the edge.
    
    indices_out_edges = x_new_coord_1>size(X_coord_lr,2); % pixel on the ouside side of the edge.
    x_new_coord_1(indices_out_edges) = size(X_coord_lr,2); % use the pixel inside the edge.
    indices_out_edges = y_new_coord_1>size(Y_coord_lr,1); % pixel on the ouside side of the edge.
    y_new_coord_1(indices_out_edges) = size(Y_coord_lr,1); % use the pixel inside the edge.
    indices_out_edges = z_new_coord_1>size(Z_coord_lr,3); % pixel on the ouside side of the edge.
    z_new_coord_1(indices_out_edges) = size(Z_coord_lr,3); % use the pixel inside the edge.
    
    % Output image:
    image_out = zeros(size(X_coord_out));
    % Fill it with linear indices:
    % I have the weights for the higher resolution matrix, if the output
    % it's in this matrix, we just fill the image. If it's the opposite
    % case, use the coordinates to index the output image:
    if (x_step_in > x_step_out) %Writting in a higher resolution matrix
        % Trilinear interpolation, repeat bilinear for two slices:
        image_out(:) = z_weight_0(:).*(x_weight_0(:).*y_weight_0(:).*image_in(sub2ind(size(image_in),y_new_coord_0(:), x_new_coord_0(:), z_new_coord_0(:))) + x_weight_0(:).*y_weight_1(:).*image_in(sub2ind(size(image_in),y_new_coord_1(:), x_new_coord_0(:), z_new_coord_0(:))) + ...
            x_weight_1(:).*y_weight_0(:) .* image_in(sub2ind(size(image_in),y_new_coord_0(:), x_new_coord_1(:), z_new_coord_0(:))) + x_weight_1(:).*y_weight_1(:) .* image_in(sub2ind(size(image_in),y_new_coord_1(:), x_new_coord_1(:), z_new_coord_0(:))));
        image_out(:) = image_out(:) + z_weight_1(:).*(x_weight_0(:).*y_weight_0(:).*image_in(sub2ind(size(image_in),y_new_coord_0(:), x_new_coord_0(:), z_new_coord_1(:))) + x_weight_0(:).*y_weight_1(:).*image_in(sub2ind(size(image_in),y_new_coord_1(:), x_new_coord_0(:), z_new_coord_1(:))) + ...
            x_weight_1(:).*y_weight_0(:) .* image_in(sub2ind(size(image_in),y_new_coord_0(:), x_new_coord_1(:), z_new_coord_1(:))) + x_weight_1(:).*y_weight_1(:) .* image_in(sub2ind(size(image_in),y_new_coord_1(:), x_new_coord_1(:), z_new_coord_1(:))));
    else
        image_out(:) = image_out(:) + accumarray(sub2ind(size(X_coord_out),y_new_coord_0(:), x_new_coord_0(:), z_new_coord_0(:)), z_weight_0(:).*x_weight_0(:).*y_weight_0(:) .* image_in(~indicesOut), [numel(X_coord_out) 1]);
        image_out(:) = image_out(:) + accumarray(sub2ind(size(X_coord_out),y_new_coord_1(:), x_new_coord_0(:), z_new_coord_0(:)), z_weight_0(:).*x_weight_0(:).*y_weight_1(:) .* image_in(~indicesOut), [numel(X_coord_out) 1]);
        image_out(:) = image_out(:) + accumarray(sub2ind(size(X_coord_out),y_new_coord_0(:), x_new_coord_1(:), z_new_coord_0(:)), z_weight_0(:).*x_weight_1(:).*y_weight_0(:) .* image_in(~indicesOut), [numel(X_coord_out) 1]);
        image_out(:) = image_out(:) + accumarray(sub2ind(size(X_coord_out),y_new_coord_1(:), x_new_coord_1(:), z_new_coord_0(:)), z_weight_0(:).*x_weight_1(:).*y_weight_1(:) .* image_in(~indicesOut), [numel(X_coord_out) 1]);
        
        image_out(:) = image_out(:) + accumarray(sub2ind(size(X_coord_out),y_new_coord_0(:), x_new_coord_0(:), z_new_coord_1(:)), z_weight_1(:).*x_weight_0(:).*y_weight_0(:) .* image_in(~indicesOut), [numel(X_coord_out) 1]);
        image_out(:) = image_out(:) + accumarray(sub2ind(size(X_coord_out),y_new_coord_1(:), x_new_coord_0(:), z_new_coord_1(:)), z_weight_1(:).*x_weight_0(:).*y_weight_1(:) .* image_in(~indicesOut), [numel(X_coord_out) 1]);
        image_out(:) = image_out(:) + accumarray(sub2ind(size(X_coord_out),y_new_coord_0(:), x_new_coord_1(:), z_new_coord_1(:)), z_weight_1(:).*x_weight_1(:).*y_weight_0(:) .* image_in(~indicesOut), [numel(X_coord_out) 1]);
        image_out(:) = image_out(:) + accumarray(sub2ind(size(X_coord_out),y_new_coord_1(:), x_new_coord_1(:), z_new_coord_1(:)), z_weight_1(:).*x_weight_1(:).*y_weight_1(:) .* image_in(~indicesOut), [numel(X_coord_out) 1]);
    end

end 