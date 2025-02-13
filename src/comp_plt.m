function c_plt = comp_plt(c)
    % Ensure the input is properly reshaped from Nz x 1 x Nc to Nz x Nc
    c = squeeze(c);  % Remove any singleton dimension
    [Nz, Nc] = size(c);
    
    % Number of horizontal points (fixed at 100)
    Nx = 500;

    % Scale mass fractions to discrete counts summing to Nx
    fractions = round(c * Nx);

    % Adjust rounding errors row-wise to ensure sum equals Nx
    diff = Nx - sum(fractions, 2);
    [~, maxIdx] = max(fractions, [], 2); % Find index of max fraction in each row
    linearIdx = sub2ind([Nz, Nc], (1:Nz)', maxIdx); % Convert to linear indices
    fractions(linearIdx) = fractions(linearIdx) + diff; % Adjust max fraction in each row

    % Construct compositional array using cumulative sum and indexing
    c_plt = zeros(Nz, Nx, 'uint8'); % Preallocate for efficiency
    boundaries = [zeros(Nz, 1), cumsum(fractions, 2)]; % Compute boundaries for each component

    % Fill the output matrix efficiently using vectorized assignments
    x = 1:Nx; % Column indices
    for ic = 1:Nc
        mask = (x > boundaries(:, ic)) & (x <= boundaries(:, ic + 1)); % Logical mask
        c_plt(mask) = ic;
    end

    % % Plot the composition as a colormap
    % imagesc(1:Nx, 1:Nz, c_plt);
    % colormap(jet(Nc));  % Use a distinct colormap
    % colorbar;
    % xlabel('Composition Distribution');
    % ylabel('Depth');
    % title('Rock Composition Along Vertical Column');
end