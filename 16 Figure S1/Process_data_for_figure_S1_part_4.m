

% Summary data
N_cells = 10;
tscore = zeros(N_cells,1);
taxes = zeros(N_cells,2);
rscore = zeros(N_cells,1);
raxes = zeros(N_cells,2);

for nc = 1:N_cells
    eval(['load Whole_field_cell_fits/Translation/cell_' num2str(nc) '.mat axes maxscore'])
    tscore(nc) = maxscore;
    taxes(nc,:) = axes
    eval(['load Whole_field_cell_fits/Rotation/cell_' num2str(nc) '.mat axes maxscore'])
    rscore(nc) = maxscore;
    raxes(nc,:) = axes
end

save whole_field_optic_flow_fits_calculated.mat tscore taxes rscore raxes