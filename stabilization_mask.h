// Stabilization mask: a periodic lattice of free wells built automatically
// from the generated initial structure (a generalization of the square-well
// network of the FsMP/kMC paper, which had to be laid out by hand for every
// polymorph). A molecule within mask_free_radius of the nearest lattice site
// feels nothing; moving further costs a penalty that rises smoothly over
// mask_ramp_width up to mask_penalty (J/mol). The penalty is multiplied by
// the local damping lambda(x), so it vanishes in the ideal gas and acts fully
// in the crystal: pores of a metastable polymorph cannot hold guest molecules
// and off-lattice phases cannot nucleate at the interface, while on-site
// molecules feel exactly nothing (zero energy, zero gradient).
//
// The mask is one unit cell of sites tiled periodically over the whole box, so
// a drift of the crystal by a whole lattice vector is degenerate. The sites
// follow every box rescaling (see pressure_balance.h), staying glued to the
// crystal. The mask enters only the Monte Carlo energies through the
// external-field slot of the molecule; like the well network in the paper it
// is kept out of the analytical pressure bookkeeping, which is justified
// because in sampled configurations molecules sit at the well bottoms.

bool mask_ready = false;
static double mask_cell_x = 0, mask_cell_y = 0;   // lattice periods, A
static vector<double> mask_site_x, mask_site_y;   // sites of one unit cell, absolute A

// Remember the site lattice from the params[0] molecules of the first generated
// unit cell; called by generate_elongated_cell once the coordinates are final
void mask_build(vector <double> &params, vector <state> &coordinates)
{
	if (!stabilization_mask) { return; }
	mask_cell_x = params[1];
	mask_cell_y = params[2];
	mask_site_x.clear();
	mask_site_y.clear();
	for (int k = 0; k < (int)params[0]; k++)
	{
		mask_site_x.push_back(coordinates[k].x);
		mask_site_y.push_back(coordinates[k].y);
	}
	mask_ready = true;
	cout << "Stabilization mask: " << mask_site_x.size() << " sites per "
	     << mask_cell_x << " x " << mask_cell_y << " A cell, free radius "
	     << mask_free_radius << " A, ramp " << mask_ramp_width << " A, penalty "
	     << mask_penalty / 1000.0 << " kJ/mol" << endl;
}

// Keep the sites glued to the crystal when the box is rescaled
void mask_rescale(double factor_x, double factor_y)
{
	if (!mask_ready) { return; }
	mask_cell_x *= factor_x;
	mask_cell_y *= factor_y;
	for (size_t k = 0; k < mask_site_x.size(); k++)
	{
		mask_site_x[k] *= factor_x;
		mask_site_y[k] *= factor_y;
	}
}

// The undamped penalty: zero within mask_free_radius of the nearest site, then
// a smoothstep rise over mask_ramp_width up to the mask_penalty plateau
double mask_penalty_at(double x, double y)
{
	double r1 = mask_free_radius + mask_ramp_width;
	double best_2 = r1 * r1;
	for (size_t k = 0; k < mask_site_x.size(); k++)
	{
		double dx = x - mask_site_x[k];
		double dy = y - mask_site_y[k];
		dx -= mask_cell_x * floor(dx / mask_cell_x + 0.5);   // nearest periodic image
		dy -= mask_cell_y * floor(dy / mask_cell_y + 0.5);
		double d_2 = dx * dx + dy * dy;
		if (d_2 < best_2) { best_2 = d_2; }
	}
	if (best_2 >= r1 * r1) { return mask_penalty; }
	if (best_2 <= mask_free_radius * mask_free_radius) { return 0.0; }
	double s = (sqrt(best_2) - mask_free_radius) / mask_ramp_width;
	return mask_penalty * s * s * (3.0 - 2.0 * s);
}

// The external field acting on a molecule at (x, y): the smooth potential
// u_ext plus the damped stabilization mask (zero in the ideal gas)
results external_field_and_mask(double x, double y, double &Lx)
{
	results field = external_field(x, Lx);
	if (mask_ready)
	{
		double sqrt_lambda = damping_field(x, Lx);
		field.energy += sqrt_lambda * sqrt_lambda * mask_penalty_at(x, y);
	}
	return field;
}
