// Read the simulation parameters from a "key = value" text file (see configs/
// for ready-to-run examples). Everything after # on a line is a comment.
//
// The parser is strict on purpose: an unknown key, a duplicated key, a value
// that does not parse, or a missing required key stops the program with an
// error naming the file and line. A production run must never fall back to a
// silent default because of a typo.
//
// The structure key holds the path to a .cell unit cell file (the reference
// cells are in samples/cells/, and the Studio saves the same format), or the
// word "calculate", which optimizes the rough cell given by the unit_cell key.
//
// Optional keys: statistics_name (default is built from the structure name),
// seed (fixes the random sequence for reproducible runs; otherwise the wall
// clock or the FSMP_RANDOM_SEED compile flag is used), constant_pressure_value
// (default 0), unit_cell (required with structure = calculate, forbidden
// otherwise), optimize_only (default false; only with structure = calculate:
// stop right after the unit cell optimization, skipping the Monte Carlo run),
// stabilization_mask with mask_free_radius, mask_ramp_width, mask_penalty
// (default false; a lattice of free wells built from the initial structure
// that keeps a metastable polymorph intact, see the key handler below).

#include <set>
#include <map>
#include <cerrno>

static bool param_seed_given = false;
static int  param_seed = 0;

// Line each key was read from, so a range error can point at it just like a
// parse error does.
static map<string, int> param_lines;

static string trim_spaces(const string & s)
{
	size_t a = s.find_first_not_of(" \t\r\n");
	if (a == string::npos) { return ""; }
	size_t b = s.find_last_not_of(" \t\r\n");
	return s.substr(a, b - a + 1);
}

static void param_error(const string & file, int line, const string & msg)
{
	cerr << "ERROR: " << file << ":" << line << ": " << msg << endl;
	exit(1);
}

static double param_double(const string & file, int line, const string & key, const string & value)
{
	char * end = 0;
	errno = 0;
	double v = strtod(value.c_str(), &end);
	if (end == value.c_str() || *end != '\0')
	{
		param_error(file, line, "key \"" + key + "\": \"" + value + "\" is not a number");
	}
	// strtod happily returns nan and inf; they would spread through every
	// energy and never trip a comparison, so they stop the run here
	if (!isfinite(v))
	{
		param_error(file, line, "key \"" + key + "\": \"" + value + "\" is not a finite number");
	}
	if (errno == ERANGE && v != 0.0)
	{
		param_error(file, line, "key \"" + key + "\": \"" + value + "\" is too large to represent");
	}
	return v;
}

static int param_int(const string & file, int line, const string & key, const string & value)
{
	char * end = 0;
	errno = 0;
	long v = strtol(value.c_str(), &end, 10);
	if (end == value.c_str() || *end != '\0')
	{
		param_error(file, line, "key \"" + key + "\": \"" + value + "\" is not an integer");
	}
	// without this, a value above 2^31 would be silently truncated, and the
	// run would use a number the file never mentioned
	if (errno == ERANGE || v < INT_MIN || v > INT_MAX)
	{
		param_error(file, line, "key \"" + key + "\": \"" + value + "\" does not fit in a 32-bit integer");
	}
	return (int)v;
}

// A value that parsed but cannot be run with. The message points at the line
// the key came from.
static void param_require(bool ok, const string & file, const string & key, const string & msg)
{
	if (!ok)
	{
		int line = param_lines.count(key) ? param_lines[key] : 0;
		param_error(file, line, "key \"" + key + "\" " + msg);
	}
}

static bool param_bool(const string & file, int line, const string & key, const string & value)
{
	if (value == "true" || value == "yes" || value == "1") { return true; }
	if (value == "false" || value == "no" || value == "0") { return false; }
	param_error(file, line, "key \"" + key + "\": \"" + value + "\" is not a boolean (use true/false)");
	return false;
}

void read_parameters(const char * path)
{
	ifstream in(path);
	if (!in)
	{
		cerr << "ERROR: cannot open parameter file \"" << path << "\"" << endl;
		exit(1);
	}
	string file = path;
	set<string> seen;
	string line;
	int lineno = 0;
	while (getline(in, line))
	{
		lineno++;
		size_t hash = line.find('#');
		if (hash != string::npos) { line.erase(hash); }
		line = trim_spaces(line);
		if (line.empty()) { continue; }
		size_t eq = line.find('=');
		if (eq == string::npos) { param_error(file, lineno, "expected \"key = value\""); }
		string key = trim_spaces(line.substr(0, eq));
		string value = trim_spaces(line.substr(eq + 1));
		if (key.empty() || value.empty()) { param_error(file, lineno, "expected \"key = value\""); }
		if (!seen.insert(key).second) { param_error(file, lineno, "duplicate key \"" + key + "\""); }
		param_lines[key] = lineno;

		if      (key == "potential")                      { p_name = value; }
		else if (key == "structure")                      { structure_name = value; }
		else if (key == "sigma_mode")
		{
			if (value != "manual" && value != "min_dist" && value != "molecule_area")
			{
				param_error(file, lineno, "key \"sigma_mode\" must be manual, min_dist or molecule_area");
			}
			sigma_mode = value;
		}
		else if (key == "sigma")                          { sigma_manual = param_double(file, lineno, key, value); }
		else if (key == "temp_from")                      { temp_from = param_double(file, lineno, key, value); }
		else if (key == "temp_to")                        { temp_to = param_double(file, lineno, key, value); }
		else if (key == "temp_step")                      { temp_step = param_double(file, lineno, key, value); }
		else if (key == "um_from")                        { um_from = param_double(file, lineno, key, value); }
		else if (key == "um_to")                          { um_to = param_double(file, lineno, key, value); }
		else if (key == "um_step")                        { um_step = param_double(file, lineno, key, value); }
		else if (key == "temperature_in_transition_zone") { temperature_in_transition_zone = param_double(file, lineno, key, value); }
		else if (key == "lambdam")                        { lambdam = param_double(file, lineno, key, value); }
		else if (key == "nSteps")                         { nSteps = param_int(file, lineno, key, value); }
		else if (key == "nStepsEq")                       { nStepsEq = param_int(file, lineno, key, value); }
		else if (key == "constant_pressure")              { constant_pressure = param_bool(file, lineno, key, value); }
		else if (key == "constant_pressure_value")        { constant_pressure_value = param_double(file, lineno, key, value); }
		else if (key == "kMC")                            { kMC = param_bool(file, lineno, key, value); }
		// Stabilization mask: a periodic lattice of free wells built automatically
		// from the initial structure (see stabilization_mask.h). Leaving the
		// mask_free_radius neighborhood of a lattice site costs up to mask_penalty
		// (J/mol) over a smooth ramp of mask_ramp_width (A). Damped by lambda(x),
		// so the ideal gas phase is unaffected. Keeps metastable porous polymorphs
		// (chicken-wire, flower phases) from decaying: pores cannot hold guest
		// molecules and off-lattice phases cannot nucleate at the interface.
		else if (key == "stabilization_mask")             { stabilization_mask = param_bool(file, lineno, key, value); }
		else if (key == "mask_free_radius")               { mask_free_radius = param_double(file, lineno, key, value); }
		else if (key == "mask_ramp_width")                { mask_ramp_width = param_double(file, lineno, key, value); }
		else if (key == "mask_penalty")                   { mask_penalty = param_double(file, lineno, key, value); }
		// Stop after the unit cell optimization: the optimizer writes its xyz
		// animation and prints the optimized cell, then the program exits without
		// entering the Monte Carlo loop. Only meaningful with structure = calculate.
		else if (key == "optimize_only")                  { optimize_only = param_bool(file, lineno, key, value); }
		else if (key == "uc_in_x")                        { uc_in_x = param_int(file, lineno, key, value); }
		else if (key == "uc_in_y")                        { uc_in_y = param_int(file, lineno, key, value); }
		else if (key == "free_space")                     { free_space = param_double(file, lineno, key, value); }
		else if (key == "molecule_model")                 { molecule_model_file = value; }
		// The ray-based visualization was replaced by atomistic models: molecules
		// are now drawn with the xyz model given by the molecule_model key.
		else if (key == "total_molecule_directions" || key == "angle_1" || key == "angle_2")
		{
			param_error(file, lineno, "key \"" + key + "\" was removed: point molecule_model at an xyz model instead (see samples/models/)");
		}
		else if (key == "delta")                          { delta = param_double(file, lineno, key, value); }
		else if (key == "delta_angle")                    { delta_angle = param_double(file, lineno, key, value); }
		else if (key == "widom_test_index")               { widom_test_index = param_bool(file, lineno, key, value); }
		else if (key == "unit_cell_name")                 { unit_cell_name = value; }
		else if (key == "xyz_name")                       { xyz_name = value; }
		else if (key == "statistics_name")                { statistics_name = value; }
		else if (key == "seed")                           { param_seed = param_int(file, lineno, key, value); param_seed_given = true; }
		else if (key == "unit_cell")
		{
			stringstream ss(value);
			string tok;
			while (ss >> tok) { unit_cell_params.push_back(param_double(file, lineno, key, tok)); }
			if (unit_cell_params.size() < 6 ||
			    (int)unit_cell_params.size() != 3 + 3 * (int)unit_cell_params[0])
			{
				param_error(file, lineno, "key \"unit_cell\" must hold: mols_per_cell cell_x cell_y, then r theta phi for every molecule (3 + 3*mols_per_cell numbers)");
			}
		}
		else { param_error(file, lineno, "unknown key \"" + key + "\""); }
	}

	static const char * required[] = {
		"potential", "structure", "sigma_mode", "temp_from", "temp_to", "temp_step",
		"um_from", "um_to", "um_step", "temperature_in_transition_zone", "lambdam",
		"nSteps", "nStepsEq", "constant_pressure", "kMC", "uc_in_x", "uc_in_y",
		"free_space", "molecule_model",
		"delta", "delta_angle", "widom_test_index", "unit_cell_name", "xyz_name" };
	string missing;
	for (size_t i = 0; i < sizeof(required) / sizeof(required[0]); i++)
	{
		if (!seen.count(required[i])) { missing += string(" ") + required[i]; }
	}
	if (!missing.empty())
	{
		cerr << "ERROR: " << file << ": missing required keys:" << missing << endl;
		exit(1);
	}
	if (sigma_mode == "manual" && !seen.count("sigma"))
	{
		cerr << "ERROR: " << file << ": sigma_mode = manual requires the sigma key (sigma in A)" << endl;
		exit(1);
	}
	if (sigma_mode != "manual" && seen.count("sigma"))
	{
		cerr << "ERROR: " << file << ": the sigma key is only used with sigma_mode = manual" << endl;
		exit(1);
	}
	if (!stabilization_mask &&
	    (seen.count("mask_free_radius") || seen.count("mask_ramp_width") || seen.count("mask_penalty")))
	{
		cerr << "ERROR: " << file << ": mask_free_radius, mask_ramp_width and mask_penalty "
		     << "are only used with stabilization_mask = true" << endl;
		exit(1);
	}
	if (stabilization_mask && (mask_free_radius <= 0 || mask_ramp_width <= 0 || mask_penalty <= 0))
	{
		cerr << "ERROR: " << file << ": mask_free_radius, mask_ramp_width and mask_penalty "
		     << "must be positive" << endl;
		exit(1);
	}
	if (structure_name == "calculate" && !seen.count("unit_cell"))
	{
		cerr << "ERROR: " << file << ": structure = calculate requires the unit_cell key" << endl;
		exit(1);
	}
	if (structure_name != "calculate" && seen.count("unit_cell"))
	{
		cerr << "ERROR: " << file << ": the unit_cell key is only used with structure = calculate "
		     << "(a cell file carries its own cell)" << endl;
		exit(1);
	}
	if (optimize_only && structure_name != "calculate")
	{
		cerr << "ERROR: " << file << ": optimize_only = true requires structure = calculate "
		     << "(a cell read from a file is used as it is)" << endl;
		exit(1);
	}

	// Values that parse but cannot be run with. Each one below either divides
	// by zero, loops forever, samples nothing, or silently produces a NaN
	// halfway through a long run, so it is refused here instead.
	param_require(temp_from > 0, file, "temp_from", "must be a positive temperature (K)");
	param_require(temp_to > 0, file, "temp_to", "must be a positive temperature (K)");
	param_require(temperature_in_transition_zone > 0, file, "temperature_in_transition_zone",
	              "must be a positive temperature (K): the damping field takes its square root");
	// the loop walks from *_from towards *_to and stops when it passes it, so
	// a zero step over a non-empty range never terminates (the sign does not
	// matter, the program takes the absolute value)
	param_require(temp_from == temp_to || temp_step != 0, file, "temp_step",
	              "must not be zero when temp_from and temp_to differ: the loop would never reach the end");
	param_require(um_from == um_to || um_step != 0, file, "um_step",
	              "must not be zero when um_from and um_to differ: the loop would never reach the end");
	param_require(lambdam >= 0, file, "lambdam", "must not be negative: it scales the damping field");
	param_require(nSteps >= 1, file, "nSteps", "must be at least one Monte Carlo step");
	param_require(nStepsEq >= 0, file, "nStepsEq", "must not be negative");
	param_require(nStepsEq <= nSteps, file, "nStepsEq",
	              "must not exceed nSteps: the run would end before anything is averaged");
	param_require(uc_in_x >= 1, file, "uc_in_x", "must be at least one unit cell");
	param_require(uc_in_y >= 1, file, "uc_in_y", "must be at least one unit cell");
	param_require(free_space >= 0 && free_space < 0.5, file, "free_space",
	              "must be in [0, 0.5): the cell length is divided by (1 - 2*free_space)");
	param_require(delta > 0, file, "delta", "must be a positive maximal displacement (A)");
	param_require(delta_angle > 0, file, "delta_angle", "must be a positive maximal rotation (deg)");
	if (seen.count("sigma"))
	{
		param_require(sigma_manual > 0, file, "sigma", "must be a positive length (A)");
	}
	if (!unit_cell_params.empty())
	{
		param_require(unit_cell_params[0] == floor(unit_cell_params[0]), file, "unit_cell",
		              "must start with a whole number of molecules per cell");
		param_require(unit_cell_params[1] > 0 && unit_cell_params[2] > 0, file, "unit_cell",
		              "must have positive cell sides");
	}

	cout << "Parameters read from " << file << endl;

	if (structure_name != "calculate")
	{
		read_cell_file(structure_name, unit_cell_params);
	}
}
