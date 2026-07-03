// Read the simulation parameters from a "key = value" text file (see configs/
// for ready-to-run examples). Everything after # on a line is a comment.
//
// The parser is strict on purpose: an unknown key, a duplicated key, a value
// that does not parse, or a missing required key stops the program with an
// error naming the file and line. A production run must never fall back to a
// silent default because of a typo.
//
// Optional keys: statistics_name (default is built from the structure name),
// seed (fixes the random sequence for reproducible runs; otherwise the wall
// clock or the FSMP_RANDOM_SEED compile flag is used), constant_pressure_value
// (default 0), unit_cell (required with structure = calculate, forbidden
// otherwise).

#include <set>

static bool param_seed_given = false;
static int  param_seed = 0;

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
	double v = strtod(value.c_str(), &end);
	if (end == value.c_str() || *end != '\0')
	{
		param_error(file, line, "key \"" + key + "\": \"" + value + "\" is not a number");
	}
	return v;
}

static int param_int(const string & file, int line, const string & key, const string & value)
{
	char * end = 0;
	long v = strtol(value.c_str(), &end, 10);
	if (end == value.c_str() || *end != '\0')
	{
		param_error(file, line, "key \"" + key + "\": \"" + value + "\" is not an integer");
	}
	return (int)v;
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

		if      (key == "potential")                      { p_name = value; }
		else if (key == "structure")                      { structure_name = value; }
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
		else if (key == "uc_in_x")                        { uc_in_x = param_int(file, lineno, key, value); }
		else if (key == "uc_in_y")                        { uc_in_y = param_int(file, lineno, key, value); }
		else if (key == "free_space")                     { free_space = param_double(file, lineno, key, value); }
		else if (key == "total_molecule_directions")      { total_molecule_directions = param_int(file, lineno, key, value); }
		else if (key == "angle_1")                        { angle_1 = param_double(file, lineno, key, value); }
		else if (key == "angle_2")                        { angle_2 = param_double(file, lineno, key, value); }
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
		"potential", "structure", "temp_from", "temp_to", "temp_step",
		"um_from", "um_to", "um_step", "temperature_in_transition_zone", "lambdam",
		"nSteps", "nStepsEq", "constant_pressure", "kMC", "uc_in_x", "uc_in_y",
		"free_space", "total_molecule_directions", "angle_1", "angle_2",
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
	if (structure_name == "calculate" && !seen.count("unit_cell"))
	{
		cerr << "ERROR: " << file << ": structure = calculate requires the unit_cell key" << endl;
		exit(1);
	}
	if (structure_name != "calculate" && seen.count("unit_cell"))
	{
		cerr << "ERROR: " << file << ": the unit_cell key is only used with structure = calculate "
		     << "(a named structure carries its own cell)" << endl;
		exit(1);
	}
	cout << "Parameters read from " << file << endl;
}
