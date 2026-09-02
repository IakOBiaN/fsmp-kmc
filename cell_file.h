#include <cerrno>

// The .cell unit cell file: the plain text format the Studio saves and the
// format of the reference cells in samples/cells/.
//
//   <n_molecules>
//   <cell_x> <cell_y> [comment]
//   x  y  phi              one line per molecule; angstroms and degrees
//
// The cell is handed on in the chained polar form the generator and the
// optimizer both use: the molecule count, the cell sides, then (r, theta, phi)
// for every molecule, each position measured from the previous one and the
// first from the cell origin.

static void cell_error(const string & file, int line, const string & msg)
{
	cerr << "ERROR: " << file << ":" << line << ": " << msg << endl;
	exit(1);
}

static double cell_number(const string & file, int line, const string & value)
{
	char * end = 0;
	errno = 0;
	double v = strtod(value.c_str(), &end);
	if (end == value.c_str() || *end != '\0' || !isfinite(v))
	{
		cell_error(file, line, "\"" + value + "\" is not a number");
	}
	return v;
}

void read_cell_file(const string & path, vector<double> & params)
{
	ifstream in(path.c_str());
	if (!in)
	{
		cerr << "ERROR: cannot open the unit cell file \"" << path << "\"" << endl;
		cerr << "The structure key takes the path to a .cell file (the reference cells "
		     << "are in samples/cells/) or \"calculate\" to optimize the rough cell "
		     << "given by the unit_cell key." << endl;
		exit(1);
	}

	vector<string> lines;
	string line;
	while (getline(in, line))
	{
		while (!line.empty() && (line[line.size() - 1] == '\r' || line[line.size() - 1] == '\n'))
		{
			line.erase(line.size() - 1);
		}
		lines.push_back(line);
	}
	if (lines.size() < 2)
	{
		cell_error(path, 1, "expected the molecule count on line 1 and \"cell_x cell_y\" on line 2");
	}

	stringstream head(lines[0]);
	string tok;
	if (!(head >> tok)) { cell_error(path, 1, "expected the number of molecules"); }
	double count_d = cell_number(path, 1, tok);
	if (count_d != floor(count_d) || count_d < 1)
	{
		cell_error(path, 1, "the number of molecules must be a positive whole number");
	}
	int count = (int)count_d;

	stringstream sides(lines[1]);
	string sx, sy;
	if (!(sides >> sx) || !(sides >> sy))
	{
		cell_error(path, 2, "expected \"cell_x cell_y\" (a comment may follow on the same line)");
	}
	double cell_x = cell_number(path, 2, sx);
	double cell_y = cell_number(path, 2, sy);
	if (cell_x <= 0 || cell_y <= 0)
	{
		cell_error(path, 2, "the cell sides must be positive");
	}
	string comment;
	if (sides >> ws) { getline(sides, comment); }

	params.clear();
	params.push_back((double)count);
	params.push_back(cell_x);
	params.push_back(cell_y);
	const double to_deg = 45.0 / atan(1.0);
	double px = 0, py = 0;
	for (int i = 0; i < count; i++)
	{
		int lineno = 3 + i;
		if ((int)lines.size() < lineno)
		{
			stringstream msg;
			msg << "expected " << count << " molecules, the file ends after " << i;
			cell_error(path, (int)lines.size(), msg.str());
		}
		stringstream mol(lines[lineno - 1]);
		string a, b, c;
		if (!(mol >> a) || !(mol >> b) || !(mol >> c))
		{
			cell_error(path, lineno, "expected \"x y phi\"");
		}
		double x = cell_number(path, lineno, a);
		double y = cell_number(path, lineno, b);
		double phi = cell_number(path, lineno, c);
		double dx = x - px, dy = y - py;
		double r = sqrt(dx * dx + dy * dy);
		params.push_back(r);
		params.push_back(r > 1e-12 ? atan2(dy, dx) * to_deg : 0.0);
		params.push_back(phi);
		px = x;
		py = y;
	}

	cout << "Unit cell read from " << path << ": " << count << " molecules in a "
	     << cell_x << " x " << cell_y << " A cell";
	if (!comment.empty()) { cout << " (" << comment << ")"; }
	cout << endl;
}

// Writes the cell described by params, with every molecule wrapped into the
// cell, in the same format read_cell_file reads.
void write_cell_file(const string & path, const vector<double> & params, const string & comment)
{
	ofstream out(path.c_str());
	if (!out)
	{
		cerr << "ERROR: cannot write the unit cell file \"" << path << "\"" << endl;
		exit(1);
	}
	int count = (int)params[0];
	double cell_x = params[1], cell_y = params[2];
	out << count << endl;
	out << fixed << setprecision(6) << cell_x << " " << cell_y;
	if (!comment.empty()) { out << "  " << comment; }
	out << endl;
	const double to_rad = atan(1.0) / 45.0;
	double x = 0, y = 0;
	for (int i = 0; i < count; i++)
	{
		int n = 3 + i * 3;
		x += params[n] * cos(params[n + 1] * to_rad);
		y += params[n] * sin(params[n + 1] * to_rad);
		double wx = fmod(fmod(x, cell_x) + cell_x, cell_x);
		double wy = fmod(fmod(y, cell_y) + cell_y, cell_y);
		double phi = fmod(fmod(params[n + 2], 360.0) + 360.0, 360.0);
		out << setprecision(6) << setw(12) << wx << " " << setw(12) << wy
		    << " " << setprecision(4) << setw(10) << phi << endl;
	}
}
