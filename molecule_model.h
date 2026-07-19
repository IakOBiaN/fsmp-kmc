// Atomistic model of the molecule used for every xyz visualization: the unit
// cell optimization animation and the trajectory. The model is a plain xyz
// file (see samples/models/): an atom count line, a comment line, then "element x y z"
// rows in A. The model must be drawn in the frame the forcefield uses: the
// molecule center at the origin and the phi = 0 orientation along the +x axis.
// The engine only rotates the model by the molecule orientation phi and
// translates it to the molecule position; z is passed through unchanged.

class model_atom {
public:
	string element;
	double x, y, z;
};

vector<model_atom> molecule_model;

void read_molecule_model(const string & path)
{
	ifstream in(path.c_str());
	if (!in)
	{
		cerr << "ERROR: cannot open molecule model \"" << path << "\"" << endl;
		exit(1);
	}
	int n = 0;
	if (!(in >> n) || n <= 0)
	{
		cerr << "ERROR: molecule model \"" << path << "\": the first line must hold a positive atom count" << endl;
		exit(1);
	}
	string line;
	getline(in, line);   // the rest of the count line
	getline(in, line);   // the comment line
	molecule_model.clear();
	for (int i = 0; i < n; i++)
	{
		model_atom a;
		if (!(in >> a.element >> a.x >> a.y >> a.z))
		{
			cerr << "ERROR: molecule model \"" << path << "\": expected " << n
			     << " \"element x y z\" rows, failed at atom " << i + 1 << endl;
			exit(1);
		}
		molecule_model.push_back(a);
	}
	cout << "Molecule model for visualization: " << path << " (" << n << " atoms)" << endl;
}
