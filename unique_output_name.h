// Return a file name that does not collide with an existing file: "name.ext"
// becomes "name_2.ext", "name_3.ext", and so on. An output of a previous run
// is never silently overwritten.
string unique_output_name(const string & name)
{
	ifstream probe(name.c_str());
	if (!probe.good()) { return name; }
	size_t dot = name.find_last_of('.');
	string base = (dot == string::npos) ? name : name.substr(0, dot);
	string ext  = (dot == string::npos) ? ""   : name.substr(dot);
	for (int n = 2; ; n++)
	{
		stringstream candidate;
		candidate << base << "_" << n << ext;
		ifstream busy(candidate.str().c_str());
		if (!busy.good())
		{
			cout << "NOTE: \"" << name << "\" already exists; writing to \""
			     << candidate.str() << "\" instead" << endl;
			return candidate.str();
		}
	}
}
