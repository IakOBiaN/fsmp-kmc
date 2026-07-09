// Append one frame to an xyz file. Every molecule is drawn with the atomistic
// model from the parameter file (the molecule_model key), rotated by the
// molecule orientation phi and translated to the molecule center. The molecule
// is kept whole: atoms of a molecule sitting near the box edge may stick
// slightly outside instead of being wrapped atom by atom.
void write_xyz_file (string name, int &nPart, double & /*density*/, double &Lx, double &Ly, double & /*temperature*/, vector <state> &coordinates, int /*frame*/, double /*distance*/, bool init)
{
	if (init) {ofstream fileOutput(name.c_str(), ios_base::trunc);fileOutput.close();}
	ofstream fileOutput(name.c_str(), ios_base::app);
	fileOutput << nPart * (int)molecule_model.size() << endl;
	fileOutput << "Lattice=\"" << Lx << " 0 0 0 " << Ly << " 0 0 0 1\"" << endl;

	for(int i = 0; i < nPart; i++)
		{
			for (size_t a = 0; a < molecule_model.size(); a++)
			{
				double x = coordinates[i].x + molecule_model[a].x * coordinates[i].cos_phi - molecule_model[a].y * coordinates[i].sin_phi;
				double y = coordinates[i].y + molecule_model[a].x * coordinates[i].sin_phi + molecule_model[a].y * coordinates[i].cos_phi;
				fileOutput << molecule_model[a].element << " " << x << " " << y << " " << molecule_model[a].z << endl;
			}
		}
	fileOutput.close();
}
