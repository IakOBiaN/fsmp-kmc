using namespace std;

void writeData (double Xdata, double Ydata, double Ydata1, double Ydata2, double Ydata3, double Ydata4, double Ydata5)
{
 stringstream name;
 name <<  "statistics.dat";
 ofstream fileOutput(name.str().c_str(), ios_base::app);
 fileOutput << Xdata << "\t" << Ydata << "\t" << Ydata1 << "\t" << Ydata2 << "\t" << Ydata3 << "\t" << Ydata4 << "\t" << Ydata5 << endl;
 fileOutput.close();
}
