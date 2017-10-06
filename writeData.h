using namespace std;

void writeData (double Xdata, double Ydata, double Ydata1, double Ydata2, double Ydata3, double Ydata4, double Ydata5, double Ydata6, double Ydata7, double Ydata8, double Ydata9, double Ydata10, double Ydata11)
{
 stringstream name;
 name <<  "statistics.dat";
 ofstream fileOutput(name.str().c_str(), ios_base::app);
 fileOutput << Xdata << "\t" << Ydata << "\t" << Ydata1 << "\t" << Ydata2 << "\t" << Ydata3 << "\t" << Ydata4 << "\t" << Ydata5 << "\t" << Ydata6 << "\t" << Ydata7 << "\t" << Ydata8 << "\t" << Ydata9 << "\t" << Ydata10 << "\t" << Ydata11 << endl;
 fileOutput.close();
}
