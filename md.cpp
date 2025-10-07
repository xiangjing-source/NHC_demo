#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <random>
#include <string>
#include <iomanip>
#include <map>
#include <cstdlib>
#include <sstream>
#include <chrono>

using namespace std;

// ------------------- Atomic structure -------------------
struct Atom {
    double x[3];
    double v[3];
    double f[3];
};

vector<Atom> atoms;
int N = 0;

// ------------------- Simulation parameters (defaults, can be overridden by input file) -------------------
double dt = 1.0;        // fs
int nsteps = 100000;
int thermo = 100;
double T_target = 180.0;
int M = 3;               // NHC chain length
double Tdamp = 5.0;      // fs
string units = "real";

// Periodic box
double Lx=21.04, Ly=21.04, Lz=21.04;

// Nose–Hoover chain variables
vector<double> xi, eta, Q;

// Lennard-Jones parameters
double epsilon = 0.234;
double sigma   = 3.504;
double rc      = 8.76;
double mass    = 39.948;

// Boltzmann constant
double kB = 0.0019872041; // real 单位

// File IO
string coords_file = "coords.data";
string output_file = "output.log";
string plot_file   = "md_energy_temp.png";

// ------------------- Periodic boundaries -------------------
inline double pbc(double x, double L) {
    while (x > 0.5*L) x -= L;
    while (x < -0.5*L) x += L;
    return x;
}

// ------------------- Read input file -------------------
void read_input(const string &filename) {
    ifstream fin(filename);
    if(!fin.is_open()) {
        cerr << "Error: cannot open input file " << filename << endl;
        exit(1);
    }
    string key;
    while(fin >> key) {
        if(key[0]=='#') { fin.ignore(1e6,'\n'); continue; }
        if(key=="units") fin >> units;
        else if(key=="dt") fin >> dt;
        else if(key=="steps") fin >> nsteps;
        else if(key=="thermo") fin >> thermo;
        else if(key=="T_target") fin >> T_target;
        else if(key=="nhc_chain") fin >> M;
        else if(key=="Tdamp") fin >> Tdamp;
        else if(key=="Lx") fin >> Lx;
        else if(key=="Ly") fin >> Ly;
        else if(key=="Lz") fin >> Lz;
        else if(key=="epsilon") fin >> epsilon;
        else if(key=="sigma") fin >> sigma;
        else if(key=="rc") fin >> rc;
        else if(key=="mass") fin >> mass;
        else if(key=="coords") fin >> coords_file;
        else if(key=="output") fin >> output_file;
        else if(key=="plot") fin >> plot_file;
        else { fin.ignore(1e6,'\n'); }
    }
    fin.close();

    // set units
    if(units=="lj") kB = 1.0;
    else kB = 0.0019872041;
}
// ------------------- Read box dimensions from coords.data -------------------
void read_box(const string &filename) {
    ifstream fin(filename);
    if(!fin.is_open()) { cerr << "Cannot open file\n"; exit(1); }

    string line;
    getline(fin, line); // first line: number of atoms, skip

    getline(fin, line); // second line: Lattice
    size_t pos1 = line.find("Lattice=\"");
        if(pos1 != string::npos){
        pos1 += 9; // 跳过 Lattice="
        size_t pos2 = line.find("\"", pos1);
        string lattice_str = line.substr(pos1, pos2 - pos1);
        double lx, m1, m2, m3, ly, m4, m5, m6, lz;
        istringstream iss(lattice_str);
        iss >> lx >> m1 >> m2 >> m3 >> ly >> m4 >> m5 >> m6 >> lz;
        Lx = lx;
        Ly = ly;
        Lz = lz;
    } else {
        cerr << "Cannot parse lattice information, using default box sizes\n";
    }

    fin.close();
}

// ------------------- Read initial structure (extended XYZ) -------------------
void read_xyz_extended(const string &filename) {
    ifstream fin(filename);
    if(!fin.is_open()) { cerr << "Cannot open coords file\n"; exit(1); }

    string line;
    getline(fin, line);
    N = stoi(line);
    getline(fin, line); // comment

    atoms.resize(N);
    string sym;
    for(int i=0;i<N;i++){
        fin >> sym >> atoms[i].x[0] >> atoms[i].x[1] >> atoms[i].x[2];
        for(int d=0; d<3; d++){ atoms[i].v[d]=0; atoms[i].f[d]=0; }
    }
    fin.close();
}

// ------------------- Initialize velocities -------------------
void init_velocities() {
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    mt19937 gen(seed);
    double stdv = sqrt(kB*T_target/mass);
    normal_distribution<double> dist(0.0,stdv);

    for(int i=0;i<N;i++)
        for(int d=0;d<3;d++)
            atoms[i].v[d] = dist(gen);

    // remove center-of-mass velocity
    double vcm[3]={0,0,0};
    for(int i=0;i<N;i++)
        for(int d=0;d<3;d++) vcm[d]+=atoms[i].v[d];
    for(int d=0;d<3;d++) vcm[d]/=N;
    for(int i=0;i<N;i++)
        for(int d=0;d<3;d++) atoms[i].v[d]-=vcm[d];
}

// ------------------- Neighbor list -------------------
struct NeighborList {
    vector<vector<int>> neigh;
    double skin;
    NeighborList(double skin_):skin(skin_) {}
    void build() {
        neigh.assign(N, vector<int>());
        for(int i=0;i<N;i++){
            for(int j=i+1;j<N;j++){
                double dx = pbc(atoms[i].x[0]-atoms[j].x[0],Lx);
                double dy = pbc(atoms[i].x[1]-atoms[j].x[1],Ly);
                double dz = pbc(atoms[i].x[2]-atoms[j].x[2],Lz);
                double r2 = dx*dx + dy*dy + dz*dz;
                if(r2 < (rc+skin)*(rc+skin)) {
                    neigh[i].push_back(j);
                    neigh[j].push_back(i);
                }
            }
        }
    }
};

NeighborList nlist(0.2); // skin=0.2Å

// ------------------- Force computation -------------------
double compute_forces() {
    for(int i=0;i<N;i++) for(int d=0;d<3;d++) atoms[i].f[d]=0;
    double epot=0;
    double shift = 4*epsilon*(pow(sigma/rc,12)-pow(sigma/rc,6));
    for(int i=0;i<N;i++){
        for(int j:nlist.neigh[i]){
            if(j>i){
                double dx = pbc(atoms[i].x[0]-atoms[j].x[0],Lx);
                double dy = pbc(atoms[i].x[1]-atoms[j].x[1],Ly);
                double dz = pbc(atoms[i].x[2]-atoms[j].x[2],Lz);
                double r2 = dx*dx + dy*dy + dz*dz;
                if(r2<rc*rc){
                    double r2i = 1.0/r2;
                    double r6i = pow(sigma*sigma*r2i,3);
                    double ff = 48*epsilon*r6i*(r6i-0.5)*r2i;
                    atoms[i].f[0]+=ff*dx; atoms[i].f[1]+=ff*dy; atoms[i].f[2]+=ff*dz;
                    atoms[j].f[0]-=ff*dx; atoms[j].f[1]-=ff*dy; atoms[j].f[2]-=ff*dz;
                    epot += 4*epsilon*r6i*(r6i-1) - shift;
                }
            }
        }
    }
    return epot;
}

// ------------------- Kinetic energy and temperature -------------------
double compute_kinetic(){
    double vcm[3] = {0.0, 0.0, 0.0};
    for(int i=0;i<N;i++) for(int d=0; d<3; d++) vcm[d] += atoms[i].v[d];
    for(int d=0; d<3; d++) vcm[d] /= (double)N;

    double ekin = 0.0;
    for(int i=0;i<N;i++){
        for(int d=0; d<3; d++){
            double vv = atoms[i].v[d] - vcm[d];
            ekin += 0.5 * mass * vv * vv;
        }
    }
    return ekin;
}
double compute_temperature(double ekin){
    return (2.0*ekin)/((3.0*N-3)*kB);
}

// ------------------- NHC half-step update -------------------
void nhc_halfstep(double &ekin){
    double dt2 = dt * 0.5;
    double dt4 = dt2 * 0.5;
    double dt8 = dt4 * 0.5;

    double G = xi[M-2]*xi[M-2]/Q[M-2] - kB*T_target;
    xi[M-1] += dt4 * G;

    for(int m = M-2; m >= 0; m--){
        double tmp = exp(-dt8 * xi[m+1]/Q[m+1]);
        if(m==0) G = 2.0*ekin - (3.0*N-3)*kB*T_target;
        else     G = xi[m-1]*xi[m-1]/Q[m-1] - kB*T_target;
        xi[m] = tmp * (tmp * xi[m] + dt4*G);
    }
    for(int m = M-1; m >=0; m--) eta[m] += dt2*xi[m]/Q[m];

    double scale = exp(-dt2*xi[0]/Q[0]);
    for(int i=0;i<N;i++) for(int d=0;d<3;d++) atoms[i].v[d] *= scale;

    for(int m =0; m<M-1; m++){
        double tmp = exp(-dt8*xi[m+1]/Q[m+1]);
        if(m==0) G = 2.0*ekin*scale*scale - (3.0*N-3)*kB*T_target;
        else     G = xi[m-1]*xi[m-1]/Q[m-1] - kB*T_target;
        xi[m] = tmp*(tmp*xi[m]+dt4*G);
    }

    G = xi[M-2]*xi[M-2]/Q[M-2] - kB*T_target;
    xi[M-1] += dt4*G;

    ekin = compute_kinetic();
}

// ------------------- Main program -------------------
int main(int argc, char* argv[]){
    if(argc<2){ cerr<<"Usage: ./md input.in"<<endl; return 1; }
    read_input(argv[1]);
    read_xyz_extended(coords_file);
    read_box(coords_file);
    init_velocities();

    xi.assign(M,0.0); eta.assign(M,0.0); Q.assign(M,0.0);
    int Nf = 3*N - 3;
    for(int i=0;i<M;i++){ Q[i] = kB*T_target*Tdamp*Tdamp; }
    Q[0] = Nf*kB*T_target*Tdamp*Tdamp;

    nlist.build();
    double epot = compute_forces();
    double ekin = compute_kinetic();

    ofstream fout(output_file);
    fout << "# Step   E_pot   E_kin   E_tot   Temp\n";

    for(int step=0; step<=nsteps; step++){
        if(step%10==0) nlist.build();
        nhc_halfstep(ekin);

        for(int i=0;i<N;i++)
            for(int d=0;d<3;d++)
                atoms[i].v[d] += 0.5*dt*atoms[i].f[d]/mass;

        for(int i=0;i<N;i++)
            for(int d=0;d<3;d++){
                atoms[i].x[d] += dt*atoms[i].v[d];
                if(atoms[i].x[d]<0) atoms[i].x[d]+= (d==0?Lx:(d==1?Ly:Lz));
                if(atoms[i].x[d]>(d==0?Lx:(d==1?Ly:Lz))) atoms[i].x[d]-=(d==0?Lx:(d==1?Ly:Lz));
            }

        epot = compute_forces();

        for(int i=0;i<N;i++)
            for(int d=0;d<3;d++)
                atoms[i].v[d] += 0.5*dt*atoms[i].f[d]/mass;

        nhc_halfstep(ekin);

        double temp = compute_temperature(ekin);
        double etot = ekin + epot;

        if(step%thermo==0){
            fout << step << "  "
                 << setw(12)<<epot<<"  "
                 << setw(12)<<ekin<<"  "
                 << setw(12)<<etot<<"  "
                 << setw(12)<<temp<<"\n";
        }
    }
    fout.close();

    // Automatically call Python plotting
    string cmd = "python3 data.py";
    int ret = system(cmd.c_str());
    if(ret != 0){
        cerr << "Warning: gnuplot command failed with code " << ret << endl;
    }


    cout << "Simulation finished. Results saved in " << output_file
         << " and " << plot_file << endl;
    return 0;
}
