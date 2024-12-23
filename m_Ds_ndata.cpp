#define CONTROL

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <complex.h>
#include <algorithm>

#include "global.hpp"
#include "resampling.hpp"
#include "resampling_new.hpp"
#include "read.hpp"
#include "m_eff.hpp"
#include "gnuplot.hpp"
#include "eigensystem.hpp"
#include "linear_fit.hpp"
#include "various_fits.hpp"
#include "mutils.hpp"
#include "functions_gm2_analysis.hpp"
#include "correlators_analysis.hpp"
#include "non_linear_fit.hpp"
#include "tower.hpp"
#include "fit_all.hpp"
#include "functions_amu.hpp"
#include "gamma_analysis.hpp"
// #include "lhs_functions.hpp"

#include <string>
#include <cstring> 
#include <string>
#include <fstream>
#include <memory>

#include <gsl/gsl_integration.h>

// equal=OS
// opposite=TM

constexpr double MDs_MeV = 1967.0;
constexpr double MDs_MeV_err = 0.4;

constexpr double Metas_MeV = 689.89;
constexpr double Metas_MeV_err = 0.49;

constexpr double Mphi_MeV = 1019.461;
constexpr double Mphi_MeV_err = 0.016;

constexpr double Mrho_MeV = 775.0;
constexpr double Mrho_MeV_err = 0.4;

// constexpr double grhopipi = 6.06;
// constexpr double grhopipi_err = 0.03;

constexpr double grhopipi = 5.95;
constexpr double grhopipi_err = 0.0001;

constexpr double Jpsi_MeV = 3.096916 * 1000;
constexpr double Jpsi_MeV_err = 0.001 * 1000;

constexpr double Metac_MeV = 2.9839 * 1000;
constexpr double Metac_MeV_err = 0.004 * 1000;

constexpr double MK_MeV = 494.6;
constexpr double MK_MeV_err = 0.3;



// constexpr double Mpi_MeV = 139;
// constexpr double Mpi_MeV_err = 0.001;

using namespace std;

struct  kinematic kinematic_2pt;

class replica_class {
public:
    int id;
    std::vector<int> i_conf;
};

class configuration_class {
public:
    std::vector<std::string> iconfs;
    std::vector<int> to_bin;//0 last, 1 alone, 2 first, 3 in the list
    std::vector<int> next_to_bin;
    int confs_after_binning;
    std::vector<int> order;
    std::vector<replica_class>   rep;



    void check_binnign() {

        confs_after_binning = 0;
        for (int i = 0; i < iconfs.size(); i++) {
            to_bin[i] = 1;// alone
            next_to_bin[i] = -1;
            confs_after_binning++;
            for (int j = i - 1; j >= 0; j--) {
                if (strcmp(iconfs[i].c_str(), iconfs[j].c_str()) == 0) {
                    if (to_bin[j] == 1) {
                        to_bin[j] = 2;// has a similar and it is first in the list
                    }
                    else if (to_bin[j] == 0) {
                        to_bin[j] = 3;// it needs to be binned with others
                    }
                    next_to_bin[j] = i;
                    to_bin[i] = 0;// last in the list of confs to bin
                    confs_after_binning--;
                    break;
                }

            }

        }
        // for (int i = 0; i < iconfs.size(); i++) {
        //     printf("%s   %d    %d\n", iconfs[i].c_str(), to_bin[i], next_to_bin[i]);
        // }

    }

    void setup_order() {

    }


    configuration_class(const char stream[NAMESIZE]) {
        printf("first construction\n");

        long int tmp;
        int s = file_head.l1 + 1 + 3;
        int count = 0;
        string line;
        ifstream file(stream);
        while (getline(file, line)) {
            count++;
            if (line.compare(0, 1, "#") == 0) {
                line.erase(8);
                line.erase(0, 1);
                iconfs.emplace_back(line);
                // printf("%s\n", line.c_str());
            }

        }
        std::cout << "lines=" << count << std::endl;
        error(count % s != 0, 1, "read_nconfs lines and T do not match", "lines=%i   T/2=%i", count, file_head.l1);
        int c = count / s;
        std::cout << "confs=" << c << std::endl;
        error(c != iconfs.size(), 1, "confs not found", "confs=%i   lines starting with \" # 0000_r0\"=%i", c, iconfs.size());

        to_bin = std::vector<int>(iconfs.size());
        next_to_bin = std::vector<int>(iconfs.size());
    }
    configuration_class(const char stream[NAMESIZE], int T) {
        printf("second construction\n");
        long int tmp;
        int s = T + 1 + 3;
        int count = 0;
        string line;
        ifstream file(stream);
        while (getline(file, line)) {
            count++;
            if (line.compare(0, 1, "#") == 0) {
                line.erase(8);
                line.erase(0, 1);
                iconfs.emplace_back(line);
                // printf("%s\n", line.c_str());
            }

        }
        std::cout << "lines=" << count << std::endl;
        error(count % s != 0, 1, "read_nconfs lines and T do not match", "lines=%i   T/2=%i", count, T);
        int c = count / s;
        std::cout << "confs=" << c << std::endl;

        to_bin = std::vector<int>(iconfs.size());
        next_to_bin = std::vector<int>(iconfs.size());
    }
    configuration_class() {

        // std::cout << "lines=" << 0 << std::endl;
        // std::cout << "confs=" << 0 << std::endl;

        to_bin = std::vector<int>(iconfs.size());
        next_to_bin = std::vector<int>(iconfs.size());
    }

};

void get_kinematic(int ik2, int r2, int ik1, int r1, int imom2, int imom1) {
    kinematic_2pt.ik2 = ik2;
    kinematic_2pt.ik1 = ik1;
    kinematic_2pt.k2 = file_head.k[ik2 + file_head.nk];
    kinematic_2pt.k1 = file_head.k[ik1 + file_head.nk];
    kinematic_2pt.r2 = 1;
    kinematic_2pt.r1 = 1;
    kinematic_2pt.mom2 = -file_head.mom[imom2][1];
    kinematic_2pt.mom1 = file_head.mom[imom1][1];

    kinematic_2pt.mom02 = file_head.mom[imom2][0];
    kinematic_2pt.mom01 = file_head.mom[imom1][0];

    kinematic_2pt.line = 1;

}


static void  write_file_head(FILE* stream) {
    int i, dsize;
    double* dstd;

    fwrite(&file_head.twist, sizeof(int), 1, stream);
    fwrite(&file_head.nf, sizeof(int), 1, stream);
    fwrite(&file_head.nsrc, sizeof(int), 1, stream);
    fwrite(&file_head.l0, sizeof(int), 1, stream);
    fwrite(&file_head.l1, sizeof(int), 1, stream);
    fwrite(&file_head.l2, sizeof(int), 1, stream);
    fwrite(&file_head.l3, sizeof(int), 1, stream);
    fwrite(&file_head.nk, sizeof(int), 1, stream);
    fwrite(&file_head.nmoms, sizeof(int), 1, stream);

    fwrite(&file_head.beta, sizeof(double), 1, stream);
    fwrite(&file_head.ksea, sizeof(double), 1, stream);
    fwrite(&file_head.musea, sizeof(double), 1, stream);
    fwrite(&file_head.csw, sizeof(double), 1, stream);

    fwrite(file_head.k, sizeof(double), 2 * file_head.nk, stream);

    for (i = 0;i < file_head.nmoms;i++)
        fwrite(file_head.mom[i], sizeof(double), 4, stream);
}


double matrix_element_GEVP(int t, double** cor, double mass) {
    double me;

    me = cor[t][0] / sqrt(exp(-mass * t) + exp(-(file_head.l0 - t) * mass));
    me *= 2 * mass;

    return  me;
}
void return_conf_binned(double* tmp_b, double** tmp, configuration_class confs, int id, int T) {
    double N = 1;
    for (int t = 0;t < T / 2 + 1;t++)
        tmp_b[t] += tmp[id][t];
    // printf("%s\t", confs.iconfs[id].c_str());

    while (confs.to_bin[id] >= 2) {
        id = confs.next_to_bin[id];
        // printf("%s\t", confs.iconfs[id].c_str());
        for (int t = 0;t < T / 2 + 1;t++)
            tmp_b[t] += tmp[id][t];
        N++;
    }
    for (int t = 0;t < T / 2 + 1;t++)
        tmp_b[t] /= N;

    // printf("\n");

}

void read_twopt(const char namefile[NAMESIZE], configuration_class& confs, int T, double**** to_write, int id, int Nb) {
    // char tmp[NAMESIZE];
    int sd = 0;
    std::fstream newfile;
    newfile.open(namefile, std::ios::in);
    double** tmp = double_malloc_2(confs.iconfs.size(), T / 2 + 2);

    confs.rep = std::vector<replica_class>(1);
    confs.rep[0].id = 0;
    printf("reading: %s\n", namefile);
    if (newfile.is_open()) { // checking whether the file is open
        std::string tp;
        for (int i = 0;i < confs.iconfs.size();i++) {
            getline(newfile, tp);
            // cout<< tp<< endl;
            getline(newfile, tp);
            // printf("%s \n", tp.c_str());
            if (tp.compare(0, 1, "#") == 0) {
                string tp1 = tp;
                int r = stoi(tp1.erase(0, 8));
                int ic = stoi(tp.erase(7).erase(0, 1).c_str());
                // printf("%s  %d   %d\n",tp.c_str(), r, ic);
                bool new_rep = true;
                for (int ir = 0; ir < confs.rep.size(); ir++) {
                    if (r == ir) {
                        confs.rep[ir].i_conf.emplace_back(i);
                        new_rep = false;
                        break;
                    }
                }
                if (new_rep) {
                    replica_class tmp;
                    tmp.id = r;
                    tmp.i_conf.emplace_back(i);
                    confs.rep.emplace_back(tmp);
                }
            }
            // cout<< tp<< endl;
            getline(newfile, tp);
            // cout<< tp<< endl;

            for (int t = 0;t < T / 2 + 1;t++) {
                getline(newfile, tp);
                // to_write[i][id][t][0] = stod(tp);
                tmp[i][t] = stod(tp);
                // printf("%.15f\n", to_write[i][id][t][0]);
            }

        }
    }
    else {
        error(0 == 0, 1, "correlators_analysis.cpp ",
            "unable to open %s", namefile);
    }

    printf("file:  %s We found %ld streams\n", namefile, confs.rep.size());
    for (auto r : confs.rep) {
        printf("  stream: %d   confs:  %ld\n", r.id, r.i_conf.size());
    }
    newfile.close();

    double** data = double_calloc_2(confs.iconfs.size(), T / 2 + 1);
    int count = 0;
    for (auto r : confs.rep) {
        if (r.id % 2 == 0) {
            for (int ic = r.i_conf.size() - 1; ic >= 0; ic--) {
                // printf("%d  %d %d\n", r.id, r.i_conf[ic], ic);
                for (int t = 0;t < T / 2 + 1;t++) {
                    data[count][t] = tmp[r.i_conf[ic]][t];
                }
                count++;
            }
        }
        if (r.id % 2 == 1) {
            for (int ic = 0; ic < r.i_conf.size(); ic++) {
                // printf("%d  %d\n", r.id, r.i_conf[ic]);

                for (int t = 0;t < T / 2 + 1;t++) {
                    data[count][t] = tmp[r.i_conf[ic]][t];
                }
                count++;
            }
        }
    }


    // double** data = double_calloc_2(confs.confs_after_binning, T / 2 + 1);
    // int count = 0;
    // for (int i = 0;i < confs.iconfs.size();i++) {
    //     if (confs.to_bin[i] == 1 || confs.to_bin[i] == 2) { // if it is alone or it is first of the list
    //         double* tmp_b = (double*)calloc((T / 2 + 2), sizeof(double));
    //         return_conf_binned(tmp_b, tmp, confs, i, T);
    //         for (int t = 0;t < T / 2 + 1;t++) {
    //             data[count][t] = tmp_b[t];
    //         }
    //         free(tmp_b);
    //         count++;
    //     }
    // }


    // int bin = confs.confs_after_binning / Nb;
    // int leftover = confs.confs_after_binning - bin * Nb;
    // printf("confs=%d Nb=%d  bin=%d  lefover=%d\n",confs.confs_after_binning ,Nb,bin,leftover);
    // for (int t = 0;t < T / 2 + 1;t++) {
    //     int binsize;
    //     int count=0;
    //     for (int l = 0;l < Nb;l++) {
    //         if (l<leftover) binsize=bin+1;
    //         else binsize=bin;
    //         for (int i = 0;i < binsize;i++){
    //             to_write[l][id][t][0] += data[count][t];
    //             count++;
    //         }
    //         to_write[l][id][t][0] /= ((double)binsize);
    //     }
    // }

    double clustSize = ((double)confs.iconfs.size()) / ((double)Nb);
    // double clustSize = ((double)confs.confs_after_binning) / ((double)Nb);


    for (size_t iClust = 0;iClust < Nb;iClust++) {
        /// Initial time of the bin
        const double binBegin = iClust * clustSize;
        /// Final time of the bin
        const double binEnd = binBegin + clustSize;
        double binPos = binBegin;
        do {
            /// Index of the configuration related to the time
            const size_t iConf = floor(binPos + 1e-10);

            ///Rectangle left point
            const double beg = binPos;

            /// Rectangle right point
            const double end = std::min(binEnd, iConf + 1.0);

            /// Rectangle horizontal size
            const double weight = end - beg;

            // Perform the operation passing the info
            for (int t = 0;t < T / 2 + 1;t++) {
                to_write[iClust][id][t][0] += weight * data[iConf][t];
            }
            // Updates the position
            binPos = end;
        } while (binEnd - binPos > 1e-10);
        for (int t = 0;t < T / 2 + 1;t++) {
            to_write[iClust][id][t][0] /= ((double)clustSize);
        }

    }

    free_2(confs.iconfs.size(), data);
    // free_2(confs.confs_after_binning, data);
    free_2(confs.iconfs.size(), tmp);

}


void zero_twopt(const char namefile[NAMESIZE], configuration_class& confs, int T, double**** to_write, int id, int Nb) {
    for (size_t iClust = 0;iClust < Nb;iClust++) {
        for (int t = 0;t < T / 2 + 1;t++) {
            to_write[iClust][id][t][0] = 0;
            to_write[iClust][id][t][1] = 0;
        }
    }
}


void bin_data(double**** to_write, int id, double**** data, int id1, int T, int Nconfs_bolla, int Nb) {


    double clustSize = ((double)Nconfs_bolla) / ((double)Nb);
    // double clustSize = ((double)confs.confs_after_binning) / ((double)Nb);

    for (size_t iClust = 0;iClust < Nb;iClust++) {
        for (int t = 0;t < T / 2 + 1;t++) {
            to_write[iClust][id][t][0] = 0;
        }
        /// Initial time of the bin
        const double binBegin = iClust * clustSize;
        /// Final time of the bin
        const double binEnd = binBegin + clustSize;
        double binPos = binBegin;
        do {
            /// Index of the configuration related to the time
            const size_t iConf = floor(binPos + 1e-10);

            ///Rectangle left point
            const double beg = binPos;

            /// Rectangle right point
            const double end = std::min(binEnd, iConf + 1.0);

            /// Rectangle horizontal size
            const double weight = end - beg;

            // Perform the operation passing the info
            for (int t = 0;t < T / 2 + 1;t++) {
                to_write[iClust][id][t][0] += weight * data[iConf][id1][t][0];
            }
            // printf("Cluster=%ld  iConf=%ld  weight=%g  size=%g  end=%g  beg=%g\n",iClust,iConf,weight,clustSize, end,beg);
            // Updates the position
            binPos = end;
        } while (binEnd - binPos > 1e-10);
        for (int t = 0;t < T / 2 + 1;t++) {
            to_write[iClust][id][t][0] /= ((double)clustSize);
        }

    }
}

// void read_twopt(FILE* stream, int iconf, double*** to_write, cluster::IO_params params, int index) {

//     int tmp = params.data.header_size;// 
//     tmp += sizeof(double) * iconf * params.data.size + sizeof(int) * (iconf + 1);


//     double* obs = (double*)malloc(params.data.size * sizeof(double));

//     fseek(stream, tmp, SEEK_SET);
//     size_t i = fread(obs, sizeof(double), params.data.size, stream);

//     for (int t = 0;t < params.data.L[0];t++) {
//         size_t  id = index + t * params.data.ncorr;
//         (*to_write)[t][0] = obs[id];

//     }
//     free(obs);


// }



void setup_single_file_jack(char* save_name, char** argv, const char* name, int Njack) {
    FILE* f;
    mysprintf(save_name, NAMESIZE, "/dev/null");
    f = fopen(save_name, "w+");
    error(f == NULL, 1, "setup_file_jack ",
        "Unable to open output file /dev/null");
    write_file_head(f);
    fwrite(&Njack, sizeof(int), 1, f);
    fclose(f);
}

void setup_single_file_jack_ASCI(char* save_name, char** argv, const char* name, int Njack) {
    FILE* f;
    mysprintf(save_name, NAMESIZE, "/dev/null");
    f = fopen(save_name, "w+");
    error(f == NULL, 1, "setup_file_jack ",
        "Unable to open output file /dev/null");
    fclose(f);
}

void setup_file_jack(char** argv, int Njack) {
    if (strcmp(argv[4], "jack") == 0) {
        setup_single_file_jack(file_jack.M_PS, argv, "/dev/null", Njack);
        setup_single_file_jack(file_jack.f_PS, argv, "/dev/null", Njack);
        setup_single_file_jack(file_jack.Zf_PS, argv, "/dev/null", Njack);

        setup_single_file_jack(file_jack.M_PS_GEVP, argv, "/dev/null", Njack);
        setup_single_file_jack(file_jack.f_PS_ls_ss, argv, "/dev/null", Njack);

    }

    if (strcmp(argv[4], "boot") == 0) {

        setup_single_file_jack(file_jack.M_PS, argv, "/dev/null", Njack);
        setup_single_file_jack(file_jack.f_PS, argv, "/dev/null", Njack);
        setup_single_file_jack(file_jack.Zf_PS, argv, "/dev/null", Njack);

        setup_single_file_jack(file_jack.M_PS_GEVP, argv, "/dev/null", Njack);
        setup_single_file_jack(file_jack.f_PS_ls_ss, argv, "/dev/null", Njack);
    }
}



void write_header_g2(FILE* stream, generic_header header) {

    fwrite(&header.T, sizeof(int), 1, stream);
    fwrite(&header.L, sizeof(int), 1, stream);
    int s = header.mus.size();
    fwrite(&s, sizeof(int), 1, stream);
    for (double mu : header.mus) {
        fwrite(&mu, sizeof(double), 1, stream);
    }
    s = header.thetas.size();
    fwrite(&s, sizeof(int), 1, stream);
    for (double theta : header.thetas) {
        fwrite(&theta, sizeof(double), 1, stream);
    }

    header.struct_size = ftell(stream);

}

void check_confs_correlated(std::vector<configuration_class> in_confs, std::vector<std::string>  correlators, int i, int j) {
    configuration_class& confs = in_confs[i];
    configuration_class& confs1 = in_confs[j];
    error(confs.iconfs.size() != confs1.iconfs.size(), 1, "check_confs_correlated", "file: %s  Nconf: %d\nfile: %s  Nconf: %d",
        correlators[i].c_str(), confs.iconfs.size(),
        correlators[j].c_str(), confs1.iconfs.size());
    error(confs.rep.size() != confs1.rep.size(), 1, "check_confs_correlated", "file: %s  reps: %d\nfile: %s  reps: %d",
        correlators[i].c_str(), confs.rep.size(),
        correlators[j].c_str(), confs1.rep.size());

    for (int k = 0; k < confs.iconfs.size();k++) {
        // if (strcmp(confs.iconfs[k].c_str(), confs1.iconfs[k].c_str())==0) {printf("strcmp say that the confs are the same\n");}
        // else  printf("strcmp say that the confs differs\n");

        // if (confs.iconfs[k]==confs1.iconfs[k]) {printf("== say that the confs are the same\n");}
        // else  printf("== say that the confs differs\n");
        // printf("%ld  %ld\n",confs.iconfs[k].length(),confs1.iconfs[k].length());
        error(strcmp(confs.iconfs[k].c_str(), confs1.iconfs[k].c_str()) != 0, 1, "check_confs_correlated", "file: %s  conf[%d]: %s\nfile: %s  conf[%d]: %s",
            correlators[i].c_str(), k, confs.iconfs[k].c_str(),
            correlators[j].c_str(), k, confs1.iconfs[k].c_str());

    }
}

double MDs_fit(int n, int Nvar, double* x, int Npar, double* P) {
    double mus = x[0];
    double muc = x[1];
    return P[0] + P[1] * mus + P[2] * muc;
}
double MDs_fit_one_ms(int n, int Nvar, double* x, int Npar, double* P) {
    double mus = x[0];
    double muc = x[1];
    return P[0] + P[1] * muc + P[2] * muc * muc;
}


double lhs_MDs(int n, int e, int j, data_all gjack, struct fit_type fit_info) {
    double r;
    r = gjack.en[e].jack[fit_info.corr_id[n]][j];

    return r;
}


int main(int argc, char** argv) {
    int size;
    int i, j, t;

    int* iconf, confs;
    double**** data, **** data_bin, ** out;
    char c;
    double* in;

    double**** M, **** vec, **** projected_O;
    double**** lambda, **** lambda0;

    double* fit, *** y, * x, * m, * me;


    double**** conf_jack, ** r, ** mt, ** met;
    int Ncorr = 1;
    int t0 = 2;

    FILE* f_ll = NULL, * f_sl = NULL, * f_ls = NULL, * f_ss = NULL;

    FILE* plateaux_masses = NULL, * plateaux_masses_GEVP = NULL;
    char namefile_plateaux[NAMESIZE];
    mysprintf(namefile_plateaux, NAMESIZE, "plateaux.txt");
    FILE* plateaux_f = NULL;
    char namefile[NAMESIZE];
    srand(1);

    int ncorr_new;

    error(argc < 12, 1, "main ",
        "usage:./g-2  blind/see/read_plateaux -p path basename -bin $bin"
        "   -L L jack/boot  -mu mus1  mus2 mul1 mul2  ");

    error(strcmp(argv[1], "blind") != 0 && strcmp(argv[1], "see") != 0 && strcmp(argv[1], "read_plateaux") != 0, 1, "main ",
        "argv[1] only options:  blind/see/read_plateaux ");




    error(strcmp(argv[5], "-bin") != 0, 1, "main", "argv[5] must be: -bin");

    error(strcmp(argv[7], "-L") != 0, 1, "main", "argv[7] must be: -L");
    error(strcmp(argv[9], "jack") != 0 && strcmp(argv[9], "boot") != 0, 1, "main",
        "argv[6] only options: jack/boot");
    error(strcmp(argv[10], "-mu") != 0, 1, "main", "argv[7] must be: -mu");

    char** option;
    option = (char**)malloc(sizeof(char*) * 7);
    option[0] = (char*)malloc(sizeof(char) * NAMESIZE);
    option[1] = (char*)malloc(sizeof(char) * NAMESIZE);
    option[2] = (char*)malloc(sizeof(char) * NAMESIZE);
    option[3] = (char*)malloc(sizeof(char) * NAMESIZE);
    option[4] = (char*)malloc(sizeof(char) * NAMESIZE);
    option[5] = (char*)malloc(sizeof(char) * NAMESIZE);
    option[6] = (char*)malloc(sizeof(char) * NAMESIZE);

    mysprintf(option[1], NAMESIZE, argv[1]); // blind/see/read_plateaux
    mysprintf(option[2], NAMESIZE, "-p"); // -p
    mysprintf(option[3], NAMESIZE, argv[3]); // path
    mysprintf(option[4], NAMESIZE, argv[9]); //resampling
    mysprintf(option[5], NAMESIZE, "no"); // pdf

    file_head.l1 = atoi(argv[8]);
    file_head.l0 = file_head.l1 * 2;
    file_head.l2 = file_head.l1;
    file_head.l3 = file_head.l1;

    double mus1 = atof(argv[11]);
    double mus2 = atof(argv[12]);

    double mul = atof(argv[13]);

    int Ndata = (argc - 10) / 2;
    std::vector<std::string> string_mus(Ndata);
    std::vector<std::string> string_muc(Ndata);
    std::vector<double> vec_mus(Ndata);
    std::vector<double> vec_muc(Ndata);

    for (size_t i = 0; i < Ndata; i++) {
        string_mus[i] = argv[11 + i * 2];
        string_muc[i] = argv[11 + i * 2 + 1];
        vec_mus[i] = atof(argv[11 + i * 2]);
        vec_muc[i] = atof(argv[11 + i * 2 + 1]);
        /* code */
    }



    generic_header header;
    header.L = file_head.l1;
    header.T = file_head.l0;
    header.mus = { vec_mus[0] };

    header.thetas = {};

    mysprintf(namefile, NAMESIZE, "%s_Ds", argv[4]);


    mysprintf(option[6], NAMESIZE, namefile); // basename

    printf("resampling %s\n", option[4]);
    char resampling[NAMESIZE];
    mysprintf(resampling, NAMESIZE, argv[9]);
    int T = file_head.l0;

    file_head.nk = 2;
    file_head.musea = mus1;
    file_head.k = (double*)malloc(sizeof(double) * file_head.nk * 2);
    file_head.k[0] = 0;file_head.k[1] = 0;
    file_head.k[2] = mus1;
    file_head.k[3] = mus1;

    file_head.nmoms = 1;
    file_head.mom = (double**)malloc(sizeof(double*) * file_head.nmoms);
    for (i = 0;i < file_head.nmoms;i++) {
        file_head.mom[i] = (double*)malloc(sizeof(double) * 4);
        file_head.mom[i][0] = 0;
        file_head.mom[i][1] = 0;
        file_head.mom[i][2] = 0;
        file_head.mom[i][3] = 0;
    }


    mysprintf(namefile, NAMESIZE, "%s/out/%s_output", argv[3], option[6]);
    printf("writing output in :\n %s \n", namefile);
    FILE* outfile = open_file(namefile, "w+");

    mysprintf(namefile, NAMESIZE, "%s/out/%s_gamma", argv[3], option[6]);
    printf("writing output in :\n %s \n", namefile);
    FILE* out_gamma = open_file(namefile, "w+");


    mysprintf(namefile, NAMESIZE, "%s/jackknife/%s_%s", argv[3], option[4], option[6]);
    FILE* jack_file = open_file(namefile, "w+");
    write_header_g2(jack_file, header);

    mysprintf(namefile, NAMESIZE, "%s/jackknife/jack_%s_Z.txt", argv[3], argv[4]);
    FILE* ASCII_Z = open_file(namefile, "w+");


    //////////////////////////////////////////////////////////////
    // declare corr to read
    //////////////////////////////////////////////////////////////
    for (auto e : string_mus) {
        std::cout << e << "\n";
    }
    std::vector<std::string>  correlators;

    for (size_t i = 0; i < Ndata; i++) {
        mysprintf(namefile, NAMESIZE, "%s/%s_r.opposite_mu.%s_mu.%s_P5P5.txt", argv[3], argv[4], string_mus[i].c_str(), string_muc[i].c_str());//0
        correlators.emplace_back(namefile);
    }




    std::vector<configuration_class> myconfs;

    // confs = myconfs[0].iconfs.size();
    int count = 0;
    for (auto name : correlators) {
        printf("reading  confs from file: %s\n", name.c_str());

        myconfs.emplace_back(name.c_str());

        myconfs[count].confs_after_binning = myconfs[count].iconfs.size();
        // }
        cout << "number of different configurations:" << myconfs[count].confs_after_binning << endl;
        count++;
    }

    //////////////////////////////////////////////////////////////
    // setup jack/boot
    //////////////////////////////////////////////////////////////

    int bin = atoi(argv[6]);
    int Neff = bin;

    // cout << "effective configurations after binning ( bin size " << confs / bin << "):  " << Neff << endl;

    int Njack;
    if (strcmp(argv[9], "jack") == 0) {
        Njack = Neff + 1;
    }
    else if (strcmp(argv[9], "boot") == 0) {
        Njack = Nbootstrap + 1;
    }
    else {
        Njack = 0;
        error(1 == 1, 1, "main", "argv[9]= %s is not jack or boot", argv[9]);
    }
    fwrite(&Njack, sizeof(int), 1, jack_file);
    header.Njack = Njack;

    setup_file_jack(option, Njack);
    get_kinematic(0, 0, 1, 0, 0, 0);

    int var = correlators.size();
    int Max_corr = var + 10;
    data = calloc_corr(bin, Max_corr /*Dmeson*/, file_head.l0);
    //////////////////////////////////////////////////////////////
    // read
    //////////////////////////////////////////////////////////////

    for (int i = 0; i < var; i++) {
        read_twopt(correlators[i].c_str(), myconfs[i], T, data, i, bin);
    }


    ncorr_new = correlators.size();
    conf_jack = create_resampling(option[4], Neff, Max_corr, file_head.l0, data);
    free_corr(Neff, Max_corr, file_head.l0, data);


    // ////////////////// symmetrization/////////////////////////////////////////////
    for (int i = 0;i < ncorr_new;i++) { symmetrise_jackboot(Njack, i, file_head.l0, conf_jack); }
    corr_counter = -1;
    //////////////////////////////////////////////////////////////
    // read lattice spacing
    //////////////////////////////////////////////////////////////
    double* a = (double*)malloc(sizeof(double) * Njack);// allocate memory 
    double* ml = (double*)malloc(sizeof(double) * Njack);// allocate memory 
    double* phys_ms = (double*)malloc(sizeof(double) * Njack);// allocate memory 
    double* phys_mc = (double*)malloc(sizeof(double) * Njack);// allocate memory 
    std::string latt;
    set_a_ml_ms_mc(argv[4], a, ml, phys_ms, phys_mc, latt);
    printf("reading a   =  %g  %g fm\n", a[Njack - 1], myres->comp_error(a));
    corr_counter = -1;
    write_jack(a, Njack, jack_file);

    check_correlatro_counter(0);
    printf("reading a   =  %g  %g fm\n", a[Njack - 1], myres->comp_error(a));
    printf("reading amul^phys=  %g  %g\n", ml[Njack - 1], myres->comp_error(ml));
    printf("reading amus^phys=  %g  %g\n", phys_ms[Njack - 1], myres->comp_error(phys_ms));
    ////////////////////////////////////////////////
    printf("################### fitting  the correlaators #################\n");
    double* zeros = (double*)calloc(Njack, sizeof(double));
    std::vector<double*> M_Ds(Ndata);
    char namefit[NAMESIZE];
    for (size_t i = 0; i < Ndata; i++) {
        mysprintf(namefit, NAMESIZE, "M_{Ds%d}", i);
        M_Ds[i] = plateau_correlator_function(option, kinematic_2pt, (char*)"P5P5", conf_jack, Njack, namefile_plateaux, outfile, i, namefit, M_eff_T, jack_file);
        check_correlatro_counter(1 + i);
    }


    //////////////////////////////////////////////////////////////
    // global fit
    //////////////////////////////////////////////////////////////
    printf("################### interpolating #################\n");

    data_all jackall;
    jackall.resampling = resampling;
    jackall.ens = Ndata;
    jackall.en = new data_single[jackall.ens];
    for (int i = 0;i < jackall.ens;i++) {
        jackall.en[i].header = header;
        jackall.en[i].Nobs = 1;
        jackall.en[i].Njack = header.Njack;
        jackall.en[i].jack = (double**)malloc(sizeof(double*) * jackall.en[i].Nobs);
        jackall.en[i].jack[0] = M_Ds[i];
    }

    fit_type fit_info;
    // fit in sigma
    fit_info.restore_default();
    fit_info.Nvar = 2;
    fit_info.Njack = header.Njack;
    fit_info.N = 1;
    int offset = 6;
    fit_info.myen = std::vector<int>(Ndata);
    for (size_t i = 0; i < Ndata; i++)fit_info.myen[i] = i;
    fit_info.entot = fit_info.myen.size() * fit_info.N;
    fit_info.malloc_x();
    count = 0;
    for (int j = 0;j < fit_info.Njack;j++) {
        for (size_t i = 0; i < Ndata; i++) {
            fit_info.x[0][i][j] = vec_mus[i]; // mul1 
            fit_info.x[1][i][j] = vec_muc[i]; // mul2
        }

    }
    count++;

    fit_info.corr_id = { 0 };
    fit_info.function = MDs_fit_one_ms;//constant_fit;
    fit_info.Npar = 3;
    bool all_mus_equal = true;
    for (int i = 0;i < string_mus.size();i++) {
        if (string_mus[i].compare(string_mus[0]) != 0) {
            fit_info.function = MDs_fit;//constant_fit;
            fit_info.Npar = 3;
            all_mus_equal = false;
        }
    }

    fit_info.linear_fit = true;
    fit_info.covariancey = false;
    fit_info.verbosity = 0;
    mysprintf(namefit, NAMESIZE, "%s/out/%s_MDs", argv[3], option[6]);
    char** temp_argv = malloc_2<char>(5, NAMESIZE);
    mysprintf(temp_argv[1], NAMESIZE, "%s", resampling);// resampling
    mysprintf(temp_argv[3], NAMESIZE, "%s/out", option[3]);// resampling

    fit_result fit_Z0_sigma = fit_all_data(temp_argv, jackall, lhs_MDs, fit_info, namefit);
    printf("chi2/dof = %g\n", fit_Z0_sigma.chi2[Njack - 1]);
    printf("Ndata = %d\n", fit_Z0_sigma.dof + fit_Z0_sigma.Npar);
    printf("Npar = %d\n", fit_Z0_sigma.Npar);
    printf("dof = %d\n", fit_Z0_sigma.dof);
    auto min = std::ranges::min_element(vec_muc.begin(), vec_muc.end());
    auto max = std::ranges::max_element(vec_muc.begin(), vec_muc.end());
    fit_info.band_range = std::vector<double>(2);
    fit_info.band_range[0] = *min - 0.02;
    fit_info.band_range[1] = *max + 0.02;
    std::vector<double> xcont = { vec_mus[0], vec_muc[0] };
    char nameband[NAMESIZE];
    for (int i = 0; i < vec_mus.size();i++) {
        xcont[0] = vec_mus[i];
        mysprintf(nameband, NAMESIZE, "mus%d_muc", i);
        print_fit_band(temp_argv, jackall, fit_info, fit_info, namefit, nameband, fit_Z0_sigma, fit_Z0_sigma, 1, 0/*en*/, 0.001, xcont);

    }
    // print_fit_band(temp_argv, jackall, fit_info, fit_info, namefit, "mus_l1", fit_Z0_sigma, fit_Z0_sigma, 1, 2/*en*/, 0.01);
    // print_fit_band(temp_argv, jackall, fit_info, fit_info, namefit, "mus_l2", fit_Z0_sigma, fit_Z0_sigma, 1, 3/*en*/, 0.01);

    // write_jack(fit_Z0_sigma.P[0], Njack, jack_file);

    double* jack_MDs_Mev = myres->create_fake(MDs_MeV, MDs_MeV_err, 1);
    double* jack_aMDs_exp = myres->create_copy(jack_MDs_Mev);
    myres->mult(jack_aMDs_exp, jack_aMDs_exp, a);
    myres->div(jack_aMDs_exp, jack_aMDs_exp, hbarc);

    double** tif = swap_indices(fit_info.Npar, Njack, fit_Z0_sigma.P);
    std::vector<double> mc_MDs(Njack);
    std::vector<double> swapped_x(fit_info.Nvar);

    // if (all_mus_equal) {
    for (int j = 0;j < Njack;j++) {

        swapped_x[0] = phys_ms[j];
        swapped_x[1] = 9999999; // it does not mater it is not used
        mc_MDs[j] = rtbis_func_eq_input(fit_info.function, 0 /*n*/, fit_info.Nvar, swapped_x.data(), fit_info.Npar, tif[j], 1/* ivar*/, jack_aMDs_exp[j], 0.1, 0.3, 1e-10, 2);
        // (jack_aMDs_exp[j] - fit_Z0_sigma.P[0][j] * phys_ms[j]) / fit_Z0_sigma.P[1][j];
    }
    // }
    // else {
    //     for (int j = 0;j < Njack;j++) {
    //         mc_MDs[j] = (jack_aMDs_exp[j] - fit_Z0_sigma.P[0][j] * phys_ms[j]) / fit_Z0_sigma.P[1][j];
    //     }
    // }

    // zero_corr(zeros, Njack, jack_file);
    // zero_corr(zeros, Njack, jack_file);


    write_jack(mc_MDs.data(), Njack, jack_file);
    // check_correlatro_counter(7);
    printf("interpolating to MDs= %g  %g  MeV\n", MDs_MeV, MDs_MeV_err);
    printf("interpolating to aMDs= %g  \n", MDs_MeV * a[Njack - 1] / hbarc);
    printf("################### result #################\n");

    const char* description = "mc(MDs)";
    fprintf(outfile, " \n\n# aMetas_exp  Zint  err\n");

    fprintf(outfile, "%.15g   %.15g   %.15g\t", mc_MDs[Njack - 1], mc_MDs[Njack - 1], error_jackboot(resampling, Njack, mc_MDs.data()));

    fprintf(outfile, "\n\n #%s fit in [%d,%d] chi2=%.5g  %.5g\n", description, 0, 0, 0.0, 0.0);
    fprintf(outfile, "   %.15g   %.15g\n", mc_MDs[Njack - 1], error_jackboot(resampling, Njack, mc_MDs.data()));
    printf("%s (%.15g) =  %.15g   %.15g\n", description, jack_aMDs_exp[Njack - 1], mc_MDs[Njack - 1], error_jackboot(resampling, Njack, mc_MDs.data()));


    mysprintf(namefile, NAMESIZE, "%s/out/mc_from_MDs_%s.txt", argv[3], latt.c_str());
    myres->write_jack_in_file(mc_MDs.data(), namefile);


    //////////////////////



}

