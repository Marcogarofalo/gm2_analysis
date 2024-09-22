#define CONTROL

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <complex.h>

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
        "   -L L jack/boot  -mu mul   ");


    error(strcmp(argv[1], "blind") != 0 && strcmp(argv[1], "see") != 0 && strcmp(argv[1], "read_plateaux") != 0, 1, "main ",
        "argv[1] only options:  blind/see/read_plateaux ");

    // cluster::IO_params params;
    // mysprintf(namefile, NAMESIZE, "%s/%s", argv[3], argv[4]);
    // FILE* infile = open_file(namefile, "r+");
    // read_header_phi4(infile, params);




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

    double mu = atof(argv[11]);
    // double mus1 = atof(argv[12]);
    // double mus2 = atof(argv[13]);

    // double muc1 = atof(argv[14]);
    // double muc2 = atof(argv[15]);
    // double muc3 = atof(argv[16]);
    double mul1 = mu;
    // if (argc > 17 && strcmp(argv[17], "three_corr") != 0) { mul1 = atof(argv[17]); }

    generic_header header;
    header.L = file_head.l1;
    header.T = file_head.l0;
    header.mus = { mu };
    // if (argc > 17 && strcmp(argv[17], "three_corr") != 0)  header.mus.emplace_back(mul1);
    header.thetas = {};

    mysprintf(namefile, NAMESIZE, "%s_mu.%f", argv[4], mu);


    mysprintf(option[6], NAMESIZE, namefile); // basename

    printf("resampling %s\n", option[4]);
    char resampling[NAMESIZE];
    mysprintf(resampling, NAMESIZE, argv[9]);
    int T = file_head.l0;

    file_head.nk = 2;
    file_head.musea = mu;
    file_head.k = (double*)malloc(sizeof(double) * file_head.nk * 2);
    file_head.k[0] = 0;file_head.k[1] = 0;
    file_head.k[2] = mu;
    file_head.k[3] = mu;

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


    std::vector<std::string>  correlators;
    // if (strcmp(argv[4], "cA.53.24") == 0 || strcmp(argv[4], "cA.40.24") == 0 || strcmp(argv[4], "cA.30.32") == 0) {
    //     printf("OS and tm inverted for light\n");
    //     // light
    //     mysprintf(namefile, NAMESIZE, "%s/%s_r.opposite_mu.%.5f_P5A0.txt", argv[3], argv[4], mu);//3
    //     correlators.emplace_back(namefile);
    //     mysprintf(namefile, NAMESIZE, "%s/%s_r.opposite_mu.%.5f_P5P5.txt", argv[3], argv[4], mu);//4
    //     correlators.emplace_back(namefile);
    //     mysprintf(namefile, NAMESIZE, "%s/%s_r.opposite_mu.%.5f_VKVK.txt", argv[3], argv[4], mu);//5
    //     correlators.emplace_back(namefile);
    //     mysprintf(namefile, NAMESIZE, "%s/%s_r.equal_mu.%.5f_P5A0.txt", argv[3], argv[4], mu);//0
    //     correlators.emplace_back(namefile);
    //     mysprintf(namefile, NAMESIZE, "%s/%s_r.equal_mu.%.5f_P5P5.txt", argv[3], argv[4], mu);//1
    //     correlators.emplace_back(namefile);
    //     mysprintf(namefile, NAMESIZE, "%s/%s_r.equal_mu.%.5f_VKVK.txt", argv[3], argv[4], mu);//2
    //     correlators.emplace_back(namefile);
    // }
    // else {
        // light
    // mysprintf(namefile, NAMESIZE, "%s/%s_r.equal_mu.%.5f_P5A0.txt", argv[3], argv[4], mu);//0
    // correlators.emplace_back(namefile);
    // mysprintf(namefile, NAMESIZE, "%s/%s_r.equal_mu.%.5f_P5P5.txt", argv[3], argv[4], mu);//1
    // correlators.emplace_back(namefile);
    // mysprintf(namefile, NAMESIZE, "%s/%s_r.equal_mu.%.5f_VKVK.txt", argv[3], argv[4], mu);//2
    // correlators.emplace_back(namefile);
    mysprintf(namefile, NAMESIZE, "%s/%s_r.opposite_mu.%.5f_P5A0.txt", argv[3], argv[4], mu);//0
    correlators.emplace_back(namefile);
    mysprintf(namefile, NAMESIZE, "%s/%s_r.opposite_mu.%.5f_P5P5.txt", argv[3], argv[4], mu);//1
    correlators.emplace_back(namefile);
    mysprintf(namefile, NAMESIZE, "%s/%s_r.opposite_mu.%.5f_VKVK.txt", argv[3], argv[4], mu);//2
    correlators.emplace_back(namefile);
    // }

    // printf("reading confs from file: %s", correlators[0].c_str());
    // auto iconfs = read_nconfs(correlators[0].c_str());
    std::vector<configuration_class> myconfs;

    // confs = myconfs[0].iconfs.size();
    int count = 0;
    for (auto name : correlators) {
        printf("reading  confs from file: %s\n", name.c_str());
        // if (argc <= 17 || strcmp(argv[17], "three_corr") == 0) {
        //     if (count >= 36) {
        //         configuration_class tmp;
        //         myconfs.emplace_back(tmp);
        //     }
        //     else
        //         myconfs.emplace_back(name.c_str());
        // }
        // else if (count == 42) {
        //     if (strcmp(argv[4], "cB.72.96") == 0) {
        //         configuration_class tmp;
        //         myconfs.emplace_back(tmp);
        //     }
        //     else {
        //         configuration_class tmp(name.c_str(), 0);
        //         myconfs.emplace_back(tmp);
        //     }
        // }
        // else if (count == 55 || count == 56) {
        //     configuration_class tmp(name.c_str(), header.T - 1 - 3 + 1);
        //     myconfs.emplace_back(tmp);
        // }
        // else
        myconfs.emplace_back(name.c_str());

        // myconfs[cout]=read_nconfs(name.c_str());
        if (strcmp(argv[4], "cB.72.96") != 0)
            myconfs[count].check_binnign();
        else {
            myconfs[count].confs_after_binning = myconfs[count].iconfs.size();
        }
        cout << "number of different configurations:" << myconfs[count].confs_after_binning << endl;
        count++;
        // printf("checking confs from file: %s\n", name.c_str());
        // // auto check_iconfs = read_nconfs(name.c_str());
        // configuration_class check_confs(name.c_str());
        // for (int i = 0;i < confs;i++) {
        //     error(myconfs.iconfs[i] != check_confs.iconfs[i], 1, "configurations id do not match", "");
        // }
        // error(confs != check_confs.iconfs.size(), 1, "reading number of configuration", "not all the files have the same confs");
    }

    // which_confs_are_the_same(std::vector<std::string> iconfs);

    // FILE* infile_equal_P5A0 = open_file(namefile, "r+");
    // int count = 0;


    //confs=confs/10;
    // compute what will be the neff after the binning 
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
    data = calloc_corr(bin, var + 3 + 12 + 2 /*Dmeson*/, file_head.l0);

    for (int i = 0; i < var; i++) {
        read_twopt(correlators[i].c_str(), myconfs[i], T, data, i, bin);
        


    }


    /////////////////////////////////////////////////////////////////////////////////////////////
    // free theory corrections
    /////////////////////////////////////////////////////////////////////////////////////////////
    // var+3 l op
    // var+4 l eq

    // var+5 s1 op
    // var+6 s1 eq
    // var+7 s2 op
    // var+8 s2 eq

    // var+9  c1 op
    // var+10 c1 eq
    // var+11 c2 op
    // var+12 c2 eq
    // var+13 c3 op
    // var+14 c3 eq



    // data_bin = binning(confs, var, file_head.l0, data, bin);
    // data_bin = binning_toNb(confs, var, file_head.l0, data, bin);
    // //if you want to do the gamma analysis you need to do before freeing the raw data

    // gamma_correlator_func(option, kinematic_2pt, (char*)"P5P5", data, Neff, namefile_plateaux, out_gamma, 5, "VK_{pi}^{TM}", identity_func_gamma);
    // // gamma_correlator_func(option, kinematic_2pt, (char*)"P5P5", data, Neff, namefile_plateaux, out_gamma, 2, "VK_{pi}^{OS}", identity_func_gamma);

    // gamma_correlator_func(option, kinematic_2pt, (char*)"P5P5", data, Neff, namefile_plateaux, out_gamma, 11, "VK_s1^{TM}", identity_func_gamma);
    // // gamma_correlator_func(option, kinematic_2pt, (char*)"P5P5", data, Neff, namefile_plateaux, out_gamma, 8, "VK_s1^{OS}", identity_func_gamma);


    // gamma_correlator_func(option, kinematic_2pt, (char*)"P5P5", data, Neff, namefile_plateaux, out_gamma, 4, "P5P5_{pi}^{TM}", identity_func_gamma);
    // // gamma_correlator_func(option, kinematic_2pt, (char*)"P5P5", data, Neff, namefile_plateaux, out_gamma, 1, "P5P5_{pi}^{OS}", identity_func_gamma);

    // gamma_correlator_func(option, kinematic_2pt, (char*)"P5P5", data, Neff, namefile_plateaux, out_gamma, 4, "M_{pi}^{TM}", mass_gamma);
    // // gamma_correlator_func(option, kinematic_2pt, (char*)"P5P5", data, Neff, namefile_plateaux, out_gamma, 1, "M_{pi}^{OS}", mass_gamma);


    // //effective_mass_phi4_gamma(  option, kinematic_2pt,   (char*) "P5P5", data,  confs ,namefile_plateaux,out_gamma,3,"M_{PS}^{ll}");
    // free_corr(bin, var, file_head.l0, data);
    ncorr_new = correlators.size();
    conf_jack = create_resampling(option[4], Neff, correlators.size(), file_head.l0, data);
    free_corr(Neff, correlators.size(), file_head.l0, data);


    // ////////////////// symmetrization/////////////////////////////////////////////
    // for (int i = 0;i <= 7;i++) { symmetrise_jackboot(Njack, i, file_head.l0, conf_jack); }

    ////////////////////////////////////////////////
    double* zeros = (double*)calloc(Njack, sizeof(double));
    corr_counter = -1;
    // double* M_PS = plateau_correlator_function(option, kinematic_2pt, (char*)"P5P5", conf_jack, Njack, namefile_plateaux, outfile, 1, "M_{PS}^{eq}", M_eff_T, jack_file);
    zero_corr(zeros, Njack, jack_file);
    check_correlatro_counter(0);

    double* M_PS_op = plateau_correlator_function(option, kinematic_2pt, (char*)"P5P5", conf_jack, Njack, namefile_plateaux, outfile, 1, "M_{PS}^{op}", M_eff_T, jack_file);
    check_correlatro_counter(1);

    for (int i = 2; i < 163;i++) {
        zero_corr(zeros, Njack, jack_file);
    }

    //////////////////////////////////////////////////////////////
    // f_pi 
    //////////////////////////////////////////////////////////////
    fit_type fit_info;

    fit_info.restore_default();
    fit_info.Nvar = 1;
    fit_info.Npar = 1;
    fit_info.N = 1;
    fit_info.Njack = Njack;
    fit_info.n_ext_P = 3;
    fit_info.malloc_ext_P();
    for (int j = 0; j < fit_info.Njack; j++) {
        fit_info.ext_P[0][j] = M_PS_op[j];
        fit_info.ext_P[1][j] = header.mus[0];
        fit_info.ext_P[2][j] = header.mus[0];
    }
    fit_info.function = constant_fit;
    fit_info.linear_fit = true;
    fit_info.T = header.T;

    fit_info.corr_id = { 1 };

    struct fit_result f_PS = fit_fun_to_fun_of_corr(
        option, kinematic_2pt, (char*)"P5P5", conf_jack, namefile_plateaux,
        outfile, lhs_function_f_PS, "f_{pi}", fit_info, jack_file);
    free_fit_result(fit_info, f_PS);
    check_correlatro_counter(163);
    // if (argc >= 19) {
    //     for (int j = 0; j < fit_info.Njack; j++) {
    //         fit_info.ext_P[0][j] = M_PS1[j];
    //         fit_info.ext_P[1][j] = mul1;
    //         fit_info.ext_P[2][j] = mul1;
    //     }
    //     fit_info.corr_id = { id_P5P5_cor_mudmu };

    //     struct fit_result f_PS1 = fit_fun_to_fun_of_corr(
    //         option, kinematic_2pt, (char*)"P5P5", conf_jack, namefile_plateaux,
    //         outfile, lhs_function_f_PS, "f_{pi}_l1", fit_info, jack_file);
    //     free_fit_result(fit_info, f_PS1);
    //     check_correlatro_counter(164);
    // }
    // else {
    //     zero_corr(zeros, Njack, jack_file);       check_correlatro_counter(164);
    // }
    zero_corr(zeros, Njack, jack_file);       check_correlatro_counter(164);
    double* mu_jack = myres->create_fake(mu, mu / 1e+6, 1);
    double* mu1_jack = myres->create_fake(mul1, mul1 / 1e+6, 1);
    write_jack(mu_jack, Njack, jack_file);        check_correlatro_counter(165);
    write_jack(mu1_jack, Njack, jack_file);        check_correlatro_counter(166);



    //////////////////////



}

