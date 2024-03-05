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

constexpr double MK_MeV = 494.2;
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

    error(argc < 17 || argc > 21, 1, "main ",
        "usage:./g-2  blind/see/read_plateaux -p path basename -bin $bin"
        "   -L L jack/boot  -mu mul    mus1 mus2     muc1 muc2 muc3   [mul'] [bolla]  [free_corr]   [three_corr]");


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
    double mus1 = atof(argv[12]);
    double mus2 = atof(argv[13]);

    double muc1 = atof(argv[14]);
    double muc2 = atof(argv[15]);
    double muc3 = atof(argv[16]);
    double mul1;
    if (argc > 17 && strcmp(argv[17], "three_corr") != 0) { mul1 = atof(argv[17]); }

    generic_header header;
    header.L = file_head.l1;
    header.T = file_head.l0;
    header.mus = { mu, mus1, mus2 ,muc1, muc2, muc3 };
    if (argc > 17 && strcmp(argv[17], "three_corr") != 0)  header.mus.emplace_back(mul1);
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
    mysprintf(namefile, NAMESIZE, "%s/%s_r.equal_mu.%.5f_P5A0.txt", argv[3], argv[4], mu);//0
    correlators.emplace_back(namefile);
    mysprintf(namefile, NAMESIZE, "%s/%s_r.equal_mu.%.5f_P5P5.txt", argv[3], argv[4], mu);//1
    correlators.emplace_back(namefile);
    mysprintf(namefile, NAMESIZE, "%s/%s_r.equal_mu.%.5f_VKVK.txt", argv[3], argv[4], mu);//2
    correlators.emplace_back(namefile);
    mysprintf(namefile, NAMESIZE, "%s/%s_r.opposite_mu.%.5f_P5A0.txt", argv[3], argv[4], mu);//3
    correlators.emplace_back(namefile);
    mysprintf(namefile, NAMESIZE, "%s/%s_r.opposite_mu.%.5f_P5P5.txt", argv[3], argv[4], mu);//4
    correlators.emplace_back(namefile);
    mysprintf(namefile, NAMESIZE, "%s/%s_r.opposite_mu.%.5f_VKVK.txt", argv[3], argv[4], mu);//5
    correlators.emplace_back(namefile);
    // }
    // stranges
    mysprintf(namefile, NAMESIZE, "%s/%s_r.equal_mu.%.3f_P5A0.txt", argv[3], argv[4], mus1);//6
    correlators.emplace_back(namefile);
    mysprintf(namefile, NAMESIZE, "%s/%s_r.equal_mu.%.3f_P5P5.txt", argv[3], argv[4], mus1);//7
    correlators.emplace_back(namefile);
    mysprintf(namefile, NAMESIZE, "%s/%s_r.equal_mu.%.3f_VKVK.txt", argv[3], argv[4], mus1);//8
    correlators.emplace_back(namefile);
    mysprintf(namefile, NAMESIZE, "%s/%s_r.opposite_mu.%.3f_P5A0.txt", argv[3], argv[4], mus1);//9
    correlators.emplace_back(namefile);
    mysprintf(namefile, NAMESIZE, "%s/%s_r.opposite_mu.%.3f_P5P5.txt", argv[3], argv[4], mus1);//10
    correlators.emplace_back(namefile);
    mysprintf(namefile, NAMESIZE, "%s/%s_r.opposite_mu.%.3f_VKVK.txt", argv[3], argv[4], mus1);//11
    correlators.emplace_back(namefile);

    mysprintf(namefile, NAMESIZE, "%s/%s_r.equal_mu.%.3f_P5A0.txt", argv[3], argv[4], mus2);//12
    correlators.emplace_back(namefile);
    mysprintf(namefile, NAMESIZE, "%s/%s_r.equal_mu.%.3f_P5P5.txt", argv[3], argv[4], mus2);//13
    correlators.emplace_back(namefile);
    mysprintf(namefile, NAMESIZE, "%s/%s_r.equal_mu.%.3f_VKVK.txt", argv[3], argv[4], mus2);//14
    correlators.emplace_back(namefile);
    mysprintf(namefile, NAMESIZE, "%s/%s_r.opposite_mu.%.3f_P5A0.txt", argv[3], argv[4], mus2);//15
    correlators.emplace_back(namefile);
    mysprintf(namefile, NAMESIZE, "%s/%s_r.opposite_mu.%.3f_P5P5.txt", argv[3], argv[4], mus2);//16
    correlators.emplace_back(namefile);
    mysprintf(namefile, NAMESIZE, "%s/%s_r.opposite_mu.%.3f_VKVK.txt", argv[3], argv[4], mus2);//17
    correlators.emplace_back(namefile);

    mysprintf(namefile, NAMESIZE, "%s/%s_r.equal_mu.%.5f_P5A0.txt", argv[3], argv[4], muc1);//21
    correlators.emplace_back(namefile);
    mysprintf(namefile, NAMESIZE, "%s/%s_r.equal_mu.%.5f_P5P5.txt", argv[3], argv[4], muc1);//22
    correlators.emplace_back(namefile);
    mysprintf(namefile, NAMESIZE, "%s/%s_r.equal_mu.%.5f_VKVK.txt", argv[3], argv[4], muc1);//23
    correlators.emplace_back(namefile);
    mysprintf(namefile, NAMESIZE, "%s/%s_r.opposite_mu.%.5f_P5A0.txt", argv[3], argv[4], muc1);//18
    correlators.emplace_back(namefile);
    mysprintf(namefile, NAMESIZE, "%s/%s_r.opposite_mu.%.5f_P5P5.txt", argv[3], argv[4], muc1);//19
    correlators.emplace_back(namefile);
    mysprintf(namefile, NAMESIZE, "%s/%s_r.opposite_mu.%.5f_VKVK.txt", argv[3], argv[4], muc1);//20
    correlators.emplace_back(namefile);

    mysprintf(namefile, NAMESIZE, "%s/%s_r.equal_mu.%.5f_P5A0.txt", argv[3], argv[4], muc2);//27
    correlators.emplace_back(namefile);
    mysprintf(namefile, NAMESIZE, "%s/%s_r.equal_mu.%.5f_P5P5.txt", argv[3], argv[4], muc2);//28
    correlators.emplace_back(namefile);
    mysprintf(namefile, NAMESIZE, "%s/%s_r.equal_mu.%.5f_VKVK.txt", argv[3], argv[4], muc2);//29
    correlators.emplace_back(namefile);
    mysprintf(namefile, NAMESIZE, "%s/%s_r.opposite_mu.%.5f_P5A0.txt", argv[3], argv[4], muc2);//24
    correlators.emplace_back(namefile);
    mysprintf(namefile, NAMESIZE, "%s/%s_r.opposite_mu.%.5f_P5P5.txt", argv[3], argv[4], muc2);//25
    correlators.emplace_back(namefile);
    mysprintf(namefile, NAMESIZE, "%s/%s_r.opposite_mu.%.5f_VKVK.txt", argv[3], argv[4], muc2);//26
    correlators.emplace_back(namefile);


    mysprintf(namefile, NAMESIZE, "%s/%s_r.equal_mu.%.5f_P5A0.txt", argv[3], argv[4], muc3);//33
    correlators.emplace_back(namefile);
    mysprintf(namefile, NAMESIZE, "%s/%s_r.equal_mu.%.5f_P5P5.txt", argv[3], argv[4], muc3);//34
    correlators.emplace_back(namefile);
    mysprintf(namefile, NAMESIZE, "%s/%s_r.equal_mu.%.5f_VKVK.txt", argv[3], argv[4], muc3);//35
    correlators.emplace_back(namefile);
    mysprintf(namefile, NAMESIZE, "%s/%s_r.opposite_mu.%.5f_P5A0.txt", argv[3], argv[4], muc3);//30
    correlators.emplace_back(namefile);
    mysprintf(namefile, NAMESIZE, "%s/%s_r.opposite_mu.%.5f_P5P5.txt", argv[3], argv[4], muc3);//31
    correlators.emplace_back(namefile);
    mysprintf(namefile, NAMESIZE, "%s/%s_r.opposite_mu.%.5f_VKVK.txt", argv[3], argv[4], muc3);//32
    correlators.emplace_back(namefile);

    if (argc > 17 && strcmp(argv[17], "three_corr") != 0) { // cott mu+dm
        mysprintf(namefile, NAMESIZE, "%s/%s_r.equal_mu.%.7f_P5A0.txt", argv[3], argv[4], mul1);//36
        correlators.emplace_back(namefile);
        mysprintf(namefile, NAMESIZE, "%s/%s_r.equal_mu.%.7f_P5P5.txt", argv[3], argv[4], mul1);//37
        correlators.emplace_back(namefile);
        mysprintf(namefile, NAMESIZE, "%s/%s_r.equal_mu.%.7f_VKVK.txt", argv[3], argv[4], mul1);//38
        correlators.emplace_back(namefile);
        mysprintf(namefile, NAMESIZE, "%s/%s_r.opposite_mu.%.7f_P5A0.txt", argv[3], argv[4], mul1);//39
        correlators.emplace_back(namefile);
        mysprintf(namefile, NAMESIZE, "%s/%s_r.opposite_mu.%.7f_P5P5.txt", argv[3], argv[4], mul1);//40
        correlators.emplace_back(namefile);
        mysprintf(namefile, NAMESIZE, "%s/%s_r.opposite_mu.%.7f_VKVK.txt", argv[3], argv[4], mul1);//41
        correlators.emplace_back(namefile);
    }
    if (argc > 18 && strcmp(argv[17], "three_corr") != 0) {
        mysprintf(namefile, NAMESIZE, "%s/%s_mu.%.5f_bolla_std.txt", argv[3], argv[4], mu);//42
        correlators.emplace_back(namefile);
    }
    if (argc > 17 && strcmp(argv[17], "three_corr") != 0) { /// corr_mu+dm
        mysprintf(namefile, NAMESIZE, "%s/%s_small_stat_r.equal_mu.%.5f_P5A0.txt", argv[3], argv[4], mu);//43
        correlators.emplace_back(namefile);
        mysprintf(namefile, NAMESIZE, "%s/%s_small_stat_r.equal_mu.%.5f_P5P5.txt", argv[3], argv[4], mu);//44
        correlators.emplace_back(namefile);
        mysprintf(namefile, NAMESIZE, "%s/%s_small_stat_r.equal_mu.%.5f_VKVK.txt", argv[3], argv[4], mu);//45
        correlators.emplace_back(namefile);
        mysprintf(namefile, NAMESIZE, "%s/%s_small_stat_r.opposite_mu.%.5f_P5A0.txt", argv[3], argv[4], mu);//46
        correlators.emplace_back(namefile);
        mysprintf(namefile, NAMESIZE, "%s/%s_small_stat_r.opposite_mu.%.5f_P5P5.txt", argv[3], argv[4], mu);//47
        correlators.emplace_back(namefile);
        mysprintf(namefile, NAMESIZE, "%s/%s_small_stat_r.opposite_mu.%.5f_VKVK.txt", argv[3], argv[4], mu);//48
        correlators.emplace_back(namefile);
    }
    if (argc > 18 && strcmp(argv[17], "three_corr") != 0) {
        mysprintf(namefile, NAMESIZE, "%s/%s_corr_bolla_r.equal_mu.%.5f_P5A0.txt", argv[3], argv[4], mu);//49
        correlators.emplace_back(namefile);
        mysprintf(namefile, NAMESIZE, "%s/%s_corr_bolla_r.equal_mu.%.5f_P5P5.txt", argv[3], argv[4], mu);//50
        correlators.emplace_back(namefile);
        mysprintf(namefile, NAMESIZE, "%s/%s_corr_bolla_r.equal_mu.%.5f_VKVK.txt", argv[3], argv[4], mu);//51
        correlators.emplace_back(namefile);
        mysprintf(namefile, NAMESIZE, "%s/%s_corr_bolla_r.opposite_mu.%.5f_P5A0.txt", argv[3], argv[4], mu);//52
        correlators.emplace_back(namefile);
        mysprintf(namefile, NAMESIZE, "%s/%s_corr_bolla_r.opposite_mu.%.5f_P5P5.txt", argv[3], argv[4], mu);//53
        correlators.emplace_back(namefile);
        mysprintf(namefile, NAMESIZE, "%s/%s_corr_bolla_r.opposite_mu.%.5f_VKVK.txt", argv[3], argv[4], mu);//54
        correlators.emplace_back(namefile);
    }


    // printf("reading confs from file: %s", correlators[0].c_str());
    // auto iconfs = read_nconfs(correlators[0].c_str());
    std::vector<configuration_class> myconfs;

    // confs = myconfs[0].iconfs.size();
    int count = 0;
    for (auto name : correlators) {
        printf("reading  confs from file: %s\n", name.c_str());
        if (count == 42) {
            if (strcmp(argv[4], "cB.72.96") == 0) {
                configuration_class tmp;
                myconfs.emplace_back(tmp);
            }
            else {
                configuration_class tmp(name.c_str(), 0);
                myconfs.emplace_back(tmp);
            }
        }
        else if (count == 55 || count == 56) {
            configuration_class tmp(name.c_str(), header.T - 1 - 3 + 1);
            myconfs.emplace_back(tmp);
        }
        else
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
        // read correlators[i] and store in data[conf][i][t][re/im]
        if (strcmp(argv[4], "cB.72.96") != 0) {
            if (i == 42) { printf("at time zero\n");    read_twopt(correlators[i].c_str(), myconfs[i], 0, data, i, bin); }
            else read_twopt(correlators[i].c_str(), myconfs[i], T, data, i, bin);
        }
        else {
            printf("here %d\n", i);
            if (i == 42)     printf("skip\n");
            else if (i >= 49 && i <= 54)     printf("skip\n");
            else read_twopt(correlators[i].c_str(), myconfs[i], T, data, i, bin);
        }


    }
    if (argc > 18 && strcmp(argv[17], "three_corr") != 0 && strcmp(argv[4], "cB.72.96") != 0) {
        int Nconf_bolla = myconfs[42].iconfs.size();
        double**** data_no_bin = calloc_corr(Nconf_bolla, 4 + 3, file_head.l0);
        // error(Nconf_bolla!= myconfs[53].iconfs.size(),1,"main","conf of P5P5 correlated bolla are not the same of bolla");
        // error(Nconf_bolla!= myconfs[51].iconfs.size(),1,"main","conf of VKVKeq correlated bolla are not the same of bolla");
        // error(Nconf_bolla!= myconfs[54].iconfs.size(),1,"main","conf of VKVKop correlated bolla are not the same of bolla");
        check_confs_correlated(myconfs, correlators, 42, 53);
        check_confs_correlated(myconfs, correlators, 42, 51);
        check_confs_correlated(myconfs, correlators, 42, 54);
        for (int i = 0;i < 6;i++)
            check_confs_correlated(myconfs, correlators, 36 + i, 43 + i);

        read_twopt(correlators[42].c_str(), myconfs[42], 0, data_no_bin, 0, Nconf_bolla);
        read_twopt(correlators[53].c_str(), myconfs[53], T, data_no_bin, 1, Nconf_bolla);
        read_twopt(correlators[51].c_str(), myconfs[51], T, data_no_bin, 2, Nconf_bolla);
        read_twopt(correlators[54].c_str(), myconfs[54], T, data_no_bin, 3, Nconf_bolla);


        for (int b = 0; b < Nconf_bolla; b++) {
            data_no_bin[b][0][0][0] = -data_no_bin[b][0][0][0];
            for (int t = 0;t < T / 2 + 1;t++) {
                data_no_bin[b][0][t][0] = data_no_bin[b][0][0][0];
                data_no_bin[b][4][t][0] = data_no_bin[b][0][0][0] * data_no_bin[b][1][t][0]; // bolla * P5P5_op // 55
                data_no_bin[b][5][t][0] = data_no_bin[b][0][0][0] * data_no_bin[b][2][t][0]; // bolla * VKVK_eq // 56
                data_no_bin[b][6][t][0] = data_no_bin[b][0][0][0] * data_no_bin[b][3][t][0]; // bolla * VKVK_op // 57
            }
        }


        bin_data(data, 42, data_no_bin, 0, T, Nconf_bolla, bin); // bolla 
        bin_data(data, var + 0, data_no_bin, 4, T, Nconf_bolla, bin); // bolla * P5P5_op // 55
        bin_data(data, var + 1, data_no_bin, 5, T, Nconf_bolla, bin); // bolla * P5P5_op // 56
        bin_data(data, var + 2, data_no_bin, 6, T, Nconf_bolla, bin); // bolla * P5P5_op // 57
        error(var + 0 != 55, 1, "main", "error index do not match var+0=%d  expected 57", var + 0);
        error(var + 1 != 56, 1, "main", "error index do not match var+1=%d  expected 57", var + 1);
        error(var + 2 != 57, 1, "main", "error index do not match var+2=%d  expected 57", var + 2);
        free_corr(Nconf_bolla, 4 + 3, file_head.l0, data_no_bin);

    }

    correlators.emplace_back("bP5P5");
    correlators.emplace_back("bVKVKeq");
    correlators.emplace_back("bVKVKop");
    myconfs.emplace_back();
    myconfs.emplace_back();
    myconfs.emplace_back();

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

    if (strcmp(argv[argc - 1], "three_corr") == 0) {
        for (int q = 0; q < 1 + 2 + 3;q++) {
            int ir = 0;
            mysprintf(namefile, NAMESIZE, "%s/Vkvk_cont/%d_m%s/SAMER", argv[3], header.L, argv[11 + q]);
            printf("reading: %s   %d\n", namefile, var + 3 + q * 2);
            correlators.emplace_back(namefile);
            FILE* SAMER = open_file(namefile, "r+");
            ir += fscanf(SAMER, "%*[^\n]\n");// skip a line
            for (int t = 0; t < header.T / 2 + 1;t++) {
                int tr;
                double a, b;
                ir += fscanf(SAMER, "%d   %lf  %lf  %lf\n", &tr, &data[0][var + 3 + q * 2][t][0], &a, &b);
                error(tr != t, 1, "read free corr", "expected t=%d   , read t=%d", tr, t);
                for (int j = 1;j < bin;j++) {
                    data[j][var + 3 + q * 2][t][0] = data[0][var + 3 + q * 2][t][0];
                }
            }

            mysprintf(namefile, NAMESIZE, "%s/Vkvk_cont/%d_m%s/OPPOR", argv[3], header.L, argv[11 + q]);
            printf("reading: %s   %d\n", namefile, var + 4 + q * 2);
            correlators.emplace_back(namefile);
            FILE* OPPOR = open_file(namefile, "r+");
            ir += fscanf(OPPOR, "%*[^\n]\n");// skip a line
            for (int t = 0; t < header.T / 2 + 1;t++) {
                int tr;
                double a, b;
                ir += fscanf(OPPOR, "%d   %lf  %lf  %lf\n", &tr, &data[0][var + 4 + q * 2][t][0], &a, &b);
                error(tr != t, 1, "read free corr", "expected t=%d   , read t=%d", tr, t);
                for (int j = 1;j < bin;j++) {
                    data[j][var + 4 + q * 2][t][0] = data[0][var + 4 + q * 2][t][0];
                }
            }
            // keep myconf of the same lenght of correlators for later usage
            myconfs.emplace_back();
            myconfs.emplace_back();
        }
        error(var + 15 != correlators.size(), 1, "main", "wrong counting");
    }
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///// adding D meson here
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // declare
    int idD = correlators.size();
    mysprintf(namefile, NAMESIZE, "%s/%s_r.opposite_mu.%.5f_mu.%.3f_P5P5.txt", argv[3], argv[4], mu, mus1);//var+15 
    correlators.emplace_back(namefile);
    mysprintf(namefile, NAMESIZE, "%s/%s_r.opposite_mu.%.5f_mu.%.3f_P5P5.txt", argv[3], argv[4], mu, mus2);//var+16
    correlators.emplace_back(namefile);

    // parse and read
    for (int i = idD; i < idD + 2;i++) {
        // check if file exist
        FILE* tmp = NULL;
        tmp = fopen(correlators[i].c_str(), "r");
        // read if file exist
        if (tmp != NULL) {
            printf("reading  confs from file: %s\n", correlators[i].c_str());
            myconfs.emplace_back(correlators[i].c_str());
            printf("%d  %ld \n", i, myconfs.size());
            myconfs[i].check_binnign();
            cout << "number of different configurations:" << myconfs[i].confs_after_binning << endl;
            read_twopt(correlators[i].c_str(), myconfs[i], T, data, i, bin);
            fclose(tmp);
        }

    }


    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // D meson reading end
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////


    // data_bin = binning(confs, var, file_head.l0, data, bin);
    // data_bin = binning_toNb(confs, var, file_head.l0, data, bin);
    // //if you want to do the gamma analysis you need to do before freeing the raw data
    // // effective_mass_phi4_gamma(option, kinematic_2pt, (char*)"P5P5", data_bin, Neff, namefile_plateaux, out_gamma, 0, "M_{PS}^{ll}");
    // // effective_mass_phi4_gamma(option, kinematic_2pt, (char*)"P5P5", data_bin, Neff, namefile_plateaux, out_gamma, 1, "M_{PS1}^{ll}");
    // //effective_mass_phi4_gamma(  option, kinematic_2pt,   (char*) "P5P5", data,  confs ,namefile_plateaux,out_gamma,3,"M_{PS}^{ll}");
    // free_corr(bin, var, file_head.l0, data);
    ncorr_new = correlators.size();
    conf_jack = create_resampling(option[4], Neff, correlators.size(), file_head.l0, data);
    free_corr(Neff, correlators.size(), file_head.l0, data);


    // ////////////////// symmetrization/////////////////////////////////////////////
    // for (int i = 0;i <= 7;i++) { symmetrise_jackboot(Njack, i, file_head.l0, conf_jack); }

    ////////////////////////////////////////////////
    corr_counter = -1;
    double* M_PS = plateau_correlator_function(option, kinematic_2pt, (char*)"P5P5", conf_jack, Njack, namefile_plateaux, outfile, 1, "M_{PS}^{eq}", M_eff_T, jack_file);
    check_correlatro_counter(0);

    double* M_PS_op = plateau_correlator_function(option, kinematic_2pt, (char*)"P5P5", conf_jack, Njack, namefile_plateaux, outfile, 4, "M_{PS}^{op}", M_eff_T, jack_file);
    check_correlatro_counter(1);

    double* M_eta = plateau_correlator_function(option, kinematic_2pt, (char*)"P5P5", conf_jack, Njack, namefile_plateaux, outfile, 1 + 6, "M_{eta}^{eq}", M_eff_T, jack_file);
    check_correlatro_counter(2);

    double* M_eta_op = plateau_correlator_function(option, kinematic_2pt, (char*)"P5P5", conf_jack, Njack, namefile_plateaux, outfile, 4 + 6, "M_{eta}^{op}", M_eff_T, jack_file);
    check_correlatro_counter(3);

    double* M_phi = plateau_correlator_function(option, kinematic_2pt, (char*)"P5P5", conf_jack, Njack, namefile_plateaux, outfile, 2 + 6, "M_{phi}^{eq}", M_eff_T, jack_file);
    check_correlatro_counter(4);

    double* M_phi_op = plateau_correlator_function(option, kinematic_2pt, (char*)"P5P5", conf_jack, Njack, namefile_plateaux, outfile, 5 + 6, "M_{phi}^{op}", M_eff_T, jack_file);
    check_correlatro_counter(5);

    double* M_eta1 = plateau_correlator_function(option, kinematic_2pt, (char*)"P5P5", conf_jack, Njack, namefile_plateaux, outfile, 1 + 12, "M_{eta1}^{eq}", M_eff_T, jack_file);
    check_correlatro_counter(6);

    double* M_eta1_op = plateau_correlator_function(option, kinematic_2pt, (char*)"P5P5", conf_jack, Njack, namefile_plateaux, outfile, 4 + 12, "M_{eta1}^{op}", M_eff_T, jack_file);
    check_correlatro_counter(7);

    double* M_phi1 = plateau_correlator_function(option, kinematic_2pt, (char*)"P5P5", conf_jack, Njack, namefile_plateaux, outfile, 2 + 12, "M_{phi1}^{eq}", M_eff_T, jack_file);
    check_correlatro_counter(8);

    double* M_phi1_op = plateau_correlator_function(option, kinematic_2pt, (char*)"P5P5", conf_jack, Njack, namefile_plateaux, outfile, 5 + 12, "M_{phi1}^{op}", M_eff_T, jack_file);
    check_correlatro_counter(9);


    fit_type fit_info;
    fit_result fit_out;

    fit_info.Nvar = 1;
    fit_info.Npar = 1;
    fit_info.N = 1;
    fit_info.Njack = Njack;
    fit_info.n_ext_P = 1;
    fit_info.ext_P = (double**)malloc(sizeof(double*) * 1);
    fit_info.function = constant_fit;

    fit_info.ext_P[0] = M_PS;
    fit_result G_PS = fit_fun_to_fun_of_corr(option, kinematic_2pt, (char*)"A1P1phi", conf_jack, namefile_plateaux, outfile, GPS_lhs<1>, "G_{PS}^{eq}", fit_info, jack_file);
    check_correlatro_counter(10);


    fit_info.ext_P[0] = M_PS_op;
    fit_result G_PS_OS = fit_fun_to_fun_of_corr(option, kinematic_2pt, (char*)"A1P1phi", conf_jack, namefile_plateaux, outfile, GPS_OS_lhs<4>, "G_{PS}^{op}", fit_info, jack_file);
    check_correlatro_counter(11);

    fit_info.ext_P[0] = M_eta;
    fit_result G_eta = fit_fun_to_fun_of_corr(option, kinematic_2pt, (char*)"A1P1phi", conf_jack, namefile_plateaux, outfile, GPS_lhs<1 + 6>, "G_{eta}^{eq}", fit_info, jack_file);
    check_correlatro_counter(12);


    fit_info.ext_P[0] = M_eta_op;
    fit_result G_eta_OS = fit_fun_to_fun_of_corr(option, kinematic_2pt, (char*)"A1P1phi", conf_jack, namefile_plateaux, outfile, GPS_OS_lhs<4 + 6>, "G_{eta}^{op}", fit_info, jack_file);
    check_correlatro_counter(13);

    fit_info.ext_P[0] = M_eta1;
    fit_result G_eta1 = fit_fun_to_fun_of_corr(option, kinematic_2pt, (char*)"A1P1phi", conf_jack, namefile_plateaux, outfile, GPS_lhs<1 + 12>, "G_{eta1}^{eq}", fit_info, jack_file);
    check_correlatro_counter(14);


    fit_info.ext_P[0] = M_eta1_op;
    fit_result G_eta1_OS = fit_fun_to_fun_of_corr(option, kinematic_2pt, (char*)"A1P1phi", conf_jack, namefile_plateaux, outfile, GPS_OS_lhs<4 + 12>, "G_{eta1}^{op}", fit_info, jack_file);
    check_correlatro_counter(15);

    //////////////////////////////////////////////////
    fit_info.n_ext_P = 0;
    fit_info.mu = mu;
    fit_result ZVl = fit_fun_to_fun_of_corr(option, kinematic_2pt, (char*)"A1P1phi", conf_jack, namefile_plateaux, outfile, ZVl_lhs<4, 3>, "Z_V(l)", fit_info, jack_file);
    check_correlatro_counter(16);
    fit_info.mu = mus1;
    fit_result ZVs = fit_fun_to_fun_of_corr(option, kinematic_2pt, (char*)"A1P1phi", conf_jack, namefile_plateaux, outfile, ZVl_lhs<4 + 6, 3 + 6>, "Z_V(s)", fit_info, jack_file);
    check_correlatro_counter(17);
    fit_info.mu = mus2;
    fit_result ZVs1 = fit_fun_to_fun_of_corr(option, kinematic_2pt, (char*)"A1P1phi", conf_jack, namefile_plateaux, outfile, ZVl_lhs<4 + 12, 3 + 12>, "Z_V(s1)", fit_info, jack_file);
    check_correlatro_counter(18);


    fit_info.ext_P[0] = nullptr;
    free(fit_info.ext_P);
    fit_info.n_ext_P = 4;
    fit_info.ext_P = (double**)malloc(sizeof(double*) * fit_info.n_ext_P);

    fit_info.ext_P[0] = M_PS_op;
    fit_info.ext_P[1] = M_PS;
    fit_info.ext_P[2] = G_PS_OS.P[0];
    fit_info.ext_P[3] = G_PS.P[0];
    fit_info.mu = mu;
    fit_result ZAl = fit_fun_to_fun_of_corr(option, kinematic_2pt, (char*)"A1P1phi", conf_jack, namefile_plateaux, outfile, ZAl_lhs<1, 0>, "Z_A(l)", fit_info, jack_file);
    check_correlatro_counter(19);


    fit_info.ext_P[0] = M_eta_op;
    fit_info.ext_P[1] = M_eta;
    fit_info.ext_P[2] = G_eta_OS.P[0];
    fit_info.ext_P[3] = G_eta.P[0];
    fit_info.mu = mus1;
    fit_result ZAs = fit_fun_to_fun_of_corr(option, kinematic_2pt, (char*)"A1P1phi", conf_jack, namefile_plateaux, outfile, ZAl_lhs<1 + 6, 0 + 6>, "Z_A(s)", fit_info, jack_file);
    check_correlatro_counter(20);

    fit_info.ext_P[0] = M_eta1_op;
    fit_info.ext_P[1] = M_eta1;
    fit_info.ext_P[2] = G_eta1_OS.P[0];
    fit_info.ext_P[3] = G_eta1.P[0];
    fit_info.mu = mus2;
    fit_result ZAs1 = fit_fun_to_fun_of_corr(option, kinematic_2pt, (char*)"A1P1phi", conf_jack, namefile_plateaux, outfile, ZAl_lhs<1 + 12, 0 + 12>, "Z_A(s1)", fit_info, jack_file);
    check_correlatro_counter(21);

    fit_info.restore_default();

    // double* amu_sd = (double*)malloc(sizeof(double) * Njack);
    double* VV = (double*)malloc(sizeof(double) * file_head.l0);
    double mean, err;
    int seed;
    line_read_param(option, "a", mean, err, seed, namefile_plateaux);
    double* a = fake_sampling(resampling, mean, err, Njack, seed);

    double* jack_Metas_MeV_exp = fake_sampling(resampling, Metas_MeV, Metas_MeV_err, Njack, 1000);

    double* jack_Mrho_MeV_exp = fake_sampling(resampling, Mrho_MeV, Mrho_MeV_err, Njack, 1001);
    double* jack_grhopipi = fake_sampling(resampling, grhopipi, grhopipi_err, Njack, 1002);
    double* jack_Mpi_MeV_exp = fake_sampling(resampling, Mpi_MeV, Mpi_MeV_err, Njack, 1003);
    double* jack_Mphi_MeV_exp = fake_sampling(resampling, Mphi_MeV, Mphi_MeV_err, Njack, 1004);
    double* jack_Jpsi_MeV_exp = fake_sampling(resampling, Jpsi_MeV, Jpsi_MeV_err, Njack, 1005);
    double* jack_Metac_MeV_exp = fake_sampling(resampling, Metac_MeV, Metac_MeV_err, Njack, 1006);

    double* jack_MK_MeV_exp = fake_sampling(resampling, MK_MeV, MK_MeV_err, Njack, 1007);

    double* zeros = (double*)calloc(Njack, sizeof(double));

    // line_read_param(option, "ZA", mean, err, seed, namefile_plateaux);
    // double* ZA = fake_sampling(resampling, mean, err, Njack, seed);
    // line_read_param(option, "ZV", mean, err, seed, namefile_plateaux);
    // double* ZV = fake_sampling(resampling, mean, err, Njack, seed);

    double* jack_aMetas_MeV_exp = (double*)malloc(sizeof(double) * Njack);
    double* jack_aMetas2_MeV_exp = (double*)malloc(sizeof(double) * Njack);
    for (int j = 0;j < Njack;j++) {
        jack_aMetas_MeV_exp[j] = jack_Metas_MeV_exp[j] * a[j] / 197.326963;
        jack_aMetas2_MeV_exp[j] = jack_Metas_MeV_exp[j] * jack_Metas_MeV_exp[j];
    }
    double* jack_aMphi_MeV_exp = (double*)malloc(sizeof(double) * Njack);
    for (int j = 0;j < Njack;j++) {
        jack_aMphi_MeV_exp[j] = jack_Mphi_MeV_exp[j] * a[j] / 197.326963;
    }
    double* jack_aJpsi_MeV_exp = (double*)malloc(sizeof(double) * Njack);
    for (int j = 0;j < Njack;j++) {
        jack_aJpsi_MeV_exp[j] = jack_Jpsi_MeV_exp[j] * a[j] / 197.326963;
    }
    double* jack_aMetac_MeV_exp = (double*)malloc(sizeof(double) * Njack);
    for (int j = 0;j < Njack;j++) {
        jack_aMetac_MeV_exp[j] = jack_Metac_MeV_exp[j] * a[j] / 197.326963;
    }
    double* jack_aMK_MeV_exp = (double*)malloc(sizeof(double) * Njack);
    for (int j = 0;j < Njack;j++) {
        jack_aMK_MeV_exp[j] = jack_MK_MeV_exp[j] * a[j] / 197.326963;
    }

    int Nstrange = 2;
    double** Meta = (double**)malloc(sizeof(double*) * Nstrange);
    Meta[0] = M_eta_op;
    Meta[1] = M_eta1_op;

    double** Mphi = (double**)malloc(sizeof(double*) * Nstrange);
    Mphi[0] = M_phi_op;
    Mphi[1] = M_phi1_op;

    double** Z = (double**)malloc(sizeof(double*) * Nstrange);
    Z[0] = ZVs.P[0];
    Z[1] = ZVs1.P[0];

    double* ZV, * ZA;

    if (strcmp("cA.53.24", argv[4]) == 0 || strcmp("cA.40.24", argv[4]) == 0 || strcmp("cA.30.32", argv[4]) == 0) {
        line_read_param(option, "ZA", mean, err, seed, namefile_plateaux);
        ZA = fake_sampling(resampling, mean, err, Njack, seed);
        line_read_param(option, "ZV", mean, err, seed, namefile_plateaux);
        ZV = fake_sampling(resampling, mean, err, Njack, seed);
    }
    else {
        Z[0] = ZVs.P[0];
        Z[1] = ZVs1.P[0];
        ZV = interpol_Z(Nstrange, Njack, Meta, Z, jack_aMetas_MeV_exp, outfile, "Z_V", resampling);


        Z[0] = ZAs.P[0];
        Z[1] = ZAs1.P[0];
        ZA = interpol_Z(Nstrange, Njack, Meta, Z, jack_aMetas_MeV_exp, outfile, "Z_A", resampling);

    }
    // if (strcmp("cD.54.96", argv[4]) == 0) {
    //     double Za_WI_strange = 0.773944;
    //     double Za_WI_strange_err = 0.00014;
    //     free(ZA);
    //     ZA = myres->create_fake(Za_WI_strange, Za_WI_strange_err, 666);
    // }

    write_jack(ZV, Njack, jack_file);
    check_correlatro_counter(22);
    write_jack(ZA, Njack, jack_file);
    check_correlatro_counter(23);
    ///
    double** musj = double_malloc_2(Nstrange, Njack);
    for (int i = 0;i < Nstrange;i++) {
        for (int j = 0;j < Njack;j++) {
            musj[i][j] = header.mus[i + 1];
        }

    }
    double* mus_phys = interpol_Z(Nstrange, Njack, Meta, musj, jack_aMetas_MeV_exp, outfile, "mu_s_phys", resampling);
    write_jack(mus_phys, Njack, jack_file);
    check_correlatro_counter(24);

    free_2(Nstrange, musj);
    for (int i = 0;i < Nstrange;i++) {
        Z[i] = nullptr;
    }
    free(Z);

    fprintf(ASCII_Z, "mu %g  %g  %g\n", header.mus[0], header.mus[1], header.mus[2]);
    fprintf(ASCII_Z, "ZAl  ZAms1  ZAms2  ZA   ZVl  ZVms1  ZVms2  ZV\n");
    for (int j = 0;j < Njack;j++) {
        fprintf(ASCII_Z, "%-16.12g%-16.12g%-16.12g%-16.12g%-16.12g%-16.12g%-16.12g%-16.12g\n"
            , ZAl.P[0][j], ZAs.P[0][j], ZAs1.P[0][j], ZA[j]
            , ZVl.P[0][j], ZVs.P[0][j], ZVs1.P[0][j], ZV[j]);
    }


    ////////////////////////////////////////////////////////////////////////////////////////////////////////
    // a_SD_l
    ///////////////////////////////////////////////////////////////////////////////////////////////////////

    double (*int_scheme)(int, int, double*);


    int_scheme = integrate_reinman;
    int isub = (strcmp(argv[argc - 1], "three_corr") == 0) ? var + 3 : -1;
    double* amu_sd = compute_amu_sd(conf_jack, 2, Njack, ZV, a, 5.0 / 9.0, int_scheme, outfile, "amu_{sd}(eq,l)", resampling, isub);
    write_jack(amu_sd, Njack, jack_file);
    check_correlatro_counter(25);
    printf("amu_sd(eq,l) = %g  %g\n", amu_sd[Njack - 1], error_jackboot(resampling, Njack, amu_sd));
    free(amu_sd);

    isub = (strcmp(argv[argc - 1], "three_corr") == 0) ? var + 4 : -1;
    amu_sd = compute_amu_sd(conf_jack, 5, Njack, ZA, a, 5.0 / 9.0, int_scheme, outfile, "amu_{sd}(op,l)", resampling, isub);
    write_jack(amu_sd, Njack, jack_file);
    check_correlatro_counter(26);
    printf("amu_sd(op,l) = %g  %g\n", amu_sd[Njack - 1], error_jackboot(resampling, Njack, amu_sd));
    free(amu_sd);


    int_scheme = integrate_simpson38;
    isub = (strcmp(argv[argc - 1], "three_corr") == 0) ? var + 3 : -1;
    amu_sd = compute_amu_sd(conf_jack, 2, Njack, ZV, a, 5.0 / 9.0, int_scheme, outfile, "amu_{sd,simpson38}(eq,l)", resampling, isub);
    write_jack(amu_sd, Njack, jack_file);
    check_correlatro_counter(27);
    printf("amu_sd_simpson38(eq,l) = %g  %g\n", amu_sd[Njack - 1], error_jackboot(resampling, Njack, amu_sd));
    free(amu_sd);

    isub = (strcmp(argv[argc - 1], "three_corr") == 0) ? var + 4 : -1;
    amu_sd = compute_amu_sd(conf_jack, 5, Njack, ZA, a, 5.0 / 9.0, int_scheme, outfile, "amu_{sd,simpson38}(op,l)", resampling, isub);
    write_jack(amu_sd, Njack, jack_file);
    check_correlatro_counter(28);
    printf("amu_sd_simpson38(op,l) = %g  %g\n", amu_sd[Njack - 1], error_jackboot(resampling, Njack, amu_sd));
    free(amu_sd);

    ////////////////////////////////////////////////////////////////////////////////////////////////////////
    // a_SD_s_eq reinman
    ///////////////////////////////////////////////////////////////////////////////////////////////////////
    constexpr double q2s = 1.0 / 9.0;
    int_scheme = integrate_reinman;
    isub = (strcmp(argv[argc - 1], "three_corr") == 0) ? var + 3 + 1 * 2 : -1;
    double* amu_sdeq_s = compute_amu_sd(conf_jack, 2 + 6, Njack, ZVs.P[0], a, q2s, int_scheme, outfile, "amu_{sd}(eq,s)", resampling, isub);
    write_jack(amu_sdeq_s, Njack, jack_file);
    check_correlatro_counter(29);
    printf("amu_sd(eq,s) = %g  %g\n", amu_sdeq_s[Njack - 1], error_jackboot(resampling, Njack, amu_sdeq_s));


    int_scheme = integrate_reinman;
    isub = (strcmp(argv[argc - 1], "three_corr") == 0) ? var + 3 + 2 * 2 : -1;
    double* amu_sdeq_s1 = compute_amu_sd(conf_jack, 2 + 12, Njack, ZVs1.P[0], a, q2s, int_scheme, outfile, "amu_{sd}(eq,s1)", resampling, isub);
    write_jack(amu_sdeq_s1, Njack, jack_file);
    check_correlatro_counter(30);
    printf("amu_sd(eq,s1) = %g  %g\n", amu_sdeq_s1[Njack - 1], error_jackboot(resampling, Njack, amu_sdeq_s1));

    double** asd_vec = (double**)malloc(sizeof(double*) * Nstrange);
    asd_vec[0] = amu_sdeq_s;
    asd_vec[1] = amu_sdeq_s1;

    double* amu_sd_sphys = interpol_Z(Nstrange, Njack, Meta, asd_vec, jack_aMetas_MeV_exp, outfile, "amu_{sd}(eq,sphys)", resampling);
    write_jack(amu_sd_sphys, Njack, jack_file);
    printf("amu_sd(eq,shys) = %g  %g\n", amu_sd_sphys[Njack - 1], error_jackboot(resampling, Njack, amu_sd_sphys));

    check_correlatro_counter(31);

    ////////////////////////////////////////////////////////////////////////////////////////////////////////
    // a_SD_s_op reinman
    ///////////////////////////////////////////////////////////////////////////////////////////////////////

    int_scheme = integrate_reinman;
    isub = (strcmp(argv[argc - 1], "three_corr") == 0) ? var + 4 + 1 * 2 : -1;
    double* amu_sdop_s = compute_amu_sd(conf_jack, 5 + 6, Njack, ZAs.P[0], a, q2s, int_scheme, outfile, "amu_{sd}(op,s)", resampling, isub);
    write_jack(amu_sdop_s, Njack, jack_file);
    check_correlatro_counter(32);
    printf("amu_sd(op,s) = %g  %g\n", amu_sdop_s[Njack - 1], error_jackboot(resampling, Njack, amu_sdop_s));


    int_scheme = integrate_reinman;
    isub = (strcmp(argv[argc - 1], "three_corr") == 0) ? var + 4 + 2 * 2 : -1;
    double* amu_sdop_s1 = compute_amu_sd(conf_jack, 5 + 12, Njack, ZAs1.P[0], a, q2s, int_scheme, outfile, "amu_{sd}(op,s1)", resampling, isub);
    write_jack(amu_sdop_s1, Njack, jack_file);
    check_correlatro_counter(33);
    printf("amu_sd(op,s1) = %g  %g\n", amu_sdop_s1[Njack - 1], error_jackboot(resampling, Njack, amu_sdop_s1));


    asd_vec[0] = amu_sdop_s;
    asd_vec[1] = amu_sdop_s1;
    amu_sd_sphys = interpol_Z(Nstrange, Njack, Meta, asd_vec, jack_aMetas_MeV_exp, outfile, "amu_{sd}(op,sphys)", resampling);
    write_jack(amu_sd_sphys, Njack, jack_file);
    printf("amu_sd(op,shys) = %g  %g\n", amu_sd_sphys[Njack - 1], error_jackboot(resampling, Njack, amu_sd_sphys));

    check_correlatro_counter(34);

    ////////////////////////////////////////////////////////////////////////////////////////////////////////
    // a_SD_s_eq simpson38
    ///////////////////////////////////////////////////////////////////////////////////////////////////////

    int_scheme = integrate_simpson38;
    isub = (strcmp(argv[argc - 1], "three_corr") == 0) ? var + 3 + 1 * 2 : -1;
    double* amu_sdeq_simp_s = compute_amu_sd(conf_jack, 2 + 6, Njack, ZVs.P[0], a, q2s, int_scheme, outfile, "amu_{sd}_simpson38(eq,s)", resampling, isub);
    write_jack(amu_sdeq_simp_s, Njack, jack_file);
    check_correlatro_counter(35);
    printf("amu_sd_simpson38(eq,s) = %g  %g\n", amu_sdeq_simp_s[Njack - 1], error_jackboot(resampling, Njack, amu_sdeq_simp_s));


    int_scheme = integrate_simpson38;
    isub = (strcmp(argv[argc - 1], "three_corr") == 0) ? var + 3 + 2 * 2 : -1;
    double* amu_sdeq_simp_s1 = compute_amu_sd(conf_jack, 2 + 12, Njack, ZVs1.P[0], a, q2s, int_scheme, outfile, "amu_{sd}_simpson38(eq,s1)", resampling, isub);
    write_jack(amu_sdeq_simp_s1, Njack, jack_file);
    check_correlatro_counter(36);
    printf("amu_sd_simpson38(eq,s1) = %g  %g\n", amu_sdeq_simp_s1[Njack - 1], error_jackboot(resampling, Njack, amu_sdeq_simp_s1));

    asd_vec[0] = amu_sdeq_simp_s;
    asd_vec[1] = amu_sdeq_simp_s1;

    amu_sd_sphys = interpol_Z(Nstrange, Njack, Meta, asd_vec, jack_aMetas_MeV_exp, outfile, "amu_{sd}_simpson38(eq,sphys)", resampling);
    write_jack(amu_sd_sphys, Njack, jack_file);
    printf("amu_sd(eq,shys) = %g  %g\n", amu_sd_sphys[Njack - 1], error_jackboot(resampling, Njack, amu_sd_sphys));
    free(amu_sd_sphys);
    check_correlatro_counter(37);


    ////////////////////////////////////////////////////////////////////////////////////////////////////////
    // a_SD_s_op simpson38
    ///////////////////////////////////////////////////////////////////////////////////////////////////////

    int_scheme = integrate_simpson38;
    isub = (strcmp(argv[argc - 1], "three_corr") == 0) ? var + 4 + 1 * 2 : -1;
    double* amu_sdop_simp_s = compute_amu_sd(conf_jack, 5 + 6, Njack, ZAs.P[0], a, q2s, int_scheme, outfile, "amu_{sd}_simpson38(op,s)", resampling, isub);
    write_jack(amu_sdop_simp_s, Njack, jack_file);
    check_correlatro_counter(38);
    printf("amu_sd_simpson38(op,s) = %g  %g\n", amu_sdop_simp_s[Njack - 1], error_jackboot(resampling, Njack, amu_sdop_simp_s));


    int_scheme = integrate_simpson38;
    isub = (strcmp(argv[argc - 1], "three_corr") == 0) ? var + 4 + 2 * 2 : -1;
    double* amu_sdop_simp_s1 = compute_amu_sd(conf_jack, 5 + 12, Njack, ZAs1.P[0], a, q2s, int_scheme, outfile, "amu_{sd}_simpson38(op,s1)", resampling, isub);
    write_jack(amu_sdop_simp_s1, Njack, jack_file);
    check_correlatro_counter(39);
    printf("amu_sd_simpson38(op,s1) = %g  %g\n", amu_sdop_simp_s1[Njack - 1], error_jackboot(resampling, Njack, amu_sdop_simp_s1));


    asd_vec[0] = amu_sdop_simp_s;
    asd_vec[1] = amu_sdop_simp_s1;
    amu_sd_sphys = interpol_Z(Nstrange, Njack, Meta, asd_vec, jack_aMetas_MeV_exp, outfile, "amu_{sd}_simpson38(op,sphys)", resampling);
    write_jack(amu_sd_sphys, Njack, jack_file);
    printf("amu_sd_simpson38(op,shys) = %g  %g\n", amu_sd_sphys[Njack - 1], error_jackboot(resampling, Njack, amu_sd_sphys));
    check_correlatro_counter(40);
    free(amu_sd_sphys);


    write_jack(a, Njack, jack_file);
    check_correlatro_counter(41);
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////
    printf("\n WINDOW \n\n");
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////


    ////////////////////////////////////////////////////////////////////////////////////////////////////////
    // a_W_l
    ///////////////////////////////////////////////////////////////////////////////////////////////////////

    int_scheme = integrate_reinman;
    double* amu_W_eq0 = compute_amu_W(conf_jack, 2, Njack, ZV, a, 5.0 / 9.0, int_scheme, outfile, "amu_{W}(eq,l)", resampling);
    write_jack(amu_W_eq0, Njack, jack_file);
    check_correlatro_counter(42);
    printf("amu_W(eq1,l) = %g  %g\n", amu_W_eq0[Njack - 1], error_jackboot(resampling, Njack, amu_W_eq0));

    double* amu_W_op0 = compute_amu_W(conf_jack, 5, Njack, ZA, a, 5.0 / 9.0, int_scheme, outfile, "amu_{W}(op,l)", resampling);
    write_jack(amu_W_op0, Njack, jack_file);
    check_correlatro_counter(43);
    printf("amu_W(op,l) = %g  %g\n", amu_W_op0[Njack - 1], error_jackboot(resampling, Njack, amu_W_op0));



    int_scheme = integrate_simpson38;
    double* amu_W = compute_amu_W(conf_jack, 2, Njack, ZV, a, 5.0 / 9.0, int_scheme, outfile, "amu_{W,simpson38}(eq,l)", resampling);
    write_jack(amu_W, Njack, jack_file);
    check_correlatro_counter(44);
    printf("amu_W_simpson38(eq,l) = %g  %g\n", amu_W[Njack - 1], error_jackboot(resampling, Njack, amu_W));
    free(amu_W);

    amu_W = compute_amu_W(conf_jack, 5, Njack, ZA, a, 5.0 / 9.0, int_scheme, outfile, "amu_{W,simpson38}(op,l)", resampling);
    write_jack(amu_W, Njack, jack_file);
    check_correlatro_counter(45);
    printf("amu_W_simpson38(op,l) = %g  %g\n", amu_W[Njack - 1], error_jackboot(resampling, Njack, amu_W));
    free(amu_W);


    ////////////////////////////////////////////////////////////////////////////////////////////////////////
    // a_W_s_eq reinman
    ///////////////////////////////////////////////////////////////////////////////////////////////////////
    int_scheme = integrate_reinman;
    double* amu_Weq_s = compute_amu_W(conf_jack, 2 + 6, Njack, ZVs.P[0], a, q2s, int_scheme, outfile, "amu_{W}(eq,s)", resampling);
    write_jack(amu_Weq_s, Njack, jack_file);
    check_correlatro_counter(46);
    printf("amu_W(eq,s) = %g  %g\n", amu_Weq_s[Njack - 1], error_jackboot(resampling, Njack, amu_Weq_s));


    int_scheme = integrate_reinman;
    double* amu_Weq_s1 = compute_amu_W(conf_jack, 2 + 12, Njack, ZVs1.P[0], a, q2s, int_scheme, outfile, "amu_{W}(eq,s1)", resampling);
    write_jack(amu_Weq_s1, Njack, jack_file);
    check_correlatro_counter(47);
    printf("amu_W(eq,s1) = %g  %g\n", amu_Weq_s1[Njack - 1], error_jackboot(resampling, Njack, amu_Weq_s1));

    asd_vec[0] = amu_Weq_s;
    asd_vec[1] = amu_Weq_s1;
    double* amu_W_sphys = interpol_Z(Nstrange, Njack, Meta, asd_vec, jack_aMetas_MeV_exp, outfile, "amu_{W}(eq,sphys)", resampling);
    write_jack(amu_W_sphys, Njack, jack_file);
    printf("amu_W(eq,shys) = %g  %g\n", amu_W_sphys[Njack - 1], error_jackboot(resampling, Njack, amu_W_sphys));

    check_correlatro_counter(48);


    ////////////////////////////////////////////////////////////////////////////////////////////////////////
    // a_W_s_op reinman
    ///////////////////////////////////////////////////////////////////////////////////////////////////////

    int_scheme = integrate_reinman;
    double* amu_Wop_s = compute_amu_W(conf_jack, 5 + 6, Njack, ZAs.P[0], a, q2s, int_scheme, outfile, "amu_{W}(op,s)", resampling);
    write_jack(amu_Wop_s, Njack, jack_file);
    check_correlatro_counter(49);
    printf("amu_W(op,l) = %g  %g\n", amu_Wop_s[Njack - 1], error_jackboot(resampling, Njack, amu_Wop_s));


    int_scheme = integrate_reinman;
    double* amu_Wop_s1 = compute_amu_W(conf_jack, 5 + 12, Njack, ZAs1.P[0], a, q2s, int_scheme, outfile, "amu_{W}(op,s1)", resampling);
    write_jack(amu_Wop_s1, Njack, jack_file);
    check_correlatro_counter(50);
    printf("amu_W(op,l) = %g  %g\n", amu_Wop_s1[Njack - 1], error_jackboot(resampling, Njack, amu_Wop_s1));


    asd_vec[0] = amu_Wop_s;
    asd_vec[1] = amu_Wop_s1;
    amu_W_sphys = interpol_Z(Nstrange, Njack, Meta, asd_vec, jack_aMetas_MeV_exp, outfile, "amu_{W}(op,sphys)", resampling);
    write_jack(amu_W_sphys, Njack, jack_file);
    printf("amu_W(op,shys) = %g  %g\n", amu_W_sphys[Njack - 1], error_jackboot(resampling, Njack, amu_W_sphys));

    check_correlatro_counter(51);


    ////////////////////////////////////////////////////////////////////////////////////////////////////////
    // a_W_s_eq simpson38
    ///////////////////////////////////////////////////////////////////////////////////////////////////////

    int_scheme = integrate_simpson38;
    double* amu_Weq_simp_s = compute_amu_W(conf_jack, 2 + 6, Njack, ZVs.P[0], a, q2s, int_scheme, outfile, "amu_{W}_simpson38(eq,s)", resampling);
    write_jack(amu_Weq_simp_s, Njack, jack_file);
    check_correlatro_counter(52);
    printf("amu_W_simpson38(eq,s) = %g  %g\n", amu_Weq_simp_s[Njack - 1], error_jackboot(resampling, Njack, amu_Weq_simp_s));


    int_scheme = integrate_simpson38;
    double* amu_Weq_simp_s1 = compute_amu_W(conf_jack, 2 + 12, Njack, ZVs1.P[0], a, q2s, int_scheme, outfile, "amu_{W}_simpson38(eq,s1)", resampling);
    write_jack(amu_Weq_simp_s1, Njack, jack_file);
    check_correlatro_counter(53);
    printf("amu_W_simpson38(eq,s1) = %g  %g\n", amu_Weq_simp_s1[Njack - 1], error_jackboot(resampling, Njack, amu_Weq_simp_s1));

    asd_vec[0] = amu_Weq_simp_s;
    asd_vec[1] = amu_Weq_simp_s1;
    amu_W_sphys = interpol_Z(Nstrange, Njack, Meta, asd_vec, jack_aMetas_MeV_exp, outfile, "amu_{W}_simpson38(eq,sphys)", resampling);
    write_jack(amu_W_sphys, Njack, jack_file);
    printf("amu_W(eq,shys) = %g  %g\n", amu_W_sphys[Njack - 1], error_jackboot(resampling, Njack, amu_W_sphys));

    check_correlatro_counter(54);

    ////////////////////////////////////////////////////////////////////////////////////////////////////////
    // a_W_s_op simpson38
    ///////////////////////////////////////////////////////////////////////////////////////////////////////

    int_scheme = integrate_simpson38;
    double* amu_Wop_simp_s = compute_amu_W(conf_jack, 5 + 6, Njack, ZAs.P[0], a, q2s, int_scheme, outfile, "amu_{W}_simpson38(op,s)", resampling);
    write_jack(amu_Wop_simp_s, Njack, jack_file);
    check_correlatro_counter(55);
    printf("amu_W_simpson38(op,l) = %g  %g\n", amu_Wop_simp_s[Njack - 1], error_jackboot(resampling, Njack, amu_Wop_simp_s));


    int_scheme = integrate_simpson38;
    double* amu_Wop_simp_s1 = compute_amu_W(conf_jack, 5 + 12, Njack, ZAs1.P[0], a, q2s, int_scheme, outfile, "amu_{W}_simpson38(op,s1)", resampling);
    write_jack(amu_Wop_simp_s1, Njack, jack_file);
    check_correlatro_counter(56);
    printf("amu_W_simpson38(op,l) = %g  %g\n", amu_Wop_simp_s1[Njack - 1], error_jackboot(resampling, Njack, amu_Wop_simp_s1));


    asd_vec[0] = amu_Wop_simp_s;
    asd_vec[1] = amu_Wop_simp_s1;
    amu_W_sphys = interpol_Z(Nstrange, Njack, Meta, asd_vec, jack_aMetas_MeV_exp, outfile, "amu_{W}_simpson38(op,sphys)", resampling);
    write_jack(amu_W_sphys, Njack, jack_file);
    printf("amu_W_simpson38(op,shys) = %g  %g\n", amu_W_sphys[Njack - 1], error_jackboot(resampling, Njack, amu_W_sphys));
    check_correlatro_counter(57);



    ////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // DV
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////

    double* MpiMev = new double[Njack];
    for (int j = 0;j < Njack;j++) {
        MpiMev[j] = M_PS_op[j] / (a[j] / 197.326963);
    }

    {   // continuum FVE GS 
        double* DV;
        if (strcmp("cA.53.24", argv[4]) == 0 || strcmp("cA.40.24", argv[4]) == 0 || strcmp("cA.30.32", argv[4]) == 0) {
            printf("\n\n                we are not computing FVE \n\n");
            DV = (double*)calloc(Njack, sizeof(double));
            write_jack(DV, Njack, jack_file);
            check_correlatro_counter(58);
        }
        else {
            DV = compute_DVt_and_integrate(header.L, Njack, MpiMev /* jack_Mpi_MeV_exp */, jack_Mrho_MeV_exp, a, jack_grhopipi, outfile, "DVt", resampling);
            write_jack(DV, Njack, jack_file);
            check_correlatro_counter(58);
        }
        printf("DV_{W}(op,l) = %g  %g\n", DV[Njack - 1], error_jackboot(resampling, Njack, DV));
        free(DV);
    }
    //////////////////////////////////////////////////////////////////////////////////////////////
    // SD s 
    ///
    double* tmp;
    asd_vec[0] = amu_sdeq_s;
    asd_vec[1] = amu_sdeq_s1;
    tmp = interpol_Z(Nstrange, Njack, Mphi, asd_vec, jack_aMphi_MeV_exp, outfile, "amu_{sd}(eq,phiphys)", resampling);
    write_jack(tmp, Njack, jack_file);
    printf("amu_sd(eq,phiphys) = %g  %g\n", tmp[Njack - 1], error_jackboot(resampling, Njack, tmp));
    check_correlatro_counter(59); free(tmp);

    asd_vec[0] = amu_sdop_s;
    asd_vec[1] = amu_sdop_s1;
    tmp = interpol_Z(Nstrange, Njack, Mphi, asd_vec, jack_aMphi_MeV_exp, outfile, "amu_{sd}(op,phiphys)", resampling);
    write_jack(tmp, Njack, jack_file);
    printf("amu_sd(op,phiphys) = %g  %g\n", tmp[Njack - 1], error_jackboot(resampling, Njack, tmp));
    check_correlatro_counter(60); free(tmp);

    asd_vec[0] = amu_sdeq_simp_s;
    asd_vec[1] = amu_sdeq_simp_s1;
    tmp = interpol_Z(Nstrange, Njack, Mphi, asd_vec, jack_aMphi_MeV_exp, outfile, "amu_{sd,simp}(eq,phiphys)", resampling);
    write_jack(tmp, Njack, jack_file);
    printf("amu_sd(eq,phiphys) = %g  %g\n", tmp[Njack - 1], error_jackboot(resampling, Njack, tmp));
    check_correlatro_counter(61); free(tmp);

    asd_vec[0] = amu_sdop_simp_s;
    asd_vec[1] = amu_sdop_simp_s1;
    tmp = interpol_Z(Nstrange, Njack, Mphi, asd_vec, jack_aMphi_MeV_exp, outfile, "amu_{sd,simp}(op,phiphys)", resampling);
    write_jack(tmp, Njack, jack_file);
    printf("amu_sd(op,phiphys) = %g  %g\n", tmp[Njack - 1], error_jackboot(resampling, Njack, tmp));
    check_correlatro_counter(62); free(tmp);

    //////////////////////////////////////////////////////////////////////////////////////////////
    // W s
    //////////
    asd_vec[0] = amu_Weq_s;
    asd_vec[1] = amu_Weq_s1;
    amu_W_sphys = interpol_Z(Nstrange, Njack, Mphi, asd_vec, jack_aMphi_MeV_exp, outfile, "amu_{W}(eq,phiphys)", resampling);
    write_jack(amu_W_sphys, Njack, jack_file);
    printf("amu_W(eq,phiphys) = %g  %g\n", amu_W_sphys[Njack - 1], error_jackboot(resampling, Njack, amu_W_sphys));
    check_correlatro_counter(63);

    asd_vec[0] = amu_Wop_s;
    asd_vec[1] = amu_Wop_s1;
    amu_W_sphys = interpol_Z(Nstrange, Njack, Mphi, asd_vec, jack_aMphi_MeV_exp, outfile, "amu_{W}(op,phiphys)", resampling);
    write_jack(amu_W_sphys, Njack, jack_file);
    printf("amu_W(op,phiphys) = %g  %g\n", amu_W_sphys[Njack - 1], error_jackboot(resampling, Njack, amu_W_sphys));
    check_correlatro_counter(64);

    asd_vec[0] = amu_Weq_simp_s;
    asd_vec[1] = amu_Weq_simp_s1;
    amu_W_sphys = interpol_Z(Nstrange, Njack, Mphi, asd_vec, jack_aMphi_MeV_exp, outfile, "amu_{W}_simpson38(eq,phiphys)", resampling);
    write_jack(amu_W_sphys, Njack, jack_file);
    printf("amu_W(eq,phiphys) = %g  %g\n", amu_W_sphys[Njack - 1], error_jackboot(resampling, Njack, amu_W_sphys));
    check_correlatro_counter(65);

    asd_vec[0] = amu_Wop_simp_s;
    asd_vec[1] = amu_Wop_simp_s1;
    amu_W_sphys = interpol_Z(Nstrange, Njack, Mphi, asd_vec, jack_aMphi_MeV_exp, outfile, "amu_{W}_simpson38(op,phiphys)", resampling);
    write_jack(amu_W_sphys, Njack, jack_file);
    printf("amu_W_simpson38(op,phiphys) = %g  %g\n", amu_W_sphys[Njack - 1], error_jackboot(resampling, Njack, amu_W_sphys));
    check_correlatro_counter(66);

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////charms
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    int Ncharm = 3;
    char name_corr[NAMESIZE];
    std::vector<std::string> name_intgr = { "rein", "simps" };
    std::vector<std::string> name_eqop = { "eq", "op" };
    std::vector<int> id_eqop = { var + 3 + 3 * 2,
                                 var + 4 + 3 * 2,
                                 var + 3 + 4 * 2,
                                 var + 4 + 4 * 2,
                                 var + 3 + 5 * 2,
                                 var + 4 + 5 * 2 };
    std::vector<std::string> name_M = { "Metac", "MJPsi" };
    int Ncharm_inter = 3;
    if (strcmp("cD.54.96", argv[4]) == 0) Ncharm_inter = 2;

    double** Metac_vec = (double**)malloc(sizeof(double*) * Ncharm);
    for (int ic = 0;ic < Ncharm;ic++) {
        mysprintf(name_corr, NAMESIZE, "M_{etac}^{op}(c%d)", ic);
        Metac_vec[ic] = plateau_correlator_function(option, kinematic_2pt, (char*)"P5P5", conf_jack, Njack, namefile_plateaux, outfile, 4 + 6 * (3 + ic), name_corr, M_eff_T, jack_file);
        check_correlatro_counter(67 + ic);
        // printf("%s = %g  %g\n", name_corr, Jpsi_vec[ic][Njack - 1], error_jackboot(resampling, Njack, Jpsi_vec[ic]));
    }
    double** Jpsi_vec = (double**)malloc(sizeof(double*) * Ncharm);
    for (int ic = 0;ic < Ncharm;ic++) {
        mysprintf(name_corr, NAMESIZE, "M_{JPsi}^{op}(c%d)", ic);
        Jpsi_vec[ic] = plateau_correlator_function(option, kinematic_2pt, (char*)"P5P5", conf_jack, Njack, namefile_plateaux, outfile, 5 + 6 * (3 + ic), name_corr, M_eff_T, jack_file);
        check_correlatro_counter(70 + ic);
        // printf("%s = %g  %g\n", name_corr, Jpsi_vec[ic][Njack - 1], error_jackboot(resampling, Njack, Jpsi_vec[ic]));
    }

    constexpr double q2c = 4.0 / 9.0;
    if (strcmp("cB.72.96", argv[4]) != 0) {
        for (int tm = 0; tm < name_eqop.size(); tm++) {// eq op
            double* Z;
            if (tm == 0) Z = ZV;
            if (tm == 1) Z = ZA;
            for (int intgr = 0; intgr < name_intgr.size(); intgr++) {// integrators
                if (intgr == 0)      int_scheme = integrate_reinman;
                else if (intgr == 1) int_scheme = integrate_simpson38;
                else exit(-1);

                double** asdc_vec = (double**)malloc(sizeof(double*) * Ncharm);
                for (int ic = 0;ic < Ncharm;ic++) {
                    mysprintf(name_corr, NAMESIZE, "amu_{sd,%s}(%s,c%d)", name_intgr[intgr].c_str(), name_eqop[tm].c_str(), ic);
                    isub = (strcmp(argv[argc - 1], "three_corr") == 0) ? id_eqop[tm + ic * 2] : -1;
                    asdc_vec[ic] = compute_amu_sd(conf_jack, 2 + 6 * (3 + ic) + 3 * tm, Njack, Z, a, q2c, int_scheme, outfile, name_corr, resampling, isub);
                    write_jack(asdc_vec[ic], Njack, jack_file);
                    check_correlatro_counter(73 + ic + intgr * (Ncharm + name_M.size()) + tm * (Ncharm + name_M.size()) * name_intgr.size());
                    printf("%s = %g  %g\n", name_corr, asdc_vec[ic][Njack - 1], error_jackboot(resampling, Njack, asdc_vec[ic]));
                }
                for (int im = 0; im < name_M.size(); im++) {// etac JPsi
                    mysprintf(name_corr, NAMESIZE, "amu_sd_%s(%s,%s)", name_intgr[intgr].c_str(), name_eqop[tm].c_str(), name_M[im].c_str());
                    double** Meson;
                    double* Meson_phys;
                    if (im == 0) { Meson = Metac_vec; Meson_phys = jack_aMetac_MeV_exp; }
                    else if (im == 1) { Meson = Jpsi_vec; Meson_phys = jack_aJpsi_MeV_exp; }

                    double* amu_sdc_sphys = interpol_Z(Ncharm_inter, Njack, Meson, asdc_vec, Meson_phys, outfile, name_corr, resampling);
                    write_jack(amu_sdc_sphys, Njack, jack_file);
                    size_t counter_id = 76 + im + intgr * (Ncharm + name_M.size()) + tm * (Ncharm + name_M.size()) * name_intgr.size();
                    check_correlatro_counter(counter_id);
                    printf("%s  id=%ld\n", name_corr, counter_id);
                    free(amu_sdc_sphys);
                }
                free_2(Ncharm, asdc_vec);
            }
        }
    }
    else if (strcmp("cB.72.96", argv[4]) == 0) {
        for (int tm = 0; tm < name_eqop.size(); tm++) {
            for (int intgr = 0; intgr < name_intgr.size(); intgr++) {
                for (int ic = 0;ic < Ncharm;ic++) {
                    zero_corr(zeros, Njack, jack_file);
                }
                for (int im = 0; im < name_M.size(); im++) {
                    zero_corr(zeros, Njack, jack_file);
                }
            }
        }
    }

    ///////////// W 


    if (strcmp("cB.72.96", argv[4]) != 0) {
        for (int tm = 0; tm < name_eqop.size(); tm++) {// eq op
            double* Z;
            if (tm == 0) Z = ZV;
            if (tm == 1) Z = ZA;
            for (int intgr = 0; intgr < name_intgr.size(); intgr++) {// integrators
                if (intgr == 0)      int_scheme = integrate_reinman;
                else if (intgr == 1) int_scheme = integrate_simpson38;
                else exit(-1);

                double** asdc_vec = (double**)malloc(sizeof(double*) * Ncharm);
                for (int ic = 0;ic < Ncharm;ic++) {
                    mysprintf(name_corr, NAMESIZE, "amu_{W,%s}(%s,c%d)", name_intgr[intgr].c_str(), name_eqop[tm].c_str(), ic);
                    asdc_vec[ic] = compute_amu_W(conf_jack, 2 + 6 * (3 + ic) + 3 * tm, Njack, Z, a, q2c, int_scheme, outfile, name_corr, resampling);
                    write_jack(asdc_vec[ic], Njack, jack_file);
                    check_correlatro_counter(93 + ic + intgr * (Ncharm + name_M.size()) + tm * (Ncharm + name_M.size()) * name_intgr.size());
                    printf("%s = %g  %g\n", name_corr, asdc_vec[ic][Njack - 1], error_jackboot(resampling, Njack, asdc_vec[ic]));
                }
                for (int im = 0; im < name_M.size(); im++) {// etac JPsi
                    mysprintf(name_corr, NAMESIZE, "amu_W_%s(%s,%s)", name_intgr[intgr].c_str(), name_eqop[tm].c_str(), name_M[im].c_str());
                    double** Meson;
                    double* Meson_phys;
                    if (im == 0) { Meson = Metac_vec; Meson_phys = jack_aMetac_MeV_exp; }
                    else if (im == 1) { Meson = Jpsi_vec; Meson_phys = jack_aJpsi_MeV_exp; }

                    double* amu_sdc_sphys = interpol_Z(Ncharm_inter, Njack, Meson, asdc_vec, Meson_phys, outfile, name_corr, resampling);
                    write_jack(amu_sdc_sphys, Njack, jack_file);
                    size_t counter_id = 96 + im + intgr * (Ncharm + name_M.size()) + tm * (Ncharm + name_M.size()) * name_intgr.size();
                    check_correlatro_counter(counter_id);
                    printf("%s  id=%ld\n", name_corr, counter_id);
                    free(amu_sdc_sphys);
                }
                free_2(Ncharm, asdc_vec);
            }
        }
    }
    else if (strcmp("cB.72.96", argv[4]) == 0) {
        for (int tm = 0; tm < name_eqop.size(); tm++) {
            for (int intgr = 0; intgr < name_intgr.size(); intgr++) {
                for (int ic = 0;ic < Ncharm;ic++) {
                    zero_corr(zeros, Njack, jack_file);
                }
                for (int im = 0; im < name_M.size(); im++) {
                    zero_corr(zeros, Njack, jack_file);
                }
            }
        }
    }



    ////////////////////////////////masses
    double** ms = (double**)malloc(sizeof(double*) * 2);
    ms[0] = fake_sampling(resampling, header.mus[1], 1e-10, Njack, 1);
    ms[1] = fake_sampling(resampling, header.mus[2], 1e-10, Njack, 1);
    double* ms_etas = interpol_Z(Nstrange, Njack, Meta, ms, jack_aMetas_MeV_exp, outfile, "ms(etas)", resampling);
    free(ms_etas);

    double* ms_phi = interpol_Z(Nstrange, Njack, Mphi, ms, jack_aMphi_MeV_exp, outfile, "ms(phi)", resampling);
    free(ms_phi);



    double** mc = (double**)malloc(sizeof(double*) * 3);
    mc[0] = fake_sampling(resampling, header.mus[3], 1e-10, Njack, 1);
    mc[1] = fake_sampling(resampling, header.mus[4], 1e-10, Njack, 1);
    mc[2] = fake_sampling(resampling, header.mus[5], 1e-10, Njack, 1);


    double* mc_etac = interpol_Z(Ncharm_inter, Njack, Metac_vec, mc, jack_aMetac_MeV_exp, outfile, "mc(etac)", resampling);
    free(mc_etac);

    double* mc_JPsi = interpol_Z(Ncharm_inter, Njack, Jpsi_vec, mc, jack_aJpsi_MeV_exp, outfile, "mc(Jpsi)", resampling);
    free(mc_JPsi);
    free_2(3, mc);

    ///////////////////////////////////// correlator
    {
        double t1 = 0.2;
        fit_type fit_info;


        fit_info.Nvar = 1;
        fit_info.N = 1;
        fit_info.Njack = Njack;
        fit_info.n_ext_P = 3;
        fit_info.ext_P = double_malloc_2(fit_info.n_ext_P, Njack);

        for (int j = 0;j < Njack; j++) {
            fit_info.ext_P[0][j] = t1;
            fit_info.ext_P[1][j] = ZV[j];
            fit_info.ext_P[2][j] = a[j];
        }
        fit_info.repeat_start = 10;
        fit_info.verbosity = 0;
        ///////////////// 4 points fits
        fit_info.codeplateaux = true;
        fit_info.tmin = ((int)(t1 / a[Njack - 1])) - 1;
        fit_info.tmax = fit_info.tmin + 3;

        fit_info.Npar = 4;
        fit_info.function = rhs_2exp;

        fit_result ct_2exp = fit_fun_to_fun_of_corr(option, kinematic_2pt, (char*)"A1P1phi", conf_jack, namefile_plateaux, outfile, lhs_ct<2>, "C_{l}^{eq}_2exp", fit_info, jack_file);
        // write_jack(ct_2exp.P[0], Njack, jack_file);
        check_correlatro_counter(113);

        free_fit_result(fit_info, ct_2exp);

        fit_info.Npar = 4;
        fit_info.function = rhs_poly;

        fit_result ct_poly3 = fit_fun_to_fun_of_corr(option, kinematic_2pt, (char*)"A1P1phi", conf_jack, namefile_plateaux, outfile, lhs_ct<2>, "C_{l}^{eq}_poly3", fit_info, jack_file);
        // write_jack(ct_poly3.P[0], Njack, jack_file);
        check_correlatro_counter(114);
        free_fit_result(fit_info, ct_poly3);

        fit_info.Npar = 3;
        fit_info.function = rhs_poly;
        fit_info.repeat_start = 1000;
        fit_info.acc = 1e-4;
        fit_info.precision_sum = 2;

        fit_result ct_poly2 = fit_fun_to_fun_of_corr(option, kinematic_2pt, (char*)"A1P1phi", conf_jack, namefile_plateaux, outfile, lhs_ct<2>, "C_{l}^{eq}_poly2", fit_info, jack_file);
        // write_jack(ct_poly2.P[0], Njack, jack_file);
        check_correlatro_counter(115);
        free_fit_result(fit_info, ct_poly2);
        //////////////////////////////////////////////////////////////////////////////////////////////////////
        /////////////////////////////////////////////////// 3 points fits
        fit_info.tmax = fit_info.tmin + 2;

        fit_info.Npar = 2;
        fit_info.function = rhs_1exp;
        fit_info.repeat_start = 1;
        fit_result ct_1exp = fit_fun_to_fun_of_corr(option, kinematic_2pt, (char*)"A1P1phi", conf_jack, namefile_plateaux, outfile, lhs_ct<2>, "C_{l}^{eq}_1exp", fit_info, jack_file);
        // write_jack(ct_1exp.P[0], Njack, jack_file);
        check_correlatro_counter(116);
        free_fit_result(fit_info, ct_1exp);



        fit_info.Npar = 3;
        fit_info.function = rhs_poly;

        fit_result ct_poly2_3pt = fit_fun_to_fun_of_corr(option, kinematic_2pt, (char*)"A1P1phi", conf_jack, namefile_plateaux, outfile, lhs_ct<2>, "C_{l}^{eq}_poly2_3pt", fit_info, jack_file);
        // write_jack(ct_poly2_3pt.P[0], Njack, jack_file);
        check_correlatro_counter(117);
        free_fit_result(fit_info, ct_poly2_3pt);

        ///// end fits
        for (int i = 0;i < fit_info.n_ext_P; i++)   free(fit_info.ext_P[i]);

        fit_info.restore_default();
    }
    {
        double t1 = 0.2;
        fit_type fit_info;

        fit_info.Nvar = 1;
        fit_info.N = 1;
        fit_info.Njack = Njack;
        fit_info.n_ext_P = 3;
        fit_info.ext_P = double_malloc_2(fit_info.n_ext_P, Njack);

        for (int j = 0;j < Njack; j++) {
            fit_info.ext_P[0][j] = t1;
            fit_info.ext_P[1][j] = ZA[j];
            fit_info.ext_P[2][j] = a[j];
        }
        fit_info.repeat_start = 10;
        fit_info.verbosity = 0;
        ///////////////// 4 points fits
        fit_info.codeplateaux = true;
        fit_info.tmin = ((int)(t1 / a[Njack - 1])) - 1;
        fit_info.tmax = fit_info.tmin + 3;

        fit_info.Npar = 4;
        fit_info.function = rhs_2exp;

        fit_result ct_2exp = fit_fun_to_fun_of_corr(option, kinematic_2pt, (char*)"A1P1phi", conf_jack, namefile_plateaux, outfile, lhs_ct<5>, "C_{l}^{op}_2exp", fit_info, jack_file);
        // write_jack(ct_2exp.P[0], Njack, jack_file);
        check_correlatro_counter(118);

        free_fit_result(fit_info, ct_2exp);
    }
    //////////////////////////////
    /////////////// Strange
    {
        double t1 = 0.2;
        fit_type fit_info;

        fit_info.Nvar = 1;
        fit_info.N = 1;
        fit_info.Njack = Njack;
        fit_info.n_ext_P = 3;
        fit_info.ext_P = double_malloc_2(fit_info.n_ext_P, Njack);

        for (int j = 0;j < Njack; j++) {
            fit_info.ext_P[0][j] = t1;
            fit_info.ext_P[1][j] = ZV[j];
            fit_info.ext_P[2][j] = a[j];
        }
        fit_info.repeat_start = 10;
        fit_info.verbosity = 0;
        ///////////////// 4 points fits
        fit_info.codeplateaux = true;
        fit_info.tmin = ((int)(t1 / a[Njack - 1])) - 1;
        fit_info.tmax = fit_info.tmin + 3;

        fit_info.Npar = 4;
        fit_info.function = rhs_2exp;

        fit_result ct_2exp = fit_fun_to_fun_of_corr(option, kinematic_2pt, (char*)"A1P1phi", conf_jack, namefile_plateaux, outfile,
            lhs_ct<8>, "C_{s1}^{eq}_2exp", fit_info, jack_file);
        check_correlatro_counter(119);


        fit_result ct_2exp1 = fit_fun_to_fun_of_corr(option, kinematic_2pt, (char*)"A1P1phi", conf_jack, namefile_plateaux, outfile,
            lhs_ct<14>, "C_{s2}^{eq}_2exp", fit_info, jack_file);
        check_correlatro_counter(120);
        printf("%.12g    %.12g\n", conf_jack[0][8][0][0], conf_jack[0][14][0][0]);
        double** cts = (double**)malloc(sizeof(double*) * Nstrange);
        cts[0] = ct_2exp.P[0];
        cts[1] = ct_2exp1.P[0];


        double* ctst0 = interpol_Z(Nstrange, Njack, Meta, cts, jack_aMetas_MeV_exp, outfile, "C_{s,Meta}^{eq}_2exp", resampling);
        write_jack(ctst0, Njack, jack_file);
        printf("amu_sd(eq,shys) = %g  %g\n", ctst0[Njack - 1], myres->comp_error(ctst0));

        check_correlatro_counter(121);
        free(ctst0);
        free_fit_result(fit_info, ct_2exp);
        free_fit_result(fit_info, ct_2exp1);
    }

    int  id_P5P5_mudmu;
    int  id_P5P5_cor_mudmu;
    int  id_VKVKeq_mudm;
    int  id_VKVKop_mudm;
    int  id_VKVKeq_SD_mudm;
    int  id_VKVKop_SD_mudm;
    int  id_sea_VKVKeq;
    int  id_sea_VKVKop;

    if (argc >= 19) {
        fit_type fit_info;
        fit_info.N = 1;
        fit_info.Njack = Njack;
        printf("Ncorrelator=%d\n", ncorr_new);
        fit_info.corr_id = { 53,40,var/*55*/ };//P5P5_corr_bolla, P5P5dmu , bolla*P5P5dmu
        fit_info.guess = { mu,mul1 };
        id_P5P5_mudmu = ncorr_new;
        add_correlators(option, ncorr_new, conf_jack, corr_plus_dm, fit_info);

        fit_info.corr_id = { 2,2,var + 1 };//VKVKeq, VKVKeqdmu , bolla*VKVKeqdmu
        fit_info.guess = { mu,mul1 };
        add_correlators(option, ncorr_new, conf_jack, corr_plus_dm, fit_info);

        fit_info.corr_id = { 5,5,var + 2 };//VKVKop, VKVKopdmu , bolla*VKVKopdmu
        fit_info.guess = { mu,mul1 };
        add_correlators(option, ncorr_new, conf_jack, corr_plus_dm, fit_info);

        fit_info.corr_id = { 4,47, 40, 53,var/*55*/ };//P5P5, P5P5_corr_dmu,  P5P5_mu+dm ,P5P5_corr_bolla, bolla*P5P5_cor_bolla
        fit_info.guess = { mu,mul1 };
        id_P5P5_cor_mudmu = ncorr_new;
        add_correlators(option, ncorr_new, conf_jack, corr_plus_dm_correlated, fit_info);

        fit_info.corr_id = { 2,45, 38, 51,var + 1/*56*/ };//VKVK, VKVK_corr_dmu,  VKVK_mu+dm ,VKVK_corr_bolla, bolla*VKVK_cor_bolla
        fit_info.guess = { mu,mul1 };
        id_VKVKeq_mudm = ncorr_new;
        add_correlators(option, ncorr_new, conf_jack, corr_plus_dm_correlated, fit_info);

        fit_info.corr_id = { 5,48, 41, 54,var + 2/*57*/ };//VKVK, VKVK_corr_dmu,  VKVK_mu+dm ,VKVK_corr_bolla, bolla*VKVK_cor_bolla //op
        fit_info.guess = { mu,mul1 };
        id_VKVKop_mudm = ncorr_new;
        add_correlators(option, ncorr_new, conf_jack, corr_plus_dm_correlated, fit_info);


        fit_info.corr_id = { 2, 45, 38, 51, var + 1/*56*/ };//VKVK, VKVK_corr_dmu,  VKVK_mu+dm ,VKVK_corr_bolla, bolla*VKVK_cor_bolla
        fit_info.guess = { mu, mul1 };
        id_sea_VKVKeq = ncorr_new;
        add_correlators(option, ncorr_new, conf_jack, mu_sea_correction, fit_info);

        fit_info.corr_id = { 5, 48, 41, 54, var + 2/*57*/ };//VKVK, VKVK_corr_dmu,  VKVK_mu+dm ,VKVK_corr_bolla, bolla*VKVK_cor_bolla //op
        fit_info.guess = { mu, mul1 };
        id_sea_VKVKop = ncorr_new;
        add_correlators(option, ncorr_new, conf_jack, mu_sea_correction, fit_info);

        printf("Ncorrelator=%d\n", ncorr_new);

        ///////////////////////////////////
        double* M_PS1 = plateau_correlator_function(option, kinematic_2pt, (char*)"P5P5", conf_jack, Njack, namefile_plateaux, outfile, id_P5P5_mudmu, "M_{PS1}^{op}", M_eff_T, jack_file);
        check_correlatro_counter(122);
        free(M_PS1);

        M_PS1 = plateau_correlator_function(option, kinematic_2pt, (char*)"P5P5", conf_jack, Njack, namefile_plateaux, outfile, id_P5P5_cor_mudmu, "M_{PS1c}^{op}", M_eff_T, jack_file);
        check_correlatro_counter(123);


        double* jack_Mpia_exp = new double[Njack];


        int Nmu = 2;
        double** Masses_pions = (double**)malloc(sizeof(double*) * Nstrange);
        Masses_pions[0] = M_PS_op;
        Masses_pions[1] = M_PS1;
        double** mul_values = double_malloc_2(2, Njack);
        for (int j = 0; j < Njack;j++) {
            jack_Mpia_exp[j] = jack_Mpi_MeV_exp[j] * (a[j] / 197.326963);
            mul_values[0][j] = mu;
            mul_values[1][j] = mul1;
        }
        double* mu_phys = interpol_Z(Nmu, Njack, Masses_pions, mul_values, jack_Mpia_exp, outfile, "mu_phys", resampling);
        printf("mu_phys = %g  %g\n", mu_phys[Njack - 1], myres->comp_error(mu_phys));
        // free(M_PS1);

        //////////////////////////// valence 
        int_scheme = integrate_reinman;
        double* amu_W_eq2 = compute_amu_W(conf_jack, 38, Njack, ZV, a, 5.0 / 9.0, int_scheme, outfile, "valence1_amu_{W}(eq,l)", resampling);
        write_jack(amu_W_eq2, Njack, jack_file);
        check_correlatro_counter(124);
        printf("valence_amu_W(eq2,l) = %g  %g\n", amu_W_eq2[Njack - 1], error_jackboot(resampling, Njack, amu_W_eq2));

        double* amu_W_eq01 = compute_amu_W(conf_jack, 45, Njack, ZV, a, 5.0 / 9.0, int_scheme, outfile, "valence3_amu_{W}(eq,l)", resampling);
        write_jack(amu_W_eq01, Njack, jack_file);
        check_correlatro_counter(125);
        printf("valence_amu_W(eq01,l) = %g  %g\n", amu_W_eq01[Njack - 1], error_jackboot(resampling, Njack, amu_W_eq01));


        double* diff_valence = (double*)malloc(sizeof(double) * Njack);
        for (int j = 0; j < Njack; j++) {
            diff_valence[j] = (amu_W_eq2[j] - amu_W_eq01[j]);
        }
        printf("diff_valence_amu_W(eq,l) = %g  %g\n", diff_valence[Njack - 1], myres->comp_error(diff_valence));


        double* amu_W_op2 = compute_amu_W(conf_jack, 41, Njack, ZA, a, 5.0 / 9.0, int_scheme, outfile, "valence_amu_{W}(op,l)", resampling);
        write_jack(amu_W_op2, Njack, jack_file);
        check_correlatro_counter(126);
        printf("valence_amu_W(op2,l) = %g  %g\n", amu_W_op2[Njack - 1], error_jackboot(resampling, Njack, amu_W_op2));

        double* amu_W_op02 = compute_amu_W(conf_jack, 48, Njack, ZA, a, 5.0 / 9.0, int_scheme, outfile, "valence_amu_{W}(op,l)", resampling);
        write_jack(amu_W_op02, Njack, jack_file);
        check_correlatro_counter(127);
        printf("valence_amu_W(op01,l) = %g  %g\n", amu_W_op02[Njack - 1], error_jackboot(resampling, Njack, amu_W_op02));
        for (int j = 0; j < Njack; j++) {
            diff_valence[j] = (amu_W_op2[j] - amu_W_op02[j]);
        }
        printf("diff_valence_amu_W(op,l) = %g  %g\n", diff_valence[Njack - 1], myres->comp_error(diff_valence));

        //////////////////////////// sea
        double* sea_eq = compute_amu_W(conf_jack, id_sea_VKVKeq, Njack, ZV, a, 5.0 / 9.0, int_scheme, outfile, "sea_diff_amu_{W}(eq,l)", resampling);
        write_jack(sea_eq, Njack, jack_file);
        printf("sea_diff_amu_{W}(eq,l) = %g  %g\n", sea_eq[Njack - 1], myres->comp_error(sea_eq));
        check_correlatro_counter(128);
        free(sea_eq);


        double* sea_op = compute_amu_W(conf_jack, id_sea_VKVKop, Njack, ZA, a, 5.0 / 9.0, int_scheme, outfile, "sea_diff_amu_{W}(op,l)", resampling);
        write_jack(sea_op, Njack, jack_file);
        printf("sea_diff_amu_{W}(op,l) = %g  %g\n", sea_op[Njack - 1], myres->comp_error(sea_op));
        check_correlatro_counter(129);
        free(sea_op);

        /// total eq
        double** amu_W_eq_mu = (double**)malloc(sizeof(double*) * 2);
        amu_W_eq_mu[0] = (double*)malloc(sizeof(double) * Njack);
        for (int j = 0;j < Njack;j++) {
            amu_W_eq_mu[0][j] = amu_W_eq0[j];
        }

        amu_W_eq_mu[1] = compute_amu_W(conf_jack, id_VKVKeq_mudm, Njack, ZV, a, 5.0 / 9.0, int_scheme, outfile, "amu_{W}(eq,l1)", resampling);
        write_jack(amu_W_eq_mu[1], Njack, jack_file);
        printf("amu_{W}(eq,l1) = %g  %g\n", amu_W_eq_mu[1][Njack - 1], myres->comp_error(amu_W_eq_mu[1]));
        check_correlatro_counter(130);

        double* amu_W_eq_phys = interpol_Z(Nmu, Njack, Masses_pions, amu_W_eq_mu, jack_Mpia_exp, outfile, "amu_{W}(eq,phys)", resampling);
        write_jack(amu_W_eq_phys, Njack, jack_file);
        check_correlatro_counter(131);

        free_2(2, amu_W_eq_mu);


        /// total op
        double** amu_W_op_mu = (double**)malloc(sizeof(double*) * 2);
        amu_W_op_mu[0] = (double*)malloc(sizeof(double) * Njack);
        for (int j = 0;j < Njack;j++) {
            amu_W_op_mu[0][j] = amu_W_op0[j];
        }
        amu_W_op_mu[1] = compute_amu_W(conf_jack, id_VKVKop_mudm, Njack, ZA, a, 5.0 / 9.0, int_scheme, outfile, "amu_{W}(op,l1)", resampling);
        printf("amu_{W}(op,l1) = %g  %g\n", amu_W_op_mu[1][Njack - 1], myres->comp_error(amu_W_op_mu[1]));
        write_jack(amu_W_op_mu[1], Njack, jack_file);
        check_correlatro_counter(132);
        double* amu_W_op_phys = interpol_Z(Nmu, Njack, Masses_pions, amu_W_op_mu, jack_Mpia_exp, outfile, "amu_{W}(op,phys)", resampling);
        write_jack(amu_W_op_phys, Njack, jack_file);
        check_correlatro_counter(133);

        free_2(2, amu_W_op_mu);

        double* tmp_b = (double*)malloc(sizeof(double) * Njack);
        for (int j = 0;j < Njack;j++) {
            tmp_b[j] = conf_jack[j][42][0][0];
        }
        printf("<B> = %.12g  %g\n", tmp_b[Njack - 1], myres->comp_error(tmp_b));


        for (int j = 0;j < Njack;j++) {
            tmp_b[j] = conf_jack[j][var + 1][0][0];
        }
        printf("<B*V_eq(t=0)> = %.12g  %g\n", tmp_b[Njack - 1], myres->comp_error(tmp_b));

        for (int j = 0;j < Njack;j++) {
            tmp_b[j] = conf_jack[j][var + 2][0][0];
        }
        printf("<B*V_op(t=0)> = %.12g  %g\n", tmp_b[Njack - 1], myres->comp_error(tmp_b));

        for (int j = 0;j < Njack;j++) {
            tmp_b[j] = conf_jack[j][51][0][0];
        }
        printf("<V_eq(t=0)> = %.12g  %g\n", tmp_b[Njack - 1], myres->comp_error(tmp_b));
        for (int j = 0;j < Njack;j++) {
            tmp_b[j] = conf_jack[j][54][0][0];
        }
        printf("<V_op(t=0)> = %.12g  %g\n", tmp_b[Njack - 1], myres->comp_error(tmp_b));

        for (int j = 0;j < Njack;j++) {
            tmp_b[j] = conf_jack[j][var + 2][0][0] - conf_jack[j][42][0][0] * conf_jack[j][54][0][0];
        }
        printf("<B*V_op(t=0)> -<B><V_op(t=0)> = %.12g  %g\n", tmp_b[Njack - 1], myres->comp_error(tmp_b));


        // double* trash = plateau_correlator_function(option, kinematic_2pt, (char*)"P5P5", conf_jack, Njack, namefile_plateaux, outfile, 5, "trash", identity, jack_file);
        // check_correlatro_counter(122);
        isub = (strcmp(argv[argc - 1], "three_corr") == 0) ? var + 3 : -1;
        double* amu_SD_eq_l1 = compute_amu_sd(conf_jack, id_VKVKeq_mudm, Njack, ZV, a, 5.0 / 9.0, int_scheme, outfile, "amu_{SD}(eq,l1)", resampling, isub);
        write_jack(amu_SD_eq_l1, Njack, jack_file);
        printf("amu_{SD}(eq,l1) = %g  %g\n", amu_SD_eq_l1[Njack - 1], myres->comp_error(amu_SD_eq_l1));
        check_correlatro_counter(134);

        isub = (strcmp(argv[argc - 1], "three_corr") == 0) ? var + 4 : -1;
        double* amu_SD_op_l1 = compute_amu_sd(conf_jack, id_VKVKop_mudm, Njack, ZA, a, 5.0 / 9.0, int_scheme, outfile, "amu_{SD}(op,l1)", resampling, isub);
        printf("amu_{SD}(op,l1) = %g  %g\n", amu_SD_op_l1[Njack - 1], myres->comp_error(amu_SD_op_l1));
        write_jack(amu_SD_op_l1, Njack, jack_file);
        check_correlatro_counter(135);

    }
    else { for (int i = 122; i <= 135;i++)  zero_corr(zeros, Njack, jack_file); }

    /////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // amu_full s
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////

    ////////////////////////////////////////////////////////////////////////////////////////////////////////
    // a_full_s_eq simpson38
    ///////////////////////////////////////////////////////////////////////////////////////////////////////
    double** afull_vec = (double**)malloc(sizeof(double*) * Nstrange);
    int_scheme = integrate_simpson38;
    isub = (strcmp(argv[argc - 1], "three_corr") == 0) ? var + 3 + 1 * 2 : -1;
    double* amu_fulleq_simp_s = compute_amu_full(conf_jack, 2 + 6, Njack, ZVs.P[0], a, q2s, int_scheme, outfile, "amu_{full}_simpson38(eq,s)", resampling, isub);
    write_jack(amu_fulleq_simp_s, Njack, jack_file);
    check_correlatro_counter(136);
    printf("amu_full_simpson38(eq,s) = %g  %g\n", amu_fulleq_simp_s[Njack - 1], error_jackboot(resampling, Njack, amu_fulleq_simp_s));


    int_scheme = integrate_simpson38;
    isub = (strcmp(argv[argc - 1], "three_corr") == 0) ? var + 3 + 2 * 2 : -1;
    double* amu_fulleq_simp_s1 = compute_amu_full(conf_jack, 2 + 12, Njack, ZVs1.P[0], a, q2s, int_scheme, outfile, "amu_{full}_simpson38(eq,s1)", resampling, isub);
    write_jack(amu_fulleq_simp_s1, Njack, jack_file);
    check_correlatro_counter(137);
    printf("amu_full_simpson38(eq,s1) = %g  %g\n", amu_fulleq_simp_s1[Njack - 1], error_jackboot(resampling, Njack, amu_fulleq_simp_s1));

    afull_vec[0] = amu_fulleq_simp_s;
    afull_vec[1] = amu_fulleq_simp_s1;

    double* amu_full_sphys = interpol_Z(Nstrange, Njack, Meta, afull_vec, jack_aMetas_MeV_exp, outfile, "amu_{full}_simpson38(eq,etaphys)", resampling);
    write_jack(amu_full_sphys, Njack, jack_file);
    free(amu_full_sphys);
    check_correlatro_counter(138);

    tmp = interpol_Z(Nstrange, Njack, Mphi, afull_vec, jack_aMphi_MeV_exp, outfile, "amu_{full,simp}(eq,phiphys)", resampling);
    write_jack(tmp, Njack, jack_file);
    check_correlatro_counter(139); free(tmp);


    ////////////////////////////////////////////////////////////////////////////////////////////////////////
    // a_full_s_op simpson38
    ///////////////////////////////////////////////////////////////////////////////////////////////////////

    int_scheme = integrate_simpson38;
    isub = (strcmp(argv[argc - 1], "three_corr") == 0) ? var + 4 + 1 * 2 : -1;
    double* amu_fullop_simp_s = compute_amu_full(conf_jack, 5 + 6, Njack, ZAs.P[0], a, q2s, int_scheme, outfile, "amu_{full}_simpson38(op,s)", resampling, isub);
    write_jack(amu_fullop_simp_s, Njack, jack_file);
    check_correlatro_counter(140);
    printf("amu_full_simpson38(op,s) = %g  %g\n", amu_fullop_simp_s[Njack - 1], error_jackboot(resampling, Njack, amu_fullop_simp_s));


    int_scheme = integrate_simpson38;
    isub = (strcmp(argv[argc - 1], "three_corr") == 0) ? var + 4 + 2 * 2 : -1;
    double* amu_fullop_simp_s1 = compute_amu_full(conf_jack, 5 + 12, Njack, ZAs1.P[0], a, q2s, int_scheme, outfile, "amu_{full}_simpson38(op,s1)", resampling, isub);
    write_jack(amu_fullop_simp_s1, Njack, jack_file);
    check_correlatro_counter(141);
    printf("amu_full_simpson38(op,s1) = %g  %g\n", amu_fullop_simp_s1[Njack - 1], error_jackboot(resampling, Njack, amu_fullop_simp_s1));


    afull_vec[0] = amu_fullop_simp_s;
    afull_vec[1] = amu_fullop_simp_s1;
    amu_full_sphys = interpol_Z(Nstrange, Njack, Meta, afull_vec, jack_aMetas_MeV_exp, outfile, "amu_{full}_simpson38(op,etaphys)", resampling);
    write_jack(amu_full_sphys, Njack, jack_file);
    check_correlatro_counter(142);
    free(amu_full_sphys);


    tmp = interpol_Z(Nstrange, Njack, Mphi, afull_vec, jack_aMphi_MeV_exp, outfile, "amu_{full,simp}(op,phiphys)", resampling);
    write_jack(tmp, Njack, jack_file);
    check_correlatro_counter(143); free(tmp);
    /////////////////
    /// K meson
    /////////
    double* M_K1_op = plateau_correlator_function(option, kinematic_2pt, (char*)"P5P5", conf_jack, Njack, namefile_plateaux, outfile, idD, "M_{K1}^{op}", M_eff_T, jack_file);
    check_correlatro_counter(144);

    double* M_K2_op = plateau_correlator_function(option, kinematic_2pt, (char*)"P5P5", conf_jack, Njack, namefile_plateaux, outfile, idD + 1, "M_{K2}^{op}", M_eff_T, jack_file);
    check_correlatro_counter(145);

    double** MK = (double**)malloc(sizeof(double*) * Nstrange);
    MK[0] = M_K1_op;
    MK[1] = M_K2_op;

    //////////////////// fit at the MK

    afull_vec[0] = amu_fulleq_simp_s;
    afull_vec[1] = amu_fulleq_simp_s1;
    tmp = interpol_Z(Nstrange, Njack, MK, afull_vec, jack_aMK_MeV_exp, outfile, "amu_{full,simp}(eq,MK)", resampling);
    write_jack(tmp, Njack, jack_file);
    check_correlatro_counter(146); free(tmp);

    afull_vec[0] = amu_fullop_simp_s;
    afull_vec[1] = amu_fullop_simp_s1;
    tmp = interpol_Z(Nstrange, Njack, MK, afull_vec, jack_aMK_MeV_exp, outfile, "amu_{full,simp}(op,MK)", resampling);
    write_jack(tmp, Njack, jack_file);
    check_correlatro_counter(147); free(tmp);

    ///  ms from MK
    double* ms_MK = interpol_Z(Nstrange, Njack, MK, ms, jack_aMK_MeV_exp, outfile, "ms(MK)", resampling);
    write_jack(tmp, Njack, jack_file);
    check_correlatro_counter(148);
    free(ms_MK);


    free(amu_fullop_simp_s);free(amu_fullop_simp_s1);
    free(amu_fulleq_simp_s);free(amu_fulleq_simp_s1);

    //////////////////////
    free_2(2, ms);
    for (int i = 0;i < Nstrange;i++) {
        Meta[i] = nullptr;
    }
    free(Meta);
    free(asd_vec);
    free(afull_vec);
    free(M_PS);free(M_PS_op);
    free_fit_result(fit_info, G_PS);free_fit_result(fit_info, G_PS_OS);
    free_fit_result(fit_info, ZVl);free_fit_result(fit_info, ZAl);
}

