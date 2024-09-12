#define CONTROL
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
// #include <time.h>
// #include <complex.h>
#include <cstring>
#include <fstream>
#include <iostream>
#include <string>

#include "global.hpp"
#include "read.hpp"
#include "resampling_new.hpp"

#include "linear_fit.hpp"
#include "mutils.hpp"
#include "various_fits.hpp"

#include "correlators_analysis.hpp"

#include "tower.hpp"
struct kinematic kinematic_2pt;


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
        std::string line;
        std::ifstream file(stream);
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

        to_bin = std::vector<int>(iconfs.size());
        next_to_bin = std::vector<int>(iconfs.size());
    }
    configuration_class(const char stream[NAMESIZE], int T) {
        printf("second construction\n");
        long int tmp;
        int s = T + 1 + 3;
        int count = 0;
        std::string line;
        std::ifstream file(stream);
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

        std::cout << "lines=" << 0 << std::endl;

        std::cout << "confs=" << 0 << std::endl;

        to_bin = std::vector<int>(iconfs.size());
        next_to_bin = std::vector<int>(iconfs.size());
    }

};



generic_header read_head(FILE* stream) {
    generic_header header;
    return header;
}
void write_header_g2(FILE* jack_file, generic_header head) {
    int fi = 0;
    fi += fwrite(&head.T, sizeof(int), 1, jack_file);
    fi += fwrite(&head.L, sizeof(int), 1, jack_file);
    int nmu = head.mus.size();
    fi += fwrite(&nmu, sizeof(int), 1, jack_file);
    for (double mu : head.mus) {
        fi += fwrite(&mu, sizeof(double), 1, jack_file);
    }
}

char** argv_to_options(char** argv) {
    char** option;
    option = (char**)malloc(sizeof(char*) * 7);
    option[0] = (char*)malloc(sizeof(char) * NAMESIZE);
    option[1] = (char*)malloc(sizeof(char) * NAMESIZE);
    option[2] = (char*)malloc(sizeof(char) * NAMESIZE);
    option[3] = (char*)malloc(sizeof(char) * NAMESIZE);
    option[4] = (char*)malloc(sizeof(char) * NAMESIZE);
    option[5] = (char*)malloc(sizeof(char) * NAMESIZE);
    option[6] = (char*)malloc(sizeof(char) * NAMESIZE);

    mysprintf(option[1], NAMESIZE, "read_plateaux"); // blind/see/read_plateaux
    mysprintf(option[2], NAMESIZE, "-p");            // -p
    mysprintf(option[3], NAMESIZE, argv[4]);         // path
    mysprintf(option[4], NAMESIZE, argv[1]);         // resampling
    mysprintf(option[5], NAMESIZE, "no");            // pdf

    std::string name(argv[4]);
    std::string en = name.substr(name.find_last_of('/') + 1);

    mysprintf(option[6], NAMESIZE, argv[5]);         // infile
    return option;
}

void init_global_head(generic_header head) {
    file_head.l1 = head.L;
    file_head.l0 = head.T;
    file_head.l2 = head.L;
    file_head.l3 = head.L;
    file_head.nk = 2;
    file_head.musea = head.mus[0];
    file_head.k = (double*)malloc(sizeof(double) * file_head.nk * 2);
    file_head.k[0] = 0;
    file_head.k[1] = 0;
    file_head.k[2] = head.mus[0];
    file_head.k[3] = head.mus[0];

    file_head.nmoms = 1;
    file_head.mom = (double**)malloc(sizeof(double*) * file_head.nmoms);
    for (int i = 0; i < file_head.nmoms; i++) {
        file_head.mom[i] = (double*)malloc(sizeof(double) * 4);
        file_head.mom[i][0] = 0;
        file_head.mom[i][1] = 0;
        file_head.mom[i][2] = 0;
        file_head.mom[i][3] = 0;
    }
}

// void read_twopt(FILE* stream, double*** to_write, generic_header head) {
//     // write your function to read the data
//     // int fi = 0;
//     // for (int k = 0; k < head.ncorr; k++) {
//     //     for (int t = 0; t < head.T;t++) {
//     //         fi += fscanf(stream, "%lf  %lf", &to_write[k][t][0],
//     //         &to_write[k][t][1]);
//     //         // printf(" corr=%d t=%d %.12g   %.12g\n", k, t, to_write[k][t][0],
//     //         to_write[k][t][1]);
//     //     }
//     // }
//     int fi = 0;
//     int id;
//     int i = fread(&id, sizeof(int), 1, stream);
//     for (int k = 0; k < head.ncorr; k++) {
//         for (int t = 0; t < head.T; t++) {
//             fi += fread(to_write[k][t], sizeof(double), 2, stream);
//         }
//     }
// }



int id_twpt(generic_header head, int ig, int is, int TMOS) {

    // return id_counter + (head.oranges.size()) * (ig
    //     + head.gammas.size() * (TMOS + head.rs.size() * (im1 + 2 * (im))));
    return  ig + is * head.gammas.size() + TMOS * head.smearing.size() * head.gammas.size();
}


double** corr_sum_VKVK(int j, double**** in, int t, struct fit_type fit_info) {
    double** r = double_calloc_2(fit_info.N, 2);

    int id0 = fit_info.corr_id[0];
    int id1 = fit_info.corr_id[1];
    int id2 = fit_info.corr_id[2];

    r[0][0] = (in[j][id0][t][0] + in[j][id1][t][0] + in[j][id2][t][0]) / 3.0;
    r[0][1] = (in[j][id0][t][1] + in[j][id1][t][1] + in[j][id2][t][1]) / 3.0;
    return r;
}

int main(int argc, char** argv) {
    error(argc < 4, 1, "nissa_mpcac ",
        "usage:./nissa_mpcac jack/boot bin -p path file   \n separate "
        "path and file please");

    char resampling[NAMESIZE];
    mysprintf(resampling, NAMESIZE, argv[1]);
    printf("resampling: %s\n", resampling);

    char** option = argv_to_options(argv);

    char namefile[NAMESIZE];
    mysprintf(namefile, NAMESIZE, "%s/%s", option[3], option[6]);

    char namefile_plateaux[NAMESIZE];
    mysprintf(namefile_plateaux, NAMESIZE, "plateaux.txt");

    FILE* infile = open_file(namefile, "r");

    //////////////////////////////////// read and setup header
    generic_header head;
    head.read_header_debug(infile);
    // head.print_header();
    init_global_head(head);
    // std::vector<std::string>  correlators;

    // std::vector<configuration_class> myconfs;
    // for (int i = 0;i < argc - 5;i++) {
    //     mysprintf(namefile, NAMESIZE, "%s/%s", argv[4], argv[5 + i]);//54
    //     printf("file: %s\n", namefile);
    //     correlators.emplace_back(namefile);
    //     configuration_class tmp(namefile, 0);
    //     myconfs.emplace_back(tmp);
    // }


    // // generic_header head = read_head(infile);
    // generic_header head;
    // int seed;
    // double tmpd;
    // head.Njack;
    // head.ncorr = 0;
    // head.mus.resize(1);

    // line_read_param(option, "LT", head.L, head.T, seed, namefile_plateaux);
    // line_read_param(option, "mu", head.mus[0], tmpd, seed, namefile_plateaux);
    // head.print_header();
    // init_global_head(head);

    //////////////////////////////////// setup jackboot and binning
    int confs = head.Njack;
    int bin = atoi(argv[3]);
    int Neff = bin;
    int Njack;
    if (strcmp(argv[1], "jack") == 0) {
        Njack = Neff + 1;
        myres = new resampling_jack(Neff);
    }
    else if (strcmp(argv[1], "boot") == 0) {
        Njack = (Neff * 2 + 1);
        myres = new resampling_boot(Neff * 2);
    }
    else {
        Njack = 0;
        error(1 == 1, 1, "main", "argv[1]= %s is not jack or boot", argv[1]);
    }
    // now Njack need to be the number of jacks
    head.Nconf = head.Njack;
    head.Njack = Njack;
    //////////////////////////////////// setup output files
    mysprintf(namefile, NAMESIZE, "%s/out/%s_output", option[3], option[6]);
    printf("writing output in :\n %s \n", namefile);
    FILE* outfile = open_file(namefile, "w+");

    mysprintf(namefile, NAMESIZE, "%s/jackknife/%s_%s", option[3], option[4],
        option[6]);
    FILE* jack_file = open_file(namefile, "w+");
    // write_header_g2(jack_file, head);
    head.write_header(jack_file);

    //////////////////////////////////// confs
    int ncorr_new = head.ncorr;
    int Max_corr = head.ncorr + head.smearing.size() * head.rs.size() + 2 * head.rs.size();
    double**** data = calloc_corr(confs, Max_corr, head.T);

    printf("confs=%d\n", confs);
    printf("ncorr=%d\n", head.ncorr);
    printf("kappa=%g\n", head.kappa);
    printf("bin=%d\n", bin);
    for (int iconf = 0; iconf < confs; iconf++) {
        read_twopt_binary(infile, data[iconf], head);
    }
    //////////////////////////////////////////////////////////////
    // normalization
    //////////////////////////////////////////////////////////////
    // for (int iconf = 0; iconf < confs; iconf++) {
    //     // for (int reim = 0;reim < 2;reim++) {
    //     for (int t = 0;t < head.T;t++) {
    //         for (int n = 0;n < head.ncorr;n++) {
    //             data[iconf][n][t][0] /= head.L * head.L * head.L;
    //             data[iconf][n][t][1] /= head.L * head.L * head.L;
    //         }
    //     }
    //     // }
    // }


    //////////////////////////////////////////////////////////////
    // convolution
    //////////////////////////////////////////////////////////////
    // for (int iconf = 0; iconf < confs; iconf++) {
    //     // for (int reim = 0;reim < 2;reim++) {
    //     for (int t = 0;t < head.T;t++) {
    //         data[iconf][1][t][0] = 0;
    //         data[iconf][1][t][1] = 0;
    //         for (int t0 = 0;t0 < head.T;t0++) {
    //             std::complex ct0(data[iconf][0][t0][0], data[iconf][0][t0][1]);
    //             std::complex ct0pt(data[iconf][0][(t0 + t) % head.T][0], data[iconf][0][(t0 + t) % head.T][1]);

    //             // data[iconf][1][t][0] += real(ct0 * ct0pt);
    //             data[iconf][1][t][0] += data[iconf][0][t0][0] * data[iconf][0][(t0 + t) % head.T][0];
    //             data[iconf][1][t][1] += imag(ct0 * ct0pt);

    //         }
    //         data[iconf][1][t][0] /= head.T;
    //         data[iconf][1][t][1] /= head.T;
    //     }
    //     // }
    // }


    double**** data_bin = binning_toNb(confs, Max_corr, head.T, data, bin);
    free_corr(confs, Max_corr, head.T, data);
    double**** conf_jack = myres->create(Neff, Max_corr, head.T, data_bin);
    free_corr(Neff, Max_corr, head.T, data_bin);

    //////////////////////////////////////////////////////////////
    // read connected
    //////////////////////////////////////////////////////////////
    // printf("////////////  read connected ////////////////\n");
    // std::string name_conn(option[6]);


    // size_t lastindex = name_conn.find_last_of(".");
    // std::string rawname = name_conn.substr(0, lastindex) + "_conn.dat";
    // mysprintf(namefile, NAMESIZE, "%s/%s", option[3], rawname.c_str());
    // printf("reading : %s \n", namefile);

    // FILE* infile_con = open_file(namefile, "r");
    // generic_header head_con;
    // head_con.read_header_debug(infile_con);
    // init_global_head(head_con);

    // double**** data_conn = calloc_corr(confs, 1, head_con.T);
    // for (int iconf = 0; iconf < head_con.Njack; iconf++) {
    //     read_twopt_binary(infile_con, data_conn[iconf], head_con);
    // }
    // double**** data_bin_conn = binning_toNb(confs, 1, head.T, data_conn, bin);
    // free_corr(confs, 1, head.T, data_conn);
    // double**** conf_jack_conn = myres->create(Neff, 1, head.T, data_bin_conn);
    // free_corr(Neff, 1, head.T, data_bin_conn);

    // printf("//////////// end read connected ////////////////\n");

    // printf("////////////  merging connected into database connected ////////////////\n");
    // for (size_t t = 0; t < head.T; t++) {
    //     for (size_t r = 0; r < 2; r++) {
    //         for (size_t j = 0; j < Njack; j++) {
    //             conf_jack[j][ncorr_new][t][r] = conf_jack_conn[j][0][t][r];
    //             conf_jack[j][ncorr_new][t][r] /= (head.L * head.L * head.L);
    //         }
    //     }
    // }
    // int id_conn_OS = ncorr_new;
    // ncorr_new++;
    // symmetrise_jackboot(Njack, ncorr_new - 1, head.T, conf_jack);
    // printf("//////////// end merging connected ////////////////\n");


    /////////////////////////////////////////////////////////////////////////////////////////////////////////
    // print all the effective masses correlators
    // set the option to not read for a plateaux
    mysprintf(namefile, NAMESIZE, "%s/out/%s_meff_correlators", option[3],
        option[6]);
    FILE* outfile_meff_corr = open_file(namefile, "w+");
    mysprintf(namefile, NAMESIZE, "%s/out/%s_raw_correlators", option[3],
        option[6]);
    FILE* outfile_raw_corr = open_file(namefile, "w+");
    mysprintf(namefile, NAMESIZE, "%s/out/%s_shifted_correlators", option[3],
        option[6]);
    FILE* outfile_shifted_corr = open_file(namefile, "w+");
    mysprintf(namefile, NAMESIZE, "%s/out/%s_log_meff_shifted", option[3],
        option[6]);
    FILE* outfile_log_meff_shifted = open_file(namefile, "w+");
    mysprintf(namefile, NAMESIZE, "%s/out/%s_gamma", option[3], option[6]);
    FILE* out_gamma = open_file(namefile, "w+");

    mysprintf(namefile, NAMESIZE, "%s/out/%s_HLT_kernel", option[3],
        option[6]);
    FILE* outfile_HLT_kernel = open_file(namefile, "w+");
    mysprintf(namefile, NAMESIZE, "%s/out/%s_HLT_AoverB", option[3],
        option[6]);
    FILE* outfile_HLT_AoverB = open_file(namefile, "w+");


    char save_option[NAMESIZE];
    sprintf(save_option, "%s", option[1]);
    sprintf(option[1], "blind");
    FILE* dev_null = open_file("/dev/null", "w");
    struct fit_type fit_info_silent;
    fit_info_silent.verbosity = -1;
    fit_info_silent.chi2_gap_jackboot = 1e+6;
    fit_info_silent.guess_per_jack = 0;
    printf("ncorr_new:  %d\n", ncorr_new);
    for (int icorr = 0; icorr < ncorr_new; icorr++) {
        // log effective mass
        double* tmp_meff_corr = plateau_correlator_function(
            option, kinematic_2pt, (char*)"P5P5", conf_jack, Njack,
            namefile_plateaux, outfile_meff_corr, icorr, "log", M_eff_log, dev_null,
            fit_info_silent);
        free(tmp_meff_corr);
        // raw correlator
        file_head.l0 = head.T * 2;
        tmp_meff_corr = plateau_correlator_function(
            option, kinematic_2pt, (char*)"P5P5", conf_jack, Njack,
            namefile_plateaux, outfile_raw_corr, icorr, "cor", identity, dev_null,
            fit_info_silent);
        free(tmp_meff_corr);
        tmp_meff_corr = plateau_correlator_function(
            option, kinematic_2pt, (char*)"P5P5", conf_jack, Njack,
            namefile_plateaux, outfile_raw_corr, icorr, "cor", identity_im,
            dev_null, fit_info_silent);
        free(tmp_meff_corr);
        file_head.l0 = head.T;
        // shifted correlator
        tmp_meff_corr = plateau_correlator_function(
            option, kinematic_2pt, (char*)"P5P5", conf_jack, Njack,
            namefile_plateaux, outfile_shifted_corr, icorr, "shift_cor", shift_corr,
            dev_null, fit_info_silent);
        free(tmp_meff_corr);
        // log_meff shifted correlator
        tmp_meff_corr = plateau_correlator_function(
            option, kinematic_2pt, (char*)"P5P5", conf_jack, Njack,
            namefile_plateaux, outfile_log_meff_shifted, icorr, "log_shift",
            M_eff_log_shift, dev_null, fit_info_silent);
        free(tmp_meff_corr);
    }
    fit_info_silent.restore_default();
    sprintf(option[1], "%s", save_option); // restore option
    corr_counter = -1;

    ////////////////////////////////////////////////////////////
    // symmetrize
    ////////////////////////////////////////////////////////////
    // symmetrise_jackboot(Njack, 0, head.T, conf_jack);
    // symmetrise_jackboot(Njack, 1, head.T, conf_jack, -1);

    ////////////////////////////////////////////////////////////
    // start fitting
    //////////////////////////////
    corr_counter = -1;

    int id_P5P5_TM_ll = id_twpt(head, 0, 0, 0);
    double* M_PS = plateau_correlator_function(
        option, kinematic_2pt, (char*)"P5P5", conf_jack, Njack,
        namefile_plateaux, outfile, 0, "M_{pi}", M_eff_T, jack_file);
    free(M_PS);
    check_correlatro_counter(0);

    //////////////////////////////////////////////////////////////
    // summing VKVK
    //////////////////////////////////////////////////////////////
    std::vector<int> id_VKVK(head.smearing.size() * head.rs.size());

    for (int is = 0;is < head.smearing.size(); is++) {
        for (int ir = 0;ir < head.rs.size(); ir++) {
            fit_type fit_info;
            fit_info.N = 1;
            fit_info.Njack = Njack;
            printf("Ncorrelator=%d\n", ncorr_new);
            int id_V1V1 = id_twpt(head, 1, is, ir);
            int id_V2V2 = id_twpt(head, 2, is, ir);
            int id_V3V3 = id_twpt(head, 3, is, ir);
            fit_info.corr_id = { id_V1V1, id_V2V2, id_V3V3 };
            add_correlators_no_alloc(option, ncorr_new, conf_jack, corr_sum_VKVK, fit_info, Max_corr);
            id_VKVK[is + ir * head.smearing.size()] = ncorr_new - 1;

        }
    }
    //////////////////////////////////////////////////////////////
    // symm sl ls
    //////////////////////////////////////////////////////////////
    for (int ir = 0;ir < head.rs.size(); ir++) {
        int a = id_VKVK[1 + ir * head.smearing.size()];
        int b = id_VKVK[2 + ir * head.smearing.size()];
        for (int t = 0;t < head.T;t++) {
            for (int j = 0;j < Njack;j++) {
                std::complex c_01 = std::complex(conf_jack[j][a][t][0], conf_jack[j][a][t][1]);
                std::complex c_10 = std::complex(conf_jack[j][b][t][0], conf_jack[j][b][t][1]);
                c_01 = std::complex(c_01.real() + c_10.real(), c_01.imag() - c_10.imag()) / 2.0;

                conf_jack[j][a][t][0] = c_01.real();
                conf_jack[j][a][t][1] = c_01.imag();

                conf_jack[j][b][t][0] = c_01.real();
                conf_jack[j][b][t][1] = -c_01.imag();

            }
        }
    }
    //////////////////////////////////////////////////////////////
    // m_eff 
    //////////////////////////////////////////////////////////////

    for (int is = 0;is < head.smearing.size(); is++) {
        for (int ir = 0;ir < head.rs.size(); ir++) {
            std::string namefit("meff_VKVK");
            if (is == 0) namefit += "_ll";
            if (is == 1) namefit += "_sl";
            if (is == 2) namefit += "_ls";
            if (is == 3) namefit += "_ss";

            if (ir == 0) namefit += "_TM";
            if (ir == 1) namefit += "_OS";

            double* M_PS = plateau_correlator_function(
                option, kinematic_2pt, (char*)"P5P5", conf_jack, Njack,
                namefile_plateaux, outfile, id_VKVK[is + ir * head.smearing.size()]
                , namefit.c_str(), M_eff_T, jack_file);
            free(M_PS);
        }
    }

    //////////////////////////////////////////////////////////////
    // GEVP
    //////////////////////////////////////////////////////////////

    for (int ir = 0;ir < head.rs.size(); ir++) {
        fit_type fit_info;
        fit_info.N = 2;
        fit_info.corr_id = { id_VKVK[0 + ir * head.smearing.size()], id_VKVK[1 + ir * head.smearing.size()] ,
                                                                     id_VKVK[3 + ir * head.smearing.size()] };
        fit_info.value_or_vector = 0; // 0= values
        fit_info.t0_GEVP = 3;
        fit_info.GEVP_ignore_warning_after_t = 10;
        fit_info.verbosity = 2;
        fit_info.Njack = Njack;
        printf(" ncorr before GEVP %d\n", ncorr_new);
        

        // add_correlators(option, ncorr_new, conf_jack, GEVP_matrix, fit_info);
        add_correlators_no_alloc(option, ncorr_new, conf_jack, GEVP_matrix, fit_info, Max_corr);
        
        printf(" ncorr after GEVP %d\n", ncorr_new);
        printf("  GEVP  l0 %d\n", ncorr_new-2);
        printf("  GEVP  l1 %d\n", ncorr_new-1);
        
        
        std::string namefit("meff_GEVP_VKVK");
        if (ir == 0) namefit += "_TM";
        if (ir == 1) namefit += "_OS";

        std::string namefit_l0 = namefit + "_l0";
        double* M_PS = plateau_correlator_function(
            option, kinematic_2pt, (char*)"P5P5", conf_jack, Njack,
            namefile_plateaux, outfile, ncorr_new - 2
            , namefit_l0.c_str(), M_eff_T, jack_file);
        free(M_PS);

        std::string namefit_l1 = namefit + "_l1";
        double* M_PS1 = plateau_correlator_function(
            option, kinematic_2pt, (char*)"P5P5", conf_jack, Njack,
            namefile_plateaux, outfile, ncorr_new - 1
            , namefit_l1.c_str(), M_eff_T, jack_file);
        free(M_PS1);

        fit_info.restore_default();
    }
    free_corr(Njack, Max_corr, head.T, conf_jack);
    for (size_t i = 0; i < 7; i++) {
        free(option[i]);
    }
    free(option);


}