bad=()
./mpi_fpi read_plateaux -p ../../g-2_new_stat/ cA.53.24  -bin 50  -L 24  jack -mu 0.00530 || bad+=("cA.53.24")
./mpi_fpi read_plateaux -p ../../g-2_new_stat/ cA.40.24  -bin 50  -L 24  jack -mu 0.00400 || bad+=("cA.40.24")
./mpi_fpi read_plateaux -p ../../g-2_new_stat/ cA.30.32  -bin 50  -L 32  jack -mu 0.00300 || bad+=("cA.30.32")
if [  $1 == "first_run" ]
then
./mpi_fpi read_plateaux -p ../../g-2_new_stat/ cB.72.64  -bin 50  -L 64  jack -mu 0.00072 || bad+=("cB.72.64")
./mpi_fpi read_plateaux -p ../../g-2_new_stat/ cB.72.96  -bin 50  -L 96  jack -mu 0.00072 || bad+=("cB.72.96")
./mpi_fpi read_plateaux -p ../../g-2_new_stat/ cC.06.80  -bin 50  -L 80  jack -mu 0.00060 || bad+=("cC.06.80")
./mpi_fpi read_plateaux -p ../../g-2_new_stat/ cC.06.112 -bin 50  -L 112 jack -mu 0.00060 || bad+=("cC.06.112")
./mpi_fpi read_plateaux -p ../../g-2_new_stat/ cD.54.96  -bin 50  -L 96  jack -mu 0.00054 || bad+=("cD.54.96") 
./mpi_fpi read_plateaux -p ../../g-2_new_stat/ cE.44.112 -bin 50  -L 112 jack -mu 0.00044 || bad+=("cE.44.112")
./mpi_fpi  read_plateaux -p ../../g-2_new_stat/ cB.14.64  -bin 50  -L 64  jack -mu 0.0014 || bad+=("mpi cB.14.64")
./mpi_fpi  read_plateaux -p ../../g-2_new_stat/ cB.25.48  -bin 50  -L 48  jack -mu 0.0025 || bad+=("mpi cB.25.48")
./mpi_fpi  read_plateaux -p ../../g-2_new_stat/ cC.20.48  -bin 50  -L 48  jack -mu 0.0020 || bad+=("mpi cB.20.48")
fi
./scale_setting jack ../../g-2_new_stat/jackknife/ ../../g-2_new_stat_try//fit_all/|| bad+=("scale setting 1")


./gm2_analysis read_plateaux -p ../../g-2_new_stat/ cA.53.24  -bin 50  -L 24  jack -mu 0.00530  0.010  0.020   0.26500   0.29000   0.30000 three_corr || bad+=("cA.53.24")
./gm2_analysis read_plateaux -p ../../g-2_new_stat/ cA.40.24  -bin 50  -L 24  jack -mu 0.00400  0.010  0.020   0.26500   0.29000   0.30000 three_corr || bad+=("cA.40.24")
./gm2_analysis read_plateaux -p ../../g-2_new_stat/ cA.30.32  -bin 50  -L 32  jack -mu 0.00300  0.010  0.020   0.26500   0.29000   0.30000 three_corr || bad+=("cA.30.32") # the strange are fake

./gm2_analysis read_plateaux -p ../../g-2_new_stat/ cB.72.64  -bin 50  -L 64  jack -mu 0.00072  0.018  0.020   0.21000  0.23000  0.25000  0.0006675 0.00072 0.00072 three_corr|| bad+=("cB.72.64")
#./gm2_analysis read_plateaux -p ../../g-2_new_stat/ cB.72.64  -bin 50  -L 64  jack -mu 0.00072  1.7500e-02  1.8500e-02  0.21000  0.23000  0.25000  0.0006675 0.00072 0.00072 three_corr || bad+=("cB.72.64")

# ./gm2_analysis read_plateaux -p ../../g-2_new_stat/ cB.72.96  -bin 50  -L 96  jack -mu 0.00072  0.018  0.019   0.210  0.230  0.250  0.0006675 0.00072 0.00072 three_corr || bad+=("cB.72.96")
./gm2_analysis read_plateaux -p ../../g-2_new_stat/ cB.72.96  -bin 50  -L 96  jack -mu 0.00072  0.018  0.019   0.21000  0.23000  0.25000  three_corr || bad+=("cB.72.96")
./gm2_analysis read_plateaux -p ../../g-2_new_stat/ cC.06.80  -bin 50  -L 80  jack -mu 0.00060  0.016  0.018   0.17500  0.19500  0.21500  0.0005850 0.00060 0.00060 three_corr || bad+=("cC.06.80")
#./gm2_analysis read_plateaux -p ../../g-2_new_stat/ cC.06.112 -bin 50  -L 112 jack -mu 0.00060  0.016  0.018   0.99999  0.99999  0.99999  three_corr|| bad+=("cC.06.112")
./gm2_analysis read_plateaux -p ../../g-2_new_stat/ cC.06.112 -bin 50  -L 112 jack -mu 0.00060  0.016  0.018   0.18000  0.19000  0.20000  three_corr|| bad+=("cC.06.112")
./gm2_analysis read_plateaux -p ../../g-2_new_stat/ cD.54.96  -bin 50  -L 96  jack -mu 0.00054  0.013  0.014   0.16500  0.17000  0.17500  0.0004964 0.00054 0.00054 three_corr || bad+=("cD.54.96")  
#./gm2_analysis read_plateaux -p ../../g-2_new_stat/ cD.54.96  -bin 50  -L 96  jack -mu 0.00054  0.013  0.014   0.150  0.160  0.170  0.0004964 0.00054 0.00054 three_corr || bad+=("cD.54.96")  # from Ds
./gm2_analysis read_plateaux -p ../../g-2_new_stat/ cE.44.112 -bin 50  -L 112 jack -mu 0.00044  0.011  0.012   0.13000  0.14000  0.15000  three_corr || bad+=("cE.44.112")

./mpi_fpi  read_plateaux -p ../../g-2_new_stat/ cB.14.64  -bin 50  -L 64  jack -mu 0.0014 || bad+=("mpi cB.14.64")
./mpi_fpi  read_plateaux -p ../../g-2_new_stat/ cB.25.48  -bin 50  -L 48  jack -mu 0.0025 || bad+=("mpi cB.25.48")
./mpi_fpi  read_plateaux -p ../../g-2_new_stat/ cC.20.48  -bin 50  -L 48  jack -mu 0.0020 || bad+=("mpi cB.20.48")

./VKVK  jack -bin 50  ../data_VKVK/ VKVK_B64_0.5_40.dat|| bad+=("VKVK_0.5_40")
./VKVK  jack -bin 50  ../data_VKVK/ VKVK_B64_0.5_80.dat|| bad+=("VKVK_0.5_80")



# ./m_K_2x2 read_plateaux -p ../../g-2_new_stat/ cB.72.64 -bin 50  -L 64 jack -mu 1.7500e-02  1.8500e-02   0.0006675 0.00072 || bad+=("mK cB.72.64") # D/Ds
# ./m_K_ndata read_plateaux -p ../../g-2_new_stat/ cB.72.64 -bin 50  -L 64 jack -mu 0.0006675 1.7500e-02 0.0006675 1.8500e-02   0.00072  1.7500e-02   0.00072 1.8500e-02|| bad+=("mK cB.72.64") # # D/Ds
# ./m_K_ndata read_plateaux -p ../../g-2_new_stat/ cB.72.64 -bin 50  -L 64 jack -mu 0.00072  0.018  0.00072 0.020   0.0006675  0.018254  0.00072  0.018254  || bad+=("mK cB.72.64") # tau
./m_K_ndata read_plateaux -p ../../g-2_new_stat/ cB.72.64 -bin 50  -L 64 jack -mu 0.00072 0.018         0.00072 0.020  \
                                                                                  0.0006675 0.018254    0.00072 0.018254  \
                                                                                  0.0006675 1.7500e-02  0.0006675 1.8500e-02 \
                                                                                  0.00072 1.7500e-02    0.00072 1.8500e-02 || bad+=("mK cB.72.64")  #all

# ./m_K_2x2 read_plateaux -p ../../g-2_new_stat/ cC.06.80 -bin 50  -L 80 jack -mu 1.6000e-02  1.7000e-02   0.0005850 0.00060 || bad+=("mK cC.06.80")
./m_K_ndata read_plateaux -p ../../g-2_new_stat/  cC.06.80 -bin 50  -L 80 jack -mu  0.00060    0.016        0.00060    0.018 \
                                                                                    0.000585   0.016067     0.00060    0.016067  \
                                                                                    0.0005850  1.6000e-02   0.0005850  1.7000e-02 \
                                                                                    0.00060    1.6000e-02   0.00060    1.7000e-02  || bad+=("mK cC.06.80")

./m_K_ndata read_plateaux -p ../../g-2_new_stat/  cD.54.96  -bin 50  -L 96  jack -mu  0.00054     0.013        0.00054    0.014 \
                                                                                      0.0004964   0.013557     0.0004964  0.013557  || bad+=("mK cD.54.96")
./m_K_ndata read_plateaux -p ../../g-2_new_stat/ cE.44.112 -bin 50  -L 112 jack -mu 0.00044 0.011    0.00044 0.012\
                                                                                    4.3100e-4 1.1759e-02  0.00044  1.1759e-02|| bad+=("mK cE.44.112")
# ./m_K_2p2 read_plateaux -p ../../g-2_new_stat/ cE.44.112 -bin 50  -L 112 jack -mu 0.011  0.012  0.00044    4.3100e-4 0.00044  1.1759e-02|| bad+=("mK cE.44.112")
# ./m_K_ndata read_plateaux -p ../../g-2_new_stat/ cE.44.112 -bin 50  -L 112 jack -mu 0.00044 0.011    0.00044 0.012     4.3100e-4 1.1759e-02  0.00044  1.1759e-02|| bad+=("mK cE.44.112")


./m_Ds_ndata read_plateaux -p ../../g-2_new_stat/ cB.72.64 -bin 50  -L 64 jack -mu 1.7500e-02 0.23000  1.7500e-02 0.24000\
                                                                                   1.8500e-02 0.23000  1.8500e-02 0.24000|| bad+=("mDs cB.72.64")
./m_Ds_ndata read_plateaux -p ../../g-2_new_stat/ cC.06.80 -bin 50  -L 80 jack -mu 1.6000e-02 0.18000    1.6000e-02 0.19000\
                                                                                   1.7000e-02 0.18000    1.7000e-02 0.19000|| bad+=("mDs cC.06.80")
./m_Ds_ndata read_plateaux -p ../../g-2_new_stat/ cD.54.96 -bin 50  -L 96 jack -mu 1.3557e-02   0.15000   1.3557e-02 0.16000 \
                                                                                   1.3557e-02 0.17000  || bad+=("mDs cD.54.96")
./m_Ds_ndata read_plateaux -p ../../g-2_new_stat/ cE.44.112 -bin 50  -L 112 jack -mu 1.1759e-02   0.13000   1.1759e-02 0.14000 \
                                                                                   1.1759e-02 0.15000  || bad+=("mDs cE.44.112")

# ./m_Ds_2x2 read_plateaux -p ../../g-2_new_stat/ cB.72.64 -bin 50  -L 64 jack -mu 1.7500e-02  1.8500e-02 0.230  0.240  || bad+=("mDs cB.72.64")
# ./m_Ds_2x2 read_plateaux -p ../../g-2_new_stat/ cC.06.80 -bin 50  -L 80 jack -mu 1.6000e-02  1.7000e-02 0.180  0.190  || bad+=("mDs cC.06.80")
# ./m_Ds read_plateaux -p ../../g-2_new_stat/ cC.06.112 -bin 50  -L 112 jack -mu 1.6067e-02 0.180  0.190  0.200 || bad+=("mDs cC.06.112")
# ./m_Ds read_plateaux -p ../../g-2_new_stat/ cD.54.96  -bin 50  -L 96  jack -mu 1.3557e-02 0.150  0.160  0.170 || bad+=("mDs cD.54.96")
# ./m_Ds read_plateaux -p ../../g-2_new_stat/ cE.44.112 -bin 50  -L 112 jack -mu 1.1759e-02 0.130  0.140  0.150 || bad+=("mDs cE.44.112")


#./g-2 read_plateaux -p ../../../g-2/ test -bin 50  -L 96 jack -mu 0.00054  0.014  0.015   0.165  0.175  0.175|| bad+=("test")
#./fit_all_gm2  jack ../../g-2_new_stat/jackknife/ ../../g-2_new_stat//fit_all/
./some_custom_fits jack ../../g-2_new_stat/jackknife/ ../../g-2_new_stat//fit_all/  || bad+=("custom fit ")

./fit_strange jack ../../g-2_new_stat/jackknife/ ../../g-2_new_stat//fit_all_strange/ ||bad+=("fit strange ")
./fit_charm jack ../../g-2_new_stat/jackknife/ ../../g-2_new_stat//fit_all_charm/ ||bad+=("fit charm ")

if [[ -n "${bad-}" ]]; then
  echo -e "\n error in:\n"
  for path in "${bad[@]}"; do
    echo "  - $path"
  done

  exit 1
fi

echo ""
echo all good
echo ""
