bad=()
./gm2_analysis read_plateaux -p ../../g-2_new_stat/ cA.53.24  -bin 50  -L 24  jack -mu 0.00530  0.010  0.020   0.265   0.290   0.300 three_corr || bad+=("cA.53.24")
./gm2_analysis read_plateaux -p ../../g-2_new_stat/ cA.40.24  -bin 50  -L 24  jack -mu 0.00400  0.010  0.020   0.265   0.290   0.300 three_corr || bad+=("cA.40.24")
./gm2_analysis read_plateaux -p ../../g-2_new_stat/ cA.30.32  -bin 50  -L 32  jack -mu 0.00300  0.010  0.020   0.265   0.290   0.300 three_corr || bad+=("cA.30.32")
./gm2_analysis read_plateaux -p ../../g-2_new_stat/ cB.72.64  -bin 50  -L 64  jack -mu 0.00072  0.018  0.020   0.210  0.230  0.250  0.0006675 0.00072 0.00072 three_corr || bad+=("cB.72.64")
# ./gm2_analysis read_plateaux -p ../../g-2_new_stat/ cB.72.96  -bin 50  -L 96  jack -mu 0.00072  0.018  0.019   0.210  0.230  0.250  0.0006675 0.00072 0.00072 three_corr || bad+=("cB.72.96")
./gm2_analysis read_plateaux -p ../../g-2_new_stat/ cB.72.96  -bin 50  -L 96  jack -mu 0.00072  0.018  0.019   0.210  0.230  0.250  three_corr || bad+=("cB.72.96")
./gm2_analysis read_plateaux -p ../../g-2_new_stat/ cC.06.80  -bin 50  -L 80  jack -mu 0.00060  0.016  0.018   0.175  0.195  0.215  0.0005850 0.00060 0.00060 three_corr || bad+=("cC.06.80")
./gm2_analysis read_plateaux -p ../../g-2_new_stat/ cC.06.112 -bin 50  -L 112 jack -mu 0.00060  0.016  0.018   0.99999  0.99999  0.99999  three_corr|| bad+=("cC.06.112")
# ./gm2_analysis read_plateaux -p ../../g-2_new_stat/ cC.06.112 -bin 50  -L 112 jack -mu 0.00060  0.016  0.018   0.180  0.190  0.200  three_corr|| bad+=("cC.06.112")
./gm2_analysis read_plateaux -p ../../g-2_new_stat/ cD.54.96  -bin 50  -L 96  jack -mu 0.00054  0.013  0.014   0.165  0.175  0.175  0.0004964 0.00054 0.00054 three_corr || bad+=("cD.54.96")  
#./gm2_analysis read_plateaux -p ../../g-2_new_stat/ cD.54.96  -bin 50  -L 96  jack -mu 0.00054  0.013  0.014   0.150  0.160  0.170  0.0004964 0.00054 0.00054 three_corr || bad+=("cD.54.96")  # from Ds
./gm2_analysis read_plateaux -p ../../g-2_new_stat/ cE.44.112 -bin 50  -L 112 jack -mu 0.00044  0.011  0.012   0.130  0.140  0.150  three_corr || bad+=("cE.44.112")

./mpi_fpi  read_plateaux -p ../../g-2_new_stat/ cB.14.64  -bin 50  -L 64  jack -mu 0.0014 || bad+=("mpi cB.14.64")
./mpi_fpi  read_plateaux -p ../../g-2_new_stat/ cB.25.48  -bin 50  -L 48  jack -mu 0.0025 || bad+=("mpi cB.25.48")
./mpi_fpi  read_plateaux -p ../../g-2_new_stat/ cC.20.48  -bin 50  -L 48  jack -mu 0.0020 || bad+=("mpi cB.20.48")

./VKVK  jack -bin 50  ../data_VKVK/ VKVK_B64_0.5_40.dat|| bad+=("VKVK_0.5_40")
./VKVK  jack -bin 50  ../data_VKVK/ VKVK_B64_0.5_80.dat|| bad+=("VKVK_0.5_80")


./m_Ds_2x2 read_plateaux -p ../../g-2_new_stat/data_strange_charm_B64 cB.72.64 -bin 50  -L 64 jack -mu 1.7500e-02  1.8500e-02 0.230  0.240  || bad+=("mDs cB.72.64")
./m_Ds_2x2 read_plateaux -p ../../g-2_new_stat/ cC.06.80 -bin 50  -L 80 jack -mu 1.6000e-02  1.7000e-02 0.180  0.190  || bad+=("mDs cC.06.80")

./m_Ds read_plateaux -p ../../g-2_new_stat/ cC.06.112 -bin 50  -L 112 jack -mu 1.6067e-02 0.180  0.190  0.200 || bad+=("mDs cC.06.112")
./m_Ds read_plateaux -p ../../g-2_new_stat/ cD.54.96  -bin 50  -L 96  jack -mu 1.3557e-02 0.150  0.160  0.170 || bad+=("mDs cD.54.96")
./m_Ds read_plateaux -p ../../g-2_new_stat/ cE.44.112 -bin 50  -L 112 jack -mu 1.1759e-02 0.130  0.140  0.150 || bad+=("mDs cE.44.112")


#./g-2 read_plateaux -p ../../../g-2/ test -bin 50  -L 96 jack -mu 0.00054  0.014  0.015   0.165  0.175  0.175|| bad+=("test")
#./fit_all_gm2  jack ../../g-2_new_stat/jackknife/ ../../g-2_new_stat//fit_all/
./some_custom_fits jack ../../g-2_new_stat/jackknife/ ../../g-2_new_stat//fit_all/  || bad+=("custom fit ")
./scale_setting jack ../../g-2_new_stat/jackknife/ ../../g-2_new_stat//fit_all/|| bad+=("scale setting ")

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
