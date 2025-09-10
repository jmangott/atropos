python3 scripts/input_generation/set_boolean_pancreatic_cancer.py -pb -r 5 10 10
./bin/hierarchical-cme -o pancreatic_pbn_r5_10_r_tau1e-2 -s 100 -t 0.01 -f 20.0 -m r

python3 scripts/input_generation/set_boolean_pancreatic_cancer.py -pb -r 5 20 20
./bin/hierarchical-cme -o pancreatic_pbn_r5_20_r_tau1e-2 -s 100 -t 0.01 -f 20.0 -m r

python3 scripts/input_generation/set_boolean_pancreatic_cancer.py -pb -r 5 30 30
./bin/hierarchical-cme -o pancreatic_pbn_r5_30_r_tau1e-2 -s 100 -t 0.01 -f 20.0 -m r

python3 scripts/input_generation/set_boolean_pancreatic_cancer.py -pb -r 10 20 20
./bin/hierarchical-cme -o pancreatic_pbn_r10_20_r_tau1e-2 -s 100 -t 0.01 -f 20.0 -m r

python3 scripts/input_generation/set_boolean_pancreatic_cancer.py -pb -r 10 30 30
./bin/hierarchical-cme -o pancreatic_pbn_r10_30_r_tau1e-2 -s 100 -t 0.01 -f 20.0 -m r

python3 scripts/input_generation/set_boolean_pancreatic_cancer.py -pb -r 10 40 40
./bin/hierarchical-cme -o pancreatic_pbn_r10_40_r_tau1e-2 -s 100 -t 0.01 -f 20.0 -m r

python3 scripts/input_generation/set_boolean_pancreatic_cancer.py -pb -r 20 30 30
./bin/hierarchical-cme -o pancreatic_pbn_r20_30_r_tau1e-2 -s 100 -t 0.01 -f 20.0 -m r

python3 scripts/input_generation/set_boolean_pancreatic_cancer.py -pb -r 20 40 40
./bin/hierarchical-cme -o pancreatic_pbn_r20_40_r_tau1e-2 -s 100 -t 0.01 -f 20.0 -m r

python3 scripts/input_generation/set_boolean_pancreatic_cancer.py -pb -r 5 5 5
./bin/hierarchical-cme -o pancreatic_pbn_r5_r_tau1e-2 -s 100 -t 0.01 -f 20.0 -m r

python3 scripts/input_generation/set_boolean_pancreatic_cancer.py -pb -r 10 10 10
./bin/hierarchical-cme -o pancreatic_pbn_r10_r_tau1e-2 -s 100 -t 0.01 -f 20.0 -m r

python3 scripts/input_generation/set_boolean_pancreatic_cancer.py -pb -r 20 20 20
./bin/hierarchical-cme -o pancreatic_pbn_r20_r_tau1e-2 -s 100 -t 0.01 -f 20.0 -m r

python3 scripts/input_generation/set_boolean_pancreatic_cancer.py -pb -r 20 50 50
./bin/hierarchical-cme -o pancreatic_pbn_r20_50_r_tau1e-2 -s 100 -t 0.01 -f 20.0 -m r


python3 scripts/input_generation/set_boolean_pancreatic_cancer.py -pb -r 5
./bin/hierarchical-cme -o pancreatic_pbn_r5_e_tau1e-2 -s 100 -t 0.01 -f 20.0 -m e

python3 scripts/input_generation/set_boolean_pancreatic_cancer.py -pb -r 10
./bin/hierarchical-cme -o pancreatic_pbn_r10_e_tau1e-2 -s 100 -t 0.01 -f 20.0 -m e

python3 scripts/input_generation/set_boolean_pancreatic_cancer.py -pb -r 20
./bin/hierarchical-cme -o pancreatic_pbn_r20_e_tau1e-2 -s 100 -t 0.01 -f 20.0 -m e


python3 scripts/input_generation/set_boolean_pancreatic_cancer.py -pr -r 5
./bin/hierarchical-cme -o pancreatic_prn_r5_e_tau1e-2 -s 100 -t 0.01 -f 20.0 -m e

python3 scripts/input_generation/set_boolean_pancreatic_cancer.py -pr -r 10
./bin/hierarchical-cme -o pancreatic_prn_r10_e_tau1e-2 -s 100 -t 0.01 -f 20.0 -m e

python3 scripts/input_generation/set_boolean_pancreatic_cancer.py -pr -r 20
./bin/hierarchical-cme -o pancreatic_prn_r20_e_tau1e-2 -s 100 -t 0.01 -f 20.0 -m e


python3 scripts/input_generation/set_boolean_pancreatic_cancer.py -pw -r 5
./bin/hierarchical-cme -o pancreatic_pwn_r5_e_tau1e-2 -s 100 -t 0.01 -f 20.0 -m e

python3 scripts/input_generation/set_boolean_pancreatic_cancer.py -pw -r 10
./bin/hierarchical-cme -o pancreatic_pwn_r10_e_tau1e-2 -s 100 -t 0.01 -f 20.0 -m e

python3 scripts/input_generation/set_boolean_pancreatic_cancer.py -pw -r 20
./bin/hierarchical-cme -o pancreatic_pwn_r20_e_tau1e-2 -s 100 -t 0.01 -f 20.0 -m e


# Reference solution
python3 scripts/input_generation/set_boolean_pancreatic_cancer.py -pb -r 50
./bin/hierarchical-cme -o pancreatic_pbn_r50_e_tau5e-3 -s 200 -t 0.005 -f 20.0 -m e

python3 scripts/input_generation/set_boolean_pancreatic_cancer.py -pb -r 60
./bin/hierarchical-cme -o pancreatic_pbn_r60_e_tau5e-3 -s 200 -t 0.005 -f 20.0 -m e