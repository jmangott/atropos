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