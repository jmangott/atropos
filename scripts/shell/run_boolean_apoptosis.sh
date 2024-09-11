python3 scripts/input_generation/set_boolean_apoptosis.py -pb -r 5
./bin/hierarchical-cme -o apoptosis_pb_r5_e_tau1e-2 -s 100 -t 0.01 -f 20.0 -m e

python3 scripts/input_generation/set_boolean_apoptosis.py -pb -r 10
./bin/hierarchical-cme -o apoptosis_pb_r10_e_tau1e-2 -s 100 -t 0.01 -f 20.0 -m e

python3 scripts/input_generation/set_boolean_apoptosis.py -pb -r 20
./bin/hierarchical-cme -o apoptosis_pb_r20_e_tau1e-2 -s 100 -t 0.01 -f 20.0 -m e


python3 scripts/input_generation/set_boolean_apoptosis.py -pr -r 5
./bin/hierarchical-cme -o apoptosis_pr_r5_e_tau1e-2 -s 100 -t 0.01 -f 20.0 -m e

python3 scripts/input_generation/set_boolean_apoptosis.py -pr -r 10
./bin/hierarchical-cme -o apoptosis_pr_r10_e_tau1e-2 -s 100 -t 0.01 -f 20.0 -m e

python3 scripts/input_generation/set_boolean_apoptosis.py -pr -r 20
./bin/hierarchical-cme -o apoptosis_pr_r20_e_tau1e-2 -s 100 -t 0.01 -f 20.0 -m e


python3 scripts/input_generation/set_boolean_apoptosis.py -pw -r 5
./bin/hierarchical-cme -o apoptosis_pw_r5_e_tau1e-2 -s 100 -t 0.01 -f 20.0 -m e

python3 scripts/input_generation/set_boolean_apoptosis.py -pw -r 10
./bin/hierarchical-cme -o apoptosis_pw_r10_e_tau1e-2 -s 100 -t 0.01 -f 20.0 -m e


# Reference solution
python3 scripts/input_generation/set_boolean_apoptosis.py -pb -r 30
./bin/hierarchical-cme -o apoptosis_pb_r30_e_tau1e-2 -s 100 -t 0.01 -f 20.0 -m e

python3 scripts/input_generation/set_boolean_apoptosis.py -pb -r 50
./bin/hierarchical-cme -o apoptosis_pb_r50_r_tau1e-2 -s 100 -t 0.01 -f 20.0 -m r