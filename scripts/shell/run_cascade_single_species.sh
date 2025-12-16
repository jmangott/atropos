# PS-TTN integrator

# Binary tree format

python3 scripts/input_generation/set_cascade.py -r 5 -p2
./bin/hierarchical-cme -o cascade_p2_r5_e_tau1e-2 -s 100 -t 0.01 -f 350.0 -m e

python3 scripts/input_generation/set_cascade.py -r 6 -p2
./bin/hierarchical-cme -o cascade_p2_r6_e_tau1e-2 -s 100 -t 0.01 -f 350.0 -m e

python3 scripts/input_generation/set_cascade.py -r 7 -p2
./bin/hierarchical-cme -o cascade_p2_r7_e_tau1e-2 -s 100 -t 0.01 -f 350.0 -m e


# Tensor train format

python3 scripts/input_generation/set_cascade.py -r 5 -p3
./bin/hierarchical-cme -o cascade_p3_r5_e_tau1e-2 -s 100 -t 0.01 -f 350.0 -m e

python3 scripts/input_generation/set_cascade.py -r 6 -p3
./bin/hierarchical-cme -o cascade_p3_r6_e_tau1e-2 -s 100 -t 0.01 -f 350.0 -m e

python3 scripts/input_generation/set_cascade.py -r 7 -p3
./bin/hierarchical-cme -o cascade_p3_r7_e_tau1e-2 -s 100 -t 0.01 -f 350.0 -m e

python3 scripts/input_generation/set_cascade.py -r 5 -p3
./bin/hierarchical-cme -o cascade_p3_r5_e_tau1e-3 -s 1000 -t 0.001 -f 350.0 -m e

python3 scripts/input_generation/set_cascade.py -r 6 -p3
./bin/hierarchical-cme -o cascade_p3_r6_e_tau1e-3 -s 1000 -t 0.001 -f 350.0 -m e

python3 scripts/input_generation/set_cascade.py -r 7 -p3
./bin/hierarchical-cme -o cascade_p3_r7_e_tau1e-3 -s 1000 -t 0.001 -f 350.0 -m e

python3 scripts/input_generation/set_cascade.py -r 10 -p3
./bin/hierarchical-cme -o cascade_p3_r10_e_tau1e-3 -s 1000 -t 0.001 -f 350.0 -m e

python3 scripts/input_generation/set_cascade.py -r 20 -p3
./bin/hierarchical-cme -o cascade_p3_r10_e_tau1e-3 -s 1000 -t 0.001 -f 350.0 -m e
