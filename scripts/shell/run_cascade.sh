# PS-TTN integrator
python3 scripts/input_generation/set_cascade.py -r 5 -p0
./bin/hierarchical-cme -o cascade_bt_r5_e_tau1e-1 -s 10 -t 0.1 -f 350.0 -m e

python3 scripts/input_generation/set_cascade.py -r 6 -p0
./bin/hierarchical-cme -o cascade_bt_r6_e_tau1e-1 -s 10 -t 0.1 -f 350.0 -m e

python3 scripts/input_generation/set_cascade.py -r 7 -p0
./bin/hierarchical-cme -o cascade_bt_r7_e_tau1e-1 -s 10 -t 0.1 -f 350.0 -m e

python3 scripts/input_generation/set_cascade.py -r 5 -p1
./bin/hierarchical-cme -o cascade_tt_r5_e_tau1e-1 -s 10 -t 0.1 -f 350.0 -m e

python3 scripts/input_generation/set_cascade.py -r 6 -p1
./bin/hierarchical-cme -o cascade_tt_r6_e_tau1e-1 -s 10 -t 0.1 -f 350.0 -m e

python3 scripts/input_generation/set_cascade.py -r 7 -p1
./bin/hierarchical-cme -o cascade_tt_r7_e_tau1e-1 -s 10 -t 0.1 -f 350.0 -m e

python3 scripts/input_generation/set_cascade.py -r 5 -p2
./bin/hierarchical-cme -o cascade_bte_r5_e_tau1e-2 -s 100 -t 0.01 -f 350.0 -m e

python3 scripts/input_generation/set_cascade.py -r 6 -p2
./bin/hierarchical-cme -o cascade_bte_r6_e_tau1e-2 -s 100 -t 0.01 -f 350.0 -m e

python3 scripts/input_generation/set_cascade.py -r 7 -p2
./bin/hierarchical-cme -o cascade_bte_r7_e_tau1e-2 -s 100 -t 0.01 -f 350.0 -m e

python3 scripts/input_generation/set_cascade.py -r 5 -p3
./bin/hierarchical-cme -o cascade_tte_r5_e_tau1e-2 -s 100 -t 0.01 -f 350.0 -m e

python3 scripts/input_generation/set_cascade.py -r 6 -p3
./bin/hierarchical-cme -o cascade_tte_r6_e_tau1e-2 -s 100 -t 0.01 -f 350.0 -m e

python3 scripts/input_generation/set_cascade.py -r 7 -p3
./bin/hierarchical-cme -o cascade_tte_r7_e_tau1e-2 -s 100 -t 0.01 -f 350.0 -m e

# SSA
python3 scripts/reference_solutions/ssa.py --model cascade --nruns 10000
python3 scripts/reference_solutions/ssa.py --model cascade --nruns 100000
python3 scripts/reference_solutions/ssa.py --model cascade --nruns 1000000
python3 scripts/reference_solutions/ssa.py --model cascade --nruns 10000000