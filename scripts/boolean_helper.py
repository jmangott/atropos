"""Helper module for Boolean reaction models."""
import bitarray
import bitarray.util
import ctypes
import inspect
import numpy as np
import regex as re
import subprocess
import tempfile
import textwrap

import scripts.reaction_class

def convertRulesToReactions(filename: str):
    """
    Converts a rule file of the Boolean CME integrator of https://bitbucket.org/mprugger/low_rank_cme to a instance of the `ReactionSystem` class.
    NOTE: This method is only capable of converting systems with a maximum of 64 species.
    """
    with open(filename) as f:
        line0 = f.readline()
        f_string = f.read()

    match = re.match(r"RULE_SET\((.*)\)", line0)
    model_info = [x.split()[0] for x in match.group(1).split(",")]

    model_name = model_info[0]
    d = int(model_info[1])
    species_names = [name[1:-1] for name in model_info[2:]]

    rules = {}
    dependencies = {}

    rule_pattern = "template<> bool {}::rule<(\d+)>\(bitset<{}> x\) \{{([\S\s]*?)\}}".format(model_name, d)
    dependency_pattern = "template<> vector<ind> {}::depends_on<(\d+)>\(\) \{{[\s]*?return \{{(.*?)\}};\s*\}}".format(model_name, d)

    rule_matches = re.finditer(rule_pattern, f_string, re.MULTILINE)
    for _, match in enumerate(rule_matches, start=1):
        i = int(match.group(1))
        rules[i] = match.group(2)

    dependency_matches = re.finditer(dependency_pattern, f_string, re.MULTILINE)
    for _, match in enumerate(dependency_matches, start=1):
        i = int(match.group(1))
        dependency = [int(d) for d in match.group(2).split(',')]
        if i not in dependency:
            dependency.append(i)
        dependencies[i] = sorted(dependency)

    with tempfile.TemporaryDirectory() as tmpdirname:
        with open(tmpdirname+"/rule_set_temp.cpp", mode="w") as f:
            str_begin = '#include <bitset>\nextern "C"\n{'
            f.write(str_begin)
            for k, v in rules.items():
                str_rule = "\n\tbool rule_{}(std::bitset<{}> x)\n\t{{{}\t}}\n"
                f.write(str_rule.format(k, d, textwrap.indent(v, "\t")))
            str_end = "}"
            f.write(str_end)

        subprocess.run(["g++",
                        "-fPIC",
                        "-shared",
                        "-o",
                        tmpdirname + "/rule_set_temp.so",
                        tmpdirname + "/rule_set_temp.cpp"])

        handle = ctypes.CDLL(tmpdirname + "/rule_set_temp.so")

    fun_x0 = lambda x: 1 - x
    fun_x1 = lambda x: x
    reactions = []

    for i in range(d):
        d_dep = len(dependencies[i])
        dx_dep = 2**d_dep
        for j in range(dx_dep):
            x = bitarray.bitarray(d, endian="little")
            x.setall(0)
            x_dep = bitarray.util.int2ba(j, length=d_dep, endian="little")
            for k, k_dep in enumerate(dependencies[i]):
                x[k_dep] = x_dep[k]
            curr_rule = handle["rule_{}".format(i)]
            curr_rule.argtypes = [ctypes.c_ulonglong] # avoid overflow for large systems
            x_i_prime = curr_rule(bitarray.util.ba2int(x))
            if x[i] != x_i_prime: # create a reaction only when output is different from input
                nu = np.zeros(d)
                nu[i] = 1 if x[i] == 0 else -1
                propensity = {}
                for k_dep in dependencies[i]:
                    propensity[k_dep] = fun_x0 if x[k_dep] == 0 else fun_x1
                reactions.append(scripts.reaction_class.Reaction(propensity, nu))

    reaction_system = scripts.reaction_class.ReactionSystem(reactions, species_names)

    return reaction_system

if __name__ == "__main__":
    import numpy as np

    reaction_system = convertRulesToReactions("scripts/models/boolean_rulefiles/pancreatic_cancer.hpp")

    for mu, reaction in enumerate(reaction_system.reactions):
        print("reaction {}, product {}".format(mu, int(np.nonzero(reaction.nu)[0])))
        for k, v in reaction.propensity.items():
            funcString = str(inspect.getsourcelines(v)[0])
            funcString = funcString.strip("['\\n']").split(" = ")[1]
            print(k, funcString)
        print(reaction.nu, "\n")