import numpy as np


def getEnergyFromFile(filename, label):
    erg = []
    with open(filename, "r") as f:
        lines = f.readlines()
    for l in lines:
        if label in l:
            erg.append(float(l.split("=")[-1]))

    if len(erg) == 1:
        erg = erg[0]
    return erg


def getLine(ffile, string):
    with open(ffile, "r") as f:
        lines = f.readlines()
    for index, l in enumerate(lines):
        if string in l:
            break
    return lines, index


def getEnergiesOpenMM(ffile, string="Single point energy"):
    lines, index = getLine(ffile, string)
    energies = {}
    for i in range(index + 1, len(lines)):
        if "#---" in lines[i]:
            break
        if "=" not in lines[i]:
            continue
        if "lambda" in lines[i]:
            continue

        s = lines[i].split("=")
        energies[s[0].replace("(kJ/mol)", "").strip()] = eval(s[1])

    return energies


def compareResults(value, expected, threshold=None):
    delta = np.abs(value - expected)
    sumAB = np.abs(value) + np.abs(expected)
    if np.abs(expected) > 0.0:
        ratio = value / expected - 1.0
    else:
        ratio = None
    if sumAB > 0.0:
        frac = 0.5 * np.abs(delta) / sumAB
    else:
        frac = None

    res = {
        "percentage": frac,
        "ratio": ratio,
        "delta": delta,
    }

    if threshold is None:
        threshold = {"percentage": 1e-5, "ratio": 1e-5, "delta": 1e-2}

    failed = False
    for k, v in threshold.items():
        if res[k] is not None and res[k] > v:
            failed = True
            print(f"  Test failed threshold {v}")
            print(f"    Computed value = {value}")
            print(f"    Expected value = {expected}")
            print(f"    {k} = {res[k]}")

    if failed:
        print("  FAILED")
    else:
        print("  PASSED")

    return failed


def run_test_openmm(cmd, key):
    import os
    import yaml

    with open("__test.yaml", "w") as f:
        yaml.dump(cmd, f, default_flow_style=False, indent=4)
    os.system("runOpenMM.py __test.yaml > __test.out")
    return getEnergyFromFile("__test.out", "NonbondedForce (kJ/mol")
