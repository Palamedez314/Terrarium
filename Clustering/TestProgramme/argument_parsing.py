from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter

parser = ArgumentParser(formatter_class=ArgumentDefaultsHelpFormatter)
parser.add_argument("datasetname", help="tbd")
parser.add_argument("-d", "--delta", type=float, default=0.05, help="Parameter, that determines the fineness of the underlying lattice")
parser.add_argument("-e", "--eps-factor", type=float, default=3.0, help="tbd")
parser.add_argument("-t", "--tau-factor", type=float, default=2.000001, help="tbd")

args = vars(parser.parse_args())

datasetname = args["datasetname"]
delta = args["delta"]
eps = args["eps_factor"]
tau = args["tau_factor"]
