import argparse
import json
from .core import canonical_orbit, label, moment_structure

parser = argparse.ArgumentParser(description="Exact (2,2,2) symmetries and observable NPA templates.")
sub = parser.add_subparsers(dest="command", required=True)
p = sub.add_parser("canonical")
p.add_argument("coordinates", nargs=8, help="rational A0 A1 B0 B1 E00 E01 E10 E11")
p = sub.add_parser("template")
p.add_argument("level", help="1, 1+AB, 2, 3, ...")
p = sub.add_parser("bound")
p.add_argument("level")
p.add_argument("coefficients", nargs=8, type=float)
p.add_argument("--probabilities", action="store_true")
p.add_argument("--solver")
args = parser.parse_args()
if args.command == "canonical":
    q, g, orbit, stabilizer = canonical_orbit(args.coordinates)
    print(json.dumps({"representative": list(map(str, q)), "generator_images": g.images,
                      "orbit_size": orbit, "stabilizer_size": stabilizer}, indent=2))
elif args.command == "template":
    basis, keys, entries = moment_structure(args.level)
    print(json.dumps({"convention": "real observable moments; no added probability constraints",
                      "word_basis": list(map(label, basis)), "moment_keys": list(map(label, keys)),
                      "matrix_indices": entries, "normalization": "y[I]=1"}, indent=2))
else:
    from .sdp import bell_bound
    print(json.dumps(bell_bound(args.coefficients, args.level,
                               probability_constraints=args.probabilities, solver=args.solver), indent=2))

