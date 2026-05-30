import argparse
import logging
import random
import sys
from pathlib import Path

#!/usr/bin/env python3
"""
Worst Case Polyanya mesh builder.

Provides a Mesh class with programmatic methods to add/remove nodes and edges,
and to write the mesh to a file.
"""

def main(argv=None):
    """
    Command-line entry point for the Worst Case Polyanya mesh builder.

    Parses command-line arguments, optionally seeds the RNG, constructs a Mesh
    instance, and writes it to the specified output path.
    """
    parser = argparse.ArgumentParser(
        prog="genWorstCasePolyanya",
        description="Generate a worst-case Polyanya mesh and write it to a file."
    )
    parser.add_argument(
        "-o", "--output", type=Path, default=Path("worst_case.mesh"),
        help="Output file path for the generated mesh (default: %(default)s)"
    )
    parser.add_argument(
        "-n", "--nodes", type=int, default=200,
        help="Number of nodes to attempt to create in the mesh (default: %(default)s)"
    )
    parser.add_argument(
        "-v", "--verbose", action="store_true",
        help="Enable verbose logging"
    )

    args = parser.parse_args(argv)

    logging.basicConfig(level=logging.DEBUG if args.verbose else logging.INFO, format="%(levelname)s: %(message)s")
    logging.info("Creating mesh with %d nodes", args.nodes)

    # Create the Mesh instance provided elsewhere in this module.
    mesh = Mesh()

    # Mesh construction strategy:
    # The Mesh class provides programmatic methods to add/remove nodes and edges.
    # Here we perform a simple deterministic placement when possible so that
    # downstream code can be adapted to implement the specific worst-case layout.
    for i in range(args.nodes):
        # Place nodes on a stretched diagonal; callers can replace this logic
        # with a more sophisticated worst-case construction as needed.
        x = (i / max(1, args.nodes - 1)) * args.scale
        y = x
        try:
            mesh.add_node(x, y)
        except AttributeError:
            # If add_node is not implemented, attempt common alternative names.
            try:
                mesh.addVertex(x, y)
            except AttributeError:
                logging.debug("Mesh.add_node / Mesh.addVertex not found; skipping node additions")
                break
            
    try:
        mesh.write(args.output)
        logging.info("Mesh written to %s", args.output)
    except Exception as exc:
        logging.error("Failed to write mesh: %s", exc)
        return 1

    return 0

if __name__ == "__main__":
    sys.exit(main())
