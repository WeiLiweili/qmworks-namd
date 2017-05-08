#! /usr/bin/env python
from os.path import join
import h5py

def main(project_name, path_hdf5) 
    with h5py.File(path_hdf5, 'r+') as f5:
        xs = list(filter(lambda x: 'coupling' in x, f5[project_name].keys()))
        paths_css = [join(project_name, x) for x in xs]
        paths_overlaps = [join(project_name,'overlaps_{}/mtx_sji_t0_corrected'.format(i)) for i in range(10000)]
        paths = paths_css + paths_overlaps
        for p in paths:
            if p in f5:
               print(p)
               del f5[p]

        del f5[join(project_name,'swaps')]

def read_cmd_line(parser):
    """
    Parse Command line options.
    """
    args = parser.parse_args()

    attributes = ['pn', 'hdf5' ]

    return [getattr(args, p) for p in attributes]

# ============<>===============
if __name__ == "__main__":

    msg = " remove_couplings -pn <ProjectName> -hdf5 <path/to/hdf5 file> "

    parser = argparse.ArgumentParser(description=msg)
    parser.add_argument('-pn', required=True, help='path to the Hamiltonian files in Pyxaid format')
    parser.add_argument('-hdf5', required=True, help='Index of the first state')

    main(*read_cmd_line(parser))




