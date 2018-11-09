import multiprocessing
import numpy as np
from collections import OrderedDict

from locora.utils.read_input import read_infile
from locora.scripts.parse_args import parse_args, run_modes

def general_worker(subworker, args, out_q, process_i):

    results = subworker(*args)
    out_q.put([process_i, results])


def do_parallel(args, option_dict, worker, worker_args):

    start  = option_dict['start'][0]
    stop   = option_dict['stop'][0]

    _start = 0
    _stop  = 0
    _steps_p_proc    = int((stop-start)/args.nproc)
    _addsteps_p_proc = int((stop-start)%args.nproc)

    processes = list()
    out_q     = multiprocessing.Queue()
    for i in range(args.nproc):

        _start = _stop
        _stop  = _start + int((stop-start)/args.nproc)

        if _addsteps_p_proc > 0:
            _stop += 1
            _addsteps_p_proc -= 1

        p_worker_args = worker_args + tuple((_start, _stop))
        p             = multiprocessing.Process(target=general_worker,
                                                args=(worker, p_worker_args, out_q, i))
        processes.append(p)
        p.start()

    out_q_list     = list()
    for p in processes:
        out_q_list.append(out_q.get())

    for p in processes:
        p.join()

    results_merged = OrderedDict()
    for i in range(args.nproc):
        for j, results in out_q_list:
            if i != j:
                continue
            for key, value in results.items():
                if i == 0:
                    results_merged[key] = value
                else:
                    results_merged[key] = np.concatenate((results_merged[key], value))

    return results_merged

def main():

    args         = parse_args()
    verbose      = args.verbose
    option_dict  = read_infile(args.input)
    is_multi_run = False
    out_q        = multiprocessing.Queue()

    for key,value in vars(args).items():
        option_dict[key] = value

    if args.nproc > 1:
        is_multi_run = True

    worker       = run_modes[args.mode][0]
    worker_out   = run_modes[args.mode][1]
    has_parallel = bool(run_modes[args.mode][2])
    if verbose:
        print "Input values:"
        print "--"
        for key, value in option_dict.items():
            s=''
            if type(value)==list:
                for s_i in value:
                    s += str(s_i)
                    if type(s_i) != str:
                        s += " "
            else:
                s=value
            print "   %s: %s" %(key, s)
        print "--"
    if is_multi_run:
        if verbose:
            print "Attempting multi-run..."
        if not has_parallel:
            raise Exception("Run mode %s cannot be executed in parallel mode. \
                             Please rerun with -nproc=1 ." %args.mode)
        worker_args = (option_dict,)
        results     = do_parallel(args, option_dict, worker, worker_args)

    else:
        if verbose:
            print "Attempting single-run..."
        worker_args = (option_dict, option_dict['start'][0], option_dict['stop'][0])
        results     = worker(*worker_args)
    if verbose:
        print "Writing output..."
    worker_out(option_dict, results)

def entry_point():

    main()

if __name__ == '__main__':

    entry_point()