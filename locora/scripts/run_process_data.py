from collections import OrderedDict

from locora.process_data.solvent_temporal_analysis import solvent_temporal_analysis


def process_data(option_dict, start, stop):

    transient = list(option_dict['transient'])
    aniso     = option_dict['anisotropic']
    prefix    = str(option_dict['prefix'])

    results=OrderedDict()

    if transient[0]<0 or transient[1]<0:
        raise ValueError("Transient re-entering time must be >=0. transient=0 means \
                         no transient re-entering allowed.")
    if transient[0]>transient[1]:
        raise ValueError("transient[0] must be smaller than transient[1].")

    signs   = ["+", "-"]
    planes  = [[0,1],[0,2],[1,2]]

    p = solvent_temporal_analysis(option_dict, start, stop)

    results['header_trans']  = "# t* "
    results['header_orient'] = "# t* "

    results['scope_trans_avg']  = ""
    results['scope_trans_std']  = ""
    results['scope_orient_avg'] = ""
    results['scope_orient_std'] = ""

    results['scope_trans_all']  = ""
    results['scope_orient_all'] = ""

    for trans in range(transient[0],transient[1]+1):

        results['scope_trans_avg']  += "%d " %trans
        results['scope_trans_std']  += "%d " %trans
        results['scope_orient_avg'] += "%d " %trans
        results['scope_orient_std'] += "%d " %trans

        p.set_transient(trans)
        if p.aniso:
            for sign in signs:
                p.set_sign(sign)
                for plane in planes:
                    p.set_plane(plane)
                    p.select()
                    p._process()
                    if p.plot and trans==transient[0]:
                        p.plot_timeaverage()
                    if trans==transient[0]:
                        results['header_trans']  += p.header_trans
                        results['header_orient'] += p.header_orient

                    results['scope_trans_avg']  += p.scope_trans_avg
                    results['scope_trans_std']  += p.scope_trans_std
                    results['scope_orient_avg'] += p.scope_orient_avg
                    results['scope_orient_std'] += p.scope_orient_std

                ### This is untested...
                for line_trans, line_orient in zip(p.all_trans, p.all_orient):
                    results['scope_trans_all']  += line_trans  + "\n"
                    results['scope_orient_all'] += line_orient + "\n"

        else:
            p.select()
            p._process()
            if p.plot and trans==transient[0]:
                p.plot_timeaverage()

            if trans==transient[0]:
                results['header_trans']  += p.header_trans
                results['header_orient'] += p.header_orient

            results['scope_trans_avg']  += p.scope_trans_avg
            results['scope_trans_std']  += p.scope_trans_std
            results['scope_orient_avg'] += p.scope_orient_avg
            results['scope_orient_std'] += p.scope_orient_std

            for line_trans, line_orient in zip(p.all_trans, p.all_orient):
                results['scope_trans_all']  += "%d " %trans
                results['scope_trans_all']  += line_trans  + "\n"
                results['scope_orient_all'] += "%d " %trans
                results['scope_orient_all'] += line_orient + "\n"

        results['scope_trans_avg']  += "\n"
        results['scope_trans_std']  += "\n"
        results['scope_orient_avg'] += "\n"
        results['scope_orient_std'] += "\n"

    return results


def process_data_out(option_dict, results):

    prefix    = str(option_dict['prefix'])

    if prefix!="":
        prefix+="_"
    trans_avg =open("%strans_avg.txt" %prefix, "w")
    trans_std =open("%strans_std.txt" %prefix, "w")
    trans_all =open("%strans_all.txt" %prefix, "w")
    orient_avg=open("%sorient_avg.txt" %prefix, "w")
    orient_std=open("%sorient_std.txt" %prefix, "w")
    orient_all=open("%sorient_all.txt" %prefix, "w")

    trans_avg.write(results['header_trans']+"\n")
    trans_std.write(results['header_trans']+"\n")
    trans_all.write(results['header_trans']+"\n")
    orient_avg.write(results['header_orient']+"\n")
    orient_std.write(results['header_orient']+"\n")
    orient_all.write(results['header_orient']+"\n")

    trans_avg.write(results['scope_trans_avg'])
    trans_std.write(results['scope_trans_std'])
    trans_all.write(results['scope_trans_all'])
    orient_avg.write(results['scope_orient_avg'])
    orient_std.write(results['scope_orient_std'])
    orient_all.write(results['scope_orient_all'])

    trans_avg.close()
    trans_std.close()
    trans_all.close()
    orient_avg.close()
    orient_std.close()
    orient_all.close()
