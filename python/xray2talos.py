import sys
import argparse
import os
import tempfile

def run_jd2_score(filename, output_prefix, bin, database):

    cleaned_pdb_filename = tempfile.mkstemp(prefix=output_prefix)[1]
    os.system("grep ATOM "+filename+" > "+cleaned_pdb_filename)
    
    #First return value of mkstemp is some kind of OS thingy
    silentfile = tempfile.mkstemp(prefix=output_prefix)[1]

    cmd = [bin]
    cmd.extend(["-database", database])
    cmd.extend(["-s", cleaned_pdb_filename])
    cmd.extend(["-out:file:silent", silentfile])
    cmd.extend(["-out:file:silent_struct_type", "protein"])
    cmd.append("-assign_ss")

    os.system(" ".join(cmd))

    sequence = ""
    ss_list  = []
    phi_psi_list = []

    with open(silentfile) as silentf:
        in_body = False
        for line in silentf:
            s_line = line.split()

            if(s_line[0] == "ANNOTATED_SEQUENCE:"):
                in_body = True
                
                #Remove everything in the annotated sequence that comes between square braces
                tmp_seq = s_line[1]
                while(tmp_seq.rfind("[") != -1):
                    tmp_seq = tmp_seq[:tmp_seq.rfind("[")]+tmp_seq[tmp_seq.rfind("]")+1:]

                sequence = tmp_seq
            elif(s_line[0] == "CHAIN_ENDINGS"):
                #Multichain proteins have a "CHAIN_ENDINGS" line after "ANNOTATED_SEQUENCE" that
                #the phi_psi parser below (in the elif) will choke on.
                continue

            elif(in_body):
                seqpos = int(s_line[0])
                ss     = s_line[1]
                ss_list.append(ss)
                phi_psi_list.append((float(s_line[2]), float(s_line[3])))

    os.remove(cleaned_pdb_filename)
    os.remove(silentfile)

    assert len(sequence) > 1
    assert len(sequence) == len(ss_list)
    assert len(phi_psi_list) == len(ss_list)
    return (sequence, "".join(ss_list), phi_psi_list) 

def write_fake_talos_predSS_file(sequence, ss_assignment, original_filename, outfilename):

    lines = ["REMARK  Neural network secondary structure prediction for "+original_filename,
             "REMARK     H-Helix    E-Strand   L-Coil ",
             "",
             "DATA FIRST_RESID 1",
             ""
             ]

    for (i, resletter) in enumerate(sequence):
        if(i % 50 == 0):
            try:
                lines.append("".join(seq_line))
                lines.append("".join(pred_line))
                lines.append("".join(conf_line))
                lines.append("")
            except UnboundLocalError:
                pass

            seq_line  = ["DATA SEQUENCE    "]
            pred_line = ["DATA PREDICTED_SS"]
            conf_line = ["DATA CONFIDENCE  "]
        if(i % 10 == 0):
            seq_line.append(" ")
            pred_line.append(" ")
            conf_line.append(" ")
        
        seq_line.append(resletter)
        pred_line.append(ss_assignment[i])
        conf_line.append("8")

    lines.append("".join(seq_line))
    lines.append("".join(pred_line))
    lines.append("".join(conf_line))
    lines.append('')

    lines.append("VARS RESID RESNAME CS_CNT CS_CNT_R2 Q_H Q_E Q_L CONFIDENCE SS_CLASS")
    lines.append("FORMAT %4d %1s %2d %2d %8.3f %8.3f %8.3f %4.2f %s")
    lines.append("")

    #    1 M  8  3    0.333    0.333    0.333 0.00 L
    for (i, resletter) in enumerate(sequence):

        prob = {"L" : 0.1,
                "H" : 0.1,
                "E" : 0.1}

        prob[ss_assignment[i]] += 0.8
        try:
            if(ss_assignment[i-1] != ss_assignment[i]):
                prob[ss_assignment[i]]   -= 0.25
                prob[ss_assignment[i-1]] += 0.25
        except IndexError:
            pass
        
        try:
            if(ss_assignment[i+1] != ss_assignment[i]):
                prob[ss_assignment[i]]   -= 0.25
                prob[ss_assignment[i+1]] += 0.25
        except IndexError:
            pass
        

        tmp_line = [str(i+1).rjust(5, " "),
                    " "+resletter,
                    " 10",
                    "  4",
                    str(round(prob["H"],3)).ljust(5, "0").rjust(9, " "),
                    str(round(prob["E"],3)).ljust(5, "0").rjust(9, " "),
                    str(round(prob["L"],3)).ljust(5, "0").rjust(9, " "),
                    " 0.85",
                    " "+ss_assignment[i]
                    ]
        
        lines.append("".join(tmp_line))
        
    with open(outfilename, 'w') as f:
        f.write("\n".join(lines))
        f.write("\n")

    return lines

def write_fake_talos_pred_file(sequence, phi_psi, original_filename, outfilename):

    lines = ["DATA FIRST_RESID 1"]

    for (i, resletter) in enumerate(sequence):
        if(i % 50 == 0):
            try:
                lines.append("".join(seq_line))
            except UnboundLocalError:
                pass

            seq_line  = ["DATA SEQUENCE    "]
        if(i % 10 == 0):
            seq_line.append(" ")
        
        seq_line.append(resletter)

    lines.append("".join(seq_line))
    lines.append('')
    lines.append("VARS   RESID RESNAME PHI PSI DPHI DPSI DIST S2 COUNT CS_COUNT CLASS ")
    lines.append("FORMAT %4d %s %8.3f %8.3f %8.3f %8.3f %8.3f %5.3f %2d %2d %s")
    lines.append("")

    for (i, resletter) in enumerate(sequence):
        phi = phi_psi[i][0]
        psi = phi_psi[i][1]

        tmp_line = [str(i+1).rjust(4),
                    " "+resletter,
                    str(phi).rjust(9),
                    str(psi).rjust(9),
                    "   10.000",
                    "   10.000",
                    "   ??.???",
                    " 0.900", " 10", " 13", " Good" ]

        lines.append("".join(tmp_line))
    
    with open(outfilename, 'w') as f:
        f.write("\n".join(lines))

    return lines

def process_command_line(argv):
    
    parser = argparse.ArgumentParser()

    parser.add_argument("pdb_file", nargs=1,
                        help="The pdb file to be converted to TALOS-style prediction files.")
    parser.add_argument("-rosetta_path", default=os.path.expanduser("~")+"/rosetta/",
                        help="Specify the path to a directory with the appropriate rosetta_source and rosetta_database directories")

    args = parser.parse_args(argv[1:])

    #I'm not really sure why this is necessary... apparently argparse returns a list of size 1.
    args.file = args.pdb_file[0]

    return args

def main(argv=None):
    
    args = process_command_line(argv)

    print args.rosetta_path 

    (sequence, ss_assignment, phi_psi) = run_jd2_score(args.file, rosetta_path=args.rosetta_path)

    predSS_outfilename = args.file[0:args.file.rfind(".")]+".predSS.tab"
    write_fake_talos_predSS_file(sequence, ss_assignment, args.file[args.file.find("/")+1:], predSS_outfilename)
    pred_outfilename = args.file[0:args.file.rfind(".")]+".pred.tab"
    write_fake_talos_pred_file(sequence, phi_psi, args.file[args.file.find("/")+1:], pred_outfilename)
   
    return 1

if __name__ == "__main__":
    status = main(sys.argv)
    sys.exit(status)
