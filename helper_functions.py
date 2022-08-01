from Bio.Seq import Seq
from Bio.SeqUtils import molecular_weight

def compute_molecular_weight_from_seq(sequence_filename):
    with open(sequence_filename, 'r') as f:
        sequence = f.read()

    # a few string operations
    sequence = str.strip(sequence)
    sequence = str.upper(sequence)

    sequence = Seq(sequence)

    return molecular_weight(sequence, double_stranded = True, circular = True)

def pre_processing(sequence_filename, sample_conc, desired_conc, sample_vol = 1):
    sample_conc = float(sample_conc)
    desired_conc = float(desired_conc)
    mol_weight = 1188128
    molecular_conc = sample_conc * 1/(mol_weight * 1.66054e-15)
    if molecular_conc < desired_conc:
        raise ValueError("Starting concentration (%.2f ng/μL) is less than desired concentration (%.2f), so dilution cannot be performed." % (sample_conc, desired_conc))
    final_volume = (sample_vol * molecular_conc)/desired_conc
    water_to_add = final_volume - sample_vol
    print("To achieve an initial concentration of", desired_conc,"molecules/μL in a final volume of",
        f"{final_volume:.2f}", "μL, add", f"{water_to_add:.2f}", "μL of water to", 
        sample_vol, "μL of sample.", "\n")

pre_processing('test.txt', 43.5, 10**10, sample_vol = 5)