def ask_cif_path():
    return input("Enter CIF file path: ").strip()

def ask_radius():
    value = input("Particle radius in A (suggested 10-50): ")
    return float(value)


def ask_shell():
    value = input("Shell thickness in A (suggested 2-10): ")
    return float(value)


def ask_disorder():
    value = input("Disorder strength in A (suggested 0.1-0.8): ")
    return float(value)


def ask_seed():
    value = input("Random seed (press Enter for random): ").strip()
    if value == "":
        return None
        
    return int(value)